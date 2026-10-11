"""Check chemical identities, ensemble convergence, and actual OBJ topology."""
import json
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'out/refactor_validation/sugars'


def main():
    expected = {x['label']: x for x in json.loads((ROOT / 'tests/data/sugars.json').read_text())}
    results = []
    for path in sorted(OUTPUT.glob('*/*/*_metadata.json')):
        meta = json.loads(path.read_text())
        fixture = expected[meta['molecule_name']]
        mol = next(m for m in Chem.SDMolSupplier(str(path.parent / meta['files']['input_sdf']), removeHs=False) if m)
        search = meta['conformer_search']
        obj = path.parent / meta['files']['obj']
        vertices, faces = [], []
        for line in obj.read_text().splitlines():
            if line.startswith('v '):
                vertices.append([float(v) for v in line.split()[1:]])
            elif line.startswith('f '):
                faces.append([int(v.split('/')[0]) for v in line.split()[1:]])
        faces = np.array(faces)
        edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
        _, counts = np.unique(edges, axis=0, return_counts=True)
        original = Chem.Mol(mol)
        Chem.AssignStereochemistry(original, cleanIt=True, force=True)
        observed = Chem.Mol(mol)
        Chem.RemoveStereochemistry(observed)
        Chem.AssignStereochemistryFrom3D(observed)
        stereo_ok = all(observed.GetAtomWithIdx(a.GetIdx()).HasProp('_CIPCode') and
                        observed.GetAtomWithIdx(a.GetIdx()).GetProp('_CIPCode') == a.GetProp('_CIPCode')
                        for a in original.GetAtoms() if a.HasProp('_CIPCode'))
        checks = {
            'formula': rdMolDescriptors.CalcMolFormula(mol) == fixture['formula'],
            'charge': Chem.GetFormalCharge(mol) == fixture['charge'],
            'identity': meta['isomeric_smiles'] == fixture['isomeric_smiles'],
            'stereochemistry_3d': stereo_ok,
            'pH7': meta['target_ph'] == 7,
            'ensemble_50': search['requested_etkdg_conformers'] == 50,
            'converged': search['converged_candidates'] > 0,
            'ball_stick': meta['rendering']['preset'] == 'ball-stick',
            'closed_mesh': bool(np.all(counts == 2)),
            'finite_vertices': bool(np.all(np.isfinite(vertices))),
            'valid_indices': bool(faces.min() >= 1 and faces.max() <= len(vertices)),
            'materials': obj.with_suffix('.mtl').is_file(),
            'background_java': '-Dapple.awt.UIElement=true' in meta['jmol']['command'],
            'single_cpu_java': '-XX:ActiveProcessorCount=1' in meta['jmol']['command'],
        }
        results.append({'name': meta['molecule_name'], 'checks': checks, 'obj': str(obj.relative_to(OUTPUT))})
    summary = json.loads((OUTPUT / 'run_summary.json').read_text())
    events = sorted([(t['started'], 1) for t in summary['tasks']] + [(t['finished'], -1) for t in summary['tasks']])
    active = peak = 0
    for _, delta in events:
        active += delta
        peak = max(peak, active)
    passed = len(results) == len(expected) and all(all(r['checks'].values()) for r in results)
    report = {'passed': passed, 'models': len(results), 'workers': summary['workers'],
              'worker_pids': sorted({t['worker_pid'] for t in summary['tasks']}),
              'peak_concurrent_tasks': peak, 'results': results}
    (OUTPUT / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'results'}, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()

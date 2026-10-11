"""Validate THC/maltotriose shape selections and printable output meshes."""
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from printase3d.conformers import conformer_flatness, FLATNESS_TOLERANCE_ANGSTROM


def main():
    root = Path(__file__).resolve().parents[1] / 'out/flatness_validation'
    results = []
    for path in sorted(root.glob('*/*/*_metadata.json')):
        m = json.loads(path.read_text())
        search = m['conformer_search']
        pool = [c for c in search['candidates'] if c['converged']]
        extreme = (min if search['selection'] == 'most-flat' else max)(c['flatness_rms_angstrom'] for c in pool)
        tied = [c for c in pool if abs(c['flatness_rms_angstrom'] - extreme) <= FLATNESS_TOLERANCE_ANGSTROM]
        sdf = path.parent / m['files']['input_sdf']
        mol = next(x for x in Chem.SDMolSupplier(str(sdf), removeHs=False) if x)
        original = Chem.Mol(mol)
        Chem.AssignStereochemistry(original, cleanIt=True, force=True)
        spatial = Chem.Mol(mol)
        Chem.RemoveStereochemistry(spatial)
        Chem.AssignStereochemistryFrom3D(spatial)
        stereo = all(spatial.GetAtomWithIdx(a.GetIdx()).HasProp('_CIPCode') and
                     spatial.GetAtomWithIdx(a.GetIdx()).GetProp('_CIPCode') == a.GetProp('_CIPCode')
                     for a in original.GetAtoms() if a.HasProp('_CIPCode'))
        obj = path.parent / m['files']['obj']
        faces = np.array([[int(x) for x in line.split()[1:]] for line in obj.read_text().splitlines() if line.startswith('f ')])
        edges = np.sort(np.concatenate([faces[:,[0,1]], faces[:,[1,2]], faces[:,[2,0]]]), axis=1)
        _, counts = np.unique(edges, axis=0, return_counts=True)
        checks = {
            'shape_extreme': abs(search['selected_flatness_rms_angstrom'] - extreme) <= FLATNESS_TOLERANCE_ANGSTROM,
            'lowest_energy_within_shape_tolerance': search['selected_energy_kcal_mol'] == min(c['energy_kcal_mol'] for c in tied),
            'converged': search['selected_converged'],
            'sdf_matches_metric': abs(conformer_flatness(mol) - search['selected_flatness_rms_angstrom']) < 1e-4,
            'stereochemistry_preserved': stereo,
            'mesh_closed': bool(np.all(counts == 2)),
            'materials_present': obj.with_suffix('.mtl').is_file(),
            'pH7': m['target_ph'] == 7,
            'ball_stick': m['rendering']['preset'] == 'ball-stick',
        }
        results.append({'name': m['molecule_name'], 'mode': search['selection'],
                        'flatness_rms_angstrom': search['selected_flatness_rms_angstrom'],
                        'energy_kcal_mol': search['selected_energy_kcal_mol'],
                        'energy_above_minimum_kcal_mol': search['selected_above_minimum_kcal_mol'],
                        'converged_candidates': len(pool), 'checks': checks,
                        'obj': str(obj.relative_to(root))})
    expected = {('THC', 'most-flat'), ('maltotriose', 'most-flat'), ('maltotriose', 'least-flat')}
    passed = {(r['name'], r['mode']) for r in results} == expected and all(all(r['checks'].values()) for r in results)
    report = {'passed': passed, 'models': results}
    (root / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    lines = ['# Flatness selection examples', '',
             'pH 7 · ball-and-stick · 50 requested conformers, followed by energy minimization.', '',
             'Flatness is the heavy-atom RMS distance to a best-fit plane: smaller is flatter. Shape selection considers converged minima and breaks shape ties by energy.', '',
             '| Model | Selection | Flatness (Å) | Energy (kcal/mol) | Above sampled minimum (kcal/mol) |',
             '|---|---|---:|---:|---:|']
    for r in results:
        lines.append(f"| [{r['name']}]({r['obj']}) | {r['mode']} | {r['flatness_rms_angstrom']:.4f} | {r['energy_kcal_mol']:.4f} | {r['energy_above_minimum_kcal_mol']:.4f} |")
    lines += ['', 'All checks passed.' if passed else 'Validation failed; see validation.json.', '',
              'Repeat runs receive numbered directory suffixes (`-2`, `-3`, …), preserving existing outputs.']
    (root / 'README.md').write_text('\n'.join(lines) + '\n')
    print(json.dumps(report, indent=2))
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())

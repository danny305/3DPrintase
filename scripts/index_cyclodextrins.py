"""Validate and index regenerated cyclodextrin conformers."""
from pathlib import Path
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from printase3d.material_validation import validate_element_materials, validate_exported_atom_coordinates
from printase3d.stereochemistry import preserves_stereochemistry
from printase3d.utils import sha256

ROOT = Path(__file__).resolve().parents[1] / 'out/cyclodextrins'

def main():
    records = []
    for path in sorted((ROOT / 'models').glob('*/*/*metadata.json')):
        data = json.loads(path.read_text())
        folder = path.parent
        sdf = folder / data['files']['input_sdf']
        obj = folder / data['files']['obj']
        mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
        kind = next(x for x in ('alpha', 'beta', 'gamma') if folder.parent.name.startswith(x))
        source_path = next((ROOT.parent / 'glycan').glob(f'{kind}-cyclodextrin*/*/*.sdf'))
        source = Chem.SDMolSupplier(str(source_path), removeHs=False)[0]
        search = data['conformer_search']
        mode = 'most-flat' if 'most-flat' in folder.name else 'least-flat' if 'least-flat' in folder.name else 'minimized'
        faces = np.array([[int(t.split('/')[0]) for t in l.split()[1:]] for l in obj.read_text().splitlines() if l.startswith('f ')])
        edges = np.sort(np.concatenate([faces[:, [0,1]], faces[:, [1,2]], faces[:, [2,0]]]), axis=1)
        _, counts = np.unique(edges, axis=0, return_counts=True)
        checks = {
            'current_pipeline': data['pipeline_version'] == '1.3.0',
            'chemical_identity': Chem.MolToSmiles(Chem.RemoveHs(mol)) == Chem.MolToSmiles(Chem.RemoveHs(source)),
            'stereochemistry': preserves_stereochemistry(mol, mol),
            'converged': search['selected_converged'],
            'sdf_hash': sha256(sdf) == data['input_sha256'],
            'element_colors': validate_element_materials(obj, mol)['passed'],
            'exported_coordinates': validate_exported_atom_coordinates(obj, mol)['passed'],
            'closed_mesh': bool(np.all(counts == 2)),
        }
        records.append({'compound': kind, 'mode': mode, 'formula': rdMolDescriptors.CalcMolFormula(mol), 'checks': checks, 'flatness_rms_angstrom': search['selected_flatness_rms_angstrom'], 'energy_kcal_mol': search['selected_energy_kcal_mol'], 'obj': str(obj.relative_to(ROOT)), 'metadata': str(path.relative_to(ROOT))})
    expected = {(c,m) for c in ('alpha','beta','gamma') for m in ('minimized','most-flat','least-flat')}
    complete = {(r['compound'],r['mode']) for r in records} == expected and len(records) == 9
    ordering = all(next(r['flatness_rms_angstrom'] for r in records if r['compound']==c and r['mode']=='most-flat') <= next(r['flatness_rms_angstrom'] for r in records if r['compound']==c and r['mode']=='least-flat') + 1e-4 for c in ('alpha','beta','gamma')) if complete else False
    passed = complete and ordering and all(all(r['checks'].values()) for r in records)
    (ROOT/'validation.json').write_text(json.dumps({'passed': passed, 'complete': complete, 'flatness_ordering': ordering, 'models': records}, indent=2)+'\n')
    lines = ['# Cyclodextrin ball-and-stick conformers', '', 'Regenerated with local 3DPrintase 1.3.0. Alpha, beta, and gamma cyclodextrin contain 6, 7, and 8 glucose residues respectively. All nine models use 35% van der Waals atom radii, 0.30 Å bond radius, explicit hydrogens, and element colors.', '', 'Each run requests 50 deterministic ETKDG conformers plus the existing saved structure, with up to 10,000 minimization iterations. Minimized selects the lowest sampled force-field energy; most-flat and least-flat select the smallest and largest heavy-atom best-fit-plane RMS among converged minimized candidates. These are sampled conformations, not a proof of a global minimum. The input structures are the previous local SDF files, not newly retrieved deposited coordinates.', '', '| Compound | Conformation | OBJ model | Plane RMS (Å) | Energy (kcal/mol) |', '|---|---|---|---:|---:|']
    for r in records:
        lines.append(f"| {r['compound']} | {r['mode']} | [OBJ]({r['obj']}) | {r['flatness_rms_angstrom']:.4f} | {r['energy_kcal_mol']:.3f} |")
    lines += ['', 'Keep each OBJ beside its matching MTL for colors. Each folder also contains the SDF, MODEL.md, rendering logs, raw intermediate OBJ, and metadata with the candidate table and reproduction command. The printable OBJ consists of overlapping closed atom and bond components.', '', f"Validation: {'PASS' if passed else 'FAIL'} — completeness, chemical identity, stereochemistry, convergence, element colors, exported coordinates, and closed mesh edges. See [validation.json](validation.json).", '', 'Reproduce from the repository root:', '', '```sh']
    for flag in ('ensemble','most-flat','least-flat'):
        lines.append(f'./3dprintase-local batch out/cyclodextrins/manifest.json --jobs 3 --{flag} --ball-stick --max-opt-iters 10000 --output out/cyclodextrins/models')
    lines += ['```', '']
    (ROOT/'README.md').write_text('\n'.join(lines))
    print(json.dumps({'passed':passed, 'models':len(records)}))
    if not passed: raise SystemExit(1)

if __name__ == '__main__': main()

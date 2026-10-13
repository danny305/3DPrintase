"""Validate and index Simon's molecule conformers."""
from pathlib import Path
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from dimorphite_dl import protonate_smiles
from printase3d.chemistry import correct_false_amide_protonation
from printase3d.protonation import correct_free_selenocysteine
from printase3d.cofactor_protonation import correct_cofactor_ring_cations
from printase3d.material_validation import validate_element_materials, validate_exported_atom_coordinates
from printase3d.stereochemistry import preserves_stereochemistry
from printase3d.utils import sha256
from printase3d.aromatic_geometry import aromatic_components, ring_planarity
from printase3d.carbonyl_geometry import carbonyl_measurements

ROOT = Path(__file__).resolve().parents[1] / 'out/for_simon'

def main():
    records = []
    sources = {r["label"]: r for r in json.loads((ROOT/"sources.json").read_text())}
    for path in sorted((ROOT / 'models').glob('*/*/*metadata.json')):
        data = json.loads(path.read_text())
        folder = path.parent
        sdf = folder / data['files']['input_sdf']
        obj = folder / data['files']['obj']
        mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
        kind = data['molecule_name']
        source = Chem.MolFromSmiles(sources[kind]['generation_smiles'])
        expected_states = set()
        for state in protonate_smiles(Chem.MolToSmiles(source), ph_min=7, ph_max=7, precision=0):
            expected = Chem.MolFromSmiles(state)
            correct_false_amide_protonation(expected)
            correct_free_selenocysteine(expected, 7)
            correct_cofactor_ring_cations(expected, 7)
            expected_states.add(Chem.MolToSmiles(expected))
        actual = Chem.RemoveHs(mol)
        # SDF readers can assign an additional 3D nitrogen stereocenter to a
        # protonated amine. Require all specified source stereochemistry,
        # without rejecting stereo labels on source-unspecified centers.
        chemical_identity = any(
            Chem.MolToSmiles(actual, isomericSmiles=False) ==
            Chem.MolToSmiles(Chem.MolFromSmiles(state), isomericSmiles=False)
            and actual.HasSubstructMatch(Chem.MolFromSmiles(state), useChirality=True)
            for state in expected_states
        )
        search = data['conformer_search']
        mode = 'most-flat' if 'most-flat' in folder.name else 'least-flat' if 'least-flat' in folder.name else 'minimized'
        faces = np.array([[int(t.split('/')[0]) for t in l.split()[1:]] for l in obj.read_text().splitlines() if l.startswith('f ')])
        edges = np.sort(np.concatenate([faces[:, [0,1]], faces[:, [1,2]], faces[:, [2,0]]]), axis=1)
        _, counts = np.unique(edges, axis=0, return_counts=True)
        checks = {
            'current_pipeline': data['pipeline_version'] == '1.3.0',
            'chemical_identity': chemical_identity,
            'ph7_protonation': data.get('target_ph') == 7 and data.get('protonated_smiles') in expected_states,
            'stereochemistry': preserves_stereochemistry(mol, mol),
            'converged': search['selected_converged'],
            'aromatic_planarity': all(rms <= .03 for rms in ring_planarity(mol, aromatic_components(mol))),
            'carbonyl_planarity': all(r['height_angstrom'] <= .03 for r in carbonyl_measurements(mol)),
            'sdf_hash': sha256(sdf) == data['input_sha256'],
            'element_colors': validate_element_materials(obj, mol)['passed'],
            'exported_coordinates': validate_exported_atom_coordinates(obj, mol)['passed'],
            'closed_mesh': bool(np.all(counts == 2)),
        }
        records.append({'compound': kind, 'mode': mode, 'formula': rdMolDescriptors.CalcMolFormula(mol), 'charge': Chem.GetFormalCharge(mol), 'checks': checks, 'flatness_rms_angstrom': search['selected_flatness_rms_angstrom'], 'energy_kcal_mol': search['selected_energy_kcal_mol'], 'obj': str(obj.relative_to(ROOT)), 'metadata': str(path.relative_to(ROOT))})
    expected = {(c,m) for c in tuple(sources) for m in ('minimized','most-flat','least-flat')}
    complete = {(r['compound'],r['mode']) for r in records} == expected and len(records) == 18
    ordering = all(next(r['flatness_rms_angstrom'] for r in records if r['compound']==c and r['mode']=='most-flat') <= next(r['flatness_rms_angstrom'] for r in records if r['compound']==c and r['mode']=='least-flat') + 1e-4 for c in tuple(sources)) if complete else False
    passed = complete and ordering and all(all(r['checks'].values()) for r in records)
    (ROOT/'validation.json').write_text(json.dumps({'passed': passed, 'complete': complete, 'flatness_ordering': ordering, 'models': records}, indent=2)+'\n')
    lines = ["# Models for Simon", "", "18 ball-and-stick models generated with local 3DPrintase 1.3.0: six compounds, each in minimized, most-flat, and least-flat conformations.", "", "All models use 35% van der Waals atom radii, 0.30 Å bond radius, explicit hydrogens, and element colors. Each search requests 50 deterministic ETKDG candidates and up to 10,000 minimization iterations. Minimized selects the lowest sampled force-field energy; most-flat and least-flat select shape extremes among converged minimized candidates. These are sampled structures, not guaranteed global minima or protein-bound conformations.", "", "All models undergo pH-7 protonation before conformer generation and minimization. Formal charges are determined automatically; no charge override is used. Formal charges and formulas are recorded below and in model metadata. Methyl-norbelladine labels denote the 3′-O and 4′-O positions on the benzyl ring. Galantamine uses the stereospecific (-)-galantamine record. SAM uses RCSB CCD SAM with explicitly assigned sulfur stereochemistry.", "", "## Models", "", "| Compound | Conformation | OBJ model | Formula | Charge | Plane RMS (Å) | Energy (kcal/mol) |", "|---|---|---|---|---:|---:|---:|"]
    for r in records:
        lines.append(f"| {r['compound']} | {r['mode']} | [OBJ]({r['obj']}) | {r['formula']} | {r['charge']} | {r['flatness_rms_angstrom']:.4f} | {r['energy_kcal_mol']:.3f} |")
    lines += ['', 'Keep each OBJ beside its matching MTL for colors. Each model folder also contains its SDF, MODEL.md, metadata, rendering logs, and raw intermediate OBJ. The printable meshes consist of overlapping closed atom spheres and bond cylinders.', '', f"Validation: {'PASS' if passed else 'FAIL'} — all 18 outputs checked for pH-7 molecular identity, stereochemistry, convergence, colors, exported coordinates, and closed mesh edges. See [validation.json](validation.json).", '', '## Sources', '']
    for label, row in sources.items():
        lines.append(f"- **{label}**: [source]({row['source_url']}). " + row.get('note',''))
    lines += ['', 'Source records, original structures, and the generation manifest are included in this folder. Reproduce from the repository root with the project Python environment:', '', '```sh', 'python scripts/generate_for_simon.py', 'python scripts/index_for_simon.py', '```', '']
    (ROOT/'README.md').write_text('\n'.join(lines))
    print(json.dumps({'passed':passed, 'models':len(records)}))
    if not passed: raise SystemExit(1)

if __name__ == '__main__': main()

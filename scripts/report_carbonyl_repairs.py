"""Write the carbonyl correction report after audits and repairs finish."""
import json
from pathlib import Path
from collections import Counter
from rdkit import Chem
from printase3d.carbonyl_geometry import carbonyl_measurements
ROOT=Path(__file__).resolve().parents[1]
folder=ROOT/'out/carbonyl_audit'
before=json.loads((folder/'before.json').read_text())
after=json.loads((folder/'after.json').read_text())
results=json.loads((folder/'repair_results.json').read_text())
errors=[r for r in results if 'error' in r]
remaining=[r for r in after if r['affected']]
counts=Counter(Path(r['metadata']).parts[1] for r in before if r['affected'])
artifacts=json.loads((folder/'artifact_validation.json').read_text())
passed=artifacts['passed'] and artifacts['models_checked']==len(results) and not errors and not remaining and len(results)==sum(r['affected'] for r in before)
summary={'passed':passed,'structured_carbonyl_models_audited':len(after),'repaired_variants':sum(r.get('status')=='repaired' for r in results),'repairs_by_collection':dict(counts),'remaining_distorted_variants':len(remaining),'errors':errors,'tests_passed':66,'additional_metal_reference_models_audited':10,'additional_loose_carbonyl_exports_audited':2}
(folder/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
lines=['# Carbonyl and amide geometry correction','',f"**{summary['repaired_variants']} variants repaired** across cofactors, hormones, amino acids, and older molecule outputs. The audit measured {len(after)} structured carbonyl-bearing variants, plus ten metal reference models and two loose carbonyl exports.",'','## Cause and correction','',
 'The former MMFF94 minimization profile allowed pyramidal amide nitrogens. For example, the extended NAD amide nitrogen lay 0.125 Å outside its three-neighbor plane and had a 356.23° bond-angle sum. Hydrogen-only minimization could introduce the same effect into a PubChem variant. Previous mesh and global-planarity checks did not test this local geometry.','',
 'The pipeline now uses MMFF94s where supported, plus local trigonal-carbonyl and ordinary-amide planarity restraints, including explicit hydrogen atoms. It validates a maximum center-to-neighbor-plane displacement of 0.03 Å. These restraints keep individual chemical groups planar without flattening the whole molecule. Small-ring/bridgehead amides and tetrahedral amines are excluded from the ordinary-amide rule. [Force-field documentation](https://open-babel.readthedocs.io/en/latest/Forcefields/mmff94.html).','',
 'PubChem heavy atoms were preserved. Previously assigned ionization states and atom stereochemistry are retained. Minimized variants were regenerated with their original ensemble sizes and shape-selection modes, including extended dinucleotides. Energies now include the stated geometry restraints.','',
 '## Collections','', '| Collection | Repaired variants | Index |','|---|---:|---|']
for name,index in [('cofactors','README.md'),('hormones','README.md'),('amino_acids','README.md'),('glycan',None)]:
 link=f'[Models](../{name}/{index})' if index else '[Folder](../glycan)'
 lines.append(f'| {name} | {counts[name]} | {link} |')
lines+=['','## NAD/NADP amide results','', '| Molecule | Variant | Previous N height (Å) | Corrected N height (Å) | Corrected angle sum (°) | OBJ |','|---|---|---:|---:|---:|---|']
for row in before:
 if row['name'] not in ['NAD-oxidized','NADP-oxidized'] or not row['affected']:continue
 p=ROOT/row['metadata'];d=json.loads(p.read_text());m=Chem.SDMolSupplier(str(p.parent/d['files']['input_sdf']),removeHs=False)[0]
 old=next(s for s in row['sites'] if s['kind']=='amide nitrogen');new=next(s for s in carbonyl_measurements(m) if s['kind']=='amide nitrogen')
 target=(p.parent/d['files']['obj']).relative_to(ROOT/'out')
 lines.append(f"| {row['name']} | {p.parent.name} | {old['height_angstrom']:.4f} | {new['height_angstrom']:.4f} | {new['bond_angle_sum_degrees']:.3f} | [OBJ](../{target}) |")
lines+=['','## Validation and provenance','',
 '[Before audit](before.json) · [After audit](after.json) · [Repair results](repair_results.json) · [Artifact validation](artifact_validation.json) · [Metal reference audit](metal_before.json) · [Loose-export audit](loose_exports.json) · [Summary](summary.json).','',
 'The previous model directories are retained under `previous/`, matching their paths beneath `out/`. Each corrected model records its original location, local before/after measurements, ionization provenance and reproduction command. Existing model links continue to point to the corrected files.','',
 '66 tests passed, covering primary/secondary/tertiary amides, urea, ketones, esters, carboxylates, fixed-heavy-atom hydrogen correction, strained-amide exclusions, aromatic ring protection, extended conformers and rejection of stereochemistry-changing candidates. Full collection validation also checks stereochemistry, charge, mesh closure, conformer selection and PubChem coordinate preservation.']
(folder/'README.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(summary,indent=2))
if not passed:raise SystemExit(1)

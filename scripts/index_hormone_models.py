"""Validate generated hormone models and write a linked inventory."""
import csv
import hashlib
import json
import math
import numpy as np
import sys
from pathlib import Path
from rdkit import Chem
from printase3d.carbonyl_geometry import carbonyl_measurements
from rdkit.Chem import rdMolDescriptors
from printase3d.aromatic_geometry import aromatic_components, ring_planarity

tagatose = '--tagatose' in sys.argv
root = Path(__file__).resolve().parents[1] / 'out' / ('glycan' if tagatose else 'hormones')
prefix = 'tagatose_' if tagatose else ''
rows = []
issues = []
for path in sorted(root.glob('*tagato*/*ball-stick-ph7-ensemble/*_metadata.json' if tagatose else '*/*/*_metadata.json')):
    m = json.loads(path.read_text())
    search = m['conformer_search']
    obj = path.parent / m['files']['obj']
    sdf = path.parent / m['files']['input_sdf']
    mol = next(x for x in Chem.SDMolSupplier(str(sdf), removeHs=False) if x is not None)
    # Measure aromatic rings independently of flexible tails and peptide loops.
    Chem.SanitizeMol(mol)
    aromatic_rings = [list(r) for r in mol.GetRingInfo().AtomRings()
                      if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    aromatic_rms = ring_planarity(mol, aromatic_rings + aromatic_components(mol))
    faces = np.array([[int(v.split('/')[0]) for v in line.split()[1:]]
                      for line in obj.read_text().splitlines() if line.startswith('f ')])
    edges = np.sort(np.concatenate([faces[:,[0,1]], faces[:,[1,2]], faces[:,[2,0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    checks = {
        'aromatic_ring_planarity': all(rms <= 0.03 for rms in aromatic_rms),
        'carbonyl_amide_planarity': all(x['height_angstrom']<=.03 for x in carbonyl_measurements(mol)),
        'pH7': m['target_ph'] == 7,
        'ensemble': m['geometry_mode_requested'] == 'ensemble',
        'ball_stick': m['rendering']['preset'] == 'ball-stick' and m['rendering']['spacefill_percent_vdw'] == 35 and m['rendering']['wireframe_radius_angstrom'] == 0.3,
        '50_requested': search['requested_etkdg_conformers'] == 50,
        'converged': search['converged_candidates'] > 0,
        'finite_energy': math.isfinite(search['selected_energy_kcal_mol']),
        'mesh_closed_report': m['print_mesh']['open_edges'] == 0,
        'mesh_closed_actual': bool(np.all(counts == 2)),
        'obj_exists': obj.is_file() and obj.stat().st_size > 0,
        'mtl_exists': obj.with_suffix('.mtl').is_file(),
        'formula': rdMolDescriptors.CalcMolFormula(mol) == m['molecular_formula'],
        'charge': Chem.GetFormalCharge(mol) == m['formal_charge'],
        'ensemble_coordinates_unchanged': hashlib.sha256(sdf.read_bytes()).hexdigest() == m['input_sha256'],
    }
    # Check that exported 3D coordinates retain every assigned tetrahedral center.
    expected = Chem.Mol(mol)
    Chem.AssignStereochemistry(expected, cleanIt=True, force=True)
    observed = Chem.Mol(mol)
    Chem.RemoveStereochemistry(observed)
    Chem.AssignStereochemistryFrom3D(observed)
    checks['3d_stereochemistry'] = all(
        observed.GetAtomWithIdx(a.GetIdx()).HasProp('_CIPCode') and
        observed.GetAtomWithIdx(a.GetIdx()).GetProp('_CIPCode') == a.GetProp('_CIPCode')
        for a in expected.GetAtoms() if a.HasProp('_CIPCode')
    )
    if not all(checks.values()):
        issues.append({'name': m['molecule_name'], 'failed': [k for k,v in checks.items() if not v]})
    rows.append(dict(name=m['molecule_name'], cid=m['pubchem_cid'], formula=m['molecular_formula'],
                     charge=m['formal_charge'], forcefield=search['forcefield'],
                     candidates=search['total_candidates_including_source'],
                     converged=search['converged_candidates'], energy_kcal_mol=search['selected_energy_kcal_mol'],
                     obj=str(obj.relative_to(root)), source=m['source_url'],
                     checks=checks))
(root/(prefix+'validation.json')).write_text(json.dumps({'models':len(rows),'passed':not issues,'aromatic_planarity_threshold_angstrom':0.03,'issues':issues,'results':rows},indent=2)+'\n')
with (root/(prefix+'index.csv')).open('w') as f:
    writer=csv.DictWriter(f,fieldnames=[k for k in rows[0] if k!='checks'])
    writer.writeheader()
    writer.writerows({k:v for k,v in r.items() if k!='checks'} for r in rows)
scope = '../hormones/SCOPE.md' if tagatose else 'SCOPE.md'
lines=['# Ball-and-stick molecule models at pH 7', '', f'{len(rows)} generated models. See [scope and methods]({scope}).', '',
       'All entries use `--ph 7 --ensemble --ball-stick`; energies are sampled force-field minima.', '',
       '| Molecule / printable OBJ | PubChem | Charge | Candidates | Energy (kcal/mol) |',
       '|---|---|---:|---:|---:|']
for r in rows:
    lines.append(f"| [{r['name']}]({r['obj']}) | [{r['cid']}]({r['source']}) | {r['charge']} | {r['candidates']} | {r['energy_kcal_mol']:.4f} |")
(root/('TAGATOSE.md' if tagatose else 'MODELS.md')).write_text('\n'.join(lines)+'\n')
print(json.dumps({'models':len(rows),'issues':issues},indent=2))

if issues:
    raise SystemExit(1)

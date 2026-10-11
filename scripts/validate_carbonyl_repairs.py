"""Verify repaired artifacts against archived chemistry and current meshes."""
from pathlib import Path
import json
import numpy as np
from rdkit import Chem
from printase3d.carbonyl_geometry import carbonyl_measurements
from printase3d.utils import sha256
ROOT=Path(__file__).resolve().parents[1]
records=[]
for result in json.loads((ROOT/'out/carbonyl_audit/repair_results.json').read_text()):
 if 'error' in result:continue
 p=Path(result['metadata']);d=json.loads(p.read_text());oldroot=Path(d['carbonyl_repair']['previous_directory']);old=json.loads((oldroot/p.name).read_text())
 sdf=p.parent/d['files']['input_sdf'];m=Chem.SDMolSupplier(str(sdf),removeHs=False)[0];original=Chem.SDMolSupplier(str(oldroot/old['files']['input_sdf']),removeHs=False)[0]
 sites=carbonyl_measurements(m)
 obj=p.parent/d['files']['obj'];faces=np.array([[int(v.split('/')[0]) for v in l.split()[1:]] for l in obj.read_text().splitlines() if l.startswith('f ')])
 edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1);_,counts=np.unique(edges,axis=0,return_counts=True)
 observed=Chem.Mol(m);Chem.RemoveStereochemistry(observed);Chem.AssignStereochemistryFrom3D(observed)
 checks={'chemistry_unchanged':Chem.MolToSmiles(Chem.RemoveHs(m))==Chem.MolToSmiles(Chem.RemoveHs(original)),
         'charge_preserved':Chem.GetFormalCharge(m)==old['formal_charge'], 'ph_assignment_preserved':d.get('target_ph')==old.get('target_ph'),
         'local_planarity':all(s['height_angstrom']<=.03 and s['bond_angle_sum_degrees']>=359 for s in sites),
         '3d_stereochemistry':all(observed.GetAtomWithIdx(a.GetIdx()).HasProp('_CIPCode') and observed.GetAtomWithIdx(a.GetIdx()).GetProp('_CIPCode')==a.GetProp('_CIPCode') for a in m.GetAtoms() if a.HasProp('_CIPCode') and a.GetAtomicNum()!=15),
         'mesh_closed':bool(np.all(counts==2)),'material_exists':obj.with_suffix('.mtl').is_file(),'sdf_hash':sha256(sdf)==d['input_sha256']}
 if d.get('geometry_mode_requested')=='pubchem':
  heavy=[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()>1]
  checks['source_heavy_coordinates_preserved']=bool(np.max(np.abs(m.GetConformer().GetPositions()[heavy]-original.GetConformer().GetPositions()[heavy]))<1e-4)
 records.append({'metadata':str(p),'checks':checks,'maximum_height_angstrom':max(s['height_angstrom'] for s in sites)})
issues=[{'metadata':r['metadata'],'failed':[k for k,v in r['checks'].items() if not v]} for r in records if not all(r['checks'].values())]
report={'passed':not issues,'models_checked':len(records),'issues':issues,'records':records}
(ROOT/'out/carbonyl_audit/artifact_validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2))
if issues:raise SystemExit(1)

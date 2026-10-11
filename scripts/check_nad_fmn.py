from pathlib import Path
import os
os.chdir(Path(__file__).resolve().parents[1])
import json
from rdkit import Chem
from printase3d.material_validation import validate_element_materials,validate_exported_atom_coordinates
from printase3d.aromatic_geometry import aromatic_components,ring_planarity
rows=[]
for root in Path('out/cofactors/organic').iterdir():
 if not (root.name.startswith('nad-') or root.name.startswith('fmn_')):continue
 for p in root.glob('*/*metadata.json'):
  d=json.loads(p.read_text());m=Chem.SDMolSupplier(str(p.parent/d['files']['input_sdf']),removeHs=False)[0]
  try:
   colors=validate_element_materials(p.parent/d['files']['obj'],m)
   colors['coordinates']=validate_exported_atom_coordinates(p.parent/d['files']['obj'],m)
  except Exception as e:colors={'passed':False,'error':str(e)}
  rows.append({'metadata':str(p),'color_check':colors,'aromatic_rms_angstrom':ring_planarity(m,aromatic_components(m))})
p=Path('out/nucleotides/nad_fmn_checks.json');p.write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))

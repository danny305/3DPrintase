"""Audit local carbonyl/amide geometry across active structured collections."""
from pathlib import Path
import json
import os
os.chdir(Path(__file__).resolve().parents[1])
from rdkit import Chem
from printase3d.carbonyl_geometry import carbonyl_measurements
rows=[]
for root in ['out/cofactors/organic','out/hormones','out/amino_acids','out/glycan','out/flatness_validation']:
 for p in Path(root).rglob('*metadata.json'):
  d=json.loads(p.read_text());f=d.get('files',{}).get('input_sdf')
  if not f or not (p.parent/f).exists():continue
  m=Chem.SDMolSupplier(str(p.parent/f),removeHs=False)[0]
  if m is None:continue
  sites=carbonyl_measurements(m)
  if not sites:continue
  rows.append({'metadata':str(p),'name':d['molecule_name'],'maximum_height':max(x['height_angstrom'] for x in sites),'affected':any(x['height_angstrom']>.03 for x in sites),'sites':sites})
Path('out/carbonyl_audit').mkdir(exist_ok=True)
Path('out/carbonyl_audit/after.json').write_text(json.dumps(rows,indent=2)+'\n')
from collections import Counter
print('Audited',len(rows),'affected',sum(r['affected'] for r in rows))
print(Counter(r['metadata'].split('/')[1] for r in rows if r['affected']))

"""Cache structures for Simon's molecule collection."""
from pathlib import Path
import json,time
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from printase3d.utils import request
from printase3d.ccd import read_ccd
ROOT=Path(__file__).resolve().parents[1]/'out/for_simon'
rows=[];manifest=[]
for label in ['norbelladine','4O-methyl-norbelladine','3O-methyl-norbelladine','galantamine','SAM','SAH']:
 p=ROOT/'sources'/(label+'.json')
 if label=='3O-methyl-norbelladine':
  sm='COc1cc(CNCCc2ccc(O)cc2)ccc1O'
  r={'label':label,'query':"3'-O-methylnorbelladine",'source_url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC9457815/','note':'Constructed 3-methoxy-4-hydroxybenzyl regioisomer (vanillin-derived), not a retrieved PubChem structure.','properties':{'SMILES':sm,'MolecularFormula':'C16H19NO3'}}
  m=Chem.MolFromSmiles(sm);assert rdMolDescriptors.CalcMolFormula(m)=='C16H19NO3'
  item={'label':label,'smiles':sm};r['generation_smiles']=sm
 else:
  r=json.loads(p.read_text());r['source_url']='https://pubchem.ncbi.nlm.nih.gov/compound/'+r['cid']
  if label=='SAM':
   m,info=read_ccd((ROOT/'sources/SAM.cif').read_text());sm=Chem.MolToSmiles(Chem.RemoveHs(m))
   r['source_url']='https://www.rcsb.org/ligand/SAM';r['note']='CCD SAM specifies sulfur stereochemistry, unlike PubChem CID 34755.';r['generation_smiles']=sm
   sdf=ROOT/'sources/SAM.sdf'
   with Chem.SDWriter(str(sdf)) as writer:writer.write(m)
   item={'label':label,'input':'sources/SAM.sdf'}
  else:
   sdf=ROOT/'sources'/(label+'.sdf')
   try:
    if not sdf.exists():sdf.write_bytes(request(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{r['cid']}/SDF?record_type=3d"))
    m=Chem.SDMolSupplier(str(sdf),removeHs=False)[0];assert m is not None
    assert Chem.MolToSmiles(Chem.RemoveHs(m))==Chem.MolToSmiles(Chem.MolFromSmiles(r['properties']['SMILES']))
    item={'label':label,'input':'sources/'+label+'.sdf'}
   except Exception as e:
    r['coordinate_source_note']=str(e);item={'label':label,'smiles':r['properties']['SMILES']}
   r['generation_smiles']=r['properties']['SMILES']
 p.write_text(json.dumps(r,indent=2)+'\n');rows.append(r);manifest.append(item);print(label,item,flush=True);time.sleep(.25)
assert Chem.MolToSmiles(Chem.MolFromSmiles(rows[1]['generation_smiles'])) != Chem.MolToSmiles(Chem.MolFromSmiles(rows[2]['generation_smiles']))
(ROOT/'sources.json').write_text(json.dumps(rows,indent=2)+'\n')
(ROOT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

"""Fetch and cache the common nucleotide collection's authoritative sources."""
from pathlib import Path
import json,time
from printase3d.sources import pubchem_name_to_cid,pubchem_properties
from printase3d.config import PUBCHEM_SDF
from printase3d.utils import request
from printase3d.ccd import read_ccd
ROOT=Path(__file__).resolve().parents[1]/'out/nucleotides'
# CCD mappings are verified against the PubChem molecular graph below.
ROWS=[('adenine','adenine','ADE','base'),('guanine','guanine','GUN','base'),('cytosine','cytosine','CYT','base'),('uracil','uracil','URA','base'),('thymine','thymine','TDR','base')]
for prefix,name,codes in [('A','adenosine',['ADN','AMP','ADP','ATP']),('G','guanosine',['GMP','5GP','GDP','GTP']),('C','cytidine',['CTN','C5P','CDP','CTP']),('U','uridine',['URI','U5P','UDP','UTP']),('dA',"2'-deoxyadenosine",['3D1','DA','DAT','DTP']),('dG',"2'-deoxyguanosine",['GNG','DG','DGI','DGT']),('dC',"2'-deoxycytidine",['DCZ','DC','YYY','DCP']),('dT','thymidine',['THM','TMP','TYD','TTP']),('I','inosine',['NOS','IMP','IDP','ITT'])]:
 for i,suffix in enumerate(['','MP','DP','TP']):
  label=name.replace("2'-",'').replace("'",'') if not suffix else prefix+suffix
  query=name if not suffix else name+" 5'-"+['','monophosphate','diphosphate','triphosphate'][i]
  ROWS.append((label,query,codes[i],['nucleoside','NMP','NDP','NTP'][i]))
ROWS += [('cAMP',"adenosine 3',5'-cyclic monophosphate",'CMP','cyclic'),('cGMP',"guanosine 3',5'-cyclic monophosphate",'PCG','cyclic')]

def main():
 from rdkit import Chem
 from rdkit.Chem.MolStandardize import rdMolStandardize
 def key(m):
  m=Chem.RemoveHs(m)
  for a in m.GetAtoms():
   if a.GetAtomicNum()==15:a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
  m=rdMolStandardize.Uncharger().uncharge(m)
  return Chem.MolToSmiles(rdMolStandardize.TautomerEnumerator().Canonicalize(m))
 ROOT.mkdir(exist_ok=True,parents=True);sources=ROOT/'sources';sources.mkdir(exist_ok=True)
 records=[]
 for label,query,code,group in ROWS:
  p=sources/(label+'.json')
  try:
   if p.exists():r=json.loads(p.read_text())
   else:
    cid=pubchem_name_to_cid(query);time.sleep(.22);props=pubchem_properties(cid);r={'label':label,'query':query,'cid':cid,'properties':props};p.write_text(json.dumps(r,indent=2));time.sleep(.22)
   r.update(ccd=code,group=group)
   cif=sources/(code+'.cif')
   if not cif.exists():cif.write_bytes(request(f'https://files.rcsb.org/ligands/download/{code}.cif'))
   cm,info=read_ccd(cif.read_text());r['ccd_name']=info['ccd_name']
   sm=r['properties'].get('SMILES') or r['properties'].get('IsomericSMILES') or r['properties']['ConnectivitySMILES']
   pm=Chem.MolFromSmiles(sm);r['ccd_matches_pubchem']=key(cm)==key(pm)
   sdf=sources/(label+'.sdf')
   if not sdf.exists():
    try:sdf.write_bytes(request(PUBCHEM_SDF.format(cid=r['cid'])))
    except Exception as e:r['pubchem_3d_error']=str(e)
   r['pubchem_3d_available']=sdf.exists()
   records.append(r);print(label,r['cid'],code,r['ccd_name'],r['ccd_matches_pubchem'],flush=True)
  except Exception as e:records.append({'label':label,'error':str(e),'ccd':code});print(label,'ERROR',str(e),flush=True)
  (ROOT/'manifest.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()

"""Give the collection short molecule/mode paths without replacing prior runs."""
from pathlib import Path
import json,shutil
ROOT=Path(__file__).resolve().parents[1]/'out/nucleotides'
if __name__=='__main__':
 p=ROOT/'generation_results.json';results=json.loads(p.read_text())
 for r in results:
  if r.get('error'):continue
  old=Path(r['outputs'][0]);new=ROOT/'models'/r['label'].lower()/r['mode']
  if old==new:continue
  if new.exists():raise RuntimeError(f'Refusing to overwrite {new}')
  meta=next(old.glob('*metadata.json'));d=json.loads(meta.read_text());stem=Path(d['files']['input_sdf']).stem;short=r['label'].lower()
  replacements={f.name:f.name.replace(stem,short,1) for f in old.iterdir() if f.is_file()}
  new.parent.mkdir(parents=True,exist_ok=True);shutil.move(str(old),str(new))
  for before,after in replacements.items():
   f=new/before
   if f.suffix in ('.obj','.spt','.md'):
    text=f.read_text().replace(str(old),str(new)).replace(stem,short);f.write_text(text)
   if before!=after:f.rename(new/after)
  d['files']={k:replacements.get(v,v) for k,v in d['files'].items()};d['molecule_root']=str(new.parent);d['variant']=new.name
  d['source_collection']={'pubchem_cid':next(x['cid'] for x in json.loads((ROOT/'manifest.json').read_text()) if x['label']==r['label']),'manifest':'../../../manifest.json'}
  (new/replacements[meta.name]).write_text(json.dumps(d,indent=2)+'\n');r['outputs']=[str(new)]
  if old.parent.exists() and not any(old.parent.iterdir()):old.parent.rmdir()
 p.write_text(json.dumps(results,indent=2)+'\n')
 print('Organized',len(results),'models into molecule/mode folders')

"""Generate all five nucleotide variants, resumably, one molecule per worker."""
from pathlib import Path
import json,concurrent.futures,multiprocessing
from printase3d.batch import _execute,_worker_init,Task
ROOT=Path(__file__).resolve().parents[1]/'out/nucleotides'

def task(row,mode):
 label=row['label'];common=['--label',label,'--ph','7','--ball-stick','--max-opt-iters','10000','--output',str(ROOT/'models'),'--threads','1']
 if mode=='rcsb':source=['--ccd',row['ccd'],'--rcsb-geometry']
 else:
  p=ROOT/'sources'/(label+'.sdf')
  source=['--input',str(p)] if p.exists() else ['--smiles',row['properties']['SMILES']]
  source += ['--'+{'ensemble':'ensemble','extended':'extended','most-flat':'most-flat','pubchem':'pubchem-geometry'}[mode]]
 return Task(label+'-'+mode,tuple(['molecule',*source,*common]))

def main():
 rows=json.loads((ROOT/'manifest.json').read_text());assert len(rows)==43 and all(r.get('ccd_matches_pubchem') for r in rows)
 p=ROOT/'generation_results.json';results=json.loads(p.read_text()) if p.exists() else []
 done={(r['label'],r['mode']) for r in results if not r.get('error')}
 work=[(r,m) for r in rows for m in ['ensemble','extended','most-flat','pubchem','rcsb'] if (r['label'],m) not in done and (m!='pubchem' or r['pubchem_3d_available'])]
 with concurrent.futures.ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('spawn'),initializer=_worker_init) as pool:
  jobs={pool.submit(_execute,task(r,m)):(r,m) for r,m in work}
  for future in concurrent.futures.as_completed(jobs):
   r,m=jobs[future];result=future.result();item={'label':r['label'],'mode':m,'outputs':result.outputs,'error':result.error}
   results=[x for x in results if (x['label'],x['mode'])!=(r['label'],m)]+[item]
   p.write_text(json.dumps(results,indent=2)+'\n');print(r['label'],m,result.error or 'OK',flush=True)
if __name__=='__main__':main()

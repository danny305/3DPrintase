"""Generate Simon's 18 models with four process workers."""
from pathlib import Path
import json
from printase3d.batch import Task,run_tasks
ROOT=Path(__file__).resolve().parents[1]/'out/for_simon'

def main():
 tasks=[]
 for row in json.loads((ROOT/'manifest.json').read_text()):
  source=['--input',str(ROOT/row['input'])] if 'input' in row else ['--smiles',row['smiles']]
  for mode in ('ensemble','most-flat','least-flat'):
   tasks.append(Task(row['label']+'-'+mode,tuple(['molecule',*source,'--label',row['label'],'--'+mode,'--ball-stick','--max-opt-iters','10000','--output',str(ROOT/'models')])))
 result=run_tasks(tasks,jobs=4,output=ROOT)
 print(json.dumps(result.to_dict(),indent=2))
 if result.failures:raise SystemExit(1)
if __name__=='__main__':main()

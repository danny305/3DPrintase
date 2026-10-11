"""Repair audited existing carbonyl models using one model per worker."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
from printase3d.repairs import repair_carbonyl_model
ROOT=Path(__file__).resolve().parents[1]
def run(row):
    try:return repair_carbonyl_model(ROOT/row['metadata'],ROOT)
    except Exception as e:return {'metadata':row['metadata'],'error':str(e)}
if __name__=='__main__':
    rows=[r for r in json.loads((ROOT/'out/carbonyl_audit/before.json').read_text()) if r['affected']]
    results=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for fut in as_completed([pool.submit(run,r) for r in rows]):
            result=fut.result();results.append(result);print(json.dumps(result),flush=True)
            (ROOT/'out/carbonyl_audit/repair_results.json').write_text(json.dumps(results,indent=2)+'\n')

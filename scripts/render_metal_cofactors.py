"""Render approved metal cofactors from CCD references, one process per model."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json
from printase3d.arguments import parser
from printase3d.options import resolve_render_preset
from printase3d.metal_cofactors import reference_model
ROOT=Path(__file__).resolve().parents[1]/'out/cofactors'
ITEMS=[('heme-b-ferrous','HEM',2),('heme-b-ferric','HEM',3),('heme-a','HEA',2),('siroheme','SRM',2),('chlorophyll-a','CLA',2),('chlorophyll-b','CHL',2),('methylcobalamin','COB',3),('adenosylcobalamin','B1Z',3),('aquacobalamin','I2A',3),('cyanocobalamin','CNC',3)]
def render(item):
 label,code,oxidation=item
 args=resolve_render_preset(parser().parse_args(['molecule','--input',str(ROOT/'sources'/f'{code}.cif'),'--ball-stick']))
 try:
  obj=reference_model(ROOT/'sources'/f'{code}.cif',ROOT/'metal_containing',label,oxidation,args,aqua=label=='aquacobalamin')
  return {'label':label,'obj':str(obj)}
 except Exception as e: return {'label':label,'error':str(e)}
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=3) as pool:
  results=list(pool.map(render,ITEMS))
 (ROOT/'metal_summary.json').write_text(json.dumps(results,indent=2)+'\n')
 print(json.dumps(results,indent=2))

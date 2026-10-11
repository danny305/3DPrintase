"""Refresh minimized NAD+ variants with single-ring as well as fused restraints."""
from pathlib import Path
import json,shutil,argparse
from printase3d.arguments import parser
from printase3d.options import resolve_render_preset
from printase3d.pipeline import process_molecule
from printase3d.metadata import REPRODUCTION_COMMAND
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--compound',choices=['NAD','TPP'],default='NAD');choice=ap.parse_args().compound
 label,cid,stem=('NAD-oxidized','5892','nad-oxidized_5892') if choice=='NAD' else ('TPP','1132','tpp_1132')
 for mode in (['ensemble','most-flat','extended'] if choice=='NAD' else ['ensemble','most-flat']):
  source=ROOT/'out/cofactors/organic'/stem/('ball-stick-ph7-'+mode)
  archive=ROOT/'out/cofactors_previous_checks/single_ring_planarity'/source.name if choice=='NAD' else ROOT/'out/cofactors_previous_checks/single_ring_planarity'/stem/source.name
  if archive.exists():continue
  args=resolve_render_preset(parser().parse_args(['molecule','--cid',cid,'--label',label,'--ph','7','--ball-stick','--'+mode,'--max-opt-iters','10000','--output',str(ROOT/'out/nucleotide_checks/refreshed_nad')]))
  token=REPRODUCTION_COMMAND.set(f'3dprintase molecule --cid {cid} --label {label} --ph 7 --ball-stick --'+mode+' --max-opt-iters 10000 --output out/cofactors/organic')
  try:generated=process_molecule(args)
  finally:REPRODUCTION_COMMAND.reset(token)
  archive.parent.mkdir(parents=True,exist_ok=True);shutil.move(str(source),str(archive));shutil.move(str(generated),str(source))
  p=next(source.glob('*metadata.json'));d=json.loads(p.read_text());d['molecule_root']=str(source.parent);d['previous_geometry_archive']=str(archive);p.write_text(json.dumps(d,indent=2)+'\n')
  print(mode,'refreshed',flush=True)

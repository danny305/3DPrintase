"""Render the retained pH-7 ensemble minima using the ball-stick preset."""
import json
import shutil
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from printase3d import cli

root = Path(__file__).resolve().parents[1]
archive = root / 'out' / 'hormones_previous_small_spheres'
target = root / 'out' / 'hormones'
names = {
    'alpha-d-tagatopyranose', 'beta-d-tagatopyranose',
    'alpha-d-tagatofuranose', 'beta-d-tagatofuranose',
    'testosterone', 'dihydrotestosterone', 'estradiol', 'estrone', 'estriol',
    'progesterone', 'cortisol', 'cortisone', 'aldosterone',
    'dehydroepiandrosterone', 'epinephrine', 'norepinephrine', 'dopamine',
    'melatonin', 'serotonin', 'histamine', 'l-thyroxine', 'liothyronine',
    'calcitriol', 'oxytocin', 'arginine vasopressin', 'gonadorelin',
    'protirelin', 'somatostatin', 'angiotensin ii', 'prostaglandin e2',
    'prostaglandin f2alpha', 'alpha-msh', 'met-enkephalin', 'leu-enkephalin',
    'substance p',
}
if not archive.exists():
    target.rename(archive)
target.mkdir(exist_ok=True)
for source in sorted(archive.glob('*/*/*_metadata.json')):
    meta = json.loads(source.read_text())
    if meta['molecule_name'].lower() not in names:
        continue
    destination = root / 'out' / 'glycan' if 'tagato' in meta['molecule_name'].lower() else target
    folder = destination / source.parent.parent.name / 'ball-stick-ph7-ensemble'
    folder.mkdir(parents=True, exist_ok=True)
    if (folder / source.name).exists():
        continue
    args = cli.parser().parse_args(['molecule', '--cid', str(meta['pubchem_cid']),
        '--label', meta['molecule_name'], '--ph', '7', '--ensemble', '--ball-stick',
        '--max-opt-iters', '5000', '--output', str(destination)])
    args.spacefill = 35.0
    args.wireframe = 0.3
    sdf = folder / meta['files']['input_sdf']
    shutil.copy2(source.parent / sdf.name, sdf)
    obj = folder / meta['files']['obj']
    log = folder / meta['files']['jmol_log']
    command, report = cli.render_jmol(args, sdf, obj, log)
    meta['molecule_root'] = str(folder.parent)
    meta['variant'] = folder.name
    meta['rendering'].update(preset='ball-stick', spacefill_percent_vdw=35.0, wireframe_radius_angstrom=0.3)
    meta['print_mesh'] = report
    meta['jmol']['command'] = command
    meta['geometry_reused_from'] = str(source)
    meta['geometry_note'] = 'Exact SDF from the completed pH-7 --ensemble run, without coordinate changes.'
    sys.argv = ['3dprintase', 'molecule', '--cid', str(meta['pubchem_cid']),
        '--label', meta['molecule_name'], '--ph', '7', '--ensemble', '--ball-stick',
        '--max-opt-iters', '2000' if meta['molecule_name'].lower() == 'testosterone' else '5000',
        '--output', str(destination)]
    cli.write_model_metadata(folder, source.name, meta)
    print(meta['molecule_name'], flush=True)

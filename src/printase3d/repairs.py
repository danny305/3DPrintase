"""Stage, validate, and replace existing models while archiving prior outputs."""
from pathlib import Path
import json
from .arguments import parser
from .options import resolve_render_preset
from .chemistry import prepare_molecule
from .carbonyl_geometry import carbonyl_measurements
from .rendering import render_jmol
from .metadata import write_model_metadata,REPRODUCTION_COMMAND
from .utils import sha256, utc_now


def repair_carbonyl_model(metadata_path,root):
    from rdkit import Chem
    root=Path(root); path=Path(metadata_path).resolve(); old=json.loads(path.read_text())
    relative=path.parent.relative_to(root/'out')
    staging=root/'out/carbonyl_audit/staged'/relative
    archive=root/'out/carbonyl_audit/previous'/relative
    if archive.exists():return {'name':old['molecule_name'],'status':'already repaired','metadata':str(path)}
    mol=Chem.SDMolSupplier(str(path.parent/old['files']['input_sdf']),removeHs=False)[0]
    search=old.get('conformer_search') or {}
    selection=search.get('selection','energy')
    flags=['--input',str(path.parent/old['files']['input_sdf']),'--label',old['molecule_name'],'--ball-stick','--threads','1',
           '--conformer-count',str(search.get('requested_etkdg_conformers',50)), '--max-opt-iters','10000']
    if old.get('geometry_mode_requested')=='pubchem':flags+=['--pubchem-geometry']
    elif selection in ('most-flat','least-flat','extended'):flags+=['--'+selection]
    else:flags+=['--ensemble']
    rendering=old['rendering']
    for key,opt in [('spacefill_percent_vdw','--spacefill'),('wireframe_radius_angstrom','--wireframe')]:
        if rendering.get(key) is not None:flags += [opt,str(rendering[key])]
    args=resolve_render_preset(parser().parse_args(['molecule',*flags]))
    # Preserve the exact already-reviewed ionization state; do not re-protonate.
    newmol,chemistry=prepare_molecule(mol,args,old['molecule_name'])
    sites=carbonyl_measurements(newmol)
    if any(s['height_angstrom']>.03 for s in sites):
        raise ValueError('carbonyl/amide planarity failed: '+str([(s['kind'],s['height_angstrom']) for s in sites if s['height_angstrom']>.03]))
    if Chem.GetFormalCharge(newmol)!=old['formal_charge']:raise ValueError('charge changed')
    staging.mkdir(parents=True,exist_ok=True)
    sdf=staging/old['files']['input_sdf'];writer=Chem.SDWriter(str(sdf));writer.write(newmol);writer.close()
    obj=staging/old['files']['obj'];log=staging/(obj.stem+'_jmol.log')
    command,mesh=render_jmol(args,sdf,obj,log)
    meta={**old,**chemistry,'created_utc':utc_now(),'print_mesh':mesh,'input_sha256':sha256(sdf)}
    for key,value in old.items():
        if key.startswith('protonation') or key=='target_ph':meta[key]=value
    meta['coordinate_method']+='; preserved previously assigned ionization state'
    meta['carbonyl_repair']={'previous_directory':str(archive),'method':'MMFF94s with local carbonyl/amide restraints','before':carbonyl_measurements(mol),'after':sites}
    meta['jmol']={**old.get('jmol',{}),'command':command}
    meta['files']={**old['files'],'jmol_raw_obj':mesh.get('jmol_raw_obj'),'jmol_raw_mtl':mesh.get('jmol_raw_mtl'),'jmol_log':log.name}
    # The archived input and full arguments make this repair reproducible.
    flags[flags.index('--input')+1]=str(archive/old['files']['input_sdf'])
    import shlex
    token=REPRODUCTION_COMMAND.set(shlex.join(['3dprintase','molecule',*flags,'--output',str(root/'out/carbonyl_audit/reproduced')]))
    write_model_metadata(staging,path.name,meta);REPRODUCTION_COMMAND.reset(token)
    archive.parent.mkdir(parents=True,exist_ok=True)
    path.parent.rename(archive);staging.rename(path.parent)
    return {'name':old['molecule_name'],'status':'repaired','metadata':str(path),'maximum_height_angstrom':max(s['height_angstrom'] for s in sites)}

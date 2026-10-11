"""Validate semantic OBJ/MTL element colors against SDF atom identities."""
import re
from pathlib import Path
from .mesh import ELEMENT_MATERIALS,BAMBU_PRESERVE_PALETTE
from .utils import PipelineError


def validate_element_materials(obj_path, mol, *, monocolor_aliphatics=False):
    obj_path=Path(obj_path);text=obj_path.read_text();mtl_name=next(l.split(None,1)[1] for l in text.splitlines() if l.startswith('mtllib '))
    mtl=(obj_path.parent/mtl_name).read_text();palette={};name=None
    for line in mtl.splitlines():
        if line.startswith('newmtl '):name=line.split()[1]
        elif line.startswith('Kd '):palette[name]=tuple(map(float,line.split()[1:4]))
    assigned={};material=None
    for line in text.splitlines():
        if line.startswith('usemtl '):material=line.split()[1]
        elif line.startswith('# component ') and line.endswith('_Atom'):
            match=re.fullmatch(r'# component ([A-Z][a-z]?)(\d+)_Atom',line)
            if not match:raise PipelineError('Unrecognized SDF atom group: '+line)
            element,index=match.groups();index=int(index)-1
            if index in assigned:raise PipelineError('Duplicate exported atom')
            if index>=mol.GetNumAtoms() or mol.GetAtomWithIdx(index).GetSymbol()!=element:raise PipelineError('OBJ atom identity differs from SDF')
            label,rgb=ELEMENT_MATERIALS.get(element,('element',(.7,.7,.7)))
            if monocolor_aliphatics and element == 'H':rgb=ELEMENT_MATERIALS['C'][1]
            if material!=f'atom_{element}_{label}' or palette.get(material)!=rgb:raise PipelineError(f'Incorrect {element} material at atom {index}')
            assigned[index]=element
    expected={a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()>1}
    if not expected<=assigned.keys():raise PipelineError('Missing heavy atom colors in OBJ')
    if BAMBU_PRESERVE_PALETTE not in mtl:raise PipelineError('Missing slicer palette-preservation marker')
    return {'passed':True,'atoms_checked':len(assigned),'oxygen_rgb':palette.get('atom_O_oxygen'),'phosphorus_rgb':palette.get('atom_P_phosphorus'),'slicer_palette_preserved':True}


def validate_exported_atom_coordinates(obj_path, mol):
    """Compare sphere centers in the finished OBJ to SDF coordinates, up to scale/rotation."""
    import numpy as np
    centers={};points=[];current=None
    def finish():
        if current is not None and points:centers[current]=np.mean(points,axis=0)
    for line in Path(obj_path).read_text().splitlines():
        if line.startswith('# component '):
            finish();points=[]
            match=re.fullmatch(r'# component [A-Z][a-z]?(\d+)_Atom',line)
            current=int(match.group(1))-1 if match else None
        elif current is not None and line.startswith('v '):points.append(list(map(float,line.split()[1:4])))
    finish()
    ids=sorted(centers)
    if len(ids)<3:raise PipelineError('Too few atom spheres to validate exported coordinates')
    source=mol.GetConformer().GetPositions()[ids];target=np.array([centers[i] for i in ids])
    source-=source.mean(axis=0);target-=target.mean(axis=0)
    u,s,vt=np.linalg.svd(source.T@target);rotation=u@vt
    scale=float(s.sum()/np.square(source).sum())
    residual=float(np.max(np.linalg.norm(source@rotation-target/scale,axis=1)))
    if residual>.002:raise PipelineError(f'Exported atom centers differ from SDF by {residual:.5f} angstrom')
    return {'passed':True,'atom_centers_checked':len(ids),'maximum_error_angstrom':residual,'export_scale':scale}

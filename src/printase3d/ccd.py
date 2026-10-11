"""RCSB CCD ideal/model conformers with explicit chemical connectivity."""
import re
from pathlib import Path
from .structures import parse_cif, loop_rows
from .utils import PipelineError, request


def read_ccd(text, coordinates='ideal'):
    from rdkit import Chem
    import numpy as np
    singles, loops = parse_cif(text)
    rows = loop_rows(loops, '_chem_comp_atom.')
    if not rows:
        raise PipelineError('CCD entry contains no atom table')
    rw = Chem.RWMol(); positions=[]; names={}
    for row in rows:
        a=lambda key: row['_chem_comp_atom.'+key]
        atom=Chem.Atom(a('type_symbol').capitalize())
        atom.SetFormalCharge(int(a('charge')) if a('charge') not in ('.','?') else 0)
        atom.SetNoImplicit(True)
        name=a('atom_id'); atom.SetProp('ccd_atom_id',name)
        names[name]=rw.AddAtom(atom)
        fields=[('pdbx_model_Cartn_'+v+'_ideal') if coordinates=='ideal' else ('model_Cartn_'+v) for v in 'xyz']
        try: positions.append([float(a(k)) for k in fields])
        except (ValueError,KeyError) as exc: raise PipelineError(f'CCD {coordinates} coordinates unavailable for {name}') from exc
    orders={'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE,'AROM':Chem.BondType.AROMATIC}
    for row in loop_rows(loops,'_chem_comp_bond.'):
        b=lambda key: row['_chem_comp_bond.'+key]
        if b('value_order') not in orders:raise PipelineError('Unsupported CCD bond order '+b('value_order'))
        rw.AddBond(names[b('atom_id_1')],names[b('atom_id_2')],orders[b('value_order')])
    mol=rw.GetMol(); Chem.SanitizeMol(mol)
    if not np.isfinite(positions).all():raise PipelineError('Nonfinite CCD coordinates')
    conf=Chem.Conformer(len(positions));conf.Set3D(True)
    for i,xyz in enumerate(positions):conf.SetAtomPosition(i,xyz)
    mol.AddConformer(conf);Chem.AssignStereochemistryFrom3D(mol)
    # CCD stereo labels are an independent check on coordinate handedness.
    for row in rows:
        code=row.get('_chem_comp_atom.pdbx_stereo_config')
        atom=mol.GetAtomWithIdx(names[row['_chem_comp_atom.atom_id']])
        if code in ('R','S') and atom.GetAtomicNum()==6 and (not atom.HasProp('_CIPCode') or atom.GetProp('_CIPCode')!=code):
            raise PipelineError('CCD coordinates disagree with declared carbon stereochemistry')
    return mol, {'ccd_id':singles.get('_chem_comp.id'), 'ccd_name':singles.get('_chem_comp.name'),
                 'coordinate_set':coordinates,'coordinate_details':singles.get('_chem_comp.pdbx_ideal_coordinates_details' if coordinates=='ideal' else '_chem_comp.pdbx_model_coordinates_details'),
                 'model_pdb_id':singles.get('_chem_comp.pdbx_model_coordinates_db_code')}


def fetch_ccd(code, coordinates='ideal'):
    code=code.upper()
    if not re.fullmatch('[A-Z0-9]{1,5}',code):raise PipelineError('CCD identifier must contain 1–5 letters/digits')
    url=f'https://files.rcsb.org/ligands/download/{code}.cif'
    raw=request(url);mol, provenance=read_ccd(raw.decode(),coordinates)
    return mol,{**provenance,'url':url},raw

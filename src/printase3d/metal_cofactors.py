"""Preserve CCD reference geometries for unsupported metal force fields.

This path performs no conformer search. Explicit CCD coordination bonds are
retained, including source base-off cobalamin conformations. It never labels
CCD coordinates as PubChem or as energy-minimized geometries.
"""
from pathlib import Path
import json
import math
from .structures import parse_cif, loop_rows
from .outputs import model_directory
from .rendering import render_jmol
from .utils import sha256


def reference_model(cif_path, output, label, oxidation, args, *, aqua=False):
    singles, loops = parse_cif(Path(cif_path).read_text())
    raw_atoms = loop_rows(loops, '_chem_comp_atom.')
    raw_bonds = loop_rows(loops, '_chem_comp_bond.')
    atoms = {}
    for row in raw_atoms:
        def a(field): return row['_chem_comp_atom.' + field]
        name=a('atom_id'); element=a('type_symbol').capitalize()
        atoms[name]={'name':name,'element':element,'charge':int(a('charge')) if a('charge') not in ('?','.') else oxidation,
                     'xyz':[float(a('pdbx_model_Cartn_'+v+'_ideal')) for v in ('x','y','z')]}
    bonds=[]
    for row in raw_bonds:
        def b(field): return row['_chem_comp_bond.'+field]
        bonds.append((b('atom_id_1'),b('atom_id_2'),{'SING':1,'DOUB':2,'TRIP':3}[b('value_order')]))
    adjacency={name:[] for name in atoms}
    for a,b,order in bonds:
        adjacency[a].append((b,order)); adjacency[b].append((a,order))
    removed=[]; corrections=[]
    # Only peripheral carboxyl OH groups. Metal donor nitrogens are untouched.
    for name,atom in list(atoms.items()):
        if atom['element']!='O': continue
        hydrogens=[n for n,o in adjacency[name] if atoms[n]['element']=='H']
        carbons=[n for n,o in adjacency[name] if atoms[n]['element']=='C' and o==1]
        if hydrogens and any(any(atoms[n]['element']=='O' and o==2 for n,o in adjacency[c]) for c in carbons):
            removed.extend(hydrogens); atom['charge']=-1; corrections.append('deprotonated carboxyl '+name)
    metal_names=[n for n,a in atoms.items() if a['element'] in ('Fe','Mg','Co')]
    assert len(metal_names)==1
    metal=metal_names[0]; atoms[metal]['charge']=oxidation
    if aqua:
        # pH-7 model uses the aqua member of the aqua/hydroxo equilibrium.
        oxygen=next(n for n,o in adjacency[metal] if atoms[n]['element']=='O')
        h=next(n for n,o in adjacency[oxygen] if atoms[n]['element']=='H')
        op=atoms[oxygen]['xyz']; hp=atoms[h]['xyz']; mp=atoms[metal]['xyz']
        u=[hp[i]-op[i] for i in range(3)]; v=[mp[i]-op[i] for i in range(3)]
        cross=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]
        norm=math.sqrt(sum(x*x for x in cross)); unorm=math.sqrt(sum(x*x for x in u))
        assert norm>1e-6
        direction=[math.cos(math.radians(104.5))*u[i]/unorm + math.sin(math.radians(104.5))*cross[i]/norm for i in range(3)]
        atoms['HAQU']={'name':'HAQU','element':'H','charge':0,'xyz':[op[i]+.96*direction[i] for i in range(3)]}
        atoms[oxygen]['charge']=0; bonds.append((oxygen,'HAQU',1)); corrections.append('axial hydroxide protonated to water; added H is constructed, heavy coordinates retained')
    for n in removed: atoms.pop(n)
    bonds=[(a,b,o) for a,b,o in bonds if a in atoms and b in atoms]
    ordered=list(atoms); index={n:i+1 for i,n in enumerate(ordered)}
    code=singles['_chem_comp.id']
    folder=model_directory(output,label,'ball-stick-ph7-rcsb-reference',unique=True,dry_run=False)
    folder.mkdir(parents=True,exist_ok=True)
    sdf=folder/(label+'.sdf')
    lines=[label,'  3DPrintase          3D','',f'{len(atoms):3d}{len(bonds):3d}  0  0  0  0            999 V2000']
    for a in atoms.values():
        x,y,z=a['xyz']; lines.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {a["element"]:<3} 0  0  0  0  0  0  0  0  0  0  0  0')
    for a,b,o in bonds: lines.append(f'{index[a]:3d}{index[b]:3d}{o:3d}  0  0  0  0')
    charges=[(index[n],a['charge']) for n,a in atoms.items() if a['charge']]
    for start in range(0,len(charges),8):
        chunk=charges[start:start+8]; lines.append(f'M  CHG{len(chunk):3d}'+''.join(f'{i:4d}{c:4d}' for i,c in chunk))
    lines+=['M  END','$$$$']; sdf.write_text('\n'.join(lines)+'\n')
    from rdkit import Chem
    from .carbonyl_geometry import carbonyl_measurements
    reference = Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=False)[0]
    reference.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(reference)
    carbonyl = carbonyl_measurements(reference)
    if any(site['height_angstrom'] > 0.03 for site in carbonyl):
        raise ValueError('Reference carbonyl geometry failed; source heavy coordinates cannot be silently changed')
    obj=folder/(label+'.obj'); log=folder/(label+'_jmol.log')
    command,mesh=render_jmol(args,sdf,obj,log)
    coordination=[{'a':a,'b':b,'distance':math.dist(atoms[a]['xyz'],atoms[b]['xyz'])} for a,b,o in bonds if metal in (a,b)]
    assert all(1.5<b['distance']<2.8 for b in coordination),coordination
    meta={'molecule_name':label,'source_url':f'https://www.rcsb.org/ligand/{code}', 'ccd':code,
          'source_sha256':sha256(cif_path),'geometry_mode':'rcsb-reference','target_ph':7,
          'ph_assignment':'peripheral carboxylates and explicit selected metal/axial state; no protein environment prediction',
          'metal_oxidation_state':oxidation,'formal_charge':sum(a['charge'] for a in atoms.values()),
          'carbonyl_geometry':carbonyl,
          'protonation_corrections':corrections,'coordinate_method':'CCD ideal reference coordinates; no force-field optimization',
          'coordinate_generator':singles.get('_chem_comp.pdbx_ideal_coordinates_details'),
          'atom_order':ordered,'atom_elements':[a['element'] for a in atoms.values()],
          'coordination_bonds':coordination,'bond_count':len(bonds),'heavy_atom_count':sum(a['element']!='H' for a in atoms.values()),
          'rendering':{'preset':'ball-stick'},'print_mesh':mesh,
          'unavailable_variants':{'pubchem':'No PubChem 3D record','lowest-energy':'No validated coordination-specific force-field parameterization','most-flat':'No validated coordination-specific force-field parameterization'},
          'files':{'obj':obj.name,'mtl':obj.with_suffix('.mtl').name,'input_sdf':sdf.name}}
    (folder/(label+'_metadata.json')).write_text(json.dumps(meta,indent=2)+'\n')
    (folder/'MODEL.md').write_text(f'# {label}\n\n[OBJ]({obj.name}) — [source CCD {code}]({meta["source_url"]}).\n\nReference geometry fallback, not PubChem or an energy minimum. Metal oxidation state {oxidation:+d}, net charge {meta["formal_charge"]:+d}. Target pH 7; peripheral carboxylates deprotonated. Protein-specific axial ligands are omitted; source intramolecular coordination is retained.\n\n'+ '\n'.join('- '+c for c in corrections)+'\n')
    return obj

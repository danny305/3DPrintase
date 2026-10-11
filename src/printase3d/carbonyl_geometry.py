"""Local carbonyl/ordinary-amide geometry, including explicit hydrogens.

Only local trigonal centers are planarized, never an entire carbonyl-bearing
molecule. Small-ring and bridgehead amides are excluded from amide restraints
because their real geometry can be pyramidal/twisted.
"""
import numpy as np


def carbonyl_sites(mol):
    from rdkit import Chem
    sites=[]
    for c,o in mol.GetSubstructMatches(Chem.MolFromSmarts('[#6X3]=[#8]')):
        neighbors=[a.GetIdx() for a in mol.GetAtomWithIdx(c).GetNeighbors()]
        if len(neighbors)==3:sites.append({'kind':'carbonyl carbon','center':c,'neighbors':neighbors})
        for n in mol.GetAtomWithIdx(c).GetNeighbors():
            if n.GetAtomicNum()!=7 or n.GetFormalCharge()!=0 or n.GetDegree()!=3:continue
            if n.IsInRingSize(3) or n.IsInRingSize(4) or sum(b.IsInRing() for b in n.GetBonds())==3:continue
            sites.append({'kind':'amide nitrogen','center':n.GetIdx(),'neighbors':[a.GetIdx() for a in n.GetNeighbors()], 'carbonyl':[c,o]})
    return sites


def carbonyl_measurements(mol,conf_id=-1):
    xyz=mol.GetConformer(conf_id).GetPositions();rows=[]
    for site in carbonyl_sites(mol):
        n=site['center']; ns=site['neighbors']; vs=xyz[ns]-xyz[n]
        normal=np.cross(xyz[ns[1]]-xyz[ns[0]],xyz[ns[2]]-xyz[ns[0]])
        h=abs(float(np.dot(xyz[n]-xyz[ns[0]],normal)/np.linalg.norm(normal)))
        angles=sum(float(np.degrees(np.arccos(np.clip(np.dot(vs[i],vs[j])/np.linalg.norm(vs[i])/np.linalg.norm(vs[j]),-1,1)))) for i,j in [(0,1),(1,2),(0,2)])
        rows.append({**site,'height_angstrom':h,'bond_angle_sum_degrees':angles})
    return rows


def add_carbonyl_restraints(ff,mol,conf_id,forcefield):
    from rdkit.Chem import rdMolTransforms
    add=ff.MMFFAddTorsionConstraint if forcefield.startswith('MMFF') else ff.UFFAddTorsionConstraint
    conf=mol.GetConformer(conf_id)
    torsions=set()
    for site in carbonyl_sites(mol):
        a,b,c=site['neighbors'];center=site['center']
        torsions.add((a,center,b,c))
        if site['kind']=='amide nitrogen':
            carbon,oxygen=site['carbonyl']
            for neighbor in site['neighbors']:
                if neighbor!=carbon:torsions.add((oxygen,carbon,center,neighbor))
    for atoms in torsions:
        angle=rdMolTransforms.GetDihedralDeg(conf,*atoms)
        target=0.0 if abs(angle)<90 else (180.0 if angle>=0 else -180.0)
        add(*atoms,False,target-.25,target+.25,100.0)
    return len(torsions)

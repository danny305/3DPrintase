"""Planarity restraints for individual and fused fully aromatic ring components.

Restrict automatic plane restraints to components of 5–14 atoms.
Larger fused systems can be intrinsically curved and are not flattened.

Perceive RDKit aromaticity on an isolated copy: MMFF atom typing mutates its
own aromaticity flags. Reduced flavin bridge nitrogens are not aromatic and
must not join its two aromatic subrings into one restrained plane.
"""
import math


def aromatic_components(mol):
    from rdkit import Chem
    copy = Chem.Mol(mol)
    Chem.SanitizeMol(copy)
    remaining = {a.GetIdx() for a in copy.GetAtoms() if a.GetIsAromatic()}
    groups = []
    while remaining:
        group = {remaining.pop()}
        pending = list(group)
        while pending:
            for atom in copy.GetAtomWithIdx(pending.pop()).GetNeighbors():
                i = atom.GetIdx()
                if i in remaining:
                    # Only ring bonds: biaryls must remain free to twist.
                    if not any(copy.GetBondBetweenAtoms(i,j).IsInRing() for j in group if copy.GetBondBetweenAtoms(i,j)):
                        continue
                    remaining.remove(i); group.add(i); pending.append(i)
        if 5 <= len(group) <= 14:
            groups.append(sorted(group))
    return groups


def add_planarity_restraints(ff, mol, conf_id, groups, forcefield):
    from rdkit.Chem import rdMolTransforms
    conf = mol.GetConformer(conf_id)
    count = 0
    for group in groups:
        allowed = set(group)
        for bond in mol.GetBonds():
            b,c = bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
            if b not in allowed or c not in allowed or not bond.IsInRing(): continue
            for aa in mol.GetAtomWithIdx(b).GetNeighbors():
                a=aa.GetIdx()
                if a==c or a not in allowed: continue
                for dd in mol.GetAtomWithIdx(c).GetNeighbors():
                    d=dd.GetIdx()
                    if d in (a,b) or d not in allowed: continue
                    angle=rdMolTransforms.GetDihedralDeg(conf,a,b,c,d)
                    target=0.0 if abs(angle)<90 else (180.0 if angle>=0 else -180.0)
                    add=ff.MMFFAddTorsionConstraint if forcefield.startswith('MMFF') else ff.UFFAddTorsionConstraint
                    add(a,b,c,d,False,target-0.25,target+0.25,100.0)
                    count+=1
    return count


def ring_planarity(mol, groups, conf_id=-1):
    import numpy as np
    xyz=np.asarray(mol.GetConformer(conf_id).GetPositions())
    return [float(np.linalg.svd(xyz[g]-xyz[g].mean(axis=0),compute_uv=False)[-1]/math.sqrt(len(g))) for g in groups]


def flavin_geometry(mol):
    """Report oxidized vs reduced isoalloxazine-like fused 6/6/6 systems.

    This is structural annotation, not a quantum aromaticity calculation.
    Reduced flavins retain aromatic subrings but their complete tricycle is
    not constrained to a plane; no unsupported universal bend angle is imposed.
    """
    from rdkit import Chem
    copy=Chem.Mol(mol); Chem.SanitizeMol(copy)
    rings=[set(r) for r in copy.GetRingInfo().AtomRings() if len(r)==6]
    reports=[]
    import itertools
    for triple in itertools.combinations(rings,3):
        core=set.union(*triple)
        if len(core)!=14 or sum(copy.GetAtomWithIdx(i).GetAtomicNum()==7 for i in core)!=4: continue
        if sum(len(a&b)==2 for a,b in itertools.combinations(triple,2))!=2: continue
        carbonyls=sum(1 for i in core for b in copy.GetAtomWithIdx(i).GetBonds()
                      if b.GetBondType()==Chem.BondType.DOUBLE and b.GetOtherAtomIdx(i) not in core and b.GetOtherAtom(copy.GetAtomWithIdx(i)).GetAtomicNum()==8)
        if carbonyls!=2: continue
        aromatic=all(copy.GetAtomWithIdx(i).GetIsAromatic() for i in core)
        reports.append({'atom_indices':sorted(core),'state':'oxidized' if aromatic else 'reduced/non-fully-aromatic',
                        'ring_rms_angstrom':ring_planarity(copy,[sorted(core)])[0],
                        'planarity_restrained':aromatic,
                        'interpretation':'oxidized conjugated core' if aromatic else 'butterfly bending permitted; aromatic subrings do not make the full tricycle aromatic'})
    return reports

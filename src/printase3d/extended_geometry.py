"""Dinucleotide extension using a terminal-ring separation restraint."""
from .utils import PipelineError


def extension_definition(mol):
    from rdkit import Chem
    import numpy as np
    bridges=[a for a in mol.GetAtoms() if a.GetAtomicNum()==8 and sum(n.GetAtomicNum()==15 for n in a.GetNeighbors())==2]
    distances=Chem.GetDistanceMatrix(mol)
    for bridge in bridges:
        ps=[n.GetIdx() for n in bridge.GetNeighbors() if n.GetAtomicNum()==15]
        groups=[]; anchors=[]
        for p in ps:
            side={p}; pending=[p]
            while pending:
                for a in mol.GetAtomWithIdx(pending.pop()).GetNeighbors():
                    i=a.GetIdx()
                    if i!=bridge.GetIdx() and i not in side:
                        side.add(i); pending.append(i)
            rings=[set(r) for r in mol.GetRingInfo().AtomRings() if set(r)<=side]
            fused=[]
            while rings:
                group=rings.pop(); changed=True
                while changed:
                    changed=False
                    for ring in rings[:]:
                        if len(ring&group)>=2:
                            group|=ring; rings.remove(ring); changed=True
                fused.append(group)
            # Require a nitrogen-containing base/chromophore on both arms.
            fused=[g for g in fused if any(mol.GetAtomWithIdx(i).GetAtomicNum()==7 for i in g)]
            if not fused:break
            group=max(fused,key=lambda g:float(np.mean([distances[p,i] for i in g])))
            groups.append(sorted(group)); anchors.append(max(group,key=lambda i:distances[p,i]))
        if len(groups)==2:
            path=Chem.GetShortestPath(mol,*anchors)
            return {'ring_groups':groups,'anchor_atoms':anchors,'minimum_anchor_distance_angstrom':(len(path)-1)*1.15,
                    'force_constant':10.0,'minimum_ring_centroid_distance_angstrom':8.0}
    # Single nucleosides/nucleotides: extend a base-to-tail heavy-atom path.
    # Rigid bases have no meaningful extendable tail; rank their span without
    # pulling a ring apart.
    heavy=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()>1]
    if len(heavy)<2:raise PipelineError('--extended requires at least two heavy atoms')
    from rdkit.Chem import Lipinski
    flexible=bool(Lipinski.NumRotatableBonds(Chem.RemoveHs(mol)))
    anchors=max(((a,b) for i,a in enumerate(heavy) for b in heavy[i+1:]),key=lambda pair:distances[pair[0],pair[1]])
    if not np.isfinite(distances[anchors[0],anchors[1]]) or distances[anchors[0],anchors[1]]>100:
        raise PipelineError('--extended requires one connected molecule')
    path=Chem.GetShortestPath(mol,*anchors)
    # Ring chords shorten the achievable length versus an all-trans chain.
    ring_steps=sum(mol.GetBondBetweenAtoms(a,b).IsInRing() for a,b in zip(path,path[1:]))
    target=(len(path)-1-ring_steps)*1.2+ring_steps*.8
    return {'ring_groups':[[anchors[0]],[anchors[1]]], 'anchor_atoms':list(anchors),
            'minimum_anchor_distance_angstrom':target, 'force_constant':10.0 if flexible else 0.0,
            'minimum_ring_centroid_distance_angstrom':target*.85 if flexible else 0.0,
            'metric':'terminal heavy-atom distance', 'rigid_no_extension':not flexible}


def add_extension_restraint(ff,definition,forcefield):
    if not definition['force_constant']:return
    add=ff.MMFFAddDistanceConstraint if forcefield.startswith('MMFF') else ff.UFFAddDistanceConstraint
    add(*definition['anchor_atoms'],False,definition['minimum_anchor_distance_angstrom'],100.0,definition['force_constant'])


def extension_measurement(mol,definition,conf_id=-1):
    import numpy as np
    xyz=np.asarray(mol.GetConformer(conf_id).GetPositions())
    a,b=definition['ring_groups']
    return float(np.linalg.norm(xyz[a].mean(axis=0)-xyz[b].mean(axis=0)))

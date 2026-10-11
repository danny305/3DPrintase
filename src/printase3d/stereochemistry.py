"""Reject minimized coordinates that invert persistent input stereocenters."""

def preserves_stereochemistry(template,candidate,conf_id=-1):
    from rdkit import Chem
    expected=Chem.Mol(template)
    Chem.AssignStereochemistry(expected,cleanIt=True,force=True)
    observed=Chem.Mol(candidate)
    Chem.RemoveStereochemistry(observed)
    Chem.AssignStereochemistryFrom3D(observed,confId=conf_id)
    Chem.AssignStereochemistry(observed,cleanIt=True,force=True)
    for atom in expected.GetAtoms():
        # Ordinary N inversion and phosphate resonance labels are not stable
        # stereocenters. Preserve carbon, sulfonium, and other assigned centers.
        if atom.GetAtomicNum() in (7,15) or not atom.HasProp('_CIPCode'):continue
        actual=observed.GetAtomWithIdx(atom.GetIdx())
        if not actual.HasProp('_CIPCode') or actual.GetProp('_CIPCode')!=atom.GetProp('_CIPCode'):
            return False
    return True

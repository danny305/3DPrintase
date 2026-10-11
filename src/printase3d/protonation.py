"""Narrow, source-backed corrections for unsupported ionization sites."""


def correct_free_selenocysteine(mol, ph):
    """Correct only the selenol of free selenocysteine (either enantiomer).

    Dimorphite-DL 2.0 lacks selenium rules. The measured free-amino-acid
    selenohydryl pKa is 5.24: https://pubmed.ncbi.nlm.nih.gov/35659493/.
    This rule deliberately does not apply to peptides or other selenols.
    """
    from rdkit import Chem
    if mol.GetNumHeavyAtoms() != 7:
        return []
    query = Chem.MolFromSmarts('[#7;D1]-[C;X4](-[CH2]-[Se;D1])-[C](=[O])-[O;D1]')
    match = mol.GetSubstructMatch(query)
    if not match:
        return []
    selenium = mol.GetAtomWithIdx(match[3])
    charge = -1 if ph >= 5.24 else 0
    if selenium.GetFormalCharge() == charge:
        return []
    selenium.SetFormalCharge(charge)
    selenium.SetNumExplicitHs(0 if charge == -1 else 1)
    selenium.SetNoImplicit(True)
    Chem.SanitizeMol(mol)
    return [{'site': 'free selenocysteine selenol', 'pka': 5.24,
             'formal_charge': charge, 'source': 'https://pubmed.ncbi.nlm.nih.gov/35659493/'}]

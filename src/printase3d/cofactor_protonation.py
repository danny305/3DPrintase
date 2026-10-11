"""Exact-structure pH-7 corrections for known cofactor ionization failures.

The registry matches the complete heavy graph, bond orders, and charges.
Stereochemistry is preserved, and no rule applies at other pH values.
"""
KNOWN_RING_CATIONS = {'NC(=O)C1=C[NH+](C2OC(COP(=O)([O-])OP(=O)([O-])OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1': {'name': 'NADH', 'source': 'https://www.ebi.ac.uk/chebi/CHEBI:57945'}, 'C[NH+]1c2c(nc(N)[nH]c2=O)NCC1CNc1ccc(C(=O)NC(CCC(=O)[O-])C(=O)[O-])cc1': {'name': '5-methyl-THF', 'source': 'https://www.ebi.ac.uk/chebi/CHEBI:189667'}, 'Nc1nc2c(c(=O)[nH]1)[NH2+]C(CNc1ccc(C(=O)NC(CCC(=O)[O-])C(=O)[O-])cc1)CN2': {'name': 'THF', 'source': 'https://www.ebi.ac.uk/chebi/CHEBI:57453'}, 'Nc1nc2c(c(=O)[nH]1)[NH+]1CN(c3ccc(C(=O)NC(CCC(=O)[O-])C(=O)[O-])cc3)CC1CN2': {'name': '5-10-methylene-THF', 'source': 'https://www.ebi.ac.uk/chebi/CHEBI:15636'}, 'NC(=O)C1=C[NH+](C2OC(COP(=O)([O-])OP(=O)([O-])OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)([O-])[O-])C3O)C(O)C2O)C=CC1': {'name': 'NADPH', 'source': 'https://www.ebi.ac.uk/chebi/CHEBI:57783'}, 'CC(O)C(O)C1CNc2nc(N)[nH]c(=O)c2[NH2+]1': {'name': 'BH4', 'source': 'https://www.ebi.ac.uk/chebi/CHEBI:15372'}}


def correct_cofactor_ring_cations(mol, ph):
    from rdkit import Chem
    if ph != 7:
        return []
    entry = KNOWN_RING_CATIONS.get(Chem.MolToSmiles(mol, isomericSmiles=False))
    if entry is None:
        return []
    changed = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1 and atom.GetTotalNumHs() and not atom.GetIsAromatic():
            atom.SetNumExplicitHs(atom.GetTotalNumHs() - 1)
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(True)
            changed.append(atom.GetIdx())
    Chem.SanitizeMol(mol)
    return [{'site': 'cofactor reduced ring nitrogen', 'compound': entry['name'],
             'atom_indices': changed, 'source': entry['source'], 'scope': 'pH 7 only'}] if changed else []

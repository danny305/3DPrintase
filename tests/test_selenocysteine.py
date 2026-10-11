import unittest
from rdkit import Chem
from printase3d.mesh import semantic_atom_material
from printase3d.protonation import correct_free_selenocysteine


class SelenocysteineTests(unittest.TestCase):
    def test_free_amino_acid_selenolate_at_ph7(self):
        mol = Chem.MolFromSmiles('[NH3+][C@@H](C[SeH])C(=O)[O-]')
        corrections = correct_free_selenocysteine(mol, 7)
        self.assertEqual(Chem.GetFormalCharge(mol), -1)
        self.assertEqual(corrections[0]['pka'], 5.24)
        self.assertEqual(correct_free_selenocysteine(mol, 7), [])
        correct_free_selenocysteine(mol, 4)
        self.assertEqual(Chem.GetFormalCharge(mol), 0)

    def test_other_selenols_and_peptides_unaffected(self):
        for smiles in ('C[SeH]', 'N[C@@H](CC[SeH])C(=O)O', 'CN[C@@H](C[SeH])C(=O)O'):
            mol = Chem.MolFromSmiles(smiles)
            self.assertEqual(correct_free_selenocysteine(mol, 7), [])

    def test_selenium_has_its_own_atom_material(self):
        for name in ('Se1_Atom', 'SE1_Atom'):
            self.assertEqual(semantic_atom_material(name)[0], 'atom_Se_selenium')
            self.assertNotEqual(semantic_atom_material(name), semantic_atom_material('S1_Atom'))


if __name__ == '__main__':
    unittest.main()

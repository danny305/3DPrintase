import unittest
from rdkit import Chem
from printase3d.cofactor_protonation import KNOWN_RING_CATIONS, correct_cofactor_ring_cations
from printase3d.mesh import semantic_atom_material

class CofactorTests(unittest.TestCase):
    def test_reduced_cofactor_charges(self):
        expected={'NADH':-2,'NADPH':-4,'THF':-2,'5-methyl-THF':-2,'5-10-methylene-THF':-2,'BH4':0}
        for smiles, entry in KNOWN_RING_CATIONS.items():
            m=Chem.MolFromSmiles(smiles)
            before=Chem.MolToSmiles(m)
            self.assertEqual(correct_cofactor_ring_cations(m,4),[])
            self.assertEqual(Chem.MolToSmiles(m),before)
            self.assertTrue(correct_cofactor_ring_cations(m,7))
            self.assertEqual(Chem.GetFormalCharge(m),expected[entry['name']])
            self.assertEqual(correct_cofactor_ring_cations(m,7),[])
    def test_unrelated_amines_remain_protonated(self):
        for s in ['C[NH2+]C','[NH3+]CC(=O)[O-]','C[N+]1=CC=CC=C1']:
            m=Chem.MolFromSmiles(s)
            self.assertEqual(correct_cofactor_ring_cations(m,7),[])
            self.assertEqual(Chem.GetFormalCharge(m),Chem.GetFormalCharge(Chem.MolFromSmiles(s)))
    def test_metal_colors(self):
        for symbol,name in [('Fe','iron'),('Mg','magnesium'),('Co','cobalt')]:
            for atom in [symbol+'1_Atom',symbol.upper()+'1_Atom']:
                self.assertEqual(semantic_atom_material(atom)[0],f'atom_{symbol}_{name}')

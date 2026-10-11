import unittest
from rdkit import Chem
from rdkit.Chem import AllChem
from printase3d.stereochemistry import preserves_stereochemistry
from printase3d.conformers import sample_lowest_energy_conformer

class StereoPreservationTests(unittest.TestCase):
    def test_inverted_coordinates_rejected(self):
        mol=Chem.AddHs(Chem.MolFromSmiles('N[C@@H](C)C(=O)O'));AllChem.EmbedMolecule(mol,randomSeed=42)
        mirrored=Chem.Mol(mol)
        for i,p in enumerate(mol.GetConformer().GetPositions()):
            mirrored.GetConformer().SetAtomPosition(i,(-p[0],p[1],p[2]))
        self.assertTrue(preserves_stereochemistry(mol,mol))
        self.assertFalse(preserves_stereochemistry(mol,mirrored))
    def test_biotin_most_flat_keeps_input_stereochemistry(self):
        mol=Chem.AddHs(Chem.MolFromSmiles('O=C([O-])CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@@H]21'))
        AllChem.EmbedMolecule(mol,randomSeed=61453)
        out,meta=sample_lowest_energy_conformer(mol,seed=61453,conformer_count=12,prune_rms=.15,max_iters=10000,selection='most-flat')
        self.assertTrue(meta['selected_converged'])
        self.assertTrue(preserves_stereochemistry(mol,out))

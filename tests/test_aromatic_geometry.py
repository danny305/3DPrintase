import unittest
from rdkit import Chem
from rdkit.Chem import AllChem
from printase3d.aromatic_geometry import aromatic_components,ring_planarity,flavin_geometry
from printase3d.conformers import sample_lowest_energy_conformer

OX='Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C)c2cc1C'
RED='Cc1cc2c(cc1C)N(C)c1[nH]c(=O)[nH]c(=O)c1N2'
class AromaticGeometryTests(unittest.TestCase):
    def test_oxidized_flavin_is_planar_after_every_selection_mode(self):
        for mode in ('energy','most-flat','least-flat'):
            mol=Chem.AddHs(Chem.MolFromSmiles(OX)); before=Chem.MolToSmiles(mol)
            AllChem.EmbedMolecule(mol,randomSeed=42)
            # Reproduce a deliberately buckled source ring.
            conf=mol.GetConformer(); p=conf.GetAtomPosition(4); p.z+=.4; conf.SetAtomPosition(4,p)
            out,meta=sample_lowest_energy_conformer(mol,seed=42,conformer_count=4,prune_rms=.1,max_iters=5000,selection=mode)
            self.assertEqual(Chem.MolToSmiles(mol),before)
            self.assertTrue(meta['selected_converged'])
            self.assertLess(max(ring_planarity(out,aromatic_components(out))),.03)
            self.assertEqual(flavin_geometry(out)[0]['state'],'oxidized')
    def test_reduced_flavin_is_not_forced_into_oxidized_plane(self):
        mol=Chem.AddHs(Chem.MolFromSmiles(RED)); AllChem.EmbedMolecule(mol,randomSeed=42)
        self.assertTrue(all(len(g)==6 for g in aromatic_components(mol)))
        self.assertFalse(flavin_geometry(mol)[0]['planarity_restrained'])
    def test_other_fused_aromatics_and_free_biaryl(self):
        for s in ['c1ccc2ccccc2c1','c1ccc2[nH]ccc2c1','Nc1ncnc2[nH]cnc12']:
            self.assertEqual(len(aromatic_components(Chem.MolFromSmiles(s))),1)
        self.assertEqual(sorted(map(len,aromatic_components(Chem.MolFromSmiles('c1ccccc1-c1ccccc1')))),[6,6])
        self.assertEqual(aromatic_components(Chem.MolFromSmiles('C1=CC=CC=CC=C1')),[])

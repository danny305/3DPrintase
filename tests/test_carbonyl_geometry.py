import unittest
from rdkit import Chem
from rdkit.Chem import AllChem
from printase3d.carbonyl_geometry import carbonyl_sites,carbonyl_measurements
from printase3d.conformers import sample_lowest_energy_conformer
from printase3d.arguments import parser
from printase3d.chemistry import prepare_molecule

class CarbonylGeometryTests(unittest.TestCase):
    def test_primary_secondary_tertiary_amides_and_other_carbonyls(self):
        for smiles in ['NC(=O)c1cccnc1','CC(=O)NC','CC(=O)N(C)C','NC(=O)N','CC(=O)C','CC(=O)OC','CC(=O)[O-]']:
            mol=Chem.AddHs(Chem.MolFromSmiles(smiles));AllChem.EmbedMolecule(mol,randomSeed=12)
            result,meta=sample_lowest_energy_conformer(mol,seed=12,conformer_count=3,prune_rms=.15,max_iters=5000)
            self.assertEqual(meta['forcefield'],'MMFF94s')
            for site in carbonyl_measurements(result):
                self.assertLess(site['height_angstrom'],.03)
                self.assertGreater(site['bond_angle_sum_degrees'],359)
    def test_preserved_heavy_atoms_with_repaired_amide_hydrogens(self):
        import numpy as np
        mol=Chem.AddHs(Chem.MolFromSmiles('CC(=O)N'));AllChem.EmbedMolecule(mol,randomSeed=12)
        AllChem.MMFFOptimizeMolecule(mol)
        before=mol.GetConformer().GetPositions()[:4].copy()
        args=parser().parse_args(['molecule','--smiles','CC(=O)N','--pubchem-geometry'])
        out,_=prepare_molecule(mol,args,'acetamide')
        self.assertTrue(np.allclose(before,out.GetConformer().GetPositions()[:4]))
        self.assertTrue(all(s['height_angstrom']<.03 for s in carbonyl_measurements(out)))
    def test_strained_amides_and_ordinary_amines_are_not_flattened(self):
        for smiles in ['CCN','O=C1CCN1','CC(=O)[NH2+]C']:
            mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
            self.assertFalse(any(s['kind']=='amide nitrogen' for s in carbonyl_sites(mol)))

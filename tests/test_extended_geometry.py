import unittest
from rdkit import Chem
from printase3d.arguments import parser
from printase3d.extended_geometry import extension_definition,extension_measurement
from printase3d.conformers import sample_lowest_energy_conformer
from printase3d.utils import PipelineError

NAD='NC(=O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)c1'
class ExtendedGeometryTests(unittest.TestCase):
    def test_dinucleotide_extension_and_selection(self):
        mol=Chem.AddHs(Chem.MolFromSmiles(NAD))
        definition=extension_definition(mol)
        self.assertEqual(sorted(map(len,definition['ring_groups'])),[6,9])
        out,meta=sample_lowest_energy_conformer(mol,seed=42,conformer_count=3,prune_rms=.15,max_iters=5000,selection='extended')
        self.assertTrue(meta['selected_converged'])
        self.assertGreater(extension_measurement(out,definition),8)
        pool=[c for c in meta['candidates'] if c['converged']]
        self.assertAlmostEqual(meta['extension']['selected_terminal_separation_angstrom'],max(c['terminal_separation_angstrom'] for c in pool),places=4)
    def test_generic_extension_and_rigid_bases(self):
        for smiles in ['CC','Nc1ncnc2[nH]cnc12','OP(=O)(O)OP(=O)(O)O']:
            definition=extension_definition(Chem.MolFromSmiles(smiles))
            self.assertEqual(len(definition['anchor_atoms']),2)
        self.assertTrue(extension_definition(Chem.MolFromSmiles('Nc1ncnc2[nH]cnc12'))['rigid_no_extension'])
        with self.assertRaises(PipelineError):extension_definition(Chem.MolFromSmiles('[Na+]'))
    def test_parser_and_api_option(self):
        from printase3d.api import molecule_tasks,Molecule
        from printase3d.options import variant_name
        args=parser().parse_args(['molecule','--smiles',NAD,'--extended','--ball-stick','--ph','7'])
        self.assertIn('extended',variant_name(args))
        self.assertIn('--extended',molecule_tasks([Molecule(smiles=NAD)],extended=True)[0].argv)

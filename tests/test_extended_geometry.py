import unittest
from rdkit import Chem
from printase3d.arguments import parser
from printase3d.extended_geometry import extension_definition,extension_measurement,peptide_backbone_atoms
from printase3d.conformers import sample_lowest_energy_conformer
from printase3d.utils import PipelineError

NAD='NC(=O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)c1'
class ExtendedGeometryTests(unittest.TestCase):
    def test_peptide_extension_uses_backbone_not_side_chains(self):
        from dimorphite_dl import protonate_smiles
        for sequence in ('RGW', 'KDE', 'NQP', 'PG', 'CC'):
            original = Chem.MolFromSequence(sequence)
            states = [original, Chem.MolFromSmiles(protonate_smiles(
                Chem.MolToSmiles(original), ph_min=7, ph_max=7, precision=0)[0])]
            for mol in states:
                with self.subTest(sequence=sequence, charged=Chem.GetFormalCharge(mol)):
                    definition = extension_definition(mol)
                    backbone = definition['backbone_atoms']
                    self.assertEqual(len(backbone), len(sequence) * 3)
                    self.assertEqual([mol.GetAtomWithIdx(i).GetSymbol() for i in backbone],
                                     ['N', 'C', 'C'] * len(sequence))
                    self.assertEqual(definition['anchor_atoms'], [backbone[0], backbone[-1]])
                    self.assertTrue(all(mol.GetBondBetweenAtoms(a, b) is not None
                                        for a, b in zip(backbone, backbone[1:])))
        self.assertIsNone(peptide_backbone_atoms(Chem.MolFromSmiles('CC(=O)NCC(=O)O')))

    def test_peptide_backbone_extension_minimizes_and_selects(self):
        mol = Chem.AddHs(Chem.MolFromSequence('RGW'))
        out, meta = sample_lowest_energy_conformer(
            mol, seed=61453, conformer_count=3, prune_rms=.15,
            max_iters=5000, selection='extended')
        self.assertTrue(meta['selected_converged'])
        extension = meta['extension']
        self.assertGreater(extension_measurement(out, extension), 9.0)
        self.assertAlmostEqual(extension_measurement(out, extension),
                               max(c['terminal_separation_angstrom'] for c in meta['candidates'] if c['converged']), places=4)

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

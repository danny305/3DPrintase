"""Peptide connectivity, stereochemistry, CLI validation, and output regression."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from printase3d.arguments import parser
from printase3d.cli import main
from printase3d.peptides import AMINO_ACIDS, build_peptide, peptide_sequence
from printase3d.stereochemistry import preserves_stereochemistry


class PeptideTests(unittest.TestCase):
    def test_rgw_has_two_peptide_bonds_and_free_termini(self):
        mol, provenance = build_peptide('rgw')
        self.assertEqual(Chem.MolToSequence(mol), 'RGW')
        self.assertEqual(rdMolDescriptors.CalcMolFormula(mol), 'C19H27N7O4')
        self.assertEqual(len(Chem.GetMolFrags(mol)), 1)
        self.assertEqual(len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CX3](=O)[NX3]'))), 2)
        self.assertEqual(len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CX3](=O)[OX2H1]'))), 1)
        self.assertEqual(len(mol.GetSubstructMatches(Chem.MolFromSmarts('[NX3H2][C@@H]C(=O)'))), 1)
        self.assertEqual([cip for _, cip in Chem.FindMolChiralCenters(mol)], ['S', 'S'])
        self.assertEqual(provenance['residue_count'], 3)
        self.assertEqual(provenance['sequence_direction'], 'N-to-C')

    def test_standard_residues_and_l_cysteine(self):
        for residue in sorted(AMINO_ACIDS):
            with self.subTest(residue=residue):
                mol, _ = build_peptide(residue)
                self.assertEqual(Chem.MolToSequence(mol), residue)
                centers = Chem.FindMolChiralCenters(mol)
                if residue == 'G':
                    self.assertFalse(centers)
                elif residue == 'C':
                    self.assertEqual(centers[0][1], 'R')
                else:
                    self.assertEqual(centers[0][1], 'S')

    def test_invalid_sequences_fail_without_output(self):
        for value in ('', 'B', 'RGX', 'RGU', 'R-G-W', 'R G W', '>RGW', 'RGW1'):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    peptide_sequence(value)
                with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as exc:
                    parser().parse_args(['peptide', value])
                self.assertEqual(exc.exception.code, 2)

    def test_peptide_cli_defaults_and_dry_run(self):
        args = parser().parse_args(['peptide', 'rgw'])
        self.assertEqual(args.sequence, 'RGW')
        self.assertEqual(args.render_preset, 'ball-stick')
        self.assertEqual(args.output, 'out/misc')
        for flag in ('--pubchem-geometry', '--rcsb-geometry'):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                parser().parse_args(['peptide', 'RGW', flag])
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            root = Path(tmp) / 'new'
            self.assertEqual(main(['peptide', 'RGW', '--dry-run', '--output', str(root)]), 0)
            self.assertFalse(root.exists())

    def test_rgw_minimized_output_and_ph_processing(self):
        template, _ = build_peptide('RGW')
        for flags, charge in (([], 0), (['--ph', '7', '--no-hydrogens', '--monocolor-aliphatics'], 1)):
            with self.subTest(flags=flags), tempfile.TemporaryDirectory() as tmp:
                argv = ['peptide', 'RGW', '--label', 'Arg Gly Trp', '--no-render',
                        '--conformer-count', '3', '--output', tmp, *flags]
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(main(argv), 0)
                sdf = next(Path(tmp).rglob('*.sdf'))
                mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
                self.assertEqual(mol.GetNumHeavyAtoms(), template.GetNumHeavyAtoms())
                self.assertGreater(mol.GetNumAtoms(), mol.GetNumHeavyAtoms())
                self.assertEqual(Chem.GetFormalCharge(mol), charge)
                self.assertTrue(mol.GetConformer().Is3D())
                if not flags:
                    self.assertTrue(preserves_stereochemistry(template, mol))
                meta = json.loads(next(Path(tmp).rglob('*_metadata.json')).read_text())
                self.assertEqual(meta['peptide']['sequence'], 'RGW')
                self.assertEqual(meta['molecule_name'], 'Arg Gly Trp')
                self.assertIn('3dprintase peptide RGW', meta['reproduce'])
                self.assertGreater(meta['conformer_search']['converged_candidates'], 0)
                self.assertTrue((sdf.parent / 'MODEL.md').is_file())

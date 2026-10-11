"""Regression checks against pre-refactor sugar chemistry and process isolation."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from printase3d import Molecule, generate_batch, generate_molecule
from printase3d import chemistry
from printase3d.batch import cpu_count
from printase3d.cli import parser, resolve_render_preset
from printase3d.rendering import jmol_command

FIXTURES = json.loads((Path(__file__).parent / 'data/sugars.json').read_text())


class LibraryTests(unittest.TestCase):
    def test_sugar_chemistry_matches_pre_refactor_baselines(self):
        from rdkit import Chem
        for sugar in FIXTURES:
            with self.subTest(sugar=sugar['label']):
                args = resolve_render_preset(parser().parse_args([
                    'molecule', '--smiles', sugar['smiles'], '--ph', '7',
                    '--ensemble', '--conformer-count', '8', '--coordinates', 'regenerate']))
                mol, meta = chemistry.prepare_molecule(Chem.MolFromSmiles(sugar['smiles']), args, sugar['label'])
                self.assertEqual(meta['molecular_formula'], sugar['formula'])
                self.assertEqual(meta['formal_charge'], sugar['charge'])
                self.assertEqual(meta['isomeric_smiles'], sugar['isomeric_smiles'])
                self.assertAlmostEqual(meta['conformer_search']['selected_energy_kcal_mol'], sugar['energy'], places=4)

    def test_process_batch_isolated_failures_and_reproducible_geometry(self):
        entries = [Molecule(smiles=x['smiles'], label=x['label']) for x in FIXTURES[:4]]
        entries.append(Molecule(smiles='invalid-smiles', label='bad-input'))
        with tempfile.TemporaryDirectory() as folder:
            result = generate_batch(entries, jobs=2, output=folder, ph=7, ensemble=True,
                                    conformer_count=8, coordinates='regenerate', no_render=True,
                                    ball_stick=True)
            self.assertEqual(len(result.outputs), 4)
            self.assertEqual(len(result.failures), 1)
            self.assertEqual(result.failures[0].label, 'bad-input')
            self.assertEqual(len({r.worker_pid for r in result.results}), 2)
            self.assertEqual([r.label for r in result.results[:4]], [x['label'] for x in FIXTURES[:4]])
            for output, expected in zip(result.outputs, FIXTURES):
                meta = json.loads(next(output.glob('*_metadata.json')).read_text())
                self.assertEqual(meta['isomeric_smiles'], expected['isomeric_smiles'])
                self.assertAlmostEqual(meta['conformer_search']['selected_energy_kcal_mol'], expected['energy'], places=4)
                self.assertIn('--ensemble', meta['reproduce'])
                self.assertIn(expected['label'], meta['reproduce'])
            summary = json.loads((Path(folder) / 'run_summary.json').read_text())
            self.assertEqual(len(summary['outputs']), 4)
            self.assertEqual(len(summary['failures']), 1)

    def test_duplicate_normalized_labels_rejected_before_writing(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(ValueError, 'Duplicate'):
                generate_batch([Molecule(smiles='CO', label='same name'),
                                Molecule(smiles='CCO', label='same-name')], jobs=2, output=folder)
            self.assertEqual(list(Path(folder).iterdir()), [])

    def test_invalid_worker_counts_and_oversubscription(self):
        entries = [Molecule(smiles='CO', label='one'), Molecule(smiles='CCO', label='two')]
        for jobs in (0, -1, 1.5, True):
            with self.assertRaises(ValueError):
                generate_batch(entries, jobs=jobs)
        with self.assertRaisesRegex(ValueError, 'threads 1'):
            generate_batch(entries, jobs=2, threads=2)

    def test_sources_and_single_molecule_api(self):
        with self.assertRaises(ValueError):
            Molecule()
        with self.assertRaises(ValueError):
            Molecule(cid=1, smiles='CO')
        with tempfile.TemporaryDirectory() as folder:
            output = generate_molecule(Molecule(smiles='CO', label='methanol'), output=folder,
                                       no_render=True, ensemble=True, conformer_count=2)
            self.assertTrue(next(output.glob('*.sdf')).is_file())

    def test_background_command_has_no_display_and_no_activation(self):
        command = jmol_command(Path('Jmol.jar'), Path('model.jmol.spt'), threads=1)
        self.assertIn('-Dapple.awt.UIElement=true', command)
        self.assertIn('-ions', command)
        self.assertIn('-XX:ActiveProcessorCount=1', command)
        self.assertNotIn('-Djava.awt.headless=true', command)

    def test_cpu_default_and_empty_batch(self):
        self.assertGreaterEqual(cpu_count(), 1)
        self.assertEqual(generate_batch([], jobs=2, dry_run=True).outputs, [])


if __name__ == '__main__':
    unittest.main()

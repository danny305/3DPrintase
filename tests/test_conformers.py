"""Shape selection, invariance, option conflicts, and parallel API coverage."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from rdkit import Chem

from printase3d import Molecule, generate_batch, generate_molecule
from printase3d.arguments import parser
from printase3d.chemistry import prepare_molecule
from printase3d.conformers import conformer_flatness, select_conformer, FLATNESS_TOLERANCE_ANGSTROM
from printase3d.options import variant_name
from printase3d.outputs import model_directory
from printase3d.utils import PipelineError


def candidate(points, energy, status=0):
    mol = Chem.MolFromSmiles('CCCC')
    conf = Chem.Conformer(4)
    for index, point in enumerate(points):
        conf.SetAtomPosition(index, tuple(float(v) for v in point))
    mol.AddConformer(conf)
    return dict(mol=mol, conf_id=0, energy=energy, status=status)


PLANE = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])
TETRA = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])


class ConformerTests(unittest.TestCase):
    def test_flatness_is_rotation_translation_invariant(self):
        q, _ = np.linalg.qr(np.array([[1., 2, 3], [4, 2, 1], [2, 3, 5]]))
        for points, expected in ((PLANE, 0), (TETRA, 1)):
            for xyz in (points, points @ q + [10, -4, 8]):
                self.assertAlmostEqual(conformer_flatness(candidate(xyz, 0)['mol']), expected, places=12)

    def test_shape_extrema_then_energy_ties(self):
        flat = candidate(PLANE, 9)
        flat_low = candidate(PLANE, 4)
        thick = candidate(TETRA, 8)
        low_energy = candidate(TETRA * .5, -10)
        candidates = [flat, flat_low, thick, low_energy]
        self.assertIs(select_conformer(candidates), low_energy)
        self.assertIs(select_conformer(candidates, 'most-flat'), flat_low)
        self.assertIs(select_conformer(candidates, 'least-flat'), thick)
        self.assertEqual(select_conformer([flat, candidate(TETRA * 1e-5, -1)], 'most-flat')['energy'], -1)

    def test_nonconverged_extreme_is_excluded(self):
        converged = candidate(TETRA, 4)
        unconverged = candidate(PLANE, -100, status=1)
        self.assertIs(select_conformer([converged, unconverged], 'most-flat'), converged)
        with self.assertRaisesRegex(PipelineError, 'converged'):
            select_conformer([unconverged], 'least-flat')
        self.assertIs(select_conformer([unconverged]), unconverged)

    def test_hydrogens_do_not_change_flatness(self):
        mol = Chem.AddHs(candidate(PLANE, 0)['mol'])
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                mol.GetConformer().SetAtomPosition(atom.GetIdx(), (5., 6., 100.))
        self.assertAlmostEqual(conformer_flatness(mol), 0)

    def test_cli_conflicts_and_separate_folders(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            parser().parse_args(['molecule', '--smiles', 'CCCC', '--most-flat', '--least-flat'])
        paths = set()
        for flags in ([], ['--most-flat'], ['--least-flat']):
            args = parser().parse_args(['molecule', '--smiles', 'CCCC', '--ensemble', *flags])
            paths.add(variant_name(args))
        self.assertEqual(len(paths), 3)
        for conflict in (['--pubchem-geometry'], ['--optimize', 'none'], ['--optimize', 'mmff']):
            args = parser().parse_args(['molecule', '--smiles', 'CCCC', '--most-flat', *conflict])
            with self.assertRaisesRegex(PipelineError, 'cannot be combined'):
                prepare_molecule(Chem.MolFromSmiles('CCCC'), args, 'butane')

    def test_unique_output_allocation_concurrent_and_dry_run(self):
        from concurrent.futures import ThreadPoolExecutor
        with tempfile.TemporaryDirectory() as folder:
            predicted = model_directory(folder, 'thc', 'most-flat', unique=True, dry_run=True)
            self.assertFalse(predicted.exists())
            with ThreadPoolExecutor(max_workers=4) as pool:
                paths = list(pool.map(lambda _: model_directory(folder, 'thc', 'most-flat', unique=True), range(8)))
            self.assertEqual(len(set(paths)), 8)
            self.assertEqual({p.name for p in paths}, {'most-flat', *{f'most-flat-{n}' for n in range(2, 9)}})

    def test_repeat_flatness_run_preserves_previous_files(self):
        with tempfile.TemporaryDirectory() as folder:
            options = dict(output=folder, no_render=True, most_flat=True, conformer_count=2)
            molecule = Molecule(smiles='CCCC', label='butane')
            first = generate_molecule(molecule, **options)
            contents = {p.name: p.read_bytes() for p in first.iterdir()}
            second = generate_molecule(molecule, **options)
            self.assertEqual(second.name, first.name + '-2')
            self.assertEqual(contents, {p.name: p.read_bytes() for p in first.iterdir()})
            meta = json.loads(next(second.glob('*_metadata.json')).read_text())
            self.assertEqual(meta['variant'], second.name)

    def test_thc_maltotriose_parallel_api_flatness_extremes(self):
        fixtures = json.loads((Path(__file__).parent / 'data/flatness_molecules.json').read_text())
        with tempfile.TemporaryDirectory() as folder:
            for flag in ('most_flat', 'least_flat'):
                selected_fixtures = fixtures if flag == 'most_flat' else fixtures[1:]
                entries = [Molecule(smiles=f['smiles'], label=f['label']) for f in selected_fixtures]
                result = generate_batch(entries, jobs=2, output=folder, no_render=True, ph=7,
                                        conformer_count=8, **{flag: True})
                self.assertEqual(result.failures, [])
                for output in result.outputs:
                    meta = json.loads(next(output.glob('*_metadata.json')).read_text())
                    search = meta['conformer_search']
                    self.assertEqual(meta['geometry_mode_requested'], 'ensemble')
                    self.assertEqual(search['selection'], flag.replace('_', '-'))
                    converged = [c for c in search['candidates'] if c['converged']]
                    extreme = (min if flag == 'most_flat' else max)(c['flatness_rms_angstrom'] for c in converged)
                    selected = next(c for c in converged if c['selected'])
                    self.assertLessEqual(abs(selected['flatness_rms_angstrom'] - extreme), FLATNESS_TOLERANCE_ANGSTROM)
                    tied = [c for c in converged if abs(c['flatness_rms_angstrom'] - extreme) <= FLATNESS_TOLERANCE_ANGSTROM]
                    self.assertEqual(selected['energy_kcal_mol'], min(c['energy_kcal_mol'] for c in tied))
                    mol = next(m for m in Chem.SDMolSupplier(str(next(output.glob('*.sdf'))), removeHs=False) if m)
                    self.assertAlmostEqual(conformer_flatness(mol), selected['flatness_rms_angstrom'], places=4)
                    self.assertIn('--' + flag.replace('_', '-'), meta['reproduce'])


if __name__ == '__main__':
    unittest.main()

"""Rendering flags must not alter the optimized molecular geometry."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from printase3d.arguments import parser
from printase3d.chemistry import prepare_molecule
from printase3d.material_validation import validate_element_materials
from printase3d.mesh import semantic_atom_material, write_print_mtl
from printase3d.options import resolve_render_preset, variant_name


class RenderOptionTests(unittest.TestCase):
    def test_explicit_small_spheres_and_monocolor_variant(self):
        args = resolve_render_preset(parser().parse_args([
            'molecule', '--smiles', 'CO', '--small-spheres', '--monocolor-aliphatics']))
        self.assertEqual((args.spacefill, args.wireframe), (45., .4))
        self.assertEqual(variant_name(args), 'small-spheres-monocolor-aliphatics')

    def test_monocolor_palette_preserves_other_elements_and_identity(self):
        mol = Chem.AddHs(Chem.MolFromSmiles('CO'))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'model.obj'
            materials = {}
            lines = ['mtllib model.mtl']
            for atom in mol.GetAtoms():
                name = f'{atom.GetSymbol()}{atom.GetIdx()+1}_Atom'
                material, rgb = semantic_atom_material(name, monocolor_aliphatics=True)
                materials[material] = rgb
                lines += [f'usemtl {material}', f'# component {name}']
                if atom.GetSymbol() != 'H':
                    self.assertEqual((material, rgb), semantic_atom_material(name))
            self.assertEqual(materials['atom_H_hydrogen'], materials['atom_C_carbon'])
            path.write_text('\n'.join(lines)+'\n')
            write_print_mtl(path.with_suffix('.mtl'), materials)
            self.assertTrue(validate_element_materials(path, mol, monocolor_aliphatics=True)['passed'])

    def test_no_hydrogens_keeps_all_atom_optimization_identical(self):
        source = Chem.AddHs(Chem.MolFromSmiles('CCCO'))
        AllChem.EmbedMolecule(source, randomSeed=123)
        modes = [[], ['--optimize', 'hydrogens'], ['--optimize', 'mmff'],
                 ['--optimize', 'uff'], ['--ensemble'], ['--most-flat'],
                 ['--least-flat'], ['--extended'], ['--coordinates', 'regenerate']]
        for mode in modes:
            with self.subTest(mode=mode):
                argv = ['molecule', '--smiles', 'CCCO', '--conformer-count', '2', *mode]
                visible, vmeta = prepare_molecule(Chem.Mol(source), parser().parse_args(argv), 'test')
                hidden, hmeta = prepare_molecule(Chem.Mol(source), parser().parse_args([*argv, '--no-hydrogens']), 'test')
                self.assertGreater(hidden.GetNumAtoms(), hidden.GetNumHeavyAtoms())
                np.testing.assert_array_equal(visible.GetConformer().GetPositions(), hidden.GetConformer().GetPositions())
                self.assertEqual(vmeta['coordinate_method'], hmeta['coordinate_method'])

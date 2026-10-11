"""Portable runtime resolution for source and wheel installations."""
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from printase3d import rendering
from printase3d.cli import main


class InstallationTests(unittest.TestCase):
    def test_explicit_environment_path_takes_precedence(self):
        with patch.dict(os.environ, {'THREEDPRINTASE_JMOL_JAR': '~/custom/Jmol.jar'}):
            self.assertEqual(rendering.resolve_jmol_jar(), str(Path.home() / 'custom/Jmol.jar'))

    def test_source_checkout_finds_local_jar_from_another_cwd(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ, {}, clear=True):
            root = Path(tmp)
            (root / 'pyproject.toml').touch()
            jar = root / 'tools/jmol/Jmol.jar'
            jar.parent.mkdir(parents=True)
            jar.touch()
            with patch.object(rendering, '__file__', str(root / 'src/printase3d/rendering.py')):
                self.assertEqual(rendering.resolve_jmol_jar(), str(jar.resolve()))

    def test_wheel_uses_user_cache_without_searching_ancestors(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ, {'XDG_CACHE_HOME': tmp}, clear=True):
            root = Path(tmp)
            (root / 'Jmol.jar').touch()
            with patch.object(rendering, '__file__', str(root / 'site-packages/printase3d/rendering.py')):
                self.assertEqual(rendering.resolve_jmol_jar(), str(root / '3dprintase/Jmol.jar'))

    def test_no_render_does_not_require_java_or_jmol(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.dict(os.environ, {'THREEDPRINTASE_JMOL_JAR': '/missing/Jmol.jar'}):
                result = main(['molecule', '--smiles', 'CCO', '--label', 'ethanol',
                               '--no-render', '--output', tmp, '--conformer-count', '2'])
            self.assertEqual(result, 0)
            self.assertEqual(len(list(Path(tmp).rglob('*.sdf'))), 1)
            self.assertEqual(len(list(Path(tmp).rglob('*.obj'))), 0)


if __name__ == '__main__':
    unittest.main()

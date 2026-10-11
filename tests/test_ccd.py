import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from rdkit import Chem
from printase3d.ccd import read_ccd,fetch_ccd
from printase3d.arguments import parser
from printase3d.chemistry import prepare_molecule,transfer_heavy_atom_coordinates
from printase3d.api import Molecule,molecule_tasks
from printase3d.extended_geometry import extension_definition,extension_measurement
from printase3d.conformers import sample_lowest_energy_conformer
from printase3d.utils import PipelineError

TEXT=(Path(__file__).parent/'fixtures/ATP.cif').read_text()
class CCDTests(unittest.TestCase):
 def test_real_ccd_ideal_and_model_coordinates(self):
  ideal,info=read_ccd(TEXT);model,_=read_ccd(TEXT,'model')
  self.assertEqual(info['ccd_id'],'ATP')
  self.assertEqual(sum(a.GetAtomicNum()==15 for a in ideal.GetAtoms()),3)
  self.assertFalse(np.allclose(ideal.GetConformer().GetPositions(),model.GetConformer().GetPositions()))
 def test_fetch_and_invalid_identifier(self):
  with patch('printase3d.ccd.request',return_value=TEXT.encode()) as request:
   _,info,raw=fetch_ccd('atp');self.assertEqual(raw,TEXT.encode());self.assertIn('/ATP.cif',info['url'])
  with self.assertRaises(PipelineError):fetch_ccd('../ATP')
 def test_ph7_retains_source_heavy_coordinates(self):
  mol,_=read_ccd(TEXT)
  args=parser().parse_args(['molecule','--ccd','ATP','--rcsb-geometry','--ph','7'])
  result,meta=prepare_molecule(mol,args,'ATP')
  self.assertEqual(meta['formal_charge'],-4)
  mapped=Chem.RemoveHs(result);ok,_=transfer_heavy_atom_coordinates(mol,mapped);self.assertTrue(ok)
  np.testing.assert_allclose(Chem.RemoveHs(result).GetConformer().GetPositions(),mapped.GetConformer().GetPositions(),atol=1e-8)
 def test_cli_and_api(self):
  t=molecule_tasks([Molecule(ccd='ATP')],rcsb_geometry=True,ccd_coordinates='model')[0]
  args=parser().parse_args(t.argv);self.assertEqual(args.geometry_mode,'rcsb');self.assertEqual(args.ccd_coordinates,'model')
  with self.assertRaises(ValueError):Molecule(ccd='ATP',cid=5957)
 def test_extended_single_nucleotide_preserves_rings(self):
  mol,_=read_ccd(TEXT);args=parser().parse_args(['molecule','--ccd','ATP','--extended','--ph','7','--conformer-count','3','--max-opt-iters','10000'])
  out,meta=prepare_molecule(mol,args,'ATP');search=meta['conformer_search']
  self.assertEqual(search['selection'],'extended');self.assertTrue(search['selected_converged'])
  self.assertLess(max(meta['aromatic_planarity']['rms_angstrom']),.03)

 def test_pipeline_saves_ccd_provenance_and_unique_outputs(self):
  import tempfile
  from printase3d.pipeline import process_molecule
  from printase3d.options import resolve_render_preset
  import json
  with tempfile.TemporaryDirectory() as tmp, patch('printase3d.ccd.request',return_value=TEXT.encode()):
   args=resolve_render_preset(parser().parse_args(['molecule','--ccd','ATP','--ph','7','--no-render','--output',tmp]))
   first=process_molecule(args);second=process_molecule(args)
   self.assertNotEqual(first,second)
   meta=json.loads(next(first.glob('*metadata.json')).read_text())
   self.assertEqual(meta['geometry_mode_requested'],'rcsb')
   self.assertTrue((first/meta['files']['ccd_source']).is_file())
   self.assertEqual(len(meta['rcsb_ccd']['source_sha256']),64)

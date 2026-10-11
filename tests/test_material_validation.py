import tempfile,unittest
from pathlib import Path
from rdkit import Chem
from printase3d.mesh import write_print_mtl,ELEMENT_MATERIALS
from printase3d.material_validation import validate_element_materials
from printase3d.utils import PipelineError
class ElementMaterialValidationTests(unittest.TestCase):
 def test_rejects_swapped_oxygen_phosphorus_materials(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'test.obj';m=Chem.MolFromSmiles('OP')
   write_print_mtl(p.with_suffix('.mtl'),{'atom_O_oxygen':ELEMENT_MATERIALS['O'][1],'atom_P_phosphorus':ELEMENT_MATERIALS['P'][1]})
   correct='mtllib test.mtl\nusemtl atom_O_oxygen\n# component O1_Atom\nusemtl atom_P_phosphorus\n# component P2_Atom\n'
   p.write_text(correct);self.assertTrue(validate_element_materials(p,m)['passed'])
   p.write_text(correct.replace('usemtl atom_O_oxygen','usemtl atom_P_phosphorus'))
   with self.assertRaises(PipelineError):validate_element_materials(p,m)
 def test_rejects_incorrect_element_name_even_with_correct_palette(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'test.obj';m=Chem.MolFromSmiles('O')
   write_print_mtl(p.with_suffix('.mtl'),{'atom_P_phosphorus':ELEMENT_MATERIALS['P'][1]})
   p.write_text('mtllib test.mtl\nusemtl atom_P_phosphorus\n# component P1_Atom\n')
   with self.assertRaises(PipelineError):validate_element_materials(p,m)

 def test_detects_distorted_export_coordinates(self):
  from printase3d.material_validation import validate_exported_atom_coordinates
  mol=Chem.MolFromSmiles('CCCC');conf=Chem.Conformer(4)
  xyz=[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
  for i,p in enumerate(xyz):conf.SetAtomPosition(i,p)
  mol.AddConformer(conf)
  with tempfile.TemporaryDirectory() as tmp:
   path=Path(tmp)/'atoms.obj'
   def write(points):path.write_text(''.join(f'# component C{i+1}_Atom\nv {p[0]} {p[1]} {p[2]}\n' for i,p in enumerate(points)))
   write(xyz);self.assertTrue(validate_exported_atom_coordinates(path,mol)['passed'])
   write([*xyz[:3],(0,0,1.5)])
   with self.assertRaises(PipelineError):validate_exported_atom_coordinates(path,mol)

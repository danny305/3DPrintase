"""Validate nucleotide chemistry, source coordinates, geometry and colored meshes."""
from pathlib import Path
import json,csv
import numpy as np
from rdkit import Chem
from printase3d.aromatic_geometry import aromatic_components,ring_planarity
from printase3d.carbonyl_geometry import carbonyl_measurements
from printase3d.chemistry import transfer_heavy_atom_coordinates
from printase3d.ccd import read_ccd
from printase3d.material_validation import validate_element_materials,validate_exported_atom_coordinates
from printase3d.stereochemistry import preserves_stereochemistry
from printase3d.utils import sha256
ROOT=Path(__file__).resolve().parents[1]/'out/nucleotides'

def chemical_key(mol):
 from rdkit.Chem.MolStandardize import rdMolStandardize
 mol=Chem.RemoveHs(mol)
 for a in mol.GetAtoms():
  if a.GetAtomicNum()==15:a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
 return Chem.MolToSmiles(rdMolStandardize.TautomerEnumerator().Canonicalize(rdMolStandardize.Uncharger().uncharge(mol)))

def validate(row,result):
 folder=Path(result['outputs'][0]);p=next(folder.glob('*metadata.json'));d=json.loads(p.read_text());sdf=folder/d['files']['input_sdf'];m=Chem.SDMolSupplier(str(sdf),removeHs=False)[0]
 source=Chem.MolFromSmiles(row['properties']['SMILES']);obj=folder/d['files']['obj']
 rms=ring_planarity(m,aromatic_components(m));carbonyl=carbonyl_measurements(m)
 expected={'base':0,'nucleoside':0,'NMP':-2,'NDP':-3,'NTP':-4,'cyclic':-1}[row['group']]
 observed=Chem.Mol(m);Chem.RemoveStereochemistry(observed);Chem.AssignStereochemistryFrom3D(observed)
 checks={'target_ph7':d['target_ph']==7,'expected_charge':d['formal_charge']==expected,'chemical_identity':chemical_key(m)==chemical_key(source),'stereochemistry':preserves_stereochemistry(m,observed),'aromatic_planarity':all(x<=.03 for x in rms),'carbonyl_planarity':all(x['height_angstrom']<=.03 for x in carbonyl),'sdf_hash':sha256(sdf)==d['input_sha256']}
 colors=validate_element_materials(obj,m);checks['element_colors']=colors['passed']
 checks['exported_coordinates']=validate_exported_atom_coordinates(obj,m)['passed']
 faces=np.array([[int(t.split('/')[0]) for t in l.split()[1:]] for l in obj.read_text().splitlines() if l.startswith('f ')])
 edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1);_,counts=np.unique(edges,axis=0,return_counts=True);checks['mesh_closed']=bool(np.all(counts==2))
 if result['mode'] in ('pubchem','rcsb'):
  src=Chem.SDMolSupplier(str(ROOT/'sources'/(row['label']+'.sdf')),removeHs=False)[0] if result['mode']=='pubchem' else read_ccd((ROOT/'sources'/(row['ccd']+'.cif')).read_text())[0]
  target=Chem.RemoveHs(m);ok,_=transfer_heavy_atom_coordinates(src,target)
  checks['source_coordinates_preserved']=bool(ok and np.max(np.abs(Chem.RemoveHs(m).GetConformer().GetPositions()-target.GetConformer().GetPositions()))<=.00011)
 else:
  search=d['conformer_search'];checks['converged']=search['selected_converged'];checks['ensemble_size']=search['requested_etkdg_conformers']==50
 return {'label':row['label'],'mode':result['mode'],'checks':checks,'charge':d['formal_charge'],'aromatic_rms_angstrom':rms,'maximum_carbonyl_height_angstrom':max((x['height_angstrom'] for x in carbonyl),default=0),'obj':str(obj.relative_to(ROOT)),'metadata':str(p.relative_to(ROOT))}

def main():
 rows=json.loads((ROOT/'manifest.json').read_text());results=json.loads((ROOT/'generation_results.json').read_text());lookup={r['label']:r for r in rows};records=[];issues=[]
 for r in results:
  if r.get('error'):issues.append({'label':r['label'],'mode':r['mode'],'error':r['error']});continue
  try:
   record=validate(lookup[r['label']],r);records.append(record)
   if not all(record['checks'].values()):issues.append({'label':r['label'],'mode':r['mode'],'failed':[k for k,v in record['checks'].items() if not v]})
  except Exception as e:issues.append({'label':r['label'],'mode':r['mode'],'error':str(e)})
 expected=sum(4+int(r['pubchem_3d_available']) for r in rows)
 if len(records)!=expected:issues.append({'expected_models':expected,'actual':len(records)})
 report={'passed':not issues,'compounds':len(rows),'models':len(records),'issues':issues,'records':records}
 (ROOT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
 links={(r['label'],r['mode']):r['obj'] for r in records}
 lines=['# Nucleic-acid building blocks at pH 7','',f"{len(rows)} compounds; {len(records)} ball-and-stick models. This collection covers free bases, nucleosides, 5′ mono-/di-/triphosphates, and cyclic nucleotides. DNA/RNA polymers require a sequence and length and are outside this monomer collection.",'','## Conformers','',
 '- **Ensemble:** lowest restrained force-field energy among 50 requested ETKDG candidates plus the supplied conformer. This is a sampled minimum, not a proven global or aqueous free-energy minimum.',
 '- **Extended:** restrained terminal separation followed by greatest end-to-end distance among converged candidates. Rigid bases have no extendable tail; their rings are not stretched.',
 '- **Most flat:** smallest heavy-atom best-fit-plane RMS among converged minimized candidates; energy breaks shape ties.',
 '- **PubChem:** deposited heavy coordinates preserved; pH-7 hydrogens are adjusted and minimized.',
 '- **RCSB CCD:** ideal heavy coordinates preserved; pH-7 hydrogens are adjusted and minimized. These are computed ideal component conformers, not necessarily a bound conformation from an experiment. The CLI also supports the CCD representative model coordinate set.',
 '', 'All variants use a representative metal-free pH-7 microstate: neutral bases/nucleosides, NMP −2, NDP −3, NTP −4, cyclic NMP −1. Ionization equilibria and metal binding can change these states. Aromatic rings remain locally planar; the sugar rings retain their three-dimensional stereochemistry.',
 '', '## Models','', '| Compound | Group | Ensemble | Extended | Most flat | PubChem | RCSB CCD |','|---|---|---|---|---|---|---|']
 for r in rows:
  cells=[]
  for mode in ['ensemble','extended','most-flat','pubchem','rcsb']:
   path=links.get((r['label'],mode));cells.append(f'[OBJ]({path})' if path else 'Unavailable')
  lines.append('| '+ ' | '.join([r['label'],r['group'],*cells])+' |')
 lines+=['','## Colors and importing','', 'Oxygen is red (RGB 1.000, 0.055, 0.055); phosphorus is orange (1.000, 0.504, 0.000). Every exported atom assignment is checked against its SDF element. Import the linked regular `.obj` with its matching `.mtl` beside it. The `_jmol_raw` files are intermediate exports. Slicer filament assignments can override OBJ materials; retain the supplied element palette.', '', '[NAD+ and FMN checks](nad_fmn_checks.json) verify the current cofactor files. NAD+ minimized variants were additionally refreshed with standalone aromatic-ring restraints. The new restraint coverage includes individual 5/6-membered aromatic rings as well as appropriate fused components.', '', '## Sources and validation','', '[Validation](validation.json) · [Source manifest and CCD mappings](manifest.json) · [Generation results](generation_results.json). Every source CCD graph was checked against its PubChem compound, including stereochemistry; full source CIF/SDF files are retained in `sources/`.', '', '[CCD documentation](https://www.wwpdb.org/data/ccd) distinguishes ideal and model coordinates. Source links use `https://www.rcsb.org/ligand/CODE` and `https://pubchem.ncbi.nlm.nih.gov/compound/CID`; exact codes/CIDs are in the manifest.', '', '```sh','3dprintase molecule --ccd ATP --rcsb-geometry --ph 7 --ball-stick','3dprintase molecule --ccd ATP --rcsb-geometry --ccd-coordinates model --ph 7 --ball-stick','3dprintase molecule --ccd ATP --extended --ph 7 --ball-stick','```']
 unavailable=[r for r in rows if not r['pubchem_3d_available']]
 if unavailable:lines+=['','### PubChem 3D unavailable','', ', '.join(r['label']+' (CID '+r['cid']+')' for r in unavailable)+'. These variants are explicitly unavailable rather than substituting a generated conformer under the PubChem label. All four other variants are provided.']
 (ROOT/'README.md').write_text('\n'.join(lines)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2))
 if issues:raise SystemExit(1)
if __name__=='__main__':main()

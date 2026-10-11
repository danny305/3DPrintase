"""Validate cannabinoid exports and write the tail comparison index."""
from pathlib import Path
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from printase3d.material_validation import validate_element_materials,validate_exported_atom_coordinates
from printase3d.stereochemistry import preserves_stereochemistry
from printase3d.aromatic_geometry import aromatic_components,ring_planarity
from printase3d.utils import sha256
ROOT=Path(__file__).resolve().parents[1]/'out/drugs'
records=[];summaries=[]
for stem in ['thc_16078','cbd_644019']:
 root=ROOT/stem;rows=json.loads((root/'tail_conformers.json').read_text());pool=json.loads((root/'ensemble_candidates.json').read_text())
 assert len(rows)==6
 assert len({r['candidate_id'] for r in rows})==6,'Selected variants include duplicate candidates'
 for row in rows:
  obj=ROOT/row['obj'];p=obj.parent;d=json.loads((p/f'{stem}_metadata.json').read_text());sdf=p/d['files']['input_sdf'];m=Chem.SDMolSupplier(str(sdf),removeHs=False)[0]
  f=np.array([[int(x.split('/')[0]) for x in line.split()[1:]] for line in obj.read_text().splitlines() if line.startswith('f ')])
  e=np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1);_,counts=np.unique(e,axis=0,return_counts=True)
  tail=d['tail_geometry']['atom_indices_0based'];xyz=m.GetConformer().GetPositions();span=float(np.linalg.norm(xyz[tail[-1]]-xyz[tail[0]]));quads=d['tail_geometry']['torsion_atom_indices_0based'];torsions=np.array([rdMolTransforms.GetDihedralDeg(m.GetConformer(),*q) for q in quads]);diff=(torsions-np.array(row['torsions'])+180)%360-180
  reference=Chem.MolFromSmiles(d['protonated_smiles'])
  checks={'identity':Chem.MolToSmiles(Chem.RemoveHs(m))==Chem.MolToSmiles(reference),'stereochemistry':preserves_stereochemistry(m,m),'ph7':d['target_ph']==7,'neutral':d['formal_charge']==0,'converged':d['conformer_search']['selected_converged'],'colors':validate_element_materials(obj,m)['passed'],'coordinates':validate_exported_atom_coordinates(obj,m)['passed'],'aromatic_planarity':all(x<=.03 for x in ring_planarity(m,aromatic_components(m))),'closed_mesh':bool(np.all(counts==2)),'tail_measurements':abs(span-row['tail_span'])<.001 and bool(np.max(np.abs(diff))<.1),'hash':sha256(sdf)==d['input_sha256']}
  values=pool['candidates']
  if row['mode']=='tail-extended':checks['selection']=abs(row['tail_span']-max(x['tail_span'] for x in values))<1e-8
  if row['mode']=='tail-folded':checks['selection']=abs(row['tail_span']-min(x['tail_span'] for x in values))<1e-8
  if row['mode']=='ensemble':checks['selection']=row['energy_above_minimum']==0
  records.append({**row,'checks':checks})
 summaries.append({'compound':rows[0]['compound'],'accepted':pool['accepted'],'requested':pool['requested'],'embedded':pool['embedded'],'rejected':pool['rejected'],'span_range':[min(x['tail_span'] for x in pool['candidates']),max(x['tail_span'] for x in pool['candidates'])]})
passed=all(all(r['checks'].values()) for r in records)
(ROOT/'cannabinoid_validation.json').write_text(json.dumps({'passed':passed,'models':records,'ensembles':summaries},indent=2)+'\n')
lines=['# THC and CBD: ball-and-stick tail conformers','','Twelve models: six each for stereospecific Δ9-THC and cannabidiol. Models retain the pH-7 neutral state, explicit hydrogens, 35% atom radii, 0.30 Å bond radii, and element colors. Keep each OBJ with its matching MTL. SDF, MODEL.md, rendering logs, and metadata are included in every model folder.','','## Selection','','Each molecule requests 250 deterministic ETKDGv3 starts plus the PubChem conformer. Candidates are MMFF94s-minimized with the current pipeline aromatic-planarity protections. Only converged, stereochemistry-preserving candidates passing aromatic geometry checks are accepted. No tail extension or folding restraints are applied.','','- **Ensemble:** lowest sampled force-field energy.','- **Most-flat / least-flat:** whole-molecule heavy-atom plane-RMS extremes.','- **Tail-extended / tail-folded:** greatest / smallest distance from the aromatic ring attachment carbon to the pentyl terminal methyl carbon.','- **Tail-alternative:** maximizes its minimum four-torsion angular difference from the five preceding choices.','','These are sampled minimized geometries; tail extremes may be higher energy. Force-field energies include planarity restraints and are not solution free energies or conformer populations. Whole-molecule flatness is separate from tail extension.','','## Models','','| Compound | Variant | OBJ | Tail span (Å) | ΔE (kcal/mol) | Tail torsions (degrees) |','|---|---|---|---:|---:|---|']
for r in records:lines.append(f"| {r['compound']} | {r['mode']} | [OBJ]({r['obj']}) | {r['tail_span']:.3f} | {r['energy_above_minimum']:.3f} | {', '.join(f'{x:.1f}' for x in r['torsions'])} |")
lines+=['','## Sampling and validation','']
for s in summaries:lines.append(f"- {s['compound']}: {s['accepted']} accepted candidates; tail span {s['span_range'][0]:.3f}–{s['span_range'][1]:.3f} Å.")
lines+=['',f"Validation: {'PASS' if passed else 'FAIL'}. [Full validation report](cannabinoid_validation.json). All six selected candidates per molecule are distinct. The molecule folders retain the full minimized ensemble SDF, candidate table, atom indices, and tail measurements.",'','Sources: [Δ9-THC, PubChem CID 16078](https://pubchem.ncbi.nlm.nih.gov/compound/16078) and [CBD, PubChem CID 644019](https://pubchem.ncbi.nlm.nih.gov/compound/644019). Original source SDF files are preserved.','','Reproduce with the project Python environment (generation expects fresh variant folders):','','```sh','python scripts/generate_cannabinoid_conformers.py','python scripts/index_cannabinoid_conformers.py','```','']
(ROOT/'CANNABINOID_CONFORMERS.md').write_text('\n'.join(lines))
print(json.dumps({'passed':passed,'models':len(records),'ensembles':summaries},indent=2));assert passed

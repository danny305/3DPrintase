"""Generate minimized cannabinoid shape and pentyl-tail variants, with provenance."""
from pathlib import Path
import json, copy
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from printase3d.arguments import parser
from printase3d.options import resolve_render_preset
from printase3d.chemistry import prepare_molecule
from printase3d.sources import pubchem_properties
from printase3d.utils import request, sha256, utc_now
from printase3d.aromatic_geometry import aromatic_components, add_planarity_restraints, ring_planarity
from printase3d.carbonyl_geometry import add_carbonyl_restraints
from printase3d.stereochemistry import preserves_stereochemistry
from printase3d.conformers import conformer_flatness
from printase3d.rendering import render_jmol
from printase3d.metadata import write_model_metadata, REPRODUCTION_COMMAND
from printase3d.config import VERSION
ROOT=Path(__file__).resolve().parents[1]/'out/drugs'
COUNT=250

def generate(label,cid):
 stem=f'{label.lower()}_{cid}';root=ROOT/stem;root.mkdir(exist_ok=True)
 source=root/'source_pubchem.sdf'
 props=pubchem_properties(str(cid))
 source.write_bytes(request(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d'))
 raw=Chem.SDMolSupplier(str(source),removeHs=False)[0]
 args=resolve_render_preset(parser().parse_args(['molecule','--cid',str(cid),'--label',label,'--ph','7','--ball-stick','--optimize','none','--max-opt-iters','10000']))
 mol,chemistry=prepare_molecule(raw,args,label)
 # The five-carbon side chain must be unique; index order is ring attachment -> terminal methyl.
 matches=mol.GetSubstructMatches(Chem.MolFromSmarts('[c]-[CH2;!R]-[CH2;!R]-[CH2;!R]-[CH2;!R]-[CH3;!R]'))
 assert len(matches)==1,matches
 tail=list(matches[0]);neighbor=min(a.GetIdx() for a in mol.GetAtomWithIdx(tail[0]).GetNeighbors() if a.GetIsAromatic())
 path=[neighbor]+tail;quads=[path[i:i+4] for i in range(4)]
 groups=aromatic_components(mol)
 sampled=Chem.Mol(mol);sampled.RemoveAllConformers();params=AllChem.ETKDGv3();params.randomSeed=61453;params.useRandomCoords=True;params.pruneRmsThresh=.1;params.numThreads=1;params.maxIterations=3000
 ids=list(AllChem.EmbedMultipleConfs(sampled,numConfs=COUNT,params=params))
 candidates=[];rejected=0
 starts=[('source',mol.GetConformer())]+[(f'ETKDG-{i}',sampled.GetConformer(i)) for i in ids]
 for origin,conf in starts:
  m=Chem.Mol(mol);m.RemoveAllConformers();m.AddConformer(Chem.Conformer(conf),assignId=True)
  work=Chem.Mol(m);ff=AllChem.MMFFGetMoleculeForceField(work,AllChem.MMFFGetMoleculeProperties(work,mmffVariant='MMFF94s'))
  add_carbonyl_restraints(ff,work,-1,'MMFF94s');add_planarity_restraints(ff,work,-1,groups,'MMFF94s');ff.Initialize();status=ff.Minimize(maxIts=10000)
  m.RemoveAllConformers();m.AddConformer(Chem.Conformer(work.GetConformer()),assignId=True)
  if status!=0 or not preserves_stereochemistry(mol,m) or any(x>.03 for x in ring_planarity(m,groups)):
   rejected+=1;continue
  xyz=m.GetConformer().GetPositions()
  candidates.append({'mol':m,'origin':origin,'energy':float(ff.CalcEnergy()),'flatness':conformer_flatness(m),'tail_span':float(np.linalg.norm(xyz[tail[-1]]-xyz[tail[0]])),'torsions':[rdMolTransforms.GetDihedralDeg(m.GetConformer(),*q) for q in quads]})
 assert len(candidates)>10
 energy=min(candidates,key=lambda c:c['energy'])
 picks=[('ensemble',energy),('most-flat',min(candidates,key=lambda c:(c['flatness'],c['energy']))),('least-flat',max(candidates,key=lambda c:c['flatness'])),('tail-extended',max(candidates,key=lambda c:c['tail_span'])),('tail-folded',min(candidates,key=lambda c:c['tail_span']))]
 def delta(a,b):
  diff=(np.array(a['torsions'])-np.array(b['torsions'])+180)%360-180
  return float(np.sqrt(np.mean(diff**2)))
 remaining=[c for c in candidates if all(c is not p for _,p in picks)]
 alternative=max(remaining,key=lambda c:min(delta(c,p) for _,p in picks))
 picks.append(('tail-alternative',alternative))
 table=[{k:v for k,v in c.items() if k!='mol'} for c in candidates]
 (root/'ensemble_candidates.json').write_text(json.dumps({'requested':COUNT,'embedded':len(ids),'accepted':len(candidates),'rejected':rejected,'tail_atom_indices_0based':tail,'torsion_atom_indices_0based':quads,'candidates':table},indent=2)+'\n')
 with Chem.SDWriter(str(root/'minimized_ensemble.sdf')) as writer:
  for i,c in enumerate(candidates):
   m=Chem.Mol(c['mol']);m.SetProp('candidate_id',str(i));m.SetProp('energy_kcal_mol',str(c['energy']));m.SetProp('tail_span_angstrom',str(c['tail_span']));writer.write(m)
 records=[]
 for mode,c in picks:
  folder=root/f'ball-stick-ph7-{mode}';folder.mkdir(exist_ok=False)
  sdf=folder/f'{stem}.sdf'
  with Chem.SDWriter(str(sdf)) as writer:writer.write(c['mol'])
  obj=folder/f'{stem}.obj';log=folder/f'{stem}_jmol.log';command,mesh=render_jmol(args,sdf,obj,log)
  meta={**copy.deepcopy(chemistry),'pipeline':'3DPrintase','pipeline_version':VERSION,'created_utc':utc_now(),'model_id':str(cid),'molecule_name':label,'pubchem_cid':cid,'molecule_root':str(root),'variant':folder.name,'source_url':f'https://pubchem.ncbi.nlm.nih.gov/compound/{cid}','source_description':f'PubChem CID {cid}; 250-candidate tail-focused ensemble','properties':props,'input_sha256':sha256(sdf),'coordinate_method':'ETKDGv3 ensemble and converged MMFF94s minimization with aromatic planarity restraints','geometry_mode_requested':'ensemble','aromatic_planarity':{'rms_angstrom':ring_planarity(c['mol'],groups)},'conformer_search':{'method':'ETKDGv3 plus source, minimized with aromatic planarity restraints','forcefield':'MMFF94s','selection':mode,'requested_etkdg_conformers':COUNT,'embedded_etkdg_conformers':len(ids),'total_candidates_including_source':len(starts),'converged_candidates':len(candidates),'selected_converged':True,'selected_origin':c['origin'],'selected_energy_kcal_mol':c['energy'],'selected_flatness_rms_angstrom':c['flatness'],'selected_above_minimum_kcal_mol':c['energy']-energy['energy'],'ranking_scope':'force-field potential including planarity restraints, not solution free energy'},'tail_geometry':{'atom_indices_0based':tail,'torsion_atom_indices_0based':quads,'ring_attachment_to_terminal_methyl_angstrom':c['tail_span'],'torsions_degrees':c['torsions'],'selection_note':'Tail span extrema are selected among freely minimized sampled conformers; no extension or folding restraint applied.'},'files':{'input_sdf':sdf.name,'obj':obj.name,'mtl':obj.with_suffix('.mtl').name,'jmol_raw_obj':mesh['jmol_raw_obj'],'jmol_raw_mtl':mesh['jmol_raw_mtl'],'jmol_log':log.name,'metadata':f'{stem}_metadata.json','model_card':'MODEL.md'},'rendering':{'preset':'ball-stick','hydrogens_visible':True,'spacefill_percent_vdw':args.spacefill,'wireframe_radius_angstrom':args.wireframe,'bond_color_mode':'uniform','multiple_bond_cylinders':False},'print_mesh':mesh,'jmol':{'command':command,'script_sha256':sha256(args.script),'scale':args.scale}}
  token=REPRODUCTION_COMMAND.set('python scripts/generate_cannabinoid_conformers.py')
  write_model_metadata(folder,f'{stem}_metadata.json',meta);REPRODUCTION_COMMAND.reset(token)
  with (folder/'MODEL.md').open('a') as f:f.write(f"\n## Pentyl tail\n\n- Ring attachment to terminal methyl: {c['tail_span']:.3f} Å\n- Four torsions, starting at the ring: {', '.join(f'{x:.1f}' for x in c['torsions'])} degrees\n- Selection: {mode}; no imposed tail restraint.\n")
  records.append({'compound':label,'mode':mode,'obj':str(obj.relative_to(ROOT)),'tail_span':c['tail_span'],'torsions':c['torsions'],'energy_above_minimum':c['energy']-energy['energy'],'candidate_id':next(i for i,x in enumerate(candidates) if x is c)})
  print(label,mode,round(c['tail_span'],3),flush=True)
 (root/'tail_conformers.json').write_text(json.dumps(records,indent=2)+'\n')
 return records

if __name__=='__main__':
 for label,cid in [('THC',16078),('CBD',644019)]:generate(label,cid)

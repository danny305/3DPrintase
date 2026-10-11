"""Audit cofactor geometry/meshes and publish a complete availability index."""
import csv,json
from pathlib import Path
import numpy as np
from rdkit import Chem
from printase3d.carbonyl_geometry import carbonyl_measurements
from printase3d.chemistry import transfer_heavy_atom_coordinates
from printase3d.conformers import conformer_flatness,FLATNESS_TOLERANCE_ANGSTROM
from printase3d.structures import parse_cif,loop_rows
from printase3d.extended_geometry import extension_measurement
from printase3d.aromatic_geometry import aromatic_components, ring_planarity, flavin_geometry
ROOT=Path(__file__).resolve().parents[1]/'out/cofactors'
CHARGES={'NAD-oxidized':-1,'NADH':-2,'NADP-oxidized':-3,'NADPH':-4,'FAD':-2,'FADH2':-2,'FMN':-2,'FMNH2':-2,'coenzyme-A':-4,'acetyl-CoA':-4,'PLP':-2,'PMP':-1,'TPP':-2,'biotin':-1,'lipoamide':0,'dihydrolipoamide':0,'THF':-2,'5-methyl-THF':-2,'5-10-methylene-THF':-2,'SAM':1,'BH4':0,'ascorbate':-1,'CoQ10':0,'ubiquinol-10':0,'vitamin-K1-hydroquinone':0,'heme-b-ferrous':-2,'heme-b-ferric':-1,'heme-a':-2,'siroheme':-8,'chlorophyll-a':0,'chlorophyll-b':0,'methylcobalamin':0,'adenosylcobalamin':0,'aquacobalamin':1,'cyanocobalamin':0}
NO_PUBCHEM={'FAD','FADH2','coenzyme-A','acetyl-CoA','CoQ10','ubiquinol-10'}
def main():
 records=[]
 for path in sorted(ROOT.glob('*/*/*/*metadata.json')):
  meta=json.loads(path.read_text()); name=meta['molecule_name']; files=meta['files']; obj=path.parent/files['obj']
  faces=np.array([[int(x.split('/')[0]) for x in l.split()[1:]] for l in obj.read_text().splitlines() if l.startswith('f ')])
  edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1)
  _,counts=np.unique(edges,axis=0,return_counts=True)
  checks={'closed_mesh':bool(np.all(counts==2)),'pH7':meta['target_ph']==7,'expected_charge':meta['formal_charge']==CHARGES[name],'ball_stick':meta['rendering']['preset']=='ball-stick','materials':obj.with_suffix('.mtl').exists()}
  metal=meta.get('geometry_mode')=='rcsb-reference'
  if metal:
   mode='rcsb-reference'; energy=None
   lines=(path.parent/files['input_sdf']).read_text().splitlines(); n=int(lines[3][:3]); nb=int(lines[3][3:6])
   xyz=np.array([[float(l[i:i+10]) for i in (0,10,20)] for l in lines[4:4+n]])
   d,loops=parse_cif((ROOT/'sources'/f"{meta['ccd']}.cif").read_text()); ats=loop_rows(loops,'_chem_comp_atom.')
   source={a['_chem_comp_atom.atom_id']:a for a in ats}
   diffs=[]
   for i,atom in enumerate(meta['atom_order']):
    if meta['atom_elements'][i]!='H':
     row=source[atom]; diffs.extend(xyz[i]-[float(row['_chem_comp_atom.pdbx_model_Cartn_'+v+'_ideal']) for v in ('x','y','z')])
   checks['reference_heavy_coordinates_preserved']=bool(np.max(np.abs(diffs))<1e-4)
   checks['all_atoms_rendered']=meta['print_mesh']['atom_parts']==n
   checks['all_bonds_rendered']=meta['print_mesh']['bond_parts']==nb
   checks['metal_distances']=all(1.5<b['distance']<2.8 for b in meta['coordination_bonds'])
   checks['metal_color']=any('atom_'+e+'_' in obj.with_suffix('.mtl').read_text() for e in ('Fe','Mg','Co'))
   flat=None
  else:
   mol=next(m for m in Chem.SDMolSupplier(str(path.parent/files['input_sdf']),removeHs=False) if m)
   checks['carbonyl_amide_planarity']=all(x['height_angstrom']<=.03 for x in carbonyl_measurements(mol))
   checks['sdf_charge']=Chem.GetFormalCharge(mol)==CHARGES[name]
   orig=Chem.Mol(mol); Chem.AssignStereochemistry(orig,cleanIt=True,force=True)
   observed=Chem.Mol(mol); Chem.RemoveStereochemistry(observed); Chem.AssignStereochemistryFrom3D(observed)
   checks['stereochemistry_3d']=all(observed.GetAtomWithIdx(a.GetIdx()).HasProp('_CIPCode') and observed.GetAtomWithIdx(a.GetIdx()).GetProp('_CIPCode')==a.GetProp('_CIPCode') for a in orig.GetAtoms() if a.HasProp('_CIPCode') and a.GetAtomicNum()!=15)
   search=meta.get('conformer_search'); mode='pubchem' if meta['geometry_mode_requested']=='pubchem' else (search['selection'] if search['selection'] in ('most-flat','extended') else 'lowest-energy')
   energy=search['selected_energy_kcal_mol'] if search else None; flat=conformer_flatness(mol)
   if mode=='pubchem':
    src=next(m for m in Chem.SDMolSupplier(str(ROOT/'sources'/f"pubchem_{meta['pubchem_cid']}_3d.sdf"),removeHs=False) if m)
    target=Chem.RemoveHs(mol); actual=np.array(target.GetConformer().GetPositions()); ok,_=transfer_heavy_atom_coordinates(src,target)
    checks['pubchem_heavy_coordinates_preserved']=ok and bool(np.max(np.abs(actual-target.GetConformer().GetPositions()))<1e-4)
   else:
    pool=[c for c in search['candidates'] if c['converged']]
    groups=aromatic_components(mol)
    checks['aromatic_ring_planarity']=all(rms<.03 for rms in ring_planarity(mol,groups))
    checks['converged']=search['selected_converged']; checks['50_candidates_requested']=search['requested_etkdg_conformers']==50
    checks['flatness_matches']=abs(flat-search['selected_flatness_rms_angstrom'])<1e-4
    if mode=='most-flat':
     extreme=min(c['flatness_rms_angstrom'] for c in pool)
     checks['flatness_extreme']=abs(search['selected_flatness_rms_angstrom']-extreme)<=FLATNESS_TOLERANCE_ANGSTROM
     pool=[c for c in pool if abs(c['flatness_rms_angstrom']-extreme)<=FLATNESS_TOLERANCE_ANGSTROM]
    if mode=='extended':
     ext=search['extension']; sep=extension_measurement(mol,ext)
     aa,bb=ext['ring_groups']; xyz=mol.GetConformer().GetPositions()
     checks['terminal_rings_not_stacked']=bool(np.min(np.linalg.norm(xyz[aa][:,None,:]-xyz[bb][None,:,:],axis=2))>6.0)
     checks['terminal_separation']=sep>=ext['minimum_ring_centroid_distance_angstrom']
     checks['separation_matches_sdf']=abs(sep-ext['selected_terminal_separation_angstrom'])<1e-4
     extreme=max(c['terminal_separation_angstrom'] for c in pool)
     checks['maximum_separation_selected']=abs(ext['selected_terminal_separation_angstrom']-extreme)<1e-4
     pool=[c for c in pool if extreme-c['terminal_separation_angstrom']<=1e-4]
    checks['energy_selection']=energy==min(c['energy_kcal_mol'] for c in pool)
  records.append({'name':name,'mode':mode,'charge':meta['formal_charge'],'obj':str(obj.relative_to(ROOT)),'energy_kcal_mol':energy,'flatness_rms_angstrom':flat,'checks':checks})
 organic=[r['label'] for r in json.loads((ROOT/'organic.json').read_text())]
 metal_names=[n for n in CHARGES if n not in organic]
 expected={(n,m) for n in organic for m in ('lowest-energy','most-flat','pubchem') if not (m=='pubchem' and n in NO_PUBCHEM)}|{(n,'rcsb-reference') for n in metal_names}
 expected |= {(r['label'],'extended') for r in json.loads((ROOT/'dinucleotides.json').read_text())}
 found={(r['name'],r['mode']) for r in records}; issues=[{'name':r['name'],'mode':r['mode'],'failed':[k for k,v in r['checks'].items() if not v]} for r in records if not all(r['checks'].values())]
 report={'passed':found==expected and not issues and len(found)==len(records),'model_count':len(records),'missing':sorted(expected-found),'issues':issues,'models':records}
 (ROOT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
 with (ROOT/'index.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=[k for k in records[0] if k!='checks']); w.writeheader(); w.writerows({k:v for k,v in r.items() if k!='checks'} for r in records)
 lines=['# Cofactor models','',f'Approved collection: **35 compounds, {len(records)} ball-and-stick OBJ models**. Target pH 7.','',
 '## Organic cofactors','', 'Lowest energy and most-flat use 50 requested ETKDG candidates, minimized with MMFF94s for corrected models (earlier unaffected models may use MMFF94), including local carbonyl/amide and aromatic planarity restraints for small fused fully aromatic ring systems. PubChem variants preserve available source heavy coordinates. These are sampled force-field minima, not guaranteed global or aqueous free-energy minima. Most-flat selects the flattest minimized candidate, with energy breaking shape ties.','',
 '| Compound | Charge | Lowest energy | PubChem | Most flat |','|---|---:|---|---|---|']
 for n in organic:
  row={r['mode']:r for r in records if r['name']==n}
  links=[f"[OBJ]({row[m]['obj']})" if m in row else 'Unavailable' for m in ('lowest-energy','pubchem','most-flat')]
  lines.append(f'| {n} | {CHARGES[n]:+d} | '+' | '.join(links)+' |')
 lines += ['', '## Extended dinucleotides', '', '[Before/after distances and model links](EXTENDED.md).', '', 'All six dinucleotides also have elongated, unstacked variants. `--extended` minimizes an ensemble with a terminal-distance restraint, selects the greatest terminal-ring separation, and uses energy to break separation ties. These are restrained extended conformers; energies include that restraint.', '', '| Dinucleotide | Extended model |', '|---|---|']
 for r in records:
  if r['mode']=='extended': lines.append(f"| {r['name']} | [OBJ]({r['obj']}) |")
 lines+=['','The six unavailable organic PubChem variants have no usable 3D record after repeat requests; saved server responses are in `sources`. Their ensemble and most-flat models use generated coordinates.','',
 '## Metal-containing cofactors','', 'These ten models use **RCSB CCD ideal reference geometry**, with all metal-coordination bonds preserved. They are not PubChem variants or energy minima. Lowest-energy and most-flat selection are unsupported because this pipeline has no validated coordination-specific force-field parameterization.','', '| Compound | Charge | Reference model |','|---|---:|---|']
 for n in metal_names:
  row=next((r for r in records if r['name']==n),None)
  link=f"[OBJ]({row['obj']})" if row else 'Unavailable'
  lines.append(f'| {n} | {CHARGES[n]:+d} | {link} |')
 lines+=['','## Chemistry and interpretation','',
 '- pH sets organic-group protonation, not the oxidation state of a metal. Heme b/a and siroheme are modeled as Fe(II), with a separate Fe(III) heme b. Chlorophyll contains Mg(II); cobalamins use Co(III). Axial protein ligands are omitted.',
 '- The aquacobalamin entry represents the predominantly water-ligated member of the approved hydroxo/aqua pair at pH 7. It retains hydroxocobalamin CCD heavy coordinates and adds one axial-water hydrogen geometrically; no relaxation is claimed. [Measured acid dissociation](https://pubs.rsc.org/en/content/articlelanding/1974/dt/dt9740002205) supports this assignment.',
 '- Methylcobalamin and cyanocobalamin retain the source CCD five-coordinate, base-off structures; adenosylcobalamin and aquacobalamin retain six-coordinate structures. These isolated conformations do not specify the protein-bound state.',
 '- Exact-structure corrections remove spurious reduced-ring nitrogen protonation in NADH/NADPH, THF derivatives and BH4. Sources and corrected sites are recorded in metadata. Examples: [NADH dianion](https://www.ebi.ac.uk/chebi/CHEBI:57945), [NADPH tetraanion](https://www.ebi.ac.uk/chebi/CHEBI:57783), [THF dianion](https://www.ebi.ac.uk/chebi/CHEBI:57453).',
 '- Oxidized FMN/FAD isoalloxazine cores are validated within 0.03 Å RMS of a plane after minimization. Fully reduced FMNH2/FADH2 bridge nitrogens are excluded from the full-core planar restraint, permitting butterfly bending. Aromatic subrings remain present; the complete reduced core is not treated as an aromatic plane. [Redox-dependent structures](https://elifesciences.org/articles/79226).',
 '- FADH2 and FMNH2 use the specified fully reduced hydroquinone representations with phosphate charges; a single model does not describe all solution tautomers or mixed protonation populations. Biotin and lipoamide are free representations of protein-attached cofactors. Folates use one glutamate; quinones and chlorophyll retain their full tails.',
 '- Preserve the adjacent MTL file for atom colors. Iron, magnesium and cobalt have distinct materials. Each atom sphere and bond cylinder is closed; intersecting parts are not a boolean-unioned solid.',
 '', '## Files and validation','', '[Carbonyl and amide correction](../carbonyl_audit/README.md) · [Aromatic geometry fix and before/after measurements](AROMATIC_FIX.md) · [Validation report](validation.json) · [CSV index](index.csv) · [Approved list](COFACTOR_REVIEW.md) · [Organic CID manifest](organic_resolved.json) · [Metal results](metal_summary.json).',
 '', 'Models are grouped under `organic/` and `metal_containing/`, with separate variant folders. Earlier rejected chemistry and failed rendering attempts are archived outside this collection.',
 '', 'Regeneration uses `3dprintase batch out/cofactors/organic_resolved.json --jobs 3 --ph 7 --ball-stick --ensemble --output NEW_OUTPUT_ROOT`; replace `--ensemble` with `--most-flat` or `--pubchem-geometry` for the other variants. The metal reference workflow is `scripts/render_metal_cofactors.py`; it calls the modular `printase3d.metal_cofactors` renderer.']
 (ROOT/'README.md').write_text('\n'.join(lines)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k!='models'},indent=2))
 return 0 if report['passed'] else 1
if __name__=='__main__': raise SystemExit(main())

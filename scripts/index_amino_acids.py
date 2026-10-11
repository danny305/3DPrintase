"""Validate and index the amino-acid collection's three variants."""
import csv
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem
from printase3d.carbonyl_geometry import carbonyl_measurements
from rdkit.Chem import rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from printase3d.chemistry import transfer_heavy_atom_coordinates
from printase3d.conformers import FLATNESS_TOLERANCE_ANGSTROM, conformer_flatness

AA = {
    'alanine': ('Ala', 'A'), 'arginine': ('Arg', 'R'), 'asparagine': ('Asn', 'N'),
    'aspartic-acid': ('Asp', 'D'), 'cysteine': ('Cys', 'C'), 'glutamic-acid': ('Glu', 'E'),
    'glutamine': ('Gln', 'Q'), 'glycine': ('Gly', 'G'), 'histidine': ('His', 'H'),
    'isoleucine': ('Ile', 'I'), 'leucine': ('Leu', 'L'), 'lysine': ('Lys', 'K'),
    'methionine': ('Met', 'M'), 'phenylalanine': ('Phe', 'F'), 'proline': ('Pro', 'P'),
    'serine': ('Ser', 'S'), 'threonine': ('Thr', 'T'), 'tryptophan': ('Trp', 'W'),
    'tyrosine': ('Tyr', 'Y'), 'valine': ('Val', 'V'),
}


def main():
    root = ROOT / 'out/amino_acids'
    sources = {int(m.GetProp('PUBCHEM_COMPOUND_CID')): m for m in
               Chem.SDMolSupplier(str(root / 'sources/pubchem_3d.sdf'), removeHs=False) if m}
    for m in Chem.SDMolSupplier(str(root / 'sources/noncanonical_pubchem_3d.sdf'), removeHs=False):
        if m:
            sources[int(m.GetProp('PUBCHEM_COMPOUND_CID'))] = m
    for m in Chem.SDMolSupplier(str(root / 'sources/pyrrolysine_pubchem_3d.sdf'), removeHs=False):
        if m:
            sources[int(m.GetProp('PUBCHEM_COMPOUND_CID'))] = m
    proposal = json.loads((root / 'noncanonical_proposal.json').read_text())
    nca = {r['compound']: r for r in proposal if r['status'].startswith('generated')}
    charges = {'L-ornithine': 1, '(2S,5R)-5-hydroxylysine': 1,
               'N6,N6,N6-trimethyl-L-lysine': 1, 'D-glutamic acid': -1,
               'L-selenocysteine': -1, 'gamma-carboxy-L-glutamic acid': -2,
               'O-phospho-L-serine': -2, 'O-phospho-L-threonine': -2, 'O-phospho-L-tyrosine': -2}
    stereo = {name: ['S'] for name in nca}
    for name in ('beta-alanine', 'GABA', 'Taurine', 'Sarcosine', '2-aminoisobutyric acid'):
        stereo[name] = []
    for name in ('D-serine', 'D-alanine', 'D-glutamic acid', 'L-selenocysteine'):
        stereo[name] = ['R']
    for name in ('trans-4-hydroxy-L-proline', '(2S,5R)-5-hydroxylysine', 'O-phospho-L-threonine'):
        stereo[name] = ['R', 'S']
    stereo['L-pyrrolysine'] = ['R', 'R', 'S']
    stereo['L-cystine'] = ['R', 'R']
    records = []
    for path in sorted(root.glob('*/*/*_metadata.json')):
        meta = json.loads(path.read_text())
        name = meta['molecule_name']
        if name not in AA and name not in nca:
            continue
        search = meta.get('conformer_search')
        mode = 'pubchem' if meta['geometry_mode_requested'] == 'pubchem' else (
            'most-flat' if search['selection'] == 'most-flat' else 'lowest-energy')
        mol = next(m for m in Chem.SDMolSupplier(str(path.parent / meta['files']['input_sdf']), removeHs=False) if m)
        obj = path.parent / meta['files']['obj']
        faces = np.array([[int(v.split('/')[0]) for v in line.split()[1:]]
                          for line in obj.read_text().splitlines() if line.startswith('f ')])
        edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
        _, counts = np.unique(edges, axis=0, return_counts=True)
        original = Chem.Mol(mol)
        Chem.AssignStereochemistry(original, cleanIt=True, force=True)
        observed = Chem.Mol(mol)
        Chem.RemoveStereochemistry(observed)
        Chem.AssignStereochemistryFrom3D(observed)
        expected_charge = {'arginine': 1, 'lysine': 1, 'aspartic-acid': -1, 'glutamic-acid': -1}.get(name, charges.get(name, 0))
        alpha_matches = original.GetSubstructMatches(Chem.MolFromSmarts('[#7]-[C;X4]-[C](=[O])-[#8]'))
        alpha_codes = [original.GetAtomWithIdx(match[1]).GetProp('_CIPCode')
                       if original.GetAtomWithIdx(match[1]).HasProp('_CIPCode') else None
                       for match in alpha_matches]
        checks = {
            'carbonyl_amide_planarity': all(x['height_angstrom']<=.03 for x in carbonyl_measurements(mol)),
            'pH7': meta['target_ph'] == 7,
            'expected_charge': meta['formal_charge'] == Chem.GetFormalCharge(mol) == expected_charge,
            'formula': meta['molecular_formula'] == rdMolDescriptors.CalcMolFormula(mol),
            'natural_alpha_stereochemistry': alpha_codes == [None if name == 'glycine' else 'R' if name == 'cysteine' else 'S'],
            'ball_stick': meta['rendering']['preset'] == 'ball-stick',
            'closed_mesh': bool(np.all(counts == 2)),
            'materials_present': obj.with_suffix('.mtl').is_file(),
            'stereochemistry_3d': all(observed.GetAtomWithIdx(a.GetIdx()).HasProp('_CIPCode') and
                observed.GetAtomWithIdx(a.GetIdx()).GetProp('_CIPCode') == a.GetProp('_CIPCode')
                for a in original.GetAtoms() if a.HasProp('_CIPCode')),
        }
        if name in nca:
            checks.pop('natural_alpha_stereochemistry')
            checks['expected_stereoisomer'] = sorted(a.GetProp('_CIPCode') for a in original.GetAtoms() if a.HasProp('_CIPCode')) == stereo[name]
        if name == 'L-selenocysteine':
            checks['selenolate'] = any(a.GetAtomicNum() == 34 and a.GetFormalCharge() == -1 for a in mol.GetAtoms())
            checks['selenium_material'] = 'atom_Se_selenium' in obj.with_suffix('.mtl').read_text()
        if mode == 'pubchem':
            target = Chem.RemoveHs(mol)
            actual = np.array(target.GetConformer().GetPositions())
            ok, detail = transfer_heavy_atom_coordinates(sources[meta['pubchem_cid']], target)
            checks['pubchem_heavy_coordinates_unchanged'] = ok and bool(np.max(np.abs(actual - target.GetConformer().GetPositions())) < 1e-4)
        else:
            pool = [c for c in search['candidates'] if c['converged']]
            checks['converged'] = search['selected_converged']
            checks['50_requested'] = search['requested_etkdg_conformers'] == 50
            checks['flatness_matches_sdf'] = abs(conformer_flatness(mol) - search['selected_flatness_rms_angstrom']) < 1e-4
            if mode == 'most-flat':
                extreme = min(c['flatness_rms_angstrom'] for c in pool)
                checks['shape_extreme'] = abs(search['selected_flatness_rms_angstrom'] - extreme) <= FLATNESS_TOLERANCE_ANGSTROM
                pool = [c for c in pool if abs(c['flatness_rms_angstrom'] - extreme) <= FLATNESS_TOLERANCE_ANGSTROM]
            checks['minimum_energy_in_selection_pool'] = search['selected_energy_kcal_mol'] == min(c['energy_kcal_mol'] for c in pool)
        records.append({'name': name, 'mode': mode, 'cid': meta['pubchem_cid'],
                        'charge': meta['formal_charge'], 'obj': str(obj.relative_to(root)),
                        'energy_kcal_mol': search['selected_energy_kcal_mol'] if search else None,
                        'flatness_rms_angstrom': conformer_flatness(mol), 'checks': checks})
    expected = {(name, mode) for name in list(AA) + list(nca) for mode in ('lowest-energy', 'pubchem', 'most-flat')}
    expected.remove(('L-selenocysteine', 'pubchem'))
    found = {(r['name'], r['mode']) for r in records}
    issues = [{'name': r['name'], 'mode': r['mode'], 'failed': [k for k,v in r['checks'].items() if not v]}
              for r in records if not all(r['checks'].values())]
    passed = found == expected and not issues
    report = {'passed': passed, 'model_count': len(records), 'unavailable': [{'name': 'L-selenocysteine', 'mode': 'pubchem', 'reason': 'PubChem CID 25076 has no 3D record (HTTP 404)'}], 'missing': sorted(expected - found), 'issues': issues, 'models': records}
    (root / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    with (root / 'index.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=[k for k in records[0] if k != 'checks'])
        writer.writeheader()
        writer.writerows({k:v for k,v in r.items() if k != 'checks'} for r in records)
    lines = ['# Amino-acid models', '',
             '20 standard + 34 noncanonical amino acids and related compounds: **161 ball-and-stick models**, prepared at **pH 7**. All three variants are included except PubChem geometry for selenocysteine, which has no deposited 3D record. Chiral standard amino acids use their natural L forms; glycine is achiral.', '',
             '## Conformer choices', '',
             '- **Lowest energy:** lowest converged force-field energy from a 50-candidate ensemble (`--ensemble`; near-duplicates are pruned).',
             '- **PubChem:** preserves deposited heavy-atom coordinates while assigning pH-7 protonation and relaxing hydrogens (`--pubchem-geometry`).',
             '- **Most flat:** energy-minimizes the ensemble, then selects the flattest converged heavy-atom structure; energy breaks shape ties (`--most-flat`).', '',
             'These are free amino acids, usually zwitterions, rather than residues embedded in a peptide. Protonation is predicted; ensemble minima are sampled force-field minima, not guaranteed global or aqueous free-energy minima.', '',
             'Each link below opens the printable OBJ. Keep its adjacent MTL for atom colors. Each variant has its own folder, SDF, metadata, and model card.', '',
             '## Standard amino acids', '',
             '| Amino acid | Code | Charge at pH 7 | Lowest energy | PubChem | Most flat |',
             '|---|---|---:|---|---|---|']
    for name, (three, one) in AA.items():
        row = {r['mode']: r for r in records if r['name'] == name}
        links = [f"[OBJ]({row[mode]['obj']})" if mode in row else 'Pending' for mode in ('lowest-energy', 'pubchem', 'most-flat')]
        charge = next(iter(row.values()))['charge'] if row else ''
        lines.append(f"| {name.replace('-', ' ').title()} | {three} / {one} | {charge:+d} | " + ' | '.join(links) + ' |')
    lines += ['', '## Noncanonical amino acids', '',
              '[Review and expand the grouped list](NONCANONICAL_REVIEW.md). Generated: 16 natural/biological examples, 7 modified residues, and 10 engineering building blocks. Pyrrolysine is also included as an additional genetically encoded amino acid.', '']
    for group in dict.fromkeys(r['group'] for r in nca.values()):
        lines += [f'### {group}', '', '| Compound | Charge | Lowest energy | PubChem | Most flat |', '|---|---:|---|---|---|']
        for name, entry in nca.items():
            if entry['group'] != group:
                continue
            row = {r['mode']: r for r in records if r['name'] == name}
            links = [f"[OBJ]({row[mode]['obj']})" if mode in row else 'Unavailable' for mode in ('lowest-energy', 'pubchem', 'most-flat')]
            lines.append(f"| {name} | {charges.get(name, 0):+d} | " + ' | '.join(links) + ' |')
        lines.append('')
    lines += ['Selenocysteine uses the selenium-deprotonated form (net −1), supported by [measured selenol pKa 5.24](https://pubmed.ncbi.nlm.nih.gov/35659493/). Its minimized structures use UFF, with selenium charge-state parameter warnings; these are approximate force-field geometries. PubChem CID 25076 returned no 3D record, so its PubChem variant is unavailable.', '',
              '## Verification and reproducibility', '',
              '[Carbonyl and amide correction](../carbonyl_audit/README.md). [Validation report](validation.json) checks pH, expected net charge, formula, 3D stereochemistry, mesh closure, conformer selection, and preservation of PubChem heavy-atom coordinates against the saved [canonical source SDF](sources/pubchem_3d.sdf) and [noncanonical source SDF](sources/noncanonical_pubchem_3d.sdf), and [pyrrolysine source SDF](sources/pyrrolysine_pubchem_3d.sdf).', '',
              'Sources: each model card links its PubChem record; [resolved generation manifest](canonical_resolved.json) and [noncanonical manifest](noncanonical_resolved.json), and [CSV index](index.csv).', '',
              'Generation commands (from the repository root):', '', '```bash',
              '3dprintase batch out/amino_acids/canonical_resolved.json --jobs 2 --ph 7 --ensemble --ball-stick --output out/amino_acids',
              '3dprintase batch out/amino_acids/canonical_resolved.json --jobs 2 --ph 7 --pubchem-geometry --ball-stick --output out/amino_acids',
              '3dprintase batch out/amino_acids/canonical_resolved.json --jobs 2 --ph 7 --most-flat --ball-stick --output out/amino_acids',
              '```', '',
              'Most-flat reruns receive numbered suffixes. To preserve all existing variants on a fresh full rerun, choose a new output root.']
    (root / 'README.md').write_text('\n'.join(lines) + '\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'models'}, indent=2))
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())

"""Cache auditable PubChem and RCSB cofactor source records; no minimization.

Existing downloads, including explicit PubChem no-record responses, are kept.
Remove a specific cached response manually when an intentional refresh is wanted.
"""
import json
import subprocess
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[1] / 'out/cofactors'
METALS = [('heme-b-ferrous', 'heme'), ('heme-b-ferric', 'ferriheme b'),
          ('heme-a', 'heme a'), ('siroheme', 'siroheme'),
          ('chlorophyll-a', 'chlorophyll a'), ('chlorophyll-b', 'chlorophyll b'),
          ('methylcobalamin', 'methylcobalamin'), ('adenosylcobalamin', 'adenosylcobalamin'),
          ('hydroxocobalamin', 'hydroxocobalamin'), ('cyanocobalamin', 'cyanocobalamin')]
CCD_IDS = ['HEM', 'HEA', 'SRM', 'CLA', 'CHL', 'COB', 'B1Z', 'I2A', 'CNC']
BASE = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound'


def fetch(url):
    return subprocess.check_output(['curl', '-L', '-sS', '--retry', '3', url])


def cached(path, url):
    if not path.exists():
        path.write_bytes(fetch(url))
    return path.read_bytes()


def main():
    folder = ROOT / 'sources'
    folder.mkdir(exist_ok=True)
    for row in json.loads((ROOT / 'organic_resolved.json').read_text()):
        cid = row['cid']
        cached(folder / f'pubchem_{cid}_3d.sdf', f'{BASE}/cid/{cid}/SDF?record_type=3d')
    for code in CCD_IDS:
        cached(folder / f'{code}.cif', f'https://files.rcsb.org/ligands/download/{code}.cif')
    manifest = ROOT / 'metal_sources.json'
    if manifest.exists():
        return
    results = []
    for label, name in METALS:
        cid = json.loads(fetch(f'{BASE}/name/{quote(name, safe="")}/cids/JSON'))['IdentifierList']['CID'][0]
        row = {'label': label, 'query': name, 'cid': cid}
        for kind in ('2d', '3d'):
            data = cached(folder / f'{label}_{kind}.sdf', f'{BASE}/cid/{cid}/SDF?record_type={kind}')
            row[kind] = b'V2000' in data or b'V3000' in data
        results.append(row)
    manifest.write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()

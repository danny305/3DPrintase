"""Sources components of the 3DPrintase pipeline."""
from __future__ import annotations

import json
import sys
import urllib.request
from .config import PRESETS, PUBCHEM_JSON, PUBCHEM_NAME, PUBCHEM_SDF, RCSB_CIF, RCSB_SEARCH
from .utils import PipelineError, request


def search_query(args):
    preset = PRESETS[args.type]
    nodes = [
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.polymer_entity_count_protein", "operator": "greater", "value": 0}},
        {"type": "group", "logical_operator": "or", "nodes": [
            {"type": "terminal", "service": "full_text", "parameters": {"value": keyword}}
            for keyword in preset["keywords"]
        ]},
    ]
    return {"query": {"type": "group", "logical_operator": "and", "nodes": nodes}, "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": args.limit}, "sort": [{"sort_by": "score", "direction": "desc"}]}}


def do_search(args):
    payload = search_query(args)
    if args.dry_run:
        print(json.dumps(payload, indent=2)); return []
    raw = request(RCSB_SEARCH, data=payload)
    if not raw.strip():
        print("No matching entries.", file=sys.stderr); return []
    result = json.loads(raw)
    ids = [x["identifier"].lower() for x in result.get("result_set", [])]
    for pdb in ids: print(pdb)
    return ids


def download_cif(pdb, path, dry_run=False):
    url = RCSB_CIF.format(pdb=pdb.upper())
    if dry_run:
        print(f"WOULD DOWNLOAD {url} -> {path}"); return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(request(url))


def pubchem_properties(cid):
    raw = request(PUBCHEM_JSON.format(cid=cid))
    return json.loads(raw)["PropertyTable"]["Properties"][0]


def pubchem_name_to_cid(name):
    url = PUBCHEM_NAME.format(name=urllib.parse.quote(name, safe=""))
    try:
        return str(json.loads(request(url))["IdentifierList"]["CID"][0])
    except (KeyError, IndexError, json.JSONDecodeError) as exc:
        raise PipelineError(f"PubChem could not resolve molecule name: {name}") from exc


def pubchem_input(cid, folder, stem, props, seed=61453):
    """Get PubChem's conformer or reproducibly build one with RDKit."""
    cid = str(int(cid))
    sdf = folder / f"{stem}.sdf"
    method = "PubChem deposited 3D conformer"
    try:
        data = request(PUBCHEM_SDF.format(cid=cid))
        if b"V2000" not in data and b"V3000" not in data:
            raise PipelineError("response was not an SDF structure")
        sdf.write_bytes(data)
    except PipelineError:
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError as exc:
            raise PipelineError("PubChem has no 3D conformer for this CID. Install fallback builder with: python3 -m pip install -r requirements-pubchem.txt") from exc
        smiles = props.get("SMILES") or props.get("ConnectivitySMILES") or props.get("CanonicalSMILES")
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        params = AllChem.ETKDGv3(); params.randomSeed = seed; params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise PipelineError(f"RDKit could not embed CID {cid}")
        from .conformers import sample_lowest_energy_conformer
        mol, search = sample_lowest_energy_conformer(
            mol, seed=seed, conformer_count=1, prune_rms=0.15, max_iters=5000)
        forcefield = search['forcefield']
        if not search['selected_converged']:
            raise PipelineError(f"{forcefield} geometry minimization did not converge for CID {cid}")
        Chem.SanitizeMol(mol)
        writer = Chem.SDWriter(str(sdf)); writer.write(mol); writer.close()
        method = f"RDKit ETKDGv3 seed={seed}, converged {forcefield} minimization with explicit hydrogens and aromatic planarity restraints"
    return sdf, props, method

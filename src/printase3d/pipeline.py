"""Pipeline components of the 3DPrintase pipeline."""
from __future__ import annotations

from pathlib import Path
import hashlib
from .chemistry import prepare_molecule, rdkit_read_file
from .config import FRUCTAN_TYPES, GLUCAN_TYPES, PUBCHEM_SDF, RCSB_CIF, VERSION, WATER
from .glycans import build_fructan, build_glucan
from .metadata import write_model_metadata
from .options import safe_name, variant_name
from .outputs import model_directory
from .rendering import render_jmol
from .sources import download_cif, pubchem_input, pubchem_name_to_cid, pubchem_properties
from .structures import clean, dna_residues, glycan_components, load_structure, write_component
from .utils import PipelineError, request, sha256, utc_now


def process_molecule(args):
    """Process peptide sequences, PubChem, line notations, or molecule files."""
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise PipelineError("The molecule command requires RDKit") from exc

    props = {}
    source_url = None
    ccd_provenance = None; ccd_raw = None
    peptide_provenance = None
    if getattr(args, "command", None) == "peptide":
        from .peptides import build_peptide
        mol, peptide_provenance = build_peptide(args.sequence)
        sequence = peptide_provenance["sequence"]
        label = args.label or sequence
        model_id = "peptide-" + sequence.lower()
        source_description = f"Linear L-peptide {sequence} (N-to-C), uncapped termini; built locally with RDKit"
    elif getattr(args, "ccd", None):
        from .ccd import fetch_ccd
        mol, ccd_provenance, ccd_raw = fetch_ccd(args.ccd, args.ccd_coordinates)
        label = args.label or args.ccd.upper()
        model_id = "ccd-" + args.ccd.upper()
        source_url = ccd_provenance['url']
        source_description = f"RCSB CCD {args.ccd.upper()} {args.ccd_coordinates} coordinates"
        if args.geometry_mode == 'pubchem':
            raise PipelineError('--ccd cannot be labeled as PubChem geometry')
        if args.geometry_mode == 'auto' and not any(getattr(args,k,False) for k in ('extended','most_flat','least_flat')):
            args.geometry_mode = 'rcsb'
    elif args.cid or args.pubchem_name:
        cid = str(int(args.cid)) if args.cid else pubchem_name_to_cid(args.pubchem_name)
        props = pubchem_properties(cid)
        label = args.label or props.get("Title") or args.pubchem_name or f"compound-{cid}"
        model_id = cid
        source_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
        try:
            raw = request(PUBCHEM_SDF.format(cid=cid))
            supplier = Chem.SDMolSupplier()
            supplier.SetData(raw.decode("utf-8", errors="replace"), removeHs=False)
            mol = next((item for item in supplier if item is not None), None)
        except PipelineError:
            mol = None
        if mol is None:
            smiles = props.get("SMILES") or props.get("IsomericSMILES") or props.get("ConnectivitySMILES") or props.get("CanonicalSMILES")
            mol = Chem.MolFromSmiles(smiles or "")
            if mol is None:
                raise PipelineError(f"PubChem supplied neither usable 3D coordinates nor SMILES for CID {cid}")
        source_description = f"PubChem CID {cid}"
    elif args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
        if mol is None:
            raise PipelineError("RDKit could not parse --smiles")
        label = args.label or "custom-molecule"
        model_id = "smiles-" + hashlib.sha256(args.smiles.encode()).hexdigest()[:10]
        source_description = "CLI SMILES"
    elif args.inchi:
        mol = Chem.MolFromInchi(args.inchi)
        if mol is None:
            raise PipelineError("RDKit could not parse --inchi")
        label = args.label or "custom-molecule"
        model_id = "inchi-" + hashlib.sha256(args.inchi.encode()).hexdigest()[:10]
        source_description = "CLI InChI"
    else:
        mol, input_path = rdkit_read_file(args.input)
        label = args.label or (mol.GetProp("_Name") if mol.HasProp("_Name") and mol.GetProp("_Name").strip() else input_path.stem)
        model_id = "file-" + sha256(input_path)[:10]
        source_description = str(input_path)

    if args.geometry_mode == "rcsb" and ccd_provenance is None:
        raise PipelineError("--rcsb-geometry requires --ccd")
    label_slug = safe_name(label)
    stem = f"{label_slug}_{model_id}"
    variant = variant_name(args, args.seed)
    folder = model_directory(args.output, stem, variant,
                             unique=bool(getattr(args, "ccd", None)) or getattr(args, "most_flat", False) or getattr(args, "least_flat", False) or getattr(args, "extended", False),
                             dry_run=args.dry_run)
    variant = folder.name
    if args.dry_run:
        print(f"WOULD PREPARE {label} ({model_id}) -> {folder}")
        return folder
    folder.mkdir(parents=True, exist_ok=True)
    if ccd_raw is not None:
        (folder / "source_ccd.cif").write_bytes(ccd_raw)
    mol, chemistry = prepare_molecule(mol, args, label)
    sdf = folder / f"{stem}.sdf"
    writer = Chem.SDWriter(str(sdf))
    writer.write(mol)
    writer.close()
    obj = folder / f"{stem}.obj"
    log = folder / f"{stem}_jmol.log"
    command = mesh_report = None
    if not args.no_render:
        command, mesh_report = render_jmol(args, sdf, obj, log)
    metadata_name = f"{stem}_metadata.json"
    meta = {
        "pipeline": "3DPrintase", "pipeline_version": VERSION, "created_utc": utc_now(),
        "model_id": model_id, "molecule_name": label, "pubchem_cid": int(model_id) if model_id.isdigit() else None,
        "molecule_root": str(folder.parent), "variant": variant, "source_url": source_url,
        "source_description": source_description, "input_sha256": sha256(sdf), "properties": props,
        **chemistry, "rcsb_ccd": {**ccd_provenance, "source_sha256": hashlib.sha256(ccd_raw).hexdigest()} if ccd_provenance else None,
        **({"peptide": peptide_provenance} if peptide_provenance else {}),
        "files": {"ccd_source": "source_ccd.cif" if ccd_raw is not None else None, "input_sdf": sdf.name, "obj": obj.name if obj.exists() else None,
                  "mtl": obj.with_suffix(".mtl").name if obj.with_suffix(".mtl").exists() else None,
                  "jmol_raw_obj": mesh_report.get("jmol_raw_obj") if mesh_report else None,
                  "jmol_raw_mtl": mesh_report.get("jmol_raw_mtl") if mesh_report else None,
                  "jmol_log": log.name if log.exists() else None, "metadata": metadata_name,
                  "model_card": "MODEL.md"},
        "rendering": {"preset": args.render_preset, "hydrogens_visible": args.hydrogens, "monocolor_aliphatics": getattr(args, "monocolor_aliphatics", False),
                      "bonds": "one centered cylinder per covalent atom pair by default; bond order retained for geometry",
                      "bond_color_mode": args.bond_colors, "multiple_bond_cylinders": args.multiple_bond_cylinders,
                      "spacefill_percent_vdw": args.spacefill,
                      "wireframe_radius_angstrom": args.wireframe},
        "print_mesh": mesh_report,
        "jmol": {"jar": str(Path(args.jmol_jar).resolve()), "script": str(Path(args.script).resolve()),
                 "script_sha256": sha256(args.script), "command": command, "style": args.style, "scale": args.scale},
    }
    write_model_metadata(folder, metadata_name, meta)
    print(obj if obj.exists() else folder)
    return folder


def process_fructan(kind, dp, args):
    if dp < 2:
        raise PipelineError("Fructan DP must be at least 2")
    definition = FRUCTAN_TYPES[kind]
    stem = f"{definition['name']}_dp{dp}"
    molecule_root = Path(args.output).resolve() / stem
    variant = variant_name(args, args.seed)
    folder = molecule_root / variant
    if args.dry_run:
        print(f"WOULD BUILD {kind} DP{dp} ({definition['linkage_label']}) -> {folder}")
        return folder
    folder.mkdir(parents=True, exist_ok=True)
    sdf = folder / f"{stem}.sdf"
    build = build_fructan(kind, dp, sdf, args.seed, args.conformation)
    obj = folder / f"{stem}.obj"
    log = folder / f"{stem}_jmol.log"
    metadata = folder / f"{stem}_metadata.json"
    command = mesh_report = None
    if not args.no_render:
        command, mesh_report = render_jmol(args, sdf, obj, log)
    meta = {
        "pipeline_version": VERSION,
        "created_utc": utc_now(),
        "model_id": f"{kind}-dp{dp}",
        "molecule_name": f"{definition['name']} DP{dp}",
        "fructan_type": kind,
        "degree_of_polymerization": dp,
        "residue_definition": "linear fructose-only beta-D-fructofuranose oligomer with a terminal reducing Fruf residue",
        "linkage": definition["linkage_label"],
        "molecule_root": str(molecule_root),
        "variant": variant,
        **build,
        "files": {
            "input_sdf": sdf.name,
            "obj": obj.name if obj.exists() else None,
            "mtl": obj.with_suffix(".mtl").name if obj.with_suffix(".mtl").exists() else None,
            "jmol_raw_obj": mesh_report.get("jmol_raw_obj") if mesh_report else None,
            "jmol_log": log.name if log.exists() else None,
            "metadata": metadata.name,
        },
        "rendering": {
            "preset": args.render_preset,
            "hydrogens_visible": args.hydrogens, "monocolor_aliphatics": getattr(args, "monocolor_aliphatics", False),
            "bonds": "covalent only; hydrogen bonds and struts disabled",
            "bond_color_mode": args.bond_colors,
            "multiple_bond_cylinders": args.multiple_bond_cylinders,
            "spacefill_percent_vdw": args.spacefill,
            "wireframe_radius_angstrom": args.wireframe,
        },
        "print_mesh": mesh_report,
        "jmol": {
            "script": str(Path(args.script).resolve()),
            "script_sha256": sha256(args.script),
            "command": command,
            "style": args.style,
            "scale": args.scale,
        },
    }
    write_model_metadata(folder, metadata.name, meta)
    return folder


def process_glucan(kind, dp, args):
    if dp < 2:
        raise PipelineError("Glucan DP must be at least 2")
    definition = GLUCAN_TYPES[kind]
    stem = f"{definition['name']}_dp{dp}"
    molecule_root = Path(args.output).resolve() / stem
    variant = variant_name(args, args.seed)
    folder = molecule_root / variant
    if args.dry_run:
        print(f"WOULD BUILD {kind} DP{dp} ({definition['linkage_label']}) -> {folder}")
        return folder
    folder.mkdir(parents=True, exist_ok=True)
    sdf = folder / f"{stem}.sdf"
    build = build_glucan(kind, dp, sdf, args.seed, args.conformation)
    obj = folder / f"{stem}.obj"
    log = folder / f"{stem}_jmol.log"
    metadata = folder / f"{stem}_metadata.json"
    command = mesh_report = None
    if not args.no_render:
        command, mesh_report = render_jmol(args, sdf, obj, log)
    meta = {
        "pipeline_version": VERSION,
        "created_utc": utc_now(),
        "model_id": f"{kind}-dp{dp}",
        "molecule_name": f"{definition['series_name']} DP{dp}",
        "glucan_type": kind,
        "degree_of_polymerization": dp,
        "residue_definition": (
            f"linear D-{definition['residue']}p oligomer with a terminal reducing "
            f"{definition['residue']}p residue"
        ),
        "linkage": definition["linkage_label"],
        "molecule_root": str(molecule_root),
        "variant": variant,
        **build,
        "files": {
            "input_sdf": sdf.name,
            "obj": obj.name if obj.exists() else None,
            "mtl": obj.with_suffix(".mtl").name if obj.with_suffix(".mtl").exists() else None,
            "jmol_raw_obj": mesh_report.get("jmol_raw_obj") if mesh_report else None,
            "jmol_log": log.name if log.exists() else None,
            "metadata": metadata.name,
        },
        "rendering": {
            "preset": args.render_preset,
            "hydrogens_visible": args.hydrogens, "monocolor_aliphatics": getattr(args, "monocolor_aliphatics", False),
            "bonds": "covalent only; hydrogen bonds and struts disabled",
            "bond_color_mode": args.bond_colors,
            "multiple_bond_cylinders": args.multiple_bond_cylinders,
            "spacefill_percent_vdw": args.spacefill,
            "wireframe_radius_angstrom": args.wireframe,
        },
        "print_mesh": mesh_report,
        "jmol": {
            "script": str(Path(args.script).resolve()),
            "script_sha256": sha256(args.script),
            "command": command,
            "style": args.style,
            "scale": args.scale,
        },
    }
    write_model_metadata(folder, metadata.name, meta)
    return folder


def process_cid(cid, args):
    cid = str(int(cid))
    variant = variant_name(args, args.seed)
    if args.dry_run:
        folder = Path(args.output).resolve() / f"<compound-name>_{cid}" / variant
        print(f"WOULD RESOLVE PubChem CID {cid} name and render under {folder}"); return folder
    props = pubchem_properties(cid)
    stem = f"{safe_name(props.get('Title', 'compound'))}_{cid}"
    molecule_root = Path(args.output).resolve() / stem
    folder = molecule_root / variant
    obj = folder / f"{stem}.obj"
    folder.mkdir(parents=True, exist_ok=True)
    sdf, props, method = pubchem_input(cid, folder, stem, props, args.seed)
    command = mesh_report = None
    log = folder / f"{stem}_jmol.log"
    metadata = folder / f"{stem}_metadata.json"
    if not args.no_render: command, mesh_report = render_jmol(args, sdf, obj, log)
    meta = {
        "pipeline_version": VERSION, "created_utc": utc_now(), "pubchem_cid": int(cid),
        "model_id": cid, "molecule_name": props.get("Title", stem),
        "molecule_root": str(molecule_root), "variant": variant,
        "source_url": PUBCHEM_SDF.format(cid=cid), "coordinate_method": method,
        "input_sha256": sha256(sdf), "properties": props,
        "files": {"input_sdf": sdf.name, "obj": obj.name if obj.exists() else None, "mtl": obj.with_suffix('.mtl').name if obj.with_suffix('.mtl').exists() else None, "jmol_raw_obj": mesh_report.get("jmol_raw_obj") if mesh_report else None, "jmol_log": log.name if log.exists() else None, "metadata": metadata.name},
        "rendering": {"preset": args.render_preset, "hydrogens_visible": args.hydrogens, "monocolor_aliphatics": getattr(args, "monocolor_aliphatics", False), "bonds": "one centered cylinder per covalent atom pair by default; hydrogen bonds and struts disabled", "bond_color_mode": args.bond_colors, "multiple_bond_cylinders": args.multiple_bond_cylinders, "spacefill_percent_vdw": args.spacefill, "wireframe_radius_angstrom": args.wireframe},
        "print_mesh": mesh_report,
        "jmol": {"script": str(Path(args.script).resolve()), "script_sha256": sha256(args.script), "command": command, "style": args.style, "scale": args.scale},
    }
    write_model_metadata(folder, metadata.name, meta)
    return folder


def process_pdb(pdb, args):
    pdb = pdb.lower()
    root = Path(args.output).resolve() / pdb
    cif = root / f"{pdb}.cif"
    if args.dry_run:
        download_cif(pdb, cif, True); return []
    root.mkdir(parents=True, exist_ok=True)
    if not cif.exists() or args.force: download_cif(pdb, cif)
    singles, loops, comp_types, atoms = load_structure(cif)
    if not any(a["group"] == "ATOM" and a["comp"] not in WATER for a in atoms):
        raise PipelineError(f"{pdb}: no protein/polymer ATOM records; not a protein–oligosaccharide complex")
    if args.extract == "dna":
        component, residue_atoms = dna_residues(atoms)
        components = [component] if component else []
    else:
        components, residue_atoms = glycan_components(loops, comp_types, atoms)
    components = [c for c in components if len(c) >= args.min_dp]
    if args.components:
        wanted = {x.upper() for x in args.components.split(",")}
        components = [c for c in components if any(r.comp in wanted for r in c)]
    if not components:
        target = "DNA" if args.extract == "dna" else "carbohydrate"
        raise PipelineError(f"{pdb}: no {target} component with DP >= {args.min_dp}")
    outputs = []
    variant = variant_name(args)
    for number, component in enumerate(components, 1):
        if args.extract == "dna":
            stem = f"p53-response-element_{pdb}"
        else:
            stem = f"glycan_{number:02d}_dp{len(component)}_{component[0].chain}"
        molecule_root = root / stem
        folder = molecule_root / variant; folder.mkdir(parents=True, exist_ok=True)
        extracted = folder / f"{stem}.pdb"
        selected = write_component(extracted, component, residue_atoms)
        obj = folder / f"{stem}.obj"
        command = mesh_report = None
        if not args.no_render: command, mesh_report = render_jmol(args, extracted, obj, folder / "jmol.log")
        meta = {
            "pipeline_version": VERSION, "created_utc": utc_now(), "pdb_id": pdb.upper(),
            "molecule_root": str(molecule_root), "variant": variant,
            "source_url": RCSB_CIF.format(pdb=pdb.upper()), "source_sha256": sha256(cif),
            "experimental_method": clean(singles.get("_exptl.method")),
            "title": clean(singles.get("_struct.title")), "coordinate_policy": "model 1; altloc blank/A; coordinates unchanged",
            "filters": {"protein": "removed", "waters": "removed", "ions": "removed", "non_target_ligands": "removed"},
            "extract_type": args.extract,
            "component": {"dp": len(component), "atom_count": len(selected), "residues": [r.__dict__ for r in component]},
            "files": {"extracted_pdb": extracted.name, "obj": obj.name if obj.exists() else None, "mtl": obj.with_suffix('.mtl').name if obj.with_suffix('.mtl').exists() else None, "jmol_raw_obj": mesh_report.get("jmol_raw_obj") if mesh_report else None},
            "rendering": {"preset": args.render_preset, "hydrogens_visible": args.hydrogens, "monocolor_aliphatics": getattr(args, "monocolor_aliphatics", False), "bonds": "one centered cylinder per covalent atom pair by default; hydrogen bonds and struts disabled", "bond_color_mode": args.bond_colors, "multiple_bond_cylinders": args.multiple_bond_cylinders, "spacefill_percent_vdw": args.spacefill, "wireframe_radius_angstrom": args.wireframe},
            "print_mesh": mesh_report,
            "jmol": {"script": str(Path(args.script).resolve()), "script_sha256": sha256(args.script), "command": command, "style": args.style, "scale": args.scale},
        }
        write_model_metadata(folder, "metadata.json", meta)
        outputs.append(folder)
    return outputs

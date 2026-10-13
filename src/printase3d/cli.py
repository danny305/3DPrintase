"""Cli components of the 3DPrintase pipeline."""
from __future__ import annotations

from pathlib import Path
import json
import os
import shutil
import sys
from .arguments import parser, add_common, add_molecule_options, positive_int
from .config import DEFAULT_OUTPUT, FRUCTAN_TYPES, GLUCAN_TYPES, PRESETS, VERSION, RCSB_SEARCH, RCSB_CIF, PUBCHEM_SDF, PUBCHEM_JSON, PUBCHEM_NAME, RENDER_PRESETS, GLYCAN_CONFORMATION_PROFILES, WATER, IONS
from .options import ph_value, resolve_render_preset, spacefill_value, wireframe_value, safe_name, setting_token, variant_name
from .pipeline import process_cid, process_fructan, process_glucan, process_molecule, process_pdb
from .rendering import resolve_jmol_jar, render_jmol
from .sources import do_search, search_query, download_cif, pubchem_properties, pubchem_name_to_cid, pubchem_input
from .utils import PipelineError, utc_now, request, sha256
from .chemistry import rdkit_read_file, correct_false_amide_protonation, transfer_heavy_atom_coordinates, sample_lowest_energy_conformer, prepare_molecule
from .glycans import _sugar_residue_topologies, infer_glycosidic_torsions, apply_glycan_conformation, glycan_backbone_metrics, embed_glycan_conformer, build_fructan, build_glucan
from .mesh import parse_jmol_obj, compact_group_mesh, CARBON_COLOR, BOND_COLOR, ELEMENT_MATERIALS, semantic_atom_material, write_print_mtl, repair_jmol_obj
from .metadata import reproduction_command, write_model_metadata
from .structures import cif_tokens, parse_cif, loop_rows, clean, Residue, load_structure, is_sugar, distance, glycan_components, pdb_atom_line, write_component, DNA_COMPONENTS, dna_residues


def doctor(args):
    checks = []
    java = shutil.which("java")
    checks.append(("Java", bool(java), java or "not found on PATH"))
    jar = Path(args.jmol_jar).expanduser().resolve()
    checks.append(("Jmol.jar", jar.is_file(), str(jar)))
    for module in ("numpy", "trimesh", "rdkit", "glyles", "dimorphite_dl"):
        try:
            imported = __import__(module)
            checks.append((module, True, getattr(imported, "__version__", "installed")))
        except ImportError:
            checks.append((module, False, "not installed"))
    output = Path(args.output).expanduser().resolve()
    parent = next((p for p in [output, *output.parents] if p.exists()), None)
    checks.append(("output", bool(parent and os.access(parent, os.W_OK)), str(output)))
    for name, ok, detail in checks:
        print(f"{'OK' if ok else 'FAIL':4}  {name:12} {detail}")
    return 0 if all(ok for _, ok, _ in checks) else 2


def task_options(args):
    """Shared parser settings, suitable for the Python API option serializer."""
    skip = {"command", "manifest", "jobs", "label", "cids", "pdb_ids", "dp", "type",
            "limit", "render_preset", "geometry_mode"}
    options = {key: value for key, value in vars(args).items() if key not in skip}
    options[args.render_preset.replace("-", "_")] = True
    geometry = getattr(args, "geometry_mode", "auto")
    if geometry != "auto":
        options[{"ensemble":"ensemble", "pubchem":"pubchem_geometry", "rcsb":"rcsb_geometry"}[geometry]] = True
    return options


def main(argv=None):
    from .api import generate_batch, _options_argv
    from .batch import Task, run_tasks
    from .metadata import REPRODUCTION_COMMAND
    import shlex

    actual_argv = list(sys.argv[1:] if argv is None else argv)
    args = resolve_render_preset(parser().parse_args(actual_argv))
    token = REPRODUCTION_COMMAND.set(shlex.join(["3dprintase", *actual_argv]))
    try:
        if args.command == "doctor":
            return doctor(args)
        if args.command == "search":
            do_search(args)
            return 0
        if args.command in ("molecule", "peptide"):
            process_molecule(args)
            return 0
        if args.command == "batch":
            manifest = Path(args.manifest).resolve()
            entries = json.loads(manifest.read_text())
            if not isinstance(entries, list) or not all(isinstance(x, dict) for x in entries):
                raise ValueError("Batch manifest must be a JSON list of molecule objects")
            for entry in entries:
                if entry.get("input"):
                    entry["input"] = str((manifest.parent / entry["input"]).resolve())
            result = generate_batch(entries, jobs=args.jobs, **task_options(args))
        else:
            common = _options_argv(task_options(args))
            if args.command == "pubchem":
                tasks = [Task(str(int(cid)), tuple(["pubchem", str(int(cid)), *common])) for cid in args.cids]
            elif args.command in ("glucan", "fructan"):
                tasks = [Task(f"{args.type}-dp{dp}", tuple([args.command, "--type", args.type,
                         "--dp", str(dp), *common])) for dp in args.dp]
            else:
                ids = args.pdb_ids if args.command == "run" else do_search(args)
                tasks = [Task(pdb.lower(), tuple(["run", pdb.lower(), *common])) for pdb in ids]
            result = run_tasks(tasks, jobs=args.jobs, output=args.output, dry_run=args.dry_run)
        for failure in result.failures:
            print(f"WARNING {failure.label}: {failure.error}", file=sys.stderr)
        print(f"Completed: {len(result.outputs)} output(s); {len(result.failures)} failure(s); "
              f"{result.workers} worker(s)", file=sys.stderr)
        return 1 if result.failures else 0
    except (PipelineError, ValueError, TypeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    finally:
        REPRODUCTION_COMMAND.reset(token)


if __name__ == "__main__":
    raise SystemExit(main())

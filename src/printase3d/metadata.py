"""Metadata components of the 3DPrintase pipeline."""
from __future__ import annotations

from pathlib import Path
from contextvars import ContextVar
import json
import shlex
import sys

REPRODUCTION_COMMAND = ContextVar("reproduction_command", default=None)


def reproduction_command():
    return REPRODUCTION_COMMAND.get() or shlex.join(["3dprintase", *sys.argv[1:]])


def write_model_metadata(folder, metadata_name, meta):
    """Write machine-readable metadata and a friendly per-variant model card."""
    folder = Path(folder)
    meta["reproduce"] = reproduction_command()
    (folder / metadata_name).write_text(json.dumps(meta, indent=2) + "\n")
    props = meta.get("properties", {})
    rendering = meta.get("rendering", {})
    conformation = meta.get("conformation_details") or {}
    backbone = meta.get("backbone_metrics") or conformation.get("candidate_selection", {}).get("backbone_metrics", {})
    conformer_search = meta.get("conformer_search") or {}
    mesh = meta.get("print_mesh") or {}
    files = meta.get("files", {})
    lines = [
        f"# {meta.get('molecule_name') or props.get('Title') or meta.get('model_id') or '3DPrintase model'}",
        "",
        "This directory is one printable rendering variant. The molecule-level parent directory groups all variants.",
        "",
        "## Identity",
        "",
        f"- ID: `{meta.get('model_id') or meta.get('pubchem_cid') or meta.get('pdb_id', 'custom')}`",
        f"- Formula: `{meta.get('molecular_formula') or props.get('MolecularFormula', 'not reported')}`",
        f"- Source: {meta.get('source_description') or meta.get('source_url', 'generated locally')}",
        f"- Coordinates: {meta.get('coordinate_method') or meta.get('coordinate_policy', 'see JSON metadata')}",
        f"- Protonation: {meta.get('protonation_method', 'source/default state')}",
        f"- Carbonyl/amide geometry: `{len(meta.get('carbonyl_geometry', []))}` local trigonal sites measured; see JSON for hydrogen-inclusive geometry",
        f"- Formal charge: `{meta.get('formal_charge', 'not reported')}`",
        f"- Geometry policy: `{meta.get('geometry_mode_requested', 'source/default')}`",
        f"- Sampled conformers: `{conformer_search.get('total_candidates_including_source', 'not sampled')}`",
        f"- Energy interpretation: {conformer_search.get('ranking_scope', 'source geometry / no ensemble')}",
        f"- Aromatic ring plane RMS: `{(meta.get('aromatic_planarity') or conformer_search.get('aromatic_planarity') or {}).get('rms_angstrom', 'not restrained')}` Å",
        f"- Selected force-field energy: `{conformer_search.get('selected_energy_kcal_mol', 'n/a')} kcal/mol`",
        f"- Conformer selection: `{conformer_search.get('selection', 'energy/source')}`",
        f"- Extended terminal separation: `{(conformer_search.get('extension') or {}).get('selected_terminal_separation_angstrom', 'not requested')}` Å",
        f"- Heavy-atom RMS distance to best-fit plane: `{conformer_search.get('selected_flatness_rms_angstrom', 'n/a')} Å`",
        f"- Energy above sampled minimum: `{conformer_search.get('selected_above_minimum_kcal_mol', 'n/a')} kcal/mol`",
        f"- Energy below supplied-coordinate start: `{conformer_search.get('selected_below_source_kcal_mol', 'n/a')} kcal/mol`",
        "",
        "## Conformation",
        "",
        f"- Mode: `{meta.get('conformation', 'source/preserved')}`",
        f"- Profile: `{conformation.get('profile', 'not applied')}`",
        f"- Intended shape: {conformation.get('shape', 'source coordinates or general molecular conformer')}",
        f"- Torsion targets (degrees): `{conformation.get('targets_degrees', 'not applied')}`",
        f"- Ring-center end-to-end distance: `{backbone.get('end_to_end_angstrom', 'n/a')} Å`",
        f"- Profile source: {conformation.get('source', 'not applicable')}",
        "",
        "## Rendering",
        "",
        f"- Preset: `{rendering.get('preset')}`",
        f"- Atom size: `{rendering.get('spacefill_percent_vdw')}%` of van der Waals radii",
        f"- Bond radius: `{rendering.get('wireframe_radius_angstrom')} Å`",
        f"- Hydrogens visible: `{rendering.get('hydrogens_visible')}`",
        f"- Hydrogen/carbon same color: `{rendering.get('monocolor_aliphatics', False)}`",
        f"- Bond colors: `{rendering.get('bond_color_mode')}`",
        f"- Multiple-bond cylinders: `{rendering.get('multiple_bond_cylinders', False)}`",
        f"- Mesh: `{mesh.get('topology', 'not rendered')}`; open edges: `{mesh.get('open_edges', 'n/a')}`",
        "",
        "## Files",
        "",
    ]
    lines.extend(f"- {key}: `{value}`" for key, value in files.items() if value)
    lines += ["", "## Reproduce", "", "```bash", meta["reproduce"], "```", ""]
    (folder / "MODEL.md").write_text("\n".join(lines))

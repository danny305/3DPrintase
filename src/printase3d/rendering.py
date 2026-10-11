"""Rendering components of the 3DPrintase pipeline."""
from __future__ import annotations

from pathlib import Path
import os
import shlex
import subprocess
from .mesh import repair_jmol_obj
from .utils import PipelineError


def render_jmol(args, input_pdb, obj_path, log_path):
    template = Path(args.script).resolve()
    jar = Path(args.jmol_jar).expanduser().resolve()
    if not template.is_file(): raise PipelineError(f"Pipeline Jmol script not found: {template}")
    if not jar.is_file():
        raise PipelineError(
            f"Jmol.jar not found: {jar}. Download the full Jmol application from "
            "https://jmol.sourceforge.net/download/ and set THREEDPRINTASE_JMOL_JAR "
            "or pass --jmol-jar /path/to/Jmol.jar. Use --no-render for SDF-only output."
        )
    runtime = obj_path.with_suffix(".jmol.spt")
    script = template.read_text()
    replacements = {
        "__INPUT_PDB__": str(input_pdb.resolve()).replace("\\", "/"),
        "__OUTPUT_OBJ__": str(obj_path.resolve()).replace("\\", "/"),
        "__STYLE__": args.style,
        "__SCALE__": str(args.scale),
        "__HYDROGENS__": "on" if args.hydrogens else "off",
        "__COLOR__": "yes" if args.color else "no",
        "__SPACEFILL__": str(args.spacefill),
        "__WIREFRAME__": str(args.wireframe),
        "__BOND_COLORS__": args.bond_colors,
        "__MULTIPLE_BONDS__": "true" if args.multiple_bond_cylinders else "false",
    }
    for old, new in replacements.items(): script = script.replace(old, new)
    runtime.write_text(script)
    # Jmol 16.3's separate '-n -s' form can exit before an external script
    # starts; its documented combined flags execute synchronously.
    cmd = jmol_command(jar, runtime, threads=getattr(args, "threads", 1))
    if args.dry_run:
        print("WOULD RUN " + shlex.join(cmd)); return None
    try:
        proc = subprocess.run(cmd, text=True, stdin=subprocess.DEVNULL,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              timeout=args.timeout)
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or b""
        log_path.write_text(output.decode(errors="replace") if isinstance(output, bytes) else output)
        raise PipelineError(f"Jmol timed out after {args.timeout}s; see {log_path}") from exc
    log_path.write_text(proc.stdout)
    if proc.returncode != 0 or not obj_path.exists() or obj_path.stat().st_size == 0:
        raise PipelineError(f"Jmol export failed (exit {proc.returncode}); see {log_path}")
    monocolor_aliphatics = getattr(args, "monocolor_aliphatics", False)
    mesh_report = repair_jmol_obj(obj_path, color=args.color, bond_colors=args.bond_colors,
                                  monocolor_aliphatics=monocolor_aliphatics)
    if args.color and Path(input_pdb).suffix.lower() in {".sdf", ".mol"} and args.hydrogens:
        from rdkit import Chem
        from .material_validation import validate_element_materials, validate_exported_atom_coordinates
        mol = Chem.SDMolSupplier(str(input_pdb), removeHs=False, sanitize=False)[0]
        mesh_report["element_colors"] = validate_element_materials(obj_path, mol, monocolor_aliphatics=monocolor_aliphatics)
        if mol.GetNumAtoms() >= 3:
            mesh_report["coordinate_validation"] = validate_exported_atom_coordinates(obj_path, mol)
    return cmd, mesh_report


def jmol_command(jar, runtime, *, threads=1):
    """Keep Jmol's exporter in no-display mode and suppress macOS activation.

    The OBJ exporter needs AWT: true headless mode either throws from Jmol.jar
    or silently omits OBJ output in JmolData.jar. UIElement makes Java a macOS
    background agent; -n (in -ions) prevents Jmol from creating its window.
    """
    return ["java", "-Dapple.awt.UIElement=true",
            f"-XX:ActiveProcessorCount={max(1, threads)}", "-Xmx2g", "-jar",
            str(jar), "-ions", str(runtime)]


def resolve_jmol_jar():
    """Find Jmol without making the installed command depend on its cwd."""
    configured = os.environ.get("THREEDPRINTASE_JMOL_JAR")
    if configured:
        return str(Path(configured).expanduser())
    here = Path(__file__).resolve()
    # An editable checkout can keep the separately downloaded renderer locally.
    project = here.parents[2]
    if (project / "pyproject.toml").is_file():
        for candidate in (project / "tools" / "jmol" / "Jmol.jar", project / "Jmol.jar"):
            if candidate.is_file():
                return str(candidate)
    # Installed wheels use a per-user cache, independent of the working directory.
    cache = Path(os.environ.get("XDG_CACHE_HOME", str(Path.home() / ".cache")))
    return str(cache / "3dprintase" / "Jmol.jar")

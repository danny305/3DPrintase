"""Mesh components of the 3DPrintase pipeline."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import re
import shutil
from .utils import PipelineError


def parse_jmol_obj(path):
    """Read Jmol's group-oriented OBJ without merging intersecting parts."""
    vertices, groups = [], []
    current = None
    mtllib = None
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("mtllib "):
            mtllib = line.split(None, 1)[1].strip()
        elif line.startswith("g "):
            current = {"name": line.split(None, 1)[1].strip(), "material": None, "faces": []}
            groups.append(current)
        elif line.startswith("usemtl ") and current is not None:
            current["material"] = line.split(None, 1)[1].strip()
        elif line.startswith("v "):
            try:
                vertices.append(tuple(float(x) for x in line.split()[1:4]))
            except (ValueError, IndexError) as exc:
                raise PipelineError(f"Invalid OBJ vertex in {path}: {line}") from exc
        elif line.startswith("f ") and current is not None:
            face = []
            for token in line.split()[1:]:
                raw = token.split("/", 1)[0]
                index = int(raw)
                if index < 0:
                    index = len(vertices) + index
                else:
                    index -= 1
                face.append(index)
            for i in range(1, len(face) - 1):
                current["faces"].append((face[0], face[i], face[i + 1]))
    if not vertices or not groups:
        raise PipelineError(f"Jmol OBJ contains no usable mesh groups: {path}")
    return vertices, groups, mtllib


def compact_group_mesh(all_vertices, faces):
    """Return a group's vertices/faces with zero-based local indices."""
    used = sorted({index for face in faces for index in face})
    remap = {old: new for new, old in enumerate(used)}
    return [all_vertices[index] for index in used], [tuple(remap[index] for index in face) for face in faces]


CARBON_COLOR = (0.180, 0.180, 0.180)


BOND_COLOR = CARBON_COLOR

# Bambu's MTL parser recognizes this exact comment (no space after '#').
# ObjColorPanel::deal_algo then uses the supplied material palette on initial
# import instead of K-means, which can merge chemically distinct O/P colors.
# This is an importer compatibility flag, not a statement of file provenance.
BAMBU_PRESERVE_PALETTE = "#FirstTimeUsingMakerLab"


ELEMENT_MATERIALS = {
    "H": ("hydrogen", (1.000, 1.000, 1.000)),
    # Keep carbon deliberately charcoal rather than mid-gray. Bambu Studio's
    # OBJ color reduction otherwise clusters carbon with white hydrogen and
    # the neutral bond material, even though the MTL names are different.
    "C": ("carbon", CARBON_COLOR),
    "N": ("nitrogen", (0.191, 0.316, 0.973)),
    "O": ("oxygen", (1.000, 0.055, 0.055)),
    "F": ("fluorine", (0.565, 0.878, 0.314)),
    "P": ("phosphorus", (1.000, 0.504, 0.000)),
    "S": ("sulfur", (1.000, 1.000, 0.188)),
    "Fe": ("iron", (0.878, 0.400, 0.200)),
    "Mg": ("magnesium", (0.541, 1.000, 0.000)),
    "Co": ("cobalt", (0.941, 0.565, 0.627)),
    "Se": ("selenium", (0.900, 0.420, 0.100)),
    "Cl": ("chlorine", (0.122, 0.941, 0.122)),
    "Br": ("bromine", (0.651, 0.161, 0.161)),
    "I": ("iodine", (0.580, 0.000, 0.580)),
}


def semantic_atom_material(group_name, *, monocolor_aliphatics=False):
    # SDF exports commonly use C1/O12-style names, while PDB nucleotides use
    # chemically meaningful atom names such as P, OP1, O5', and C1'.  In both
    # conventions the leading element symbol is authoritative for the elements
    # currently supported by the print palette.
    match = re.match(r"\d*(Cl|Br|Se|Fe|Mg|Co|CL|BR|SE|FE|MG|CO|[A-Z])", group_name)
    element = match.group(1) if match else "X"
    element = {"CL": "Cl", "BR": "Br", "SE": "Se", "FE": "Fe", "MG": "Mg", "CO": "Co"}.get(element, element)
    label, rgb = ELEMENT_MATERIALS.get(element, ("element", (0.700, 0.700, 0.700)))
    if monocolor_aliphatics and element == "H":
        rgb = CARBON_COLOR
    return f"atom_{element}_{label}", rgb


def write_print_mtl(path, materials):
    lines = ["# 3DPrintase semantic element materials for reliable slicer import",
             BAMBU_PRESERVE_PALETTE, ""]
    for name, rgb in materials.items():
        lines += [
            f"newmtl {name}", "Ns 163", "Tr 0", "Ni 1.0", "illum 2",
            "d 1.000", f"Ka {rgb[0]:.3f} {rgb[1]:.3f} {rgb[2]:.3f}",
            f"Kd {rgb[0]:.3f} {rgb[1]:.3f} {rgb[2]:.3f}",
            "Ks 0.000 0.000 0.000", "",
        ]
    path.write_text("\n".join(lines))


def repair_jmol_obj(obj_path, *, color=True, bond_colors="uniform", monocolor_aliphatics=False):
    """Replace Jmol's open three-piece bonds with closed printable cylinders.

    Jmol writes every bond as an open cylinder plus two independent spheres.
    Those pieces look correct in a viewer but force slicers to repair thousands
    of open edges and overlapping shells.  Atom spheres are already closed, so
    retain them; discard the bond-cap spheres and reconstruct each CylinderN as
    one capped, watertight cylinder with the same axis, radius, and material.
    """
    try:
        import numpy as np
        import trimesh
    except ImportError as exc:
        raise PipelineError(
            "Printable OBJ repair requires numpy and trimesh. Recreate the "
            "3D-glycan Conda environment from environment.yml."
        ) from exc

    obj_path = Path(obj_path)
    raw_path = obj_path.with_name(f"{obj_path.stem}_jmol_raw.obj")
    shutil.copy2(obj_path, raw_path)
    all_vertices, groups, mtllib = parse_jmol_obj(raw_path)
    mtl_path = obj_path.parent / (mtllib or f"{obj_path.stem}.mtl")
    raw_mtl_path = obj_path.with_name(f"{obj_path.stem}_jmol_raw.mtl")
    if mtl_path.is_file():
        shutil.copy2(mtl_path, raw_mtl_path)
        raw_text = raw_path.read_text()
        if mtllib:
            raw_text = raw_text.replace(f"mtllib {mtllib}", f"mtllib {raw_mtl_path.name}", 1)
        raw_path.write_text(raw_text)
    meshes = []
    print_materials = {}
    atom_count = bond_count = 0

    for group in groups:
        name = group["name"]
        if name.endswith("_Atom"):
            vertices, faces = compact_group_mesh(all_vertices, group["faces"])
            mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
            if not mesh.is_watertight:
                raise PipelineError(f"Jmol atom mesh is unexpectedly open: {name}")
            if color:
                material, rgb = semantic_atom_material(name, monocolor_aliphatics=monocolor_aliphatics)
            else:
                material, rgb = "model_monochrome", (0.500, 0.500, 0.500)
            print_materials[material] = rgb
            meshes.append((name, material, mesh))
            atom_count += 1
        elif re.fullmatch(r"Cylinder\d+", name):
            vertices, _ = compact_group_mesh(all_vertices, group["faces"])
            points = np.asarray(vertices, dtype=float)
            center = points.mean(axis=0)
            covariance = np.cov((points - center).T)
            values, vectors = np.linalg.eigh(covariance)
            axis = vectors[:, int(np.argmax(values))]
            axial = (points - center) @ axis
            start = center + axis * axial.min()
            end = center + axis * axial.max()
            radial = points - center - np.outer(axial, axis)
            radius = float(np.median(np.linalg.norm(radial, axis=1)))
            if radius <= 0 or np.linalg.norm(end - start) <= 0:
                raise PipelineError(f"Could not recover printable cylinder geometry: {name}")
            mesh = trimesh.creation.cylinder(radius=radius, segment=[start, end], sections=36)
            if not mesh.is_watertight:
                raise PipelineError(f"Rebuilt bond cylinder is unexpectedly open: {name}")
            if color and bond_colors == "uniform":
                material, rgb = "bond_uniform", BOND_COLOR
            elif color:
                material, rgb = f"bond_{group['material'] or 'split'}", (0.504, 0.504, 0.504)
            else:
                material, rgb = "model_monochrome", (0.500, 0.500, 0.500)
            print_materials[material] = rgb
            meshes.append((name, material, mesh))
            bond_count += 1

    if atom_count == 0:
        raise PipelineError(f"No atom geometry found while repairing {raw_path}")

    lines = [
        "# Watertight printable OBJ rebuilt from Jmol geometry",
        "# Top-level OBJ objects are semantic element parts for slicer import",
    ]
    lines.append(f"mtllib {mtl_path.name}")
    vertex_offset = 1
    face_count = 0
    by_material = defaultdict(list)
    for name, material, mesh in meshes:
        by_material[material or "model_unassigned"].append((name, mesh))
    for material, material_meshes in by_material.items():
        lines.append(f"o part_{material}")
        lines.append(f"g part_{material}")
        lines.append(f"usemtl {material}")
        for name, mesh in material_meshes:
            lines.append(f"# component {name}")
            lines.extend(f"v {v[0]:.7f} {v[1]:.7f} {v[2]:.7f}" for v in mesh.vertices)
            for face in mesh.faces:
                a, b, c = (int(index) + vertex_offset for index in face)
                lines.append(f"f {a} {b} {c}")
                face_count += 1
            vertex_offset += len(mesh.vertices)
    obj_path.write_text("\n".join(lines) + "\n")
    write_print_mtl(mtl_path, print_materials)

    return {
        "topology": "closed atom spheres plus capped bond cylinders" if bond_count else "closed atom spheres",
        "atom_parts": atom_count,
        "bond_parts": bond_count,
        "parts": len(meshes),
        "open_edges": 0,
        "faces": face_count,
        "jmol_raw_obj": raw_path.name,
        "jmol_raw_mtl": raw_mtl_path.name if raw_mtl_path.exists() else None,
        "material_policy": "top-level semantic element parts plus a separate uniform bond part",
        "materials": sorted(print_materials),
    }

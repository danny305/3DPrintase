"""Structures components of the 3DPrintase pipeline."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import math
import re
import shlex
from .config import PRESETS


def cif_tokens(text):
    """Small mmCIF tokenizer supporting quotes and semicolon text fields."""
    lines = iter(text.splitlines())
    for line in lines:
        if line.startswith(";"):
            value = [line[1:]]
            for line in lines:
                if line.startswith(";"):
                    break
                value.append(line)
            yield "\n".join(value)
            continue
        lex = shlex.shlex(line, posix=True)
        lex.whitespace_split = True
        lex.commenters = "#"
        yield from lex


def parse_cif(text):
    tokens = list(cif_tokens(text))
    singles, loops = {}, []
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token == "loop_":
            i += 1
            names = []
            while i < len(tokens) and tokens[i].startswith("_"):
                names.append(tokens[i]); i += 1
            rows = []
            while i < len(tokens) and tokens[i] != "loop_" and not tokens[i].startswith("_") and not tokens[i].startswith("data_"):
                if i + len(names) > len(tokens):
                    break
                rows.append(dict(zip(names, tokens[i:i + len(names)])))
                i += len(names)
            loops.append((names, rows))
        elif token.startswith("_") and i + 1 < len(tokens):
            singles[token] = tokens[i + 1]; i += 2
        else:
            i += 1
    return singles, loops


def loop_rows(loops, prefix):
    for names, rows in loops:
        if names and names[0].startswith(prefix):
            return rows
    return []


def clean(v, default=""):
    return default if v in (None, ".", "?") else v


@dataclass(frozen=True)
class Residue:
    chain: str
    seq: str
    comp: str
    entity: str


def load_structure(cif_path):
    singles, loops = parse_cif(cif_path.read_text(errors="replace"))
    comp_types = {clean(r.get("_chem_comp.id")): clean(r.get("_chem_comp.type")).upper() for r in loop_rows(loops, "_chem_comp.")}
    atoms = []
    for r in loop_rows(loops, "_atom_site."):
        model = clean(r.get("_atom_site.pdbx_PDB_model_num"), "1")
        alt = clean(r.get("_atom_site.label_alt_id"))
        if model != "1" or alt not in ("", "A"):
            continue
        try:
            xyz = tuple(float(r[f"_atom_site.Cartn_{axis}"]) for axis in "xyz")
        except (KeyError, ValueError):
            continue
        comp = clean(r.get("_atom_site.label_comp_id"), clean(r.get("_atom_site.auth_comp_id")))
        atom = {
            "group": clean(r.get("_atom_site.group_PDB"), "HETATM"),
            "id": clean(r.get("_atom_site.id")), "element": clean(r.get("_atom_site.type_symbol"), "C"),
            "name": clean(r.get("_atom_site.auth_atom_id"), clean(r.get("_atom_site.label_atom_id"))),
            "comp": comp, "chain": clean(r.get("_atom_site.auth_asym_id"), clean(r.get("_atom_site.label_asym_id"), "A")),
            "label_chain": clean(r.get("_atom_site.label_asym_id")),
            "seq": clean(r.get("_atom_site.auth_seq_id"), clean(r.get("_atom_site.label_seq_id"), "1")),
            "entity": clean(r.get("_atom_site.label_entity_id")), "xyz": xyz,
            "occ": clean(r.get("_atom_site.occupancy"), "1.00"), "b": clean(r.get("_atom_site.B_iso_or_equiv"), "0.00"),
        }
        atoms.append(atom)
    return singles, loops, comp_types, atoms


def is_sugar(comp, comp_types):
    t = comp_types.get(comp, "")
    return "SACCHARIDE" in t or comp in {c for p in PRESETS.values() for c in p["components"]}


def distance(a, b):
    return math.dist(a["xyz"], b["xyz"])


def glycan_components(loops, comp_types, atoms):
    residues = defaultdict(list)
    for a in atoms:
        if is_sugar(a["comp"], comp_types):
            residues[Residue(a["chain"], a["seq"], a["comp"], a["entity"])].append(a)
    graph = {r: set() for r in residues}
    # Experiment-deposited covalent connectivity is preferred.
    by_label = {(a["label_chain"], a["seq"], a["comp"]): r for r, aa in residues.items() for a in aa[:1]}
    for row in loop_rows(loops, "_struct_conn."):
        if "covale" not in clean(row.get("_struct_conn.conn_type_id")).lower():
            continue
        keys = []
        for n in ("1", "2"):
            keys.append((
                clean(row.get(f"_struct_conn.ptnr{n}_label_asym_id")),
                clean(row.get(f"_struct_conn.ptnr{n}_auth_seq_id")),
                clean(row.get(f"_struct_conn.ptnr{n}_label_comp_id")),
            ))
        a, b = (by_label.get(k) for k in keys)
        if a and b and a != b:
            graph[a].add(b); graph[b].add(a)
    # Fill missing glycosidic links conservatively using close heavy O/C atoms.
    rs = list(residues)
    for i, a in enumerate(rs):
        for b in rs[i + 1:]:
            if b in graph[a] or a.entity != b.entity:
                continue
            if any(x["element"] != "H" and y["element"] != "H" and distance(x, y) <= 1.9 for x in residues[a] for y in residues[b]):
                graph[a].add(b); graph[b].add(a)
    out, seen = [], set()
    for start in rs:
        if start in seen:
            continue
        stack, component = [start], []
        while stack:
            r = stack.pop()
            if r in seen: continue
            seen.add(r); component.append(r); stack.extend(graph[r] - seen)
        out.append(sorted(component, key=lambda r: (r.chain, int(r.seq) if r.seq.lstrip("-").isdigit() else 999999, r.seq)))
    return sorted(out, key=lambda c: (-len(c), c[0].chain, c[0].seq)), residues


def pdb_atom_line(serial, atom):
    name = atom["name"][:4]
    comp = atom["comp"][:3]
    chain = (atom["chain"] or "A")[:1]
    seq_match = re.search(r"-?\d+", atom["seq"])
    seq = int(seq_match.group()) if seq_match else serial
    x, y, z = atom["xyz"]
    try: occ = float(atom["occ"])
    except ValueError: occ = 1.0
    try: bf = float(atom["b"])
    except ValueError: bf = 0.0
    record = "ATOM  " if atom.get("group") == "ATOM" else "HETATM"
    return f"{record}{serial:5d} {name:>4s} {comp:>3s} {chain}{seq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bf:6.2f}          {atom['element'][:2]:>2s}"


def write_component(path, component, residue_atoms):
    selected = [a for r in component for a in residue_atoms[r]]
    lines = ["REMARK 900 EXPERIMENTAL COORDINATES EXTRACTED FROM PDB MMCIF"]
    for i, atom in enumerate(selected, 1):
        lines.append(pdb_atom_line(i, atom))
    # Jmol infers standard bonds; explicit close contacts preserve cross-residue links.
    for i, a in enumerate(selected):
        bonded = []
        for j, b in enumerate(selected):
            if j <= i or a["element"] == "H" or b["element"] == "H": continue
            cutoff = 1.9 if "H" not in (a["element"], b["element"]) else 1.25
            if distance(a, b) <= cutoff: bonded.append(j + 1)
        if bonded: lines.append(f"CONECT{i + 1:5d}" + "".join(f"{x:5d}" for x in bonded))
    lines += ["END", ""]
    path.write_text("\n".join(lines))
    return selected


DNA_COMPONENTS = {"DA", "DC", "DG", "DT", "DI"}


def dna_residues(atoms):
    """Return deposited DNA residues grouped by strand, preserving coordinates."""
    residues = defaultdict(list)
    for atom in atoms:
        if atom["comp"].upper() in DNA_COMPONENTS:
            residue = Residue(atom["chain"], atom["seq"], atom["comp"], atom["entity"])
            residues[residue].append(atom)
    strands = defaultdict(list)
    for residue in residues:
        strands[(residue.entity, residue.chain)].append(residue)
    ordered = []
    for strand in strands.values():
        ordered.extend(sorted(strand, key=lambda r: (int(r.seq) if r.seq.lstrip("-").isdigit() else 999999, r.seq)))
    return ordered, residues

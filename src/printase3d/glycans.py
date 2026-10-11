"""Glycans components of the 3DPrintase pipeline."""
from __future__ import annotations

import math
from .config import FRUCTAN_TYPES, GLUCAN_TYPES, GLYCAN_CONFORMATION_PROFILES
from .utils import PipelineError


def _sugar_residue_topologies(mol):
    """Infer pyranose/furanose carbon numbering from an RDKit graph.

    This deliberately uses connectivity rather than atom indices or residue
    names, so the same code works for Glc, Gal, Man, and related aldoses.
    """
    residues = []
    for ring_tuple in mol.GetRingInfo().AtomRings():
        ring = set(ring_tuple)
        oxygens = [i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 8]
        carbons = [i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 6]
        if len(oxygens) != 1 or (len(ring) == 6 and len(carbons) != 5) or (len(ring) == 5 and len(carbons) != 4):
            continue
        if len(ring) not in (5, 6):
            continue
        ring_o = oxygens[0]
        o_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(ring_o).GetNeighbors() if n.GetIdx() in ring]
        if len(o_neighbors) != 2:
            continue

        def exocyclic_carbons(atom_idx):
            return [
                n.GetIdx() for n in mol.GetAtomWithIdx(atom_idx).GetNeighbors()
                if n.GetAtomicNum() == 6 and n.GetIdx() not in ring
            ]

        if len(ring) == 6:
            # C5 bears C6; the other ring-O neighbor is the anomeric C1.
            with_exo = [i for i in o_neighbors if exocyclic_carbons(i)]
            if len(with_exo) != 1:
                continue
            c5 = with_exo[0]
            c1 = next(i for i in o_neighbors if i != c5)
            path = [c1]
            previous = ring_o
            while path[-1] != c5:
                candidates = [
                    n.GetIdx() for n in mol.GetAtomWithIdx(path[-1]).GetNeighbors()
                    if n.GetIdx() in ring and n.GetIdx() != previous and n.GetAtomicNum() == 6
                ]
                if len(candidates) != 1:
                    path = []
                    break
                previous, nxt = path[-1], candidates[0]
                path.append(nxt)
            if len(path) != 5:
                continue
            positions = {number: atom_idx for number, atom_idx in enumerate(path, 1)}
            positions[6] = exocyclic_carbons(c5)[0]
            ring_kind, anomeric = "pyranose", 1
        else:
            # In a ketohexofuranose C2 has more exocyclic heavy-atom
            # substitution than C5. Walk C2-C3-C4-C5 around the ring.
            scores = {
                i: sum(1 for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() not in ring and n.GetAtomicNum() > 1)
                for i in o_neighbors
            }
            c2 = max(o_neighbors, key=lambda i: scores[i])
            c5 = next(i for i in o_neighbors if i != c2)
            if scores[c2] == scores[c5]:
                continue
            path = [c2]
            previous = ring_o
            while path[-1] != c5:
                candidates = [
                    n.GetIdx() for n in mol.GetAtomWithIdx(path[-1]).GetNeighbors()
                    if n.GetIdx() in ring and n.GetIdx() != previous and n.GetAtomicNum() == 6
                ]
                if len(candidates) != 1:
                    path = []
                    break
                previous, nxt = path[-1], candidates[0]
                path.append(nxt)
            if len(path) != 4:
                continue
            positions = {number: atom_idx for number, atom_idx in zip((2, 3, 4, 5), path)}
            c1s, c6s = exocyclic_carbons(c2), exocyclic_carbons(c5)
            if len(c1s) != 1 or len(c6s) != 1:
                continue
            positions[1], positions[6] = c1s[0], c6s[0]
            ring_kind, anomeric = "furanose", 2
        residues.append({
            "ring": ring,
            "ring_kind": ring_kind,
            "ring_oxygen": ring_o,
            "anomeric_position": anomeric,
            "positions": positions,
            "atom_to_position": {atom_idx: number for number, atom_idx in positions.items()},
        })
    return residues


def infer_glycosidic_torsions(mol):
    """Return standard heavy-atom phi/psi/(omega) torsions for each link."""
    residues = _sugar_residue_topologies(mol)
    atom_residue = {
        atom_idx: residue
        for residue in residues
        for atom_idx in residue["atom_to_position"]
    }
    links = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.IsInRing() or atom.GetDegree() != 2:
            continue
        endpoints = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(endpoints) != 2 or any(i not in atom_residue for i in endpoints):
            continue
        ra, rb = (atom_residue[i] for i in endpoints)
        if ra is rb:
            continue
        pa, pb = ra["atom_to_position"][endpoints[0]], rb["atom_to_position"][endpoints[1]]
        candidates = []
        if pa == ra["anomeric_position"]:
            candidates.append((ra, endpoints[0], rb, endpoints[1], pb))
        if pb == rb["anomeric_position"]:
            candidates.append((rb, endpoints[1], ra, endpoints[0], pa))
        if len(candidates) != 1:
            continue
        donor, donor_c, acceptor, acceptor_c, acceptor_position = candidates[0]
        toward_ring_position = acceptor_position - 1
        if acceptor_position == 1:
            toward_ring_position = acceptor["anomeric_position"]
        if toward_ring_position not in acceptor["positions"]:
            continue
        toward_ring = acceptor["positions"][toward_ring_position]
        bridge_o = atom.GetIdx()
        link = {
            "donor_ring_kind": donor["ring_kind"],
            "donor_position": donor["anomeric_position"],
            "acceptor_position": acceptor_position,
            "bridge_oxygen": bridge_o,
            "phi": (donor["ring_oxygen"], donor_c, bridge_o, acceptor_c),
            "psi": (donor_c, bridge_o, acceptor_c, toward_ring),
        }
        if acceptor_position in (1, 6):
            # Standard exocyclic definitions are O1-C1-C2-C3 for 2->1
            # fructans and O6-C6-C5-C4 for 2->6 fructans.
            further_position = toward_ring_position + (1 if acceptor_position == 1 else -1)
            if further_position in acceptor["positions"]:
                link["omega"] = (
                    bridge_o, acceptor_c, toward_ring,
                    acceptor["positions"][further_position],
                )
        links.append(link)
    return links


def apply_glycan_conformation(mol, profile_name, *, max_iters=5000):
    """Apply a canonical linkage profile and perform restrained MMFF relaxation."""
    from rdkit.Chem import AllChem, rdMolTransforms

    if profile_name not in GLYCAN_CONFORMATION_PROFILES:
        raise PipelineError(f"Unknown glycan conformation profile: {profile_name}")
    profile = GLYCAN_CONFORMATION_PROFILES[profile_name]
    links = infer_glycosidic_torsions(mol)
    if not links:
        raise PipelineError("No glycosidic torsions could be inferred from the generated structure")
    conf = mol.GetConformer()
    # Two passes make this independent of atom ordering: rotating a downstream
    # linkage on pass one cannot leave an upstream linkage off target.
    for _ in range(2):
        for link in links:
            for angle_name in ("phi", "psi", "omega"):
                if angle_name in link and angle_name in profile:
                    rdMolTransforms.SetDihedralDeg(conf, *link[angle_name], float(profile[angle_name]))

    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    if props is None:
        raise PipelineError("MMFF94s parameters are incomplete for restrained glycan minimization")
    forcefield = AllChem.MMFFGetMoleculeForceField(mol, props)
    constrained = 0
    for link in links:
        for angle_name in ("phi", "psi", "omega"):
            if angle_name in link and angle_name in profile:
                forcefield.MMFFAddTorsionConstraint(*link[angle_name], True, -3.0, 3.0, 100.0)
                constrained += 1
    from .carbonyl_geometry import add_carbonyl_restraints
    add_carbonyl_restraints(forcefield, mol, -1, "MMFF94s")
    forcefield.Initialize()
    status = forcefield.Minimize(maxIts=int(max_iters), forceTol=1e-4, energyTol=1e-6)
    measured = []
    for link in links:
        values = {
            key: round(rdMolTransforms.GetDihedralDeg(conf, *link[key]), 3)
            for key in ("phi", "psi", "omega") if key in link
        }
        measured.append(values)
    return {
        "profile": profile_name,
        "shape": profile["shape"],
        "source": profile["source"],
        "targets_degrees": {k: profile[k] for k in ("phi", "psi", "omega") if k in profile},
        "measured_degrees": measured,
        "link_count": len(links),
        "constrained_torsion_count": constrained,
        "minimization": "MMFF94s with local carbonyl/amide restraints and glycosidic torsions restrained to +/-3 degrees",
        "minimization_converged": status == 0,
        "final_mmff_energy": round(forcefield.CalcEnergy(), 6),
    }


def glycan_backbone_metrics(mol):
    """Measure the ordered sugar-ring centerline for QA and model selection."""
    import numpy as np

    residues = _sugar_residue_topologies(mol)
    conf = mol.GetConformer()
    centers = {
        id(residue): np.mean([
            [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
            for i in residue["ring"]
        ], axis=0)
        for residue in residues
    }
    atom_residue = {i: id(r) for r in residues for i in r["atom_to_position"]}
    graph = {key: [] for key in centers}
    for atom in mol.GetAtoms():
        endpoints = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if atom.GetAtomicNum() != 8 or atom.IsInRing() or len(endpoints) != 2:
            continue
        if not all(i in atom_residue for i in endpoints):
            continue
        a, b = (atom_residue[i] for i in endpoints)
        if a != b and b not in graph[a]:
            graph[a].append(b); graph[b].append(a)
    terminals = [key for key, neighbors in graph.items() if len(neighbors) == 1]
    if len(residues) < 2 or len(terminals) != 2 or any(len(v) > 2 for v in graph.values()):
        return {"ring_count": len(residues), "linear_centerline_found": False}
    order, previous, current = [], None, terminals[0]
    while current is not None:
        order.append(current)
        following = [key for key in graph[current] if key != previous]
        previous, current = current, (following[0] if following else None)
    points = np.asarray([centers[key] for key in order])
    steps = np.diff(points, axis=0)
    lengths = np.linalg.norm(steps, axis=1)
    turns = [
        float(np.dot(steps[i - 1], steps[i]) / (lengths[i - 1] * lengths[i]))
        for i in range(1, len(steps)) if lengths[i - 1] and lengths[i]
    ]
    contour = float(lengths.sum())
    end_to_end = float(np.linalg.norm(points[-1] - points[0]))
    return {
        "ring_count": len(residues),
        "linear_centerline_found": True,
        "end_to_end_angstrom": round(end_to_end, 4),
        "contour_angstrom": round(contour, 4),
        "end_to_contour_ratio": round(end_to_end / contour, 6) if contour else 0.0,
        "step_length_stddev": round(float(np.std(lengths)), 6),
        "turn_cosine_stddev": round(float(np.std(turns)), 6) if turns else 0.0,
    }


def embed_glycan_conformer(template, *, seed, profile_name=None):
    """Embed candidates and select a regular canonical carbohydrate backbone."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    profile = GLYCAN_CONFORMATION_PROFILES.get(profile_name) if profile_name else None
    candidate_seeds = [int(seed)] if not profile else [int(seed), 2, 8, 270911]
    candidates, errors = [], []
    for attempt, candidate_seed in enumerate(dict.fromkeys(candidate_seeds)):
        candidate = Chem.Mol(template)
        candidate.RemoveAllConformers()
        params = AllChem.ETKDGv3()
        params.randomSeed = candidate_seed
        params.useRandomCoords = False
        params.maxIterations = 2000
        params.numThreads = 1
        if AllChem.EmbedMolecule(candidate, params) != 0:
            errors.append(f"seed {candidate_seed}: embedding failed")
            continue
        if not AllChem.MMFFHasAllMoleculeParams(candidate):
            raise PipelineError("MMFF94s parameters are incomplete for this glycan")
        from .carbonyl_geometry import add_carbonyl_restraints
        properties = AllChem.MMFFGetMoleculeProperties(candidate, mmffVariant="MMFF94s")
        ff = AllChem.MMFFGetMoleculeForceField(candidate, properties)
        add_carbonyl_restraints(ff, candidate, -1, "MMFF94s")
        ff.Initialize()
        status = ff.Minimize(maxIts=5000)
        if status != 0:
            errors.append(f"seed {candidate_seed}: initial MMFF94s minimization did not converge")
            continue
        details = apply_glycan_conformation(candidate, profile_name) if profile else None
        metrics = glycan_backbone_metrics(candidate)
        if not metrics.get("linear_centerline_found"):
            errors.append(f"seed {candidate_seed}: could not trace the sugar-ring centerline")
            continue
        regularity_penalty = metrics["step_length_stddev"] + 3.0 * metrics["turn_cosine_stddev"]
        if profile and profile["shape"].startswith("extended"):
            score = metrics["end_to_end_angstrom"] - regularity_penalty
        elif profile and "helix_residues_per_turn" in profile:
            steps = max(1, metrics["ring_count"] - 1)
            residues_per_turn = profile["helix_residues_per_turn"]
            pitch = profile["helix_pitch_angstrom"]
            radius = profile["helix_radius_angstrom"]
            axial_rise = pitch * steps / residues_per_turn
            chord = 2.0 * radius * abs(math.sin(math.pi * steps / residues_per_turn))
            expected_end = math.hypot(axial_rise, chord)
            score = -abs(metrics["end_to_end_angstrom"] - expected_end) - 0.15 * regularity_penalty
        else:
            # A repeated torsion profile should produce a regular screw/helix.
            # Prefer geometric repeat regularity, then the restrained energy.
            energy = details["final_mmff_energy"] if details else 0.0
            score = -regularity_penalty - 1e-6 * energy
        candidates.append((score, candidate_seed, candidate, details, metrics))
    if not candidates:
        raise PipelineError("Could not generate a usable glycan conformer; " + "; ".join(errors))
    _, selected_seed, selected, details, metrics = max(candidates, key=lambda item: item[0])
    selection = {
        "candidate_count": len(candidates),
        "selected_embedding_seed": selected_seed,
        "selection_policy": (
            "maximum extended centerline with repeat-regularity penalty"
            if profile and profile["shape"].startswith("extended")
            else "closest canonical helix dimensions with repeat-regularity penalty"
            if profile and "helix_residues_per_turn" in profile
            else "minimum repeat-irregularity with restrained-energy tie-breaker"
        ),
        "backbone_metrics": metrics,
    }
    if details:
        details["candidate_selection"] = selection
    return selected, details, selection


def build_fructan(kind, dp, sdf, seed=61453, conformation="natural"):
    """Build a deterministic linear fructose-only fructan oligomer."""
    if conformation not in {"natural", "generic"}:
        raise PipelineError("conformation must be 'natural' or 'generic'")
    definition = FRUCTAN_TYPES[kind]
    sequence = "Fruf" + "".join(f"({definition['linkage']})Fruf" for _ in range(dp - 1))
    try:
        from glyles import convert
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolDescriptors
    except ImportError as exc:
        raise PipelineError(
            "Fructan generation requires GlyLES and RDKit. Recreate the "
            "3D-glycan Conda environment from environment.yml."
        ) from exc
    converted = convert(glycan=sequence)
    if not converted or not converted[0][1]:
        raise PipelineError(f"GlyLES could not convert fructan sequence: {sequence}")
    smiles = converted[0][1]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise PipelineError(f"RDKit could not parse GlyLES output for {kind} DP{dp}")
    mol = Chem.AddHs(mol)
    profile_name = definition["conformation_profile"] if conformation == "natural" else None
    mol, conformation_details, selection = embed_glycan_conformer(
        mol, seed=seed, profile_name=profile_name
    )
    mol.SetProp("_Name", f"{definition['name']} DP{dp}")
    mol.SetProp("IUPAC_condensed", sequence)
    mol.SetProp("linkage", definition["linkage_label"])
    writer = Chem.SDWriter(str(sdf))
    writer.write(mol)
    writer.close()
    return {
        "iupac_condensed": sequence,
        "isomeric_smiles": Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True),
        "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
        "atom_count_with_hydrogens": mol.GetNumAtoms(),
        "coordinate_method": (
            f"GlyLES sequence conversion; RDKit ETKDGv3 seed={seed}; "
            + (f"{conformation_details['minimization']}; canonical {conformation_details['shape']}"
               if conformation_details else "MMFF94s with local carbonyl/amide restraints (generic legacy mode)")
        ),
        "conformation": conformation,
        "conformation_details": conformation_details,
        "backbone_metrics": selection["backbone_metrics"],
    }


def build_glucan(kind, dp, sdf, seed=61453, conformation="natural"):
    """Build a deterministic linear (1→4)-linked pyranose oligomer."""
    if conformation not in {"natural", "generic"}:
        raise PipelineError("conformation must be 'natural' or 'generic'")
    definition = GLUCAN_TYPES[kind]
    residue = definition["residue"]
    sequence = residue + "".join(f"({definition['linkage']}){residue}" for _ in range(dp - 1))
    try:
        from glyles import convert
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolDescriptors
    except ImportError as exc:
        raise PipelineError(
            "Glucan generation requires GlyLES and RDKit. Recreate the "
            "3D-glycan Conda environment from environment.yml."
        ) from exc
    converted = convert(glycan=sequence)
    if not converted or not converted[0][1]:
        raise PipelineError(f"GlyLES could not convert glucan sequence: {sequence}")
    smiles = converted[0][1]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise PipelineError(f"RDKit could not parse GlyLES output for {kind} DP{dp}")
    mol = Chem.AddHs(mol)
    profile_name = definition["conformation_profile"] if conformation == "natural" else None
    mol, conformation_details, selection = embed_glycan_conformer(
        mol, seed=seed, profile_name=profile_name
    )
    mol.SetProp("_Name", f"{definition['series_name']} DP{dp}")
    mol.SetProp("IUPAC_condensed", sequence)
    mol.SetProp("linkage", definition["linkage_label"])
    writer = Chem.SDWriter(str(sdf))
    writer.write(mol)
    writer.close()
    return {
        "iupac_condensed": sequence,
        "isomeric_smiles": Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True),
        "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
        "atom_count_with_hydrogens": mol.GetNumAtoms(),
        "coordinate_method": (
            f"GlyLES sequence conversion; RDKit ETKDGv3 seed={seed}; "
            + (f"{conformation_details['minimization']}; canonical {conformation_details['shape']}"
               if conformation_details else "MMFF94s with local carbonyl/amide restraints (generic legacy mode)")
        ),
        "conformation": conformation,
        "conformation_details": conformation_details,
        "backbone_metrics": selection["backbone_metrics"],
    }

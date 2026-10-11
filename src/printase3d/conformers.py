"""Conformer sampling, minimization, shape measurement, and selection."""
from __future__ import annotations

from .utils import PipelineError
from .stereochemistry import preserves_stereochemistry
from .carbonyl_geometry import add_carbonyl_restraints, carbonyl_measurements
from .extended_geometry import extension_definition, add_extension_restraint, extension_measurement
from .aromatic_geometry import aromatic_components, add_planarity_restraints, ring_planarity, flavin_geometry


FLATNESS_TOLERANCE_ANGSTROM = 1e-4


def conformer_flatness(mol, conf_id=-1):
    """Heavy-atom RMS distance from a least-squares plane, in angstroms.

    Translation/rotation invariant; smaller means flatter. One to three atoms
    always fit a plane. Hydrogen-only structures use all atoms instead.
    This measures existing coordinates; it never projects or deforms them.
    """
    import numpy as np
    indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    if not indices:
        indices = list(range(mol.GetNumAtoms()))
    xyz = np.asarray(mol.GetConformer(conf_id).GetPositions())[indices]
    if not len(xyz) or not np.isfinite(xyz).all():
        raise PipelineError("Cannot measure flatness of missing/nonfinite coordinates")
    if len(xyz) <= 3:
        return 0.0
    singular = np.linalg.svd(xyz - xyz.mean(axis=0), compute_uv=False)
    return float(singular[-1] / np.sqrt(len(xyz)))


def select_conformer(candidates, selection="energy"):
    """Select a minimized candidate by energy or the requested shape extreme.

    Flatness modes require converged minima. Within 0.0001 angstrom of the
    shape extreme, choose the lowest energy to ignore numerical shape noise.
    Energy mode retains the historical unconverged fallback.
    """
    import math
    if selection not in {"energy", "most-flat", "least-flat", "extended"}:
        raise PipelineError(f"Unknown conformer selection: {selection}")
    finite = [c for c in candidates if math.isfinite(c['energy'])]
    for candidate in finite:
        candidate['flatness_rms_angstrom'] = conformer_flatness(candidate['mol'], candidate['conf_id'])
    converged = [c for c in finite if c['status'] == 0]
    if selection != "energy" and not converged:
        raise PipelineError("Shape selection requires converged minimized conformers; increase --max-opt-iters")
    pool = converged or finite
    if not pool:
        raise PipelineError("No finite-energy conformer candidates")
    if selection == "extended":
        extreme = max(c['terminal_separation_angstrom'] for c in pool)
        pool = [c for c in pool if extreme-c['terminal_separation_angstrom'] <= 1e-4]
    elif selection != "energy":
        extreme = (min if selection == "most-flat" else max)(c['flatness_rms_angstrom'] for c in pool)
        pool = [c for c in pool if abs(c['flatness_rms_angstrom'] - extreme) <= FLATNESS_TOLERANCE_ANGSTROM]
    return min(pool, key=lambda c: c['energy'])


def sample_lowest_energy_conformer(
    mol, *, seed, conformer_count, prune_rms, max_iters,
    current_origin="supplied-coordinate start", threads=1, selection="energy",
):
    """Minimize a source-started candidate plus a deterministic 3D ensemble."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    extension = extension_definition(mol) if selection == "extended" else None
    aromatic_groups = aromatic_components(mol)
    use_mmff = AllChem.MMFFHasAllMoleculeParams(Chem.Mol(mol))
    forcefield = "MMFF94s" if use_mmff else "UFF"
    candidates = []
    rejected_stereochemistry = []

    def optimize_single(candidate, origin):
        if use_mmff:
            props = AllChem.MMFFGetMoleculeProperties(candidate, mmffVariant="MMFF94s")
            ff = AllChem.MMFFGetMoleculeForceField(candidate, props)
        else:
            ff = AllChem.UFFGetMoleculeForceField(candidate)
        add_carbonyl_restraints(ff, candidate, -1, forcefield)
        if aromatic_groups:
            add_planarity_restraints(ff, candidate, -1, aromatic_groups, forcefield)
        if extension:
            add_extension_restraint(ff, extension, forcefield)
        ff.Initialize()
        status = ff.Minimize(maxIts=max_iters)
        if not preserves_stereochemistry(mol, candidate):
            rejected_stereochemistry.append(origin)
            return
        separation = extension_measurement(candidate, extension) if extension else None
        candidates.append({
            "terminal_separation_angstrom": separation,
            "energy": float(ff.CalcEnergy()), "status": int(status),
            "origin": origin, "mol": candidate, "conf_id": candidate.GetConformer().GetId(),
        })

    if mol.GetNumConformers() and mol.GetConformer().Is3D():
        optimize_single(Chem.Mol(mol), current_origin)

    sampled = Chem.Mol(mol)
    sampled.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useRandomCoords = True
    params.pruneRmsThresh = float(prune_rms)
    params.maxIterations = 3000
    params.numThreads = threads
    conformer_ids = list(AllChem.EmbedMultipleConfs(
        sampled, numConfs=int(conformer_count), params=params
    ))
    if not conformer_ids and not candidates:
        raise PipelineError("RDKit could not generate any conformers for energy ranking")
    if conformer_ids:
        for conf in list(sampled.GetConformers()):
            candidate = Chem.Mol(sampled)
            candidate.RemoveAllConformers()
            candidate.AddConformer(Chem.Conformer(conf), assignId=True)
            optimize_single(candidate, "ETKDG ensemble")

    converged = [candidate for candidate in candidates if candidate["status"] == 0]
    pool = converged or candidates
    selected = select_conformer(candidates, selection)
    output = Chem.Mol(mol)
    output.RemoveAllConformers()
    output.AddConformer(
        Chem.Conformer(selected["mol"].GetConformer(selected["conf_id"])), assignId=True
    )
    if extension and extension_measurement(output, extension) < extension["minimum_ring_centroid_distance_angstrom"]:
        raise PipelineError("Extended conformer failed terminal-ring separation validation")
    carbonyl = carbonyl_measurements(output)
    if any(site["height_angstrom"] > 0.03 for site in carbonyl):
        raise PipelineError("Local carbonyl/amide planarity validation failed")
    planarity = ring_planarity(output, aromatic_groups)
    if any(value > 0.03 for value in planarity):
        raise PipelineError("Aromatic fused-ring planarity validation failed after minimization")
    ranked = sorted(candidates, key=lambda candidate: candidate["energy"])
    source_candidate = next(
        (candidate for candidate in candidates if candidate["origin"] != "ETKDG ensemble"), None
    )
    return output, {
        "method": "deterministic ETKDGv3 conformer ensemble followed by force-field minimization",
        "forcefield": forcefield,
        "rejected_stereochemistry_candidates": len(rejected_stereochemistry),
        "carbonyl_geometry": carbonyl,
        "extension": {**extension, "selected_terminal_separation_angstrom": extension_measurement(output, extension)} if extension else None,
        "flavin_geometry": flavin_geometry(output),
        "aromatic_planarity": {"atom_groups": aromatic_groups, "rms_angstrom": planarity,
            "maximum_allowed_rms_angstrom": 0.03, "torsion_window_degrees": 0.25,
            "force_constant": 100.0, "scope": "5–14 atom aromatic ring components; excludes reduced flavin bridges"},
        "requested_etkdg_conformers": int(conformer_count),
        "embedded_etkdg_conformers": len(conformer_ids),
        "total_candidates_including_source": len(candidates),
        "converged_candidates": len(converged),
        "prune_rms_threshold_angstrom": float(prune_rms),
        "selected_origin": selected["origin"],
        "selected_energy_kcal_mol": round(selected["energy"], 6),
        "selection": selection,
        "flatness_metric": "heavy-atom RMS distance to best-fit plane (angstrom); all atoms if hydrogen-only",
        "flatness_tie_tolerance_angstrom": FLATNESS_TOLERANCE_ANGSTROM,
        "selected_flatness_rms_angstrom": selected["flatness_rms_angstrom"],
        "selected_converged": selected["status"] == 0,
        "minimum_energy_kcal_mol": round(min(c["energy"] for c in pool), 6),
        "selected_above_minimum_kcal_mol": round(selected["energy"] - min(c["energy"] for c in pool), 6),
        "candidates": [
            {"terminal_separation_angstrom": c.get("terminal_separation_angstrom"), "origin": c["origin"], "converged": c["status"] == 0,
             "energy_kcal_mol": round(c["energy"], 6),
             "flatness_rms_angstrom": c.get("flatness_rms_angstrom"),
             "selected": c is selected}
            for c in candidates
        ],
        "source_start_energy_kcal_mol": (
            round(source_candidate["energy"], 6) if source_candidate else None
        ),
        "selected_below_source_kcal_mol": (
            round(source_candidate["energy"] - selected["energy"], 6)
            if source_candidate else None
        ),
        "next_energy_kcal_mol": round(ranked[1]["energy"], 6) if len(ranked) > 1 else None,
        "ranking_scope": "force-field potential energy including carbonyl/amide, aromatic planarity and any extension restraints; not aqueous free energy",
    }

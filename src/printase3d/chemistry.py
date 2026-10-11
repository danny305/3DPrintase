"""Chemistry components of the 3DPrintase pipeline."""
from __future__ import annotations

from pathlib import Path
from .utils import PipelineError
from .stereochemistry import preserves_stereochemistry
from .carbonyl_geometry import add_carbonyl_restraints, carbonyl_measurements
from .conformers import sample_lowest_energy_conformer
from .aromatic_geometry import aromatic_components, add_planarity_restraints, flavin_geometry, ring_planarity
from .cofactor_protonation import correct_cofactor_ring_cations
from .protonation import correct_free_selenocysteine


def rdkit_read_file(path):
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise PipelineError("General molecule input requires RDKit; reinstall the 3DPrintase environment") from exc
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise PipelineError(f"Molecule input file not found: {path}")
    suffix = path.suffix.lower()
    if suffix in {".sdf", ".sd"}:
        mol = next((m for m in Chem.SDMolSupplier(str(path), removeHs=False) if m is not None), None)
    elif suffix in {".mol", ".mdl"}:
        mol = Chem.MolFromMolFile(str(path), removeHs=False)
    elif suffix == ".mol2":
        mol = Chem.MolFromMol2File(str(path), removeHs=False)
    elif suffix in {".pdb", ".ent"}:
        mol = Chem.MolFromPDBFile(str(path), removeHs=False)
    else:
        raise PipelineError("Supported molecule files: .sdf, .sd, .mol, .mdl, .mol2, .pdb, .ent")
    if mol is None:
        raise PipelineError(f"RDKit could not parse molecule file: {path}")
    return mol, path


def correct_false_amide_protonation(mol):
    """Neutralize protonated amide nitrogens from a known Dimorphite edge case."""
    from rdkit import Chem

    pattern = Chem.MolFromSmarts("[N+;!$(N~[!#6])]-[C,c](=O)")
    corrected = 0
    for match in mol.GetSubstructMatches(pattern):
        nitrogen = mol.GetAtomWithIdx(match[0])
        if nitrogen.GetFormalCharge() != 1:
            continue
        nitrogen.SetFormalCharge(0)
        explicit_h = nitrogen.GetNumExplicitHs()
        if explicit_h:
            nitrogen.SetNumExplicitHs(explicit_h - 1)
        nitrogen.SetNoImplicit(False)
        corrected += 1
    if corrected:
        Chem.SanitizeMol(mol)
    return corrected


def transfer_heavy_atom_coordinates(source, target):
    """Copy a stereochemistry-aware full heavy-atom coordinate match.

    Protonation changes formal charge and hydrogen count, but normally leave
    the covalent heavy-atom graph untouched. An element/bond MCS lets that
    graph retain deposited coordinates rather than being re-embedded.
    """
    from rdkit import Chem
    from rdkit.Chem import rdFMCS

    if not source.GetNumConformers() or not source.GetConformer().Is3D():
        return False, "source has no 3D conformer"
    source_heavy = Chem.RemoveHs(source)
    target_heavy = Chem.RemoveHs(target)
    if source_heavy.GetNumAtoms() != target_heavy.GetNumAtoms():
        return False, "protonation changed the heavy-atom count"
    result = rdFMCS.FindMCS(
        [source_heavy, target_heavy],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        # Equivalent absolute stereochemistry may carry opposite RDKit CW/CCW
        # tags after SMILES atom reordering. Connectivity mapping is followed
        # by unchanged source coordinates, so strict tag-text matching here
        # would incorrectly reject a valid full graph match.
        matchChiralTag=False,
        timeout=10,
    )
    if result.canceled or result.numAtoms != target_heavy.GetNumAtoms():
        return False, f"only {result.numAtoms}/{target_heavy.GetNumAtoms()} heavy atoms matched"
    query = Chem.MolFromSmarts(result.smartsString)
    source_match = source_heavy.GetSubstructMatch(query, useChirality=False)
    target_match = target_heavy.GetSubstructMatch(query, useChirality=False)
    if len(source_match) != result.numAtoms or len(target_match) != result.numAtoms:
        return False, "stereochemistry-aware heavy-atom mapping failed"
    source_conf = source_heavy.GetConformer()
    target.RemoveAllConformers()
    target_conf = Chem.Conformer(target.GetNumAtoms())
    target_conf.Set3D(True)
    for source_idx, target_idx in zip(source_match, target_match):
        target_conf.SetAtomPosition(target_idx, source_conf.GetAtomPosition(source_idx))
    target.AddConformer(target_conf, assignId=True)
    return True, f"transferred {result.numAtoms} deposited heavy-atom coordinates"


def prepare_molecule(mol, args, label):
    """Standardize fragments, add H, establish 3D coordinates, and optimize."""
    most_flat = getattr(args, "most_flat", False)
    least_flat = getattr(args, "least_flat", False)
    extended = getattr(args, "extended", False)
    if sum((most_flat, least_flat, extended)) > 1:
        raise PipelineError("--most-flat, --least-flat, and --extended are mutually exclusive")
    selection = "extended" if extended else "most-flat" if most_flat else "least-flat" if least_flat else "energy"
    if selection != "energy":
        if args.geometry_mode in {"pubchem", "rcsb"} or args.optimize != "auto":
            raise PipelineError("Shape-selection flags imply ensemble minimization and cannot be combined with --pubchem-geometry or explicit --optimize")
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolDescriptors
    except ImportError as exc:
        raise PipelineError("General molecule preparation requires RDKit") from exc
    if args.conformer_count < 1:
        raise PipelineError("--conformer-count must be at least 1")
    if args.conformer_prune_rms < 0:
        raise PipelineError("--conformer-prune-rms cannot be negative")

    if not args.keep_all_fragments:
        fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if fragments:
            mol = max(fragments, key=lambda item: item.GetNumHeavyAtoms())
    coordinate_template = Chem.Mol(mol)
    coordinate_transfer = {
        "attempted": False,
        "succeeded": False,
        "detail": "not required",
    }
    protonation = {
        "target_ph": None,
        "protonation_method": "source/default protonation state",
        "protonation_variants": 1,
        "protonation_variant_selected": 1,
    }
    if args.ph is not None:
        try:
            from dimorphite_dl import protonate_smiles
        except ImportError as exc:
            raise PipelineError("--ph requires Dimorphite-DL; reinstall the 3DPrintase environment") from exc
        source_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
        generated = protonate_smiles(
            source_smiles, ph_min=args.ph, ph_max=args.ph,
            precision=args.ph_precision, max_variants=args.max_protonation_variants,
        )
        canonical = []
        corrections_by_smiles = {}
        amide_corrections = 0
        for candidate in generated:
            candidate_mol = Chem.MolFromSmiles(candidate)
            if candidate_mol is not None:
                amide_corrections += correct_false_amide_protonation(candidate_mol)
                corrections = correct_free_selenocysteine(candidate_mol, args.ph)
                corrections += correct_cofactor_ring_cations(candidate_mol, args.ph)
                value = Chem.MolToSmiles(candidate_mol, isomericSmiles=True)
                corrections_by_smiles[value] = corrections
                if value not in canonical:
                    canonical.append(value)
        canonical.sort()
        if not canonical:
            raise PipelineError(f"Dimorphite-DL produced no valid protonation state for {label} at pH {args.ph:g}")
        if not 1 <= args.protonation_variant <= len(canonical):
            raise PipelineError(
                f"--protonation-variant {args.protonation_variant} is out of range; "
                f"{len(canonical)} state(s) were generated at pH {args.ph:g}"
            )
        selected_smiles = canonical[args.protonation_variant - 1]
        mol = Chem.MolFromSmiles(selected_smiles)
        coordinate_transfer["attempted"] = True
        transferred, transfer_detail = transfer_heavy_atom_coordinates(coordinate_template, mol)
        coordinate_transfer.update({"succeeded": transferred, "detail": transfer_detail})
        protonation = {
            "target_ph": args.ph,
            "protonation_method": (
                f"Dimorphite-DL exact-pH ionization (pH={args.ph:g}, precision={args.ph_precision:g}); "
                f"variant {args.protonation_variant} of {len(canonical)}"
                + ("; curated site-specific ionization correction (see provenance)"
                   if corrections_by_smiles.get(selected_smiles) else "")
            ),
            "protonation_source_smiles": source_smiles,
            "protonated_smiles": selected_smiles,
            "protonation_variants": len(canonical),
            "protonation_variant_selected": args.protonation_variant,
            "protonation_amide_corrections": amide_corrections,
            "protonation_curated_corrections": corrections_by_smiles.get(selected_smiles, []),
        }
        # When the heavy graph is unchanged, deposited coordinates have been
        # transferred above. Only a genuinely changed/unmappable graph falls
        # back to fresh embedding.
        if not transferred:
            mol.RemoveAllConformers()
    formal_charge = int(Chem.GetFormalCharge(mol))
    if args.charge is not None and formal_charge != args.charge:
        raise PipelineError(
            f"Requested formal charge {args.charge}, but the supplied structure has charge {formal_charge}. "
            "Provide the desired protonation/ionization state as SMILES or InChI."
        )

    had_3d = mol.GetNumConformers() > 0 and mol.GetConformer().Is3D()
    if args.coordinates == "regenerate" and mol.GetNumConformers():
        mol.RemoveAllConformers()
        had_3d = False
    # Hydrogen visibility is a rendering choice: all coordinate preparation,
    # optimization, and validation must use explicit H, even for --no-hydrogens.
    # Keep them in the SDF; the Jmol export removes them after preparation.
    mol = Chem.AddHs(mol, addCoords=had_3d)
    embedded = False
    if not had_3d:
        params = AllChem.ETKDGv3()
        params.randomSeed = int(args.seed)
        params.useRandomCoords = True
        params.maxIterations = 3000
        params.numThreads = getattr(args, "threads", 1)
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise PipelineError(f"RDKit could not generate a 3D conformer for {label}")
        embedded = True

    method = "preserved supplied 3D coordinates"
    if coordinate_transfer["succeeded"]:
        method = "preserved supplied heavy-atom 3D coordinates across protonation change"
    if embedded:
        method = f"RDKit ETKDGv3 conformer, seed={args.seed}"
    stereochemistry_template = Chem.Mol(mol)
    aromatic_groups = aromatic_components(mol)
    optimize = args.optimize
    conformer_search = None
    geometry_mode = "ensemble" if selection != "energy" else args.geometry_mode
    if geometry_mode in {"pubchem", "rcsb"}:
        if args.coordinates == "regenerate":
            raise PipelineError("Source-coordinate preservation conflicts with --coordinates regenerate")
        if not had_3d:
            raise PipelineError("Source-coordinate preservation requires usable supplied/PubChem 3D coordinates")
        if optimize not in {"auto", "hydrogens", "none"}:
            raise PipelineError("Source-coordinate preservation supports --optimize auto, hydrogens, or none")
        if optimize == "auto":
            optimize = "hydrogens"
    elif geometry_mode == "ensemble":
        if optimize not in {"auto"}:
            raise PipelineError("--ensemble controls optimization; do not combine it with an explicit --optimize mode")
        optimize = "ensemble"
    if optimize != "none":
        chosen = optimize
        if chosen == "auto":
            # A changed charge state or a regenerated structure requires a
            # global search. Unchanged deposited structures retain their
            # experimental/PubChem heavy-atom conformer by default.
            chosen = "ensemble" if args.ph is not None or not had_3d or args.coordinates == "regenerate" else "hydrogens"
        if chosen == "ensemble":
            current_origin = (
                "transferred supplied-coordinate start" if coordinate_transfer["succeeded"]
                else "initial ETKDG start" if embedded
                else "supplied-coordinate start"
            )
            mol, conformer_search = sample_lowest_energy_conformer(
                mol, seed=args.seed, conformer_count=args.conformer_count,
                prune_rms=args.conformer_prune_rms, max_iters=args.max_opt_iters,
                threads=getattr(args, "threads", 1),
                current_origin=current_origin,
                selection=selection,
            )
            method += (
                f"; selected {'lowest-energy' if selection == 'energy' else selection} of {conformer_search['total_candidates_including_source']} "
                f"{conformer_search['forcefield']} candidates"
            )
            status = 0 if conformer_search["converged_candidates"] else 1
            forcefield = conformer_search["forcefield"]
        elif chosen == "hydrogens":
            if AllChem.MMFFHasAllMoleculeParams(mol):
                properties = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
                ff = AllChem.MMFFGetMoleculeForceField(mol, properties)
                forcefield = "fixed-heavy-atom MMFF94s"
            else:
                ff = AllChem.UFFGetMoleculeForceField(mol)
                forcefield = "fixed-heavy-atom UFF"
            add_carbonyl_restraints(ff, mol, -1, "MMFF94s" if "MMFF" in forcefield else "UFF")
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() > 1:
                    ff.AddFixedPoint(atom.GetIdx())
            ff.Initialize()
            status = ff.Minimize(maxIts=args.max_opt_iters)
        elif chosen == "mmff":
            if not AllChem.MMFFHasAllMoleculeParams(mol):
                raise PipelineError(f"MMFF94 parameters are incomplete for {label}; use --optimize uff or none")
            properties = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
            ff = AllChem.MMFFGetMoleculeForceField(mol, properties)
            add_carbonyl_restraints(ff, mol, -1, "MMFF94s" if chosen=="mmff" else "UFF")
            add_planarity_restraints(ff, mol, -1, aromatic_groups, "MMFF94")
            ff.Initialize()
            status = ff.Minimize(maxIts=args.max_opt_iters)
            forcefield = "MMFF94s with carbonyl/amide and aromatic planarity restraints"
        else:
            ff = AllChem.UFFGetMoleculeForceField(mol)
            add_carbonyl_restraints(ff, mol, -1, "MMFF94s" if chosen=="mmff" else "UFF")
            add_planarity_restraints(ff, mol, -1, aromatic_groups, "UFF")
            ff.Initialize()
            status = ff.Minimize(maxIts=args.max_opt_iters)
            forcefield = "UFF with aromatic planarity restraints"
        if chosen != "ensemble":
            convergence = "converged" if status == 0 else f"reached iteration limit ({args.max_opt_iters})"
            method += f"; {convergence} {forcefield} minimization"
    else:
        method += "; no force-field minimization"

    if optimize != "none" and any(s['height_angstrom'] > 0.03 for s in carbonyl_measurements(mol)):
        raise PipelineError("Local carbonyl/amide planarity validation failed; preserved heavy coordinates may require a separate regenerated variant")
    if not preserves_stereochemistry(stereochemistry_template, mol):
        raise PipelineError("Minimization changed input stereochemistry")
    if any(x > .03 for x in ring_planarity(mol, aromatic_groups)):
        raise PipelineError("Aromatic ring planarity validation failed; source coordinates may require a regenerated variant")
    # MMFF typing changes aromaticity flags; restore chemical serialization.
    Chem.SanitizeMol(mol)
    mol.SetProp("_Name", label)
    return mol, {
        "aromatic_planarity": {"atom_groups": aromatic_groups, "rms_angstrom": ring_planarity(mol, aromatic_groups), "maximum_allowed_rms_angstrom": .03},
        "carbonyl_geometry": carbonyl_measurements(mol),
        "flavin_geometry": flavin_geometry(mol),
        "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
        "isomeric_smiles": Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True),
        "formal_charge": int(Chem.GetFormalCharge(mol)),
        "atom_count_with_hydrogens": mol.GetNumAtoms(),
        "coordinate_method": method,
        "geometry_mode_requested": geometry_mode,
        "heavy_atom_coordinate_transfer": coordinate_transfer,
        "conformer_search": conformer_search,
        **protonation,
    }

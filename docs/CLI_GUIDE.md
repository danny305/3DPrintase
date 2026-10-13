# CLI and Python API guide

3DPrintase converts general molecules, generated oligosaccharides, and experimentally observed PDB glycans into watertight OBJ/MTL models for 3D printing. It combines PubChem/RCSB retrieval, RDKit chemistry, a copied Jmol rendering script, and mesh repair that replaces fragile open bond shells with capped cylinders.

Read [PIPELINE_OVERVIEW.md](PIPELINE_OVERVIEW.md) for the complete shareable explanation and [BAMBU_STUDIO_PRINTING.md](BAMBU_STUDIO_PRINTING.md) for slicer guidance.

## Install

Follow [installation and troubleshooting](INSTALLATION.md), or the [quick start](../README.md). Commands below are run from the repository root unless otherwise noted.

## Background rendering

Jmol exports run in no-display mode with standard input disconnected. On macOS,
Java runs with `-Dapple.awt.UIElement=true`, preventing normal app/Dock activation
while you work in other applications. No extra flag is needed. The pipeline still
waits for the export and checks its result. This uses the full `Jmol.jar`: the
installed JmolData/headless configuration cannot produce these OBJ exports.

## Short linear peptides

```bash
# Arg–Gly–Trp, written from the N terminus to the C terminus.
3dprintase peptide RGW

# Select pH-7 protonation and an extended conformer for printing.
3dprintase peptide RGW --ph 7 --extended --no-hydrogens

# Generate only the optimized SDF and metadata, without Java/Jmol.
3dprintase peptide rgw --no-render --label Arg-Gly-Trp --output out/peptides
```

The required positional `SEQUENCE` uses the 20 standard one-letter amino-acid
codes `ACDEFGHIKLMNPQRSTVWY`; lowercase is accepted as the same L sequence.
Sequences run N-to-C. The builder creates a single covalent chain with peptide
bonds, L stereochemistry (glycine is achiral), and free, uncapped amino and
carboxyl termini. Ambiguous/modified residues, D-residue notation, separators,
FASTA headers, terminal caps, and cyclization are not supported.

Without `--ph`, the RDKit builder's neutral formal-charge representation is
retained. `--ph PH` uses the same protonation processing as `molecule`.
Coordinates are generated locally, then explicit-hydrogen conformers are
minimized and ranked by the existing molecular pipeline. Linear connectivity
does not imply a straight or biologically bound conformation. `--extended`,
`--most-flat`, `--least-flat`, `--ensemble`, `--conformer-count`, `--seed`, and
`--optimize` are available; deposited-coordinate flags are not applicable.

For uncapped linear peptides, `--extended` applies its separation restraint
between the backbone N-terminus and C-terminal carbonyl carbon, excluding side
chains from the endpoint selection. It selects the most extended converged
candidate while retaining bond geometry and stereochemistry; it does not force
the atoms onto a geometrically straight line.

All rendering controls apply, including `--monocolor-aliphatics` and
`--no-hydrogens`. Hydrogen removal happens only when rendering, after
optimization; the SDF retains the full molecule. Default output for `RGW` is
`out/misc/rgw_peptide-rgw/ball-stick/`, containing SDF, OBJ/MTL (unless
`--no-render`), JSON provenance, and `MODEL.md`. Metadata records the normalized
sequence, residue count, termini, stereochemistry, and reproduction command.
This subcommand accepts one sequence per invocation; JSON molecule manifests
retain their existing source types.

## Parallel molecule batches

Create a JSON list of sources and readable labels, for example `molecules.json`:

```json
[
  {"cid": 79025, "label": "alpha-D-glucose"},
  {"cid": 64689, "label": "beta-D-glucose"},
  {"cid": 439709, "label": "beta-D-fructofuranose"}
]
```

```bash
3dprintase batch molecules.json --ph 7 --ensemble --ball-stick --output out/sugars
# Limit the batch to four workers, leaving other cores available:
3dprintase batch molecules.json --jobs 4 --ph 7 --ensemble --ball-stick --output out/sugars
```

Each entry accepts exactly one of `cid`, `pubchem_name`, `smiles`, `inchi`, or
`input`, plus an optional `label`. Input file paths are relative to the manifest.
Batch settings apply to every entry. Normalized labels must be unique to prevent
output collisions. The included `tests/data/sugar_batch.json` runs without network
access because it contains explicit stereochemical SMILES.

Batch commands default to one worker per available **logical CPU**, capped at the
number of inputs. Each worker processes one molecule at a time; RDKit, Java, and
native math pools use one computational thread per worker. `--jobs N` controls
concurrency. `--threads N` allows more RDKit/JVM threads for a single molecule;
parallel batches require `--threads 1`. Java may create auxiliary service threads,
and workers are not pinned to particular cores. Each exporter has a 2 GB maximum
Java heap; reduce `--jobs` for large molecules or limited memory.

`pubchem`, `glucan`, `fructan`, `run`, and `auto` also accept `--jobs`. PDB jobs are
scheduled per PDB entry and render its extracted components sequentially.
Failures are recorded per input without stopping other molecules. The parent
writes `run_summary.json`, including worker PIDs and timings; any failed input
makes the batch exit with status 1. `--dry-run` does not create output files.

## Python library

```python
from printase3d import Molecule, generate_molecule, generate_batch

if __name__ == "__main__":
    result = generate_batch(
        [Molecule(cid=79025, label="alpha-D-glucose"),
         Molecule(cid=64689, label="beta-D-glucose")],
        jobs=4, ph=7, ensemble=True, ball_stick=True, output="out/sugars",
    )
    print(result.outputs)   # successful model folders, in input order
    print(result.failures)  # per-molecule errors
```

Use the `__main__` guard in executable scripts because workers use the portable
`spawn` process method. `generate_molecule(Molecule(...), **options)` returns one
model folder and raises `PipelineError` on failure. Keyword options correspond
to CLI flags, using underscores, such as `conformer_count=50` and `no_render=True`.

### Module layout

| Module | Responsibility |
|---|---|
| `api.py`, `batch.py` | Public library interface and process workers |
| `arguments.py`, `cli.py` | Shared option definitions and command dispatch |
| `chemistry.py`, `conformers.py` | Protonation, coordinate transfer, ensemble minimization, shape selection |
| `glycans.py` | Sugar builders and glycosidic conformations |
| `structures.py`, `sources.py` | mmCIF/PDB handling, PubChem and RCSB retrieval |
| `rendering.py`, `mesh.py` | Background Jmol export, OBJ repair, materials |
| `pipeline.py`, `metadata.py` | End-to-end processing and provenance |
| `config.py`, `options.py`, `utils.py` | Presets, naming, validation, shared utilities |

Existing `printase3d.cli` function imports remain available for compatibility.

## Any molecule

```bash
3dprintase molecule --pubchem-name cocaine --ball-stick
3dprintase molecule --cid 71485 --label sucralose
3dprintase molecule --smiles 'CC(=O)OC1=CC=CC=C1C(=O)O' --label aspirin
3dprintase molecule --input ./structure.sdf --label my-molecule
```

Useful controls:

- `--small-spheres` selects the original style: 45% van der Waals atom radii and 0.40 Å bond radius.
- `--ball-stick` is the default and uses 35% atom radii and 0.30 Å bond radius.
- `--spacefill` and `--wireframe` override radii without changing atom positions or chemical bond lengths.
- `--monocolor-aliphatics` gives all hydrogen atoms the same charcoal color as carbon while preserving other element colors. `--no-color` takes precedence.
- `--no-hydrogens` removes hydrogens only during rendering: embedding, conformer selection, minimization, and coordinate validation retain explicit hydrogens, as does the saved SDF.
- `--no-hydrogens`, `--no-color`, and `--bond-colors split` change visible parts/materials.
- Double, triple, and aromatic bonds render as one centered cylinder per atom pair by default. Their bond order still controls RDKit geometry; `--multiple-bond-cylinders` opts back into parallel cylinders when a chemical diagram style is specifically wanted.
- `--coordinates regenerate` makes a deterministic ETKDGv3 conformer.
- `--pubchem-geometry` explicitly preserves/maps PubChem heavy-atom coordinates and relaxes only hydrogens. This is the right choice when the deposited conformer itself is the desired model.
- `--ensemble` explicitly includes the PubChem conformer as one candidate, generates a deterministic ETKDG ensemble, minimizes every candidate, and selects the lowest restrained MMFF94s (or UFF fallback) potential energy. `--conformer-count` defaults to 50 and `--conformer-prune-rms` defaults to 0.15 Å.
- With neither geometry flag, `auto` preserves an unchanged supplied conformer but performs the ensemble search after protonation or coordinate regeneration. The sampled force-field minimum is not a proof of the global minimum and is not an aqueous free-energy calculation.
- `--optimize auto|hydrogens|mmff|uff|none` remains a lower-level control. `hydrogens` fixes heavy atoms, explicit `mmff`/`uff` performs one local all-atom minimization, and `none` performs no minimization.
- `--largest-fragment` removes salts/counterions by default; `--keep-all-fragments` preserves disconnected components.
- `--charge N` validates, but never silently changes, formal charge.
- `--ph 7` generates an explicit pH-7 ionization state with Dimorphite-DL, preserves mapped heavy-atom coordinates, and places/relaxes the new hydrogens. The result is kept separately in `ball-stick-ph7`; `--charge N` can enforce an expected net charge.
- When `--ph` changes only protonation, 3DPrintase now maps the complete heavy-atom graph back onto the deposited coordinates before adding hydrogens. It regenerates coordinates only when that full graph mapping is impossible, and records the decision in JSON metadata.

Geometry-policy examples:

```bash
# Reproduce PubChem's conformer, including across a pH-state change.
3dprintase molecule --cid 441477 --ph 7 --pubchem-geometry --ball-stick

# Sample charge-specific conformers; PubChem is included as one candidate.
3dprintase molecule --cid 441477 --ph 7 --ensemble --ball-stick

# Exclude PubChem coordinates and search entirely from regenerated starts.
3dprintase molecule --cid 441477 --ph 7 --coordinates regenerate --ensemble --ball-stick
```

### Flattest and least-flat minimized conformers

```bash
3dprintase molecule --cid 16078 --label THC --ph 7 --most-flat --ball-stick --output out/models
3dprintase molecule --cid 192826 --label maltotriose --ph 7 --most-flat --ball-stick --output out/models
3dprintase molecule --cid 192826 --label maltotriose --ph 7 --least-flat --ball-stick --output out/models
3dprintase batch molecules.json --jobs 4 --ph 7 --most-flat --ball-stick --output out/sugars
```

`--most-flat` and `--least-flat` automatically enable ensemble sampling and
minimization; adding `--ensemble` explicitly is also supported. Both are available
on `molecule` and `batch`, and through the Python API as `most_flat=True` or
`least_flat=True`. They are mutually exclusive and conflict with
`--pubchem-geometry` or an explicit `--optimize` mode.

Every candidate is first energy-minimized. Among **converged** candidates,
`--most-flat` selects the smallest heavy-atom RMS distance to a best-fit plane;
`--least-flat` selects the largest. This rotation-independent metric is measured
in angstroms. Candidates within **0.0001 Å** of the flatness extreme are treated
as tied, and the lowest-energy candidate wins. Hydrogens are excluded unless the
molecule contains only hydrogen; molecules with at most three measured atoms
always fit a plane and are therefore selected by energy.

These flags choose the shape extreme among minimized structures, which can have
higher energy than the ordinary `--ensemble` winner. They do not flatten atoms
onto a plane, force a distorted geometry, or guarantee a global energy minimum.
If no candidate converges, the run fails with a suggestion to increase
`--max-opt-iters`. Increase `--conformer-count` to explore more candidates.

Short folders such as `ball-stick-ph7-most-flat` and `ball-stick-ph7-least-flat`
preserve each selection. Repeated flatness runs automatically create `-2`, `-3`,
and so on, using exclusive directory creation to prevent overwriting existing
outputs even across concurrent runs. JSON metadata records the full candidate table,
selected flatness, convergence, and energy above the sampled minimum; `MODEL.md`
summarizes the selection. Core sampling and shape-selection code lives in
`conformers.py`, with protonation and preparation in `chemistry.py`.
- `--seed`, `--max-opt-iters`, `--scale`, `--timeout`, `--dry-run`, and `--no-render` support automation.

## Glycans and carbohydrate builders

```bash
3dprintase search --type cellulose --limit 20
3dprintase run 3ACI 1RP8 --min-dp 5
3dprintase run 3ACI --min-dp 5 --components BGC --dry-run
3dprintase run 3TS8 --extract dna --min-dp 20 --no-hydrogens
3dprintase glucan --type maltose --dp 3 4 5 6 7 8 9 10
3dprintase glucan --type cellulose --dp 3 4 5 6 7 8 9 10
3dprintase glucan --type galactan --dp 3 4 5 6 7 8 9 10
3dprintase glucan --type mannan --dp 3 4 5 6 7 8 9 10
3dprintase fructan --type inulin --dp 3 4 5 6 7 8 9 10
3dprintase fructan --type levan --dp 3 4 5 6 7 8 9 10
```

Generated oligosaccharides now use `--conformation natural` by default. The builder detects sugar rings and glycosidic bridge oxygens from the molecular graph, applies linkage/residue-specific heavy-atom φ/ψ/ω targets, evaluates a small deterministic conformer set, and performs restrained MMFF94s relaxation. This produces an extended cellulose/mannan centerline and canonical helical models for amylose-like maltose oligomers, β-(1→4)-galactan, and inulin; levan uses its documented low-energy β-(2→6) repeat. These are representative single conformers for visualization, not complete solution ensembles.

Use `--conformation generic` to reproduce the old unconstrained ETKDG/MMFF behavior. It writes to a separate `*-generic-conformer` style folder. Every model's JSON and `MODEL.md` record the profile, source, target and measured torsions, selected seed, centerline metrics, and minimization method.

PDB extraction preserves model-1 Cartesian coordinates and blank/A alternate locations. The complete downloaded mmCIF remains alongside extracted components. `--extract glycan` is the default; `--extract dna` keeps all deposited DNA strands together as a duplex while removing protein, waters, ions, and non-target ligands. X-ray structures normally contain no deposited hydrogen coordinates, so `--no-hydrogens` accurately describes those DNA exports.

## Output

The default root is `out/misc`, relative to the current working directory.

```text
cocaine_446220/
└── ball-stick/
    ├── cocaine_446220.obj
    ├── cocaine_446220.mtl
    ├── cocaine_446220.sdf
    ├── cocaine_446220_jmol_raw.obj
    ├── cocaine_446220_metadata.json
    └── MODEL.md
```

Each molecule has one short subfolder per style. Presets produce `small-spheres` or `ball-stick`; only non-default controls are appended, for example `ball-stick-no-h` or `ball-stick-wf0p4`. `MODEL.md` and JSON hold the full settings and reproduction command.

A pH-specific state receives its own folder, for example `ball-stick-ph7`, so it never overwrites the source/default protonation state. Protonation software predicts a discrete structural state; it does not imply that only one microstate exists in solution.

Colored exports use semantic MTL names such as `atom_O_oxygen`, `atom_P_phosphorus`, and `bond_uniform`. This prevents slicers from confusing opaque Jmol color identifiers and makes element parts easier to identify and recolor.

For Bambu Studio, the printable OBJ also writes one top-level `o part_atom_<element>` object per chemical element and `o part_bond_uniform` for bonds. Carbon and bonds use the exact same CPK-style charcoal by default, while hydrogen remains white; the bond part stays separately named so it can still be reassigned manually. In Bambu's Multi-Color dialog, choose at least the number of chemical-element colors, click **Apply**, and disable **Auto-merge same filament** only if distinct element colors are being mapped to the same loaded spool.

The semantic palette and part mapping cover H, C, N, O, S, F, Cl, Br, I, and P. Two-letter halogens are recognized in both title-case (`Cl`, `Br`) and uppercase (`CL`, `BR`) atom-name conventions. Every non-carbon element has a distinct default color; uniform bonds intentionally share carbon's color but remain a separately named part.

The primary OBJ contains closed atom spheres plus capped 36-sided covalent bond cylinders with zero open edges. Hydrogen bonds and Jmol struts are disabled. The raw Jmol OBJ is retained for provenance, and the MTL stays adjacent for element materials.

## Jmol attribution

The adapted template is `src/printase3d/jmol/3DP-Jmol.pipeline.spt`. See [third-party notices](../THIRD_PARTY_NOTICES.md) for the original source checksum, attribution, and upstream license.

## Validate

```bash
3dprintase doctor
python -m unittest discover -s tests -v
```

### Fused aromatic ring geometry

Ensemble and most/least-flat minimization preserve near-planarity of fully aromatic fused components containing 7–14 atoms using torsion restraints (±0.25°, force constant 100). The selected component must have plane RMS below 0.03 Å. Larger potentially curved aromatic systems and freely rotating biaryls are excluded from this automatic rule. MMFF typing runs on private copies and chemical aromaticity is restored before serialization. Energies include the restraint potential and are not unconstrained quantum or aqueous minima.

Oxidized FMN/FAD cores receive these restraints. Reduced FMNH2/FADH2 retain their distinct bond orders and nonaromatic bridge nitrogens, allowing butterfly bending without imposing a universal bend angle. “Most flat” measures the entire molecule and must not change redox state or flatten a reduced ring by projection. Geometry metadata records flavin state and core RMS. PubChem source heavy coordinates remain unchanged.

### Extended dinucleotides

Use `--extended` with `molecule`, `batch`, or `extended=True` in the Python API to generate an elongated dinucleotide. This identifies two nitrogen-containing terminal ring systems across a diphosphate bridge, applies a lower-bound anchor-distance restraint during the ensemble minimization, and selects the converged candidate with the greatest ring-centroid separation (energy breaks ties). Aromatic geometry protections remain active. Non-dinucleotides are rejected. The option is mutually exclusive with most/least-flat and source-geometry preservation.

Example: `3dprintase molecule --cid 5892 --extended --ph 7 --ball-stick --output out/cofactors/organic`. Extended outputs receive unique `ball-stick-ph7-extended` variant directories. These are restraint-minimized presentation conformers, not claims of an unconstrained global energy minimum.

### Carbonyl and amide geometry

Molecular minimization now uses the static MMFF94s parameter set where supported. The former MMFF94 choice could leave ordinary amide nitrogens pyramidal, including their explicit hydrogens. Local restraints maintain trigonal carbonyl carbon geometry and the conjugated planes of ordinary neutral amides (including ureas and carbamates). The acceptance threshold is 0.03 Å displacement of each center from the plane of its three bonded neighbors. This does not flatten entire molecules, tetrahedral amines, or ordinary sp3 carbons. Small-ring and bridgehead amides are excluded because actual strain can break normal amide planarity.

The same protections apply to ensemble, most/least-flat, extended, and hydrogen-only minimization. PubChem mode continues to fix all source heavy atoms; if those atoms prevent a valid local geometry, it fails explicitly instead of changing the source coordinates silently. Energies include the documented geometry restraints.

## RCSB CCD molecule conformers

Use `--ccd CODE` as a molecule source through the CLI, batch manifests, or `Molecule(ccd="ATP")`. By default it preserves CCD ideal heavy coordinates, assigns the requested pH state and minimizes hydrogens. `--ccd-coordinates model` selects the representative model coordinates instead; neither source is silently substituted if unavailable. Source CIF and coordinate provenance are saved.

```sh
./3dprintase-local molecule --ccd ATP --rcsb-geometry --ph 7 --ball-stick
./3dprintase-local molecule --ccd ATP --rcsb-geometry --ccd-coordinates model --ph 7 --ball-stick
./3dprintase-local molecule --ccd ATP --ensemble --ph 7 --ball-stick
./3dprintase-local molecule --ccd ATP --extended --ph 7 --ball-stick
```

`--extended` now also supports nucleosides and individual nucleotides using terminal heavy-atom separation. Rigid bases are ranked by span without extension forces. The earlier dinucleotide terminal-ring metric is retained. Aromatic planarity restraints now include individual 5/6-membered rings, and source-preserving variants reject distorted source rings instead of quietly changing their heavy coordinates. Colored SDF exports validate each atom's OBJ material against its chemical element; oxygen is red and phosphorus orange.

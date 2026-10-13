# 3DPrintase

Turn molecular structures into colored, printable ball-and-stick models from the command line. 3DPrintase combines RDKit geometry and conformer selection, Jmol rendering, and mesh repair to produce OBJ/MTL files with closed atom spheres and capped bond cylinders.

Inputs include PubChem names/CIDs, SMILES, InChI, local structures, RCSB chemical components, and PDB glycans or DNA. Generated linear peptides, glucan chains, and fructan chains are also supported.

## Quick start

Requirements: Git and Conda (or Miniforge). The environment installs Python 3.12, Java 17, and the Python dependencies. Jmol is downloaded separately.

```sh
# Replace OWNER with the GitHub account hosting your copy.
git clone https://github.com/OWNER/3DPrintase.git
cd 3DPrintase
conda env create -f environment.yml
conda activate 3D-glycan
3dprintase --help
```

Download the **full desktop Jmol** binary archive from the [official downloads page](https://jmol.sourceforge.net/download/), extract it, and copy `Jmol.jar` into `tools/jmol/`. The renderer was validated with Jmol 16.3.47; use the full `Jmol.jar`, not `JmolData.jar` or JSmol. Jmol is not bundled with this repository.

```sh
cp /path/to/extracted/Jmol.jar tools/jmol/Jmol.jar
3dprintase doctor

# First printable model: works without a network connection after installation.
3dprintase molecule --smiles 'CCO' --label ethanol --ball-stick --output out/demo
```

Open the resulting `.obj` in your slicer, keeping its `.mtl` beside it for colors. Each model folder also includes the `.sdf`, `MODEL.md`, provenance/geometry JSON, and raw Jmol export. Output defaults to `out/misc`, relative to the current directory.

**Try the chemistry without Jmol or Java:**

```sh
3dprintase molecule --smiles 'CCO' --label ethanol --no-render --output out/chemistry
```

For pip/venv installation, Windows, custom Jmol paths, and headless Linux, see [installation and troubleshooting](docs/INSTALLATION.md).

## Common commands

```sh
# Resolve a PubChem molecule, select its pH-7 state, and minimize an ensemble.
3dprintase molecule --pubchem-name caffeine --ph 7 --ensemble --ball-stick

# Select shape extremes among converged minimized candidates.
3dprintase molecule --cid 192826 --most-flat --ball-stick
3dprintase molecule --cid 192826 --least-flat --ball-stick

# Batch processing from an included offline example.
3dprintase batch examples/sugars.json --jobs 2 --ensemble --ball-stick --output out/sugars

# Generate carbohydrate chains or extract a deposited glycan.
3dprintase glucan --type cellulose --dp 3 4 5 --output out/cellulose
3dprintase run 3ACI --min-dp 5 --output out/pdb

# Use RCSB Chemical Component Dictionary coordinates.
3dprintase molecule --ccd ATP --rcsb-geometry --ph 7 --ball-stick

# Build a short linear peptide from its N-to-C one-letter sequence.
3dprintase peptide RGW
3dprintase peptide RGW --ph 7 --extended --no-hydrogens
```

PubChem/RCSB commands need network access. SMILES, local files, peptide builders, and carbohydrate builders can run offline after installation. Use `3dprintase COMMAND --help` for all options.

`peptide` accepts the 20 standard amino-acid one-letter codes (case-insensitive), builds L residues with uncapped amino/carboxyl termini, and samples and minimizes 3D conformers. `RGW` means Arg–Gly–Trp. Use `--ph 7` for pH-dependent protonation; without `--ph`, the builder's neutral formal-charge representation is retained. Linear describes the chain connectivity, not a guaranteed straight conformation; `--extended` selects an extended conformer using the existing distance-restraint search. Cyclic peptides, terminal caps, D residues, and modified residues are not supported by this subcommand.

The ball-and-stick preset shows explicit hydrogens, 35% van der Waals atom radii, and 0.30 Å bond radii. Conformer energies are sampled force-field estimates, not guarantees of a global minimum or a biological bound state. Printed meshes contain overlapping closed components rather than a single boolean-unioned solid.

## Documentation

- [Installation and troubleshooting](docs/INSTALLATION.md)
- [CLI and Python API reference](docs/CLI_GUIDE.md)
- [Pipeline and geometry](docs/PIPELINE_OVERVIEW.md)
- [Bambu Studio printing](docs/BAMBU_STUDIO_PRINTING.md)
- [Development, tests, and packaging](CONTRIBUTING.md)
- [Third-party attribution and license status](THIRD_PARTY_NOTICES.md)

## Repository layout

| Directory | Contents |
|---|---|
| `src/printase3d/` | Installable Python package and Jmol template |
| `examples/` | Small, reproducible CLI input manifests |
| `tests/` | Offline unit/regression tests and fixtures |
| `docs/` | Setup, usage, geometry, and slicer guides |
| `scripts/` | Optional collection-generation and maintenance workflows |
| `tools/` | Optional Blender preview utility and local Jmol location |
| `out/` | Generated models and reports; ignored by Git |

Generated model collections and downloaded binaries are not required to clone or install the project and are not tracked. Some maintenance scripts depend on local collections; see [scripts/README.md](scripts/README.md).

## License status

A license for the original 3DPrintase code has not yet been selected. The adapted 3DP-Jmol template retains its upstream MIT license and attribution; see [third-party notices](THIRD_PARTY_NOTICES.md).

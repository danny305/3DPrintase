# How 3DPrintase works

3DPrintase is a reproducible command-line bridge from chemical structure data to a watertight OBJ that can be imported into Bambu Studio. It supports ordinary small molecules, generated oligosaccharides, and experimentally observed carbohydrate coordinates extracted from protein–glycan PDB entries.

## Pipeline at a glance

1. **Resolve the molecule.** `molecule` accepts a PubChem CID or name, SMILES, InChI, or a local SDF/MOL/MOL2/PDB file. PubChem names are resolved to stable CIDs. The PDB commands download mmCIF coordinates; the glucan/fructan commands construct a defined sequence.
2. **Choose the geometry policy.** `--pubchem-geometry` preserves the complete deposited heavy-atom conformer. `--ensemble` includes that conformer alongside deterministic ETKDG candidates, minimizes each, and selects the lowest sampled MMFF94s/UFF potential energy. With neither flag, protonation/regeneration triggers ensemble sampling while an unchanged deposited conformer is preserved. Energy ranking is a force-field potential-energy approximation, not an aqueous free-energy calculation.
3. **Prepare chemistry.** RDKit validates connectivity and stereochemistry, removes salts by retaining the largest covalent component by default, adds explicit hydrogens, and minimizes geometry with MMFF94s or UFF. When `--ph` is supplied, Dimorphite-DL constructs the requested ionization state and 3DPrintase maps the unchanged heavy-atom graph onto source coordinates before applying the selected geometry policy.
4. **Render with Jmol.** A pipeline-owned copy of the original Jmol script applies element radii, bond cylinders, materials, visibility, and global scale. Hydrogen bonds and structural struts are disabled. Double, triple, and aromatic bond orders remain in the chemical structure and constrain optimized geometry, but the print-oriented default represents every bonded atom pair with one centered cylinder instead of parallel cylinders.
5. **Repair for printing.** Jmol's visually closed but topologically open bond shells are replaced by 36-sided, capped cylinders. Atom spheres remain separate parts, and opaque Jmol color codes are replaced by semantic materials such as `atom_O_oxygen`, `atom_P_phosphorus`, and `bond_uniform`. The primary OBJ therefore has zero open edges, while the untouched `*_jmol_raw.obj` and MTL remain for provenance.
6. **Record everything.** Every style folder includes JSON metadata and `MODEL.md` with identity, source, coordinate method, exact render settings, files, topology, and the reproduction command.

## Output organization

The default output root, relative to the current directory, is:

`out/misc`

Molecule roots use `<lay-name>_<stable-id>`. Render variants have short names such as `small-spheres`, `ball-stick`, or `ball-stick-no-h`. Only non-preset overrides are appended, such as `ball-stick-wf0p4`.

```text
out/misc/cocaine_446220/
└── ball-stick/
    ├── cocaine_446220.obj
    ├── cocaine_446220.mtl
    ├── cocaine_446220.sdf
    ├── cocaine_446220_jmol_raw.obj
    ├── cocaine_446220_metadata.json
    └── MODEL.md
```

## Coordinate accuracy versus printable appearance

`--spacefill` changes displayed atom radii, `--wireframe` changes bond-cylinder radius, and `--scale` changes the whole model uniformly. None changes chemical bond lengths. The optional `--small-spheres` preset is 45% van der Waals radii with 0.40 Å bonds. `--ball-stick` is 35% with 0.30 Å bonds. Explicit hydrogens are visible by default.

Use `--multiple-bond-cylinders` only when parallel double/triple-bond cylinders are desired for diagrammatic emphasis. It is off by default because a single centered cylinder produces cleaner intersections, stronger prints, and a more consistent ball-and-stick appearance.

PubChem conformers and supplied 3D files are used as minimization starting points. With the default `--optimize auto`, MMFF94s is used when completely parameterized and UFF otherwise. Generated conformations are useful physical models, not experimental claims about a molecule's conformation in a crystal, solution, receptor, or living system.

## Main calls in the code

- `main()` parses the command and dispatches it.
- `process_molecule()` resolves general inputs and creates the output bundle.
- `prepare_molecule()` standardizes fragments, adds hydrogens, embeds missing 3D coordinates, and minimizes geometry.
- `process_pdb()` downloads and extracts experimentally observed glycan components without moving their atoms.
- `build_glucan()` and `build_fructan()` define carbohydrate connectivity through GlyLES and prepare coordinates with RDKit. Their default `natural` mode identifies pyranose/furanose rings and glycosidic bridges from connectivity, applies a documented linkage-specific φ/ψ/ω profile, selects a deterministic candidate with the expected extended or helical global shape, and relaxes bond lengths and unconstrained degrees of freedom with restrained MMFF94s. `--conformation generic` retains the legacy unconstrained embedding as an explicit fallback.
- `render_jmol()` fills the copied Jmol script template and invokes Jmol through Java.
- `repair_jmol_obj()` replaces open Jmol bond shells with capped printable cylinders.
- `write_model_metadata()` writes the JSON record and readable `MODEL.md`.

## Installation and sharing

Jmol is Java software, not a Python wheel, so `pyproject.toml` declares the Python libraries while `environment.yml` declares Java 17 and installs 3DPrintase itself. The Jmol JAR is located automatically in this checkout or can be supplied with `--jmol-jar` or `THREEDPRINTASE_JMOL_JAR`.

Follow [installation and troubleshooting](INSTALLATION.md) for the complete setup. After installation, `3dprintase` works from any directory while its Python environment is active.

## Examples

```bash
3dprintase molecule --pubchem-name sucralose --ball-stick
3dprintase molecule --cid 71485 --label sucralose --wireframe 0.4
3dprintase molecule --smiles 'CC(=O)OC1=CC=CC=C1C(=O)O' --label aspirin
3dprintase molecule --input ./molecule.sdf --label my-molecule --coordinates preserve
3dprintase molecule --input ./molecule.mol2 --label my-molecule --coordinates regenerate --optimize uff
3dprintase glucan --type maltose --dp 3 4 5 6 7 8 9 10
3dprintase glucan --type cellulose --dp 3 4 5 6 7 8 9 10
3dprintase glucan --type galactan --dp 3 4 5 6 7 8 9 10
3dprintase glucan --type mannan --dp 3 4 5 6 7 8 9 10
3dprintase fructan --type inulin --dp 3 4 5 6 7 8 9 10
3dprintase fructan --type levan --dp 3 4 5 6 7 8 9 10
3dprintase run 3ACI 1RP8 --min-dp 5
```

Run `3dprintase molecule --help` for all general-molecule controls and `3dprintase doctor` to check a shared installation.

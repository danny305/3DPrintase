# Installation and troubleshooting

## Conda (recommended)

From the repository root:

```sh
conda env create -f environment.yml
conda activate 3D-glycan
3dprintase --version
```

The environment uses Python 3.12 and OpenJDK 17 and installs this checkout in editable mode. To update an existing environment, use `conda env update -f environment.yml`. Activate it in each terminal before running the CLI. On Windows, use a Conda-enabled terminal.

## pip and venv

Use Python 3.12 for the tested setup. The package declares Python 3.10 or newer; wheels for its native dependencies must also be available for your Python/platform combination. Java 17 must be installed separately for rendering.

```sh
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
3dprintase --help
```

On Windows, create the environment with `py -3.12 -m venv .venv` and activate it with `.venv\Scripts\Activate.ps1` in PowerShell. Conda is an alternative if RDKit or other native wheels are unavailable.

`python -m printase3d` is equivalent to `3dprintase`. The optional `./3dprintase-local` launcher uses the active `python3`; it does not activate an environment or install dependencies.

## Install Jmol

Download a full desktop binary archive from [Jmol's official download page](https://jmol.sourceforge.net/download/) and extract `Jmol.jar`. This project has been render-tested with version **16.3.47**. Jmol is an external runtime, not a Python dependency, and is not redistributed here.

For a checkout, place the file at `tools/jmol/Jmol.jar`. For an installed wheel or a shared installation, set an absolute path:

```sh
export THREEDPRINTASE_JMOL_JAR=/absolute/path/to/Jmol.jar
3dprintase doctor
```

PowerShell:

```powershell
$env:THREEDPRINTASE_JMOL_JAR = 'C:\tools\Jmol.jar'
3dprintase doctor
```

Lookup order:

1. The command's explicit `--jmol-jar` argument.
2. `THREEDPRINTASE_JMOL_JAR` (an invalid explicit path is reported, not silently replaced).
3. `tools/jmol/Jmol.jar`, then `Jmol.jar`, in the source checkout.
4. `$XDG_CACHE_HOME/3dprintase/Jmol.jar`, or `~/.cache/3dprintase/Jmol.jar` when `XDG_CACHE_HOME` is unset.

Run `3dprintase doctor` to inspect the resolved path, Java, Python dependencies, and output access. `doctor` checks dependency availability, not a complete render. Confirm rendering with:

```sh
3dprintase molecule --smiles 'CCO' --label ethanol --ball-stick --output out/demo
```

## Rendering on a server

OBJ export needs the **full Jmol application** and AWT support. `JmolData.jar` and Java's `-Djava.awt.headless=true` mode do not support this export path. On Linux without a display, install Xvfb and run the renderer inside it:

```sh
xvfb-run -a 3dprintase molecule --smiles 'CCO' --label ethanol --ball-stick
```

On macOS the pipeline uses `-Dapple.awt.UIElement=true` to avoid normal Java Dock activation. The locally verified platform is macOS with Java 17; CI checks Python/chemistry portability on Linux and macOS without rendering.

## Common problems

| Symptom | Action |
|---|---|
| `3dprintase: command not found` | Activate the environment, run `python -m pip install -e .`, or try `python -m printase3d --help`. |
| Missing RDKit, glyles, or Dimorphite-DL | Install the entire project, not individual scripts; use `3dprintase doctor`. |
| `Jmol.jar not found` | Set the environment variable or pass `--jmol-jar` to the subcommand. |
| Java/AWT error or missing OBJ | Check `java -version`, use full `Jmol.jar`, inspect the model's `_jmol.log`, and use Xvfb on displayless Linux. |
| PubChem/RCSB request fails | Check connectivity, retry, or use an explicit stereochemical SMILES/local SDF. |
| Memory pressure during a batch | Lower `--jobs`; each Java exporter can use up to 2 GB of heap. |
| Only chemistry is needed | Add `--no-render`; no Java/Jmol execution is required. |

Outputs are relative to your current directory unless `--output` is absolute. Keep OBJ and MTL files together. Never use an output directory containing valuable unrelated files as a scratch directory for collection-repair scripts.

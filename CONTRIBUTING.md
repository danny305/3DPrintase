# Development

Create the environment described in [installation](docs/INSTALLATION.md), then install packaging tools:

```sh
python -m pip install -e '.[dev]'
python -m unittest discover -s tests -v
```

The test suite uses local fixtures and mocks for remote services. It does not need Java, Jmol, or a network connection after dependencies have been installed. A separate rendering smoke test is documented in the installation guide.

## Build and check a distribution

```sh
python -m build
python -m twine check dist/*
```

The wheel must contain `printase3d/jmol/3DP-Jmol.pipeline.spt`. It does not contain Jmol binaries, generated model collections, or local environments. Install the wheel into a fresh environment, change out of the checkout, and verify:

```sh
3dprintase --help
3dprintase molecule --smiles 'CCO' --label ethanol --no-render --output smoke-output
```

GitHub Actions tests editable installs, builds the sdist/wheel, and exercises the installed wheel from outside the checkout. Java rendering is checked separately because it needs the full Jmol binary and a display/AWT runtime.

## Working conventions

- Use a feature/review branch; keep generated models in `out/`.
- Keep dependency requirements in `pyproject.toml`; `requirements.txt` delegates to it. `environment.yml` also installs Java and the recommended Python version.
- Preserve chemical identity, stereochemistry, and provenance when changing conformer code.
- Run the tests relevant to your changes and the full suite before preparing a release.
- Keep source/data fixtures small. Do not commit credentials, local absolute paths, environments, JAR files, or generated mesh collections.

## Publishing this initial repository

The local review branch is `review/github-ready`. No GitHub repository or remote has been created. After review, create an empty GitHub repository and configure its URL as `origin`. Update the `OWNER` placeholder in the README clone example. Push the review branch when ready; choosing a default branch and merging are separate review decisions.

The original-code license remains undecided at the owner's request. Preserve `THIRD_PARTY_NOTICES.md` and the upstream license bundled with the Jmol template.

# Third-party notices

## 3DPrintase

No license has yet been selected for the original 3DPrintase code. This notice does not grant a license for that code.

## Adapted 3DP-Jmol script

`src/printase3d/jmol/3DP-Jmol.pipeline.spt` is adapted from [3DP-Jmol](https://github.com/mariusmihasan/3DP-Jmol), by Marius Mihășan and Angel Herráez. The upstream MIT notice is reproduced in [licenses/3DP-Jmol-MIT.txt](licenses/3DP-Jmol-MIT.txt), with a packaged copy alongside the template.

Upstream copyright: Copyright (c) 2025 Mihasan Marius and Angel Herráez.

The original local source was `3DP-Jmol.v.alfa1.public`, SHA-256 `abb2079198f0044afa260d7fb9850566d4a1df3f7e11e40b414de18a12b1508c`. Pipeline adaptations add parameter substitution, rendering controls, and export integration. The original source file is not distributed here. Existing author credits in the adapted script are retained.

The upstream license text was obtained from https://raw.githubusercontent.com/mariusmihasan/3DP-Jmol/main/LICENSE.md on 2026-10-10.

## External dependencies and data

Jmol is installed separately under its own license; see the [official distribution](https://jmol.sourceforge.net/download/). Python dependencies retain their respective licenses. The repository does not bundle their distributions.

PubChem and RCSB structures and identifiers retain their source provenance. `tests/fixtures/ATP.cif` is an RCSB Chemical Component Dictionary fixture; source: https://files.rcsb.org/ligands/download/ATP.cif. Generated models are excluded from the source repository.

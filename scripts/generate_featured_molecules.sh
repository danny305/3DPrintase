#!/usr/bin/env bash
set -euo pipefail

# PubChem resolves each unambiguous name to a CID; --label controls the short,
# layperson-friendly output name. The default output root is used.
3dprintase molecule --pubchem-name "dextroamphetamine" --label dextroamphetamine --ball-stick
3dprintase molecule --pubchem-name "hydrocodone" --label hydrocodone --ball-stick
3dprintase molecule --pubchem-name "alprazolam" --label alprazolam --ball-stick
3dprintase molecule --pubchem-name "lisdexamfetamine" --label lisdexamfetamine-vyvanse --ball-stick
3dprintase molecule --pubchem-name "delta-9-tetrahydrocannabinol" --label thc --ball-stick
3dprintase molecule --pubchem-name "lysergic acid diethylamide" --label lsd --ball-stick
3dprintase molecule --pubchem-name "psilocybin" --label psilocybin --ball-stick
3dprintase molecule --pubchem-name "cocaine" --label cocaine --ball-stick

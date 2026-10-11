#!/usr/bin/env bash
set -euo pipefail

# Exact-pH Dimorphite-DL states, regenerated and minimized after ionization.
# Expected net charges provide a fail-fast chemistry check for this curated set.
3dprintase molecule --pubchem-name "dextroamphetamine" --label dextroamphetamine --ball-stick --ph 7 --charge 1
3dprintase molecule --pubchem-name "hydrocodone" --label hydrocodone --ball-stick --ph 7 --charge 1
3dprintase molecule --pubchem-name "alprazolam" --label alprazolam --ball-stick --ph 7 --charge 0
3dprintase molecule --pubchem-name "lisdexamfetamine" --label lisdexamfetamine-vyvanse --ball-stick --ph 7 --charge 2
3dprintase molecule --pubchem-name "delta-9-tetrahydrocannabinol" --label thc --ball-stick --ph 7 --charge 0
3dprintase molecule --pubchem-name "lysergic acid diethylamide" --label lsd --ball-stick --ph 7 --charge 1
3dprintase molecule --pubchem-name "psilocybin" --label psilocybin --ball-stick --ph 7 --charge -1
3dprintase molecule --pubchem-name "cocaine" --label cocaine --ball-stick --ph 7 --charge 1
3dprintase molecule --pubchem-name "esketamine" --label s-ketamine-esketamine --ball-stick --ph 7 --charge 1
3dprintase molecule --pubchem-name "arketamine" --label r-ketamine-arketamine --ball-stick --ph 7 --charge 1

"""Shared command and library option definitions."""
from pathlib import Path
import argparse
from .config import DEFAULT_OUTPUT, FRUCTAN_TYPES, GLUCAN_TYPES, PRESETS, VERSION
from .options import ph_value, spacefill_value, wireframe_value
from .rendering import resolve_jmol_jar


def positive_int(value):
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number



def add_common(p, *, pdb_options=False):
    base = Path(__file__).resolve().parent
    p.add_argument("--jobs", type=positive_int, default=None,
                   help="parallel molecule workers (default: available logical CPUs)")
    p.add_argument("--threads", type=positive_int, default=1,
                   help="RDKit/JVM threads per molecule (default: 1; parallel batches require 1)")
    p.add_argument("--output", default=DEFAULT_OUTPUT, help=f"output root (default: {DEFAULT_OUTPUT})")
    if pdb_options:
        p.add_argument("--extract", choices=["glycan", "dna"], default="glycan",
                       help="molecular component to extract from each PDB entry (default: glycan)")
        p.add_argument("--min-dp", type=int, default=5)
        p.add_argument("--components", help="comma-separated CCD component filter, e.g. GLC,BGC")
        p.add_argument("--force", action="store_true", help="redownload cached mmCIF files")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--no-render", action="store_true", help="extract and write metadata without Jmol")
    p.add_argument("--jmol-jar", default=resolve_jmol_jar(), help="Jmol.jar path (or set THREEDPRINTASE_JMOL_JAR)")
    p.add_argument("--script", default=str(base / "jmol" / "3DP-Jmol.pipeline.spt"))
    p.add_argument("--style", choices=["ballandstick", "spacefill"], default="ballandstick")
    p.add_argument("--scale", type=float, default=0.3)
    visual = p.add_mutually_exclusive_group()
    visual.add_argument("--small-spheres", dest="render_preset", action="store_const", const="small-spheres", help="use the original small-spheres preset (spacefill 45, wireframe 0.4 Å)")
    visual.add_argument("--ball-stick", dest="render_preset", action="store_const", const="ball-stick", help="use the default ball-and-stick preset (spacefill 35, wireframe 0.3 Å)")
    p.set_defaults(render_preset="ball-stick")
    p.add_argument("--spacefill", type=spacefill_value, default=None, metavar="PERCENT", help="override the selected preset's atom radius percentage")
    p.add_argument("--wireframe", type=wireframe_value, default=None, metavar="ANGSTROMS", help="override the selected preset's covalent bond radius in Å")
    p.add_argument("--bond-colors", choices=["uniform", "split"], default="uniform", help="uniform makes each bond one print-safe cylinder (default); split uses two element-colored half-cylinders")
    p.add_argument("--multiple-bond-cylinders", action="store_true",
                   help="opt in to parallel cylinders for double/triple bonds; default uses one centered cylinder per bonded atom pair")
    p.add_argument("--hydrogens", action=argparse.BooleanOptionalAction, default=True, help="show explicit hydrogens (default); use --no-hydrogens to remove them only at rendering, after all coordinate optimization")
    p.add_argument("--monocolor-aliphatics", action="store_true",
                   help="color all hydrogen atoms the same charcoal as carbon")
    p.add_argument("--color", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--timeout", type=int, default=600)



def add_molecule_options(m, *, label=True):
    if label:
        m.add_argument("--label", help="readable molecule name used in output paths")
    m.add_argument("--seed", type=int, default=61453, help="deterministic conformer seed")
    geometry = m.add_mutually_exclusive_group()
    geometry.add_argument("--pubchem-geometry", dest="geometry_mode", action="store_const", const="pubchem",
                          help="preserve/map PubChem heavy-atom coordinates and relax only hydrogens")
    geometry.add_argument("--ensemble", dest="geometry_mode", action="store_const", const="ensemble",
                          help="rank PubChem plus sampled charge-specific conformers and select the lowest force-field energy")
    geometry.add_argument("--rcsb-geometry", dest="geometry_mode", action="store_const", const="rcsb", help="preserve RCSB CCD heavy-atom coordinates and relax hydrogens")
    m.add_argument("--ccd-coordinates", choices=["ideal", "model"], default="ideal", help="CCD coordinate set (default: ideal; model is a representative deposited geometry)")
    m.set_defaults(geometry_mode="auto")
    flatness = m.add_mutually_exclusive_group()
    flatness.add_argument("--extended", action="store_true",
                          help="minimize with terminal separation restraints; select the most extended conformer")
    flatness.add_argument("--most-flat", action="store_true",
                          help="minimize an ensemble, then select the flattest converged conformer (energy breaks shape ties)")
    flatness.add_argument("--least-flat", action="store_true",
                          help="minimize an ensemble, then select the least planar converged conformer (energy breaks shape ties)")
    m.add_argument("--coordinates", choices=["preserve", "regenerate"], default="preserve",
                   help="use supplied 3D coordinates as the minimization start (default), or regenerate")
    m.add_argument("--optimize", choices=["auto", "hydrogens", "mmff", "uff", "none"], default="auto",
                   help="low-level minimizer; geometry flags normally select this automatically")
    m.add_argument("--conformer-count", type=int, default=50, metavar="N",
                   help="ETKDG candidates requested by auto/ensemble global search (default: 50)")
    m.add_argument("--conformer-prune-rms", type=float, default=0.15, metavar="ANGSTROMS",
                   help="discard near-duplicate ETKDG candidates below this RMSD (default: 0.15 Å)")
    m.add_argument("--max-opt-iters", type=int, default=5000, metavar="N")
    fragments = m.add_mutually_exclusive_group()
    fragments.add_argument("--largest-fragment", dest="keep_all_fragments", action="store_false",
                           help="discard salts/counterions and retain the largest covalent component (default)")
    fragments.add_argument("--keep-all-fragments", action="store_true", help="retain salts and disconnected components")
    m.set_defaults(keep_all_fragments=False)
    m.add_argument("--charge", type=int, help="validate the expected total formal charge")
    m.add_argument("--ph", type=ph_value, metavar="PH", help="generate an ionization/protonation state at this pH")
    m.add_argument("--ph-precision", type=float, default=0.0, metavar="SIGMA",
                   help="Dimorphite-DL pKa uncertainty multiplier (default: 0 for one exact-pH state)")
    m.add_argument("--protonation-variant", type=int, default=1, metavar="N",
                   help="select a deterministic 1-based state if pH processing enumerates alternatives")
    m.add_argument("--max-protonation-variants", type=int, default=32, metavar="N")



def parser():
    p = argparse.ArgumentParser(prog="3dprintase", description=__doc__)
    p.add_argument("--version", action="version", version=VERSION)
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("search", help="list candidate experimental protein–glycan PDB entries")
    s.add_argument("--type", choices=sorted(PRESETS), default="glycan")
    s.add_argument("--limit", type=int, default=20); s.add_argument("--dry-run", action="store_true")
    r = sub.add_parser("run", help="download, extract, and render one or more PDB IDs")
    r.add_argument("pdb_ids", nargs="+"); add_common(r, pdb_options=True)
    a = sub.add_parser("auto", help="search then process candidates")
    a.add_argument("--type", choices=sorted(PRESETS), default="glycan"); a.add_argument("--limit", type=int, default=20)
    add_common(a, pdb_options=True)
    c = sub.add_parser("pubchem", help="fetch or build PubChem CID 3D models and render them")
    c.add_argument("cids", nargs="+"); c.add_argument("--seed", type=int, default=61453)
    add_common(c)
    f = sub.add_parser("fructan", help="build inulin- or levan-type fructooligosaccharides")
    f.add_argument("--type", choices=sorted(FRUCTAN_TYPES), required=True)
    f.add_argument("--dp", nargs="+", type=int, required=True, metavar="N")
    f.add_argument("--seed", type=int, default=61453)
    f.add_argument("--conformation", choices=["natural", "generic"], default="natural",
                   help="use a linkage-specific canonical backbone (default), or legacy unconstrained embedding")
    add_common(f)
    g = sub.add_parser("glucan", help="build canonical Glc-, Gal-, or Man-(1→4) oligosaccharides")
    g.add_argument("--type", choices=sorted(GLUCAN_TYPES), required=True)
    g.add_argument("--dp", nargs="+", type=int, required=True, metavar="N")
    g.add_argument("--seed", type=int, default=61453)
    g.add_argument("--conformation", choices=["natural", "generic"], default="natural",
                   help="use a linkage/residue-specific canonical backbone (default), or legacy unconstrained embedding")
    add_common(g)
    m = sub.add_parser("molecule", help="render any molecule from PubChem, SMILES, InChI, or a structure file")
    source = m.add_mutually_exclusive_group(required=True)
    source.add_argument("--ccd", help="RCSB chemical component identifier; defaults to preserved ideal coordinates")
    source.add_argument("--cid", help="PubChem compound ID")
    source.add_argument("--pubchem-name", metavar="NAME", help="name or synonym resolved by PubChem")
    source.add_argument("--smiles", help="SMILES string")
    source.add_argument("--inchi", help="InChI string")
    source.add_argument("--input", metavar="FILE", help="SDF, MOL, MOL2, or PDB file")
    add_molecule_options(m)
    add_common(m)
    b = sub.add_parser("batch", help="generate molecules from a JSON manifest, one per CPU worker")
    b.add_argument("manifest", help="JSON list of molecule sources and labels")
    add_molecule_options(b, label=False)
    add_common(b)
    d = sub.add_parser("doctor", help="check Java, Jmol, Python dependencies, and output access")
    d.add_argument("--jmol-jar", default=resolve_jmol_jar())
    d.add_argument("--output", default=DEFAULT_OUTPUT)
    return p

"""Options components of the 3DPrintase pipeline."""
from __future__ import annotations

import argparse
import re
from .config import RENDER_PRESETS


def spacefill_value(value):
    value = float(value)
    if not 1.0 <= value <= 100.0:
        raise argparse.ArgumentTypeError("spacefill must be between 1 and 100 percent")
    return value


def wireframe_value(value):
    value = float(value)
    if not 0.05 <= value <= 2.0:
        raise argparse.ArgumentTypeError("wireframe radius must be between 0.05 and 2.0 Å")
    return value


def ph_value(value):
    value = float(value)
    if not 0.0 <= value <= 14.0:
        raise argparse.ArgumentTypeError("pH must be between 0 and 14")
    return value


def resolve_render_preset(args):
    """Apply a named visual preset while allowing numeric CLI overrides."""
    if not hasattr(args, "render_preset"):
        return args
    settings = RENDER_PRESETS[args.render_preset]
    if args.spacefill is None:
        args.spacefill = settings["spacefill"]
    if args.wireframe is None:
        args.wireframe = settings["wireframe"]
    return args


def safe_name(value):
    """Return a readable, portable filename fragment."""
    value = value.strip().lower().replace("α", "alpha").replace("β", "beta").replace("γ", "gamma")
    return re.sub(r"[^a-z0-9]+", "-", value).strip("-") or "compound"


def setting_token(value):
    """Compact filesystem-safe representation of a numeric CLI setting."""
    return f"{float(value):g}".replace("-", "m").replace(".", "p")


def variant_name(args, seed=None):
    """Readable folder name; append only settings that differ from a preset."""
    resolve_render_preset(args)
    name = args.render_preset if args.style == "ballandstick" else "spacefill"
    defaults = RENDER_PRESETS[args.render_preset]
    modifiers = []
    if args.spacefill != defaults["spacefill"]:
        modifiers.append(f"sf{setting_token(args.spacefill)}")
    if args.wireframe != defaults["wireframe"]:
        modifiers.append(f"wf{setting_token(args.wireframe)}")
    if args.scale != 0.3:
        modifiers.append(f"scale{setting_token(args.scale)}")
    if not args.hydrogens:
        modifiers.append("no-h")
    if getattr(args, "monocolor_aliphatics", False):
        modifiers.append("monocolor-aliphatics")
    if not args.color:
        modifiers.append("mono")
    if args.bond_colors != "uniform":
        modifiers.append("split-bonds")
    if args.multiple_bond_cylinders:
        modifiers.append("multiple-bonds")
    if getattr(args, "ph", None) is not None:
        modifiers.append(f"ph{setting_token(args.ph)}")
    if getattr(args, "geometry_mode", "auto") == "rcsb":
        modifiers.append("rcsb-" + getattr(args, "ccd_coordinates", "ideal"))
    elif getattr(args, "geometry_mode", "auto") == "pubchem":
        modifiers.append("pubchem-geometry")
    elif (getattr(args, "geometry_mode", "auto") == "ensemble"
          and not (getattr(args, "most_flat", False) or getattr(args, "least_flat", False) or getattr(args, "extended", False))):
        modifiers.append("ensemble")
    if getattr(args, "extended", False):
        modifiers.append("extended")
    if getattr(args, "most_flat", False):
        modifiers.append("most-flat")
    if getattr(args, "least_flat", False):
        modifiers.append("least-flat")
    if getattr(args, "coordinates", "preserve") == "regenerate":
        modifiers.append("regenerated")
    if getattr(args, "optimize", "auto") != "auto":
        modifiers.append(f"opt-{safe_name(args.optimize)}")
    if getattr(args, "conformer_count", 50) != 50:
        modifiers.append(f"conf{int(args.conformer_count)}")
    if getattr(args, "conformation", "natural") == "generic":
        modifiers.append("generic-conformer")
    if seed is not None and int(seed) != 61453:
        modifiers.append(f"seed{int(seed)}")
    return "-".join([name, *modifiers])

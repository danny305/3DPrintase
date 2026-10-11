"""Output directory allocation without overwriting shape-selection runs."""
from pathlib import Path


def model_directory(root, stem, variant, *, unique=False, dry_run=False):
    """Allocate a variant folder, using -2, -3, ... when unique is requested.

    Exclusive mkdir makes allocation safe between concurrent processes. Dry runs
    only predict a path; the real run reserves its own directory atomically.
    """
    parent = Path(root).resolve() / stem
    if not dry_run:
        parent.mkdir(parents=True, exist_ok=True)
    number = 1
    while True:
        folder = parent / (variant if number == 1 else f"{variant}-{number}")
        if dry_run:
            if not unique or not folder.exists():
                return folder
        else:
            try:
                folder.mkdir(exist_ok=not unique)
                return folder
            except FileExistsError:
                if not unique:
                    raise
        number += 1

"""Public Python API for molecule generation and parallel batches."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .batch import Task, run_tasks
from .utils import PipelineError


@dataclass(frozen=True)
class Molecule:
    """Exactly one structure source, with an optional readable output label."""
    cid: int | str | None = None
    pubchem_name: str | None = None
    smiles: str | None = None
    inchi: str | None = None
    input: str | Path | None = None
    label: str | None = None
    ccd: str | None = None

    def __post_init__(self):
        if sum(getattr(self, key) is not None for key in
               ("cid", "pubchem_name", "smiles", "inchi", "input", "ccd")) != 1:
            raise ValueError("Specify exactly one molecule source")


def _options_argv(options):
    """Map Python keyword options to the same validated options as the CLI."""
    args = []
    for key, value in options.items():
        flag = "--" + key.replace("_", "-")
        if value is None or value is False:
            if value is False and key in ("hydrogens", "color"):
                args.append("--no-" + key)
            continue
        if value is True:
            args.append(flag)
        else:
            args.extend((flag, str(value)))
    return args


def molecule_tasks(molecules, **options):
    tasks = []
    for index, molecule in enumerate(molecules, 1):
        if isinstance(molecule, dict):
            molecule = Molecule(**molecule)
        label = molecule.label or molecule.pubchem_name or (
            f"cid-{molecule.cid}" if molecule.cid is not None else
            Path(molecule.input).stem if molecule.input else
            f"ccd-{molecule.ccd}" if molecule.ccd else f"molecule-{index:04d}")
        argv = ["molecule", "--label", label]
        for key in ("cid", "pubchem_name", "smiles", "inchi", "input", "ccd"):
            value = getattr(molecule, key)
            if value is not None:
                argv.extend(("--" + key.replace("_", "-"), str(value)))
        argv.extend(_options_argv(options))
        tasks.append(Task(label, tuple(argv)))
    return tasks


def generate_molecule(molecule, **options):
    """Generate one molecule and return its output folder; raise on failure."""
    result = run_tasks(molecule_tasks([molecule], **options), jobs=1,
                       dry_run=options.get("dry_run", False))
    if result.failures:
        raise PipelineError(result.failures[0].error)
    return result.outputs[0]


def generate_batch(molecules, *, jobs=None, **options):
    """Generate one molecule per worker (default: available logical CPU count).

    Use an ``if __name__ == '__main__':`` guard in executable Python scripts.
    The returned BatchResult includes successful folders and per-molecule errors.
    """
    from .config import DEFAULT_OUTPUT
    return run_tasks(molecule_tasks(molecules, **options), jobs=jobs,
                     output=options.get("output", DEFAULT_OUTPUT),
                     dry_run=options.get("dry_run", False))

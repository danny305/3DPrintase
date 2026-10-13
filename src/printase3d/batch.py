"""Process-isolated batch execution, with one molecule per CPU worker."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
import multiprocessing
import os
from pathlib import Path
import shlex
import time

from .metadata import REPRODUCTION_COMMAND
from .options import safe_name
from .utils import PipelineError, utc_now


def cpu_count():
    """Number of logical CPUs available to this process."""
    if hasattr(os, "process_cpu_count"):
        return os.process_cpu_count() or 1
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return os.cpu_count() or 1


@dataclass(frozen=True)
class Task:
    """One ordinary CLI invocation, without the executable name."""
    label: str
    argv: tuple[str, ...]


@dataclass
class TaskResult:
    label: str
    outputs: list[str] = field(default_factory=list)
    error: str | None = None
    worker_pid: int = 0
    started: float = 0
    finished: float = 0


@dataclass
class BatchResult:
    results: list[TaskResult]
    workers: int

    @property
    def outputs(self):
        return [Path(path) for result in self.results for path in result.outputs]

    @property
    def failures(self):
        return [result for result in self.results if result.error is not None]

    def to_dict(self):
        return {"created_utc": utc_now(), "workers": self.workers,
                "outputs": [str(path) for path in self.outputs],
                "failures": [asdict(result) for result in self.failures],
                "tasks": [asdict(result) for result in self.results]}


def _worker_init():
    # Spawn imports the lightweight package before any RDKit/NumPy work.
    # Limit native math pools as well as RDKit/JVM's explicit thread controls.
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"


def _execute(task):
    from .arguments import parser
    from .options import resolve_render_preset
    from .pipeline import process_molecule, process_cid, process_pdb, process_glucan, process_fructan

    result = TaskResult(task.label, worker_pid=os.getpid(), started=time.time())
    token = REPRODUCTION_COMMAND.set(shlex.join(["3dprintase", *task.argv]))
    try:
        args = resolve_render_preset(parser().parse_args(task.argv))
        if args.command in ("molecule", "peptide"):
            paths = [process_molecule(args)]
        elif args.command == "pubchem":
            paths = [process_cid(args.cids[0], args)]
        elif args.command == "run":
            paths = process_pdb(args.pdb_ids[0], args)
        elif args.command in ("glucan", "fructan"):
            build = process_glucan if args.command == "glucan" else process_fructan
            paths = [build(args.type, args.dp[0], args)]
        else:
            raise PipelineError(f"Unsupported batch task: {args.command}")
        result.outputs = [str(path) for path in paths]
    except (Exception, SystemExit) as exc:
        result.error = str(exc)
    finally:
        REPRODUCTION_COMMAND.reset(token)
        result.finished = time.time()
    return result


def run_tasks(tasks, *, jobs=None, output=None, dry_run=False):
    """Run independent tasks in spawn workers; collect failures without stopping siblings.

    Results preserve input order. Parallel tasks must request one native thread
    each. The parent alone writes the summary, so workers cannot overwrite it.
    """
    import json
    tasks = list(tasks)
    jobs = cpu_count() if jobs is None else jobs
    if isinstance(jobs, bool) or not isinstance(jobs, int) or jobs < 1:
        raise ValueError("jobs must be a positive integer")
    if len({safe_name(task.label) for task in tasks}) != len(tasks):
        raise ValueError("Duplicate task labels would risk overwriting output")
    workers = min(jobs, len(tasks))
    if workers > 1:
        for task in tasks:
            if "--threads" in task.argv and task.argv[task.argv.index("--threads") + 1] != "1":
                raise ValueError("Parallel batches require --threads 1 to avoid CPU oversubscription")
    if workers <= 1 or dry_run:
        results = [_execute(task) for task in tasks]
    else:
        results = [None] * len(tasks)
        with ProcessPoolExecutor(max_workers=workers,
                                 mp_context=multiprocessing.get_context("spawn"),
                                 initializer=_worker_init) as pool:
            futures = {pool.submit(_execute, task): index for index, task in enumerate(tasks)}
            for future in as_completed(futures):
                index = futures[future]
                try:
                    results[index] = future.result()
                except Exception as exc:
                    results[index] = TaskResult(tasks[index].label, error=str(exc))
    result = BatchResult(results, workers)
    if output is not None and not dry_run:
        root = Path(output).resolve()
        root.mkdir(parents=True, exist_ok=True)
        temporary = root / f".run_summary.{os.getpid()}.json"
        temporary.write_text(json.dumps(result.to_dict(), indent=2) + "\n")
        temporary.replace(root / "run_summary.json")
    return result

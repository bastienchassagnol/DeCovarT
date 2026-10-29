"""Shared scenario loop of the competitor runners.

Scenarios are independent fits (each method fits the bulks of one scenario
jointly or one by one, never across scenarios), so they run in separate
processes. Each finished scenario prints one line, so the progress reaches the
log of the calling R script. Linear-algebra libraries are kept to one thread
per process; the caller sets OPENBLAS_NUM_THREADS and OMP_NUM_THREADS.
"""

import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd


def _run_one(fit_and_write, scenario_id, in_dir, out_dir, options):
    scenario_out = out_dir / scenario_id
    scenario_out.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    try:
        fit_and_write(in_dir / scenario_id, scenario_out, **options)
        return (scenario_id, "ok", "", time.perf_counter() - start)
    except Exception as error:  # record the failure, no fallback fit
        traceback.print_exc()
        message = f"{type(error).__name__}: {error}"
        return (scenario_id, "error", message, time.perf_counter() - start)


def run_batch(label, fit_and_write, in_dir: Path, out_dir: Path, workers=1, **options):
    """Fit every scenario of `in_dir/manifest.csv`; write `out_dir/status.csv`."""
    ids = pd.read_csv(in_dir / "manifest.csv")["scenario_id"].astype(str).tolist()
    workers = max(1, min(int(workers), len(ids)))
    print(f"{label}: {len(ids)} scenarios on {workers} process(es).", flush=True)
    status = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs = [
            pool.submit(_run_one, fit_and_write, sid, in_dir, out_dir, options)
            for sid in ids
        ]
        for done, job in enumerate(as_completed(jobs), start=1):
            row = job.result()
            status[row[0]] = row
            detail = row[1] if row[1] == "ok" else f"error ({row[2]})"
            print(
                f"{label} [{done}/{len(ids)}] {row[0]}: {detail}, {row[3]:.1f} s",
                flush=True,
            )
    pd.DataFrame(
        [status[sid] for sid in ids],
        columns=["scenario_id", "status", "message", "elapsed_sec"],
    ).to_csv(out_dir / "status.csv", index=False)

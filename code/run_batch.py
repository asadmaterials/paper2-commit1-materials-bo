"""
run_batch.py — launch the registered Paper 2 runs (pilot or confirmatory), v1.1
===============================================================================
Registration: https://doi.org/10.5281/zenodo.23118101 (PREREGISTRATION.md v1.1, Section 9)

  pilot         primary universe; beta = 0 only (no prior); Y1 and Y2; all 20
                splits; seeds 10-11 (80 runs). --pilot-rerun uses seeds 12-13
                (the pre-declared repeat after an engine fix).
  confirmatory  seeds 100-102 on 20 splits; 7,020 runs (3,000 + 2,640 + 1,380).
                Refuses to start unless (a) the SHA-256 of analyze_runs.py in
                this directory is written in the decision log given by
                --decision-log (it is logged there before the first
                confirmatory run), and (b) --pilot-check points to a
                pilot_check.json whose engine check passed.

Each run is one call of run_search.py v1.1 with single-threaded BLAS. Existing
records are kept, so an interrupted batch can be resumed; a record that cannot
be read stops the batch (nothing is deleted or overwritten). A manifest with
every run hash is written at the end.

Usage
-----
  python run_batch.py --stage pilot --data-dir data/derived --out-dir runs/pilot --workers 2
"""

import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import analyze_runs as ar

HERE = Path(__file__).resolve().parent
PILOT_SEEDS, PILOT_RERUN_SEEDS = (10, 11), (12, 13)


def pilot_runs(seeds=PILOT_SEEDS):
    return [("primary", sp, obj, "P0", 0.0, seed)
            for sp in ar.all_splits() for seed in seeds for obj in ("Y1", "Y2")]


def record_name(key):
    u, sp, obj, prior, beta, seed = key
    return f"{u}_{sp}_D_{obj}_{prior.replace(':', '-')}_b{beta}_s{seed}.json"


def command(key, data_dir, out_dir, registered):
    u, sp, obj, prior, beta, seed = key
    d = Path(data_dir)
    return [sys.executable, str(HERE / "run_search.py"),
            "--features", str(d / ar.FILES["features"].format(u=u)), "--pool", str(d / ar.FILES["pool"]),
            "--splits", str(d / ar.FILES["splits"]), "--priors", str(d / ar.FILES["priors"].format(u=u)),
            "--universe", u, "--split", sp, "--objective", obj, "--prior", prior,
            "--beta", str(beta), "--seed", str(seed), "--side", "D",
            "--registered", registered, "--out-dir", str(out_dir)]


def run_one(key, data_dir, out_dir, registered):
    f = Path(out_dir) / record_name(key)
    if f.exists():
        try:
            return key, json.loads(f.read_text())["run_hash"], "kept"
        except (ValueError, KeyError):
            return key, None, "existing record cannot be read; inspect and remove it by hand"
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    p = subprocess.run(command(key, data_dir, out_dir, registered), env=env, capture_output=True, text=True)
    if p.returncode != 0:
        return key, None, (p.stderr or p.stdout).strip()[-400:]
    return key, json.loads(f.read_text())["run_hash"], "ran"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["pilot", "confirmatory"])
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--pilot-rerun", action="store_true")
    ap.add_argument("--decision-log", default=None)
    ap.add_argument("--pilot-check", default=None)
    a = ap.parse_args(argv)

    if ar.sha256_file(HERE / "run_search.py") != ar.ENGINE_SHA256:
        raise SystemExit("run_search.py is not the registered engine (hash mismatch)")
    if a.stage == "pilot":
        keys = pilot_runs(PILOT_RERUN_SEEDS if a.pilot_rerun else PILOT_SEEDS)
    else:
        analysis_sha = ar.sha256_file(HERE / "analyze_runs.py")
        try:
            logged = analysis_sha in Path(a.decision_log).read_text(encoding="utf-8")
        except (TypeError, OSError):
            logged = False
        if not logged:
            raise SystemExit("Refusing the confirmatory stage: the SHA-256 of analyze_runs.py "
                             f"({analysis_sha}) is not written in the decision log given by --decision-log")
        try:
            ok = json.loads(Path(a.pilot_check).read_text())["check3_engine"]["pass"] is True
        except (TypeError, OSError, ValueError, KeyError):
            ok = False
        if not ok:
            raise SystemExit("Refusing the confirmatory stage: --pilot-check must point to a "
                             "pilot_check.json whose engine check passed")
        keys = ar.expected_runs()
    Path(a.out_dir).mkdir(parents=True, exist_ok=True)
    print(f"{a.stage}: {len(keys)} runs, {a.workers} workers", flush=True)
    done, failed = [], []
    with ThreadPoolExecutor(a.workers) as ex:
        for i, (key, h, status) in enumerate(ex.map(
                lambda k: run_one(k, a.data_dir, a.out_dir, ar.REGISTRATION), keys), 1):
            (done if h else failed).append((key, h, status))
            if i % 50 == 0 or i == len(keys):
                print(f"  {i}/{len(keys)}", flush=True)
    manifest = {"stage": a.stage, "registration": ar.REGISTRATION, "n_runs": len(done),
                "engine_sha256": ar.ENGINE_SHA256, "analysis_sha256": ar.sha256_file(HERE / "analyze_runs.py"),
                "runs": {record_name(k): h for k, h, _ in sorted(done)}}
    (Path(a.out_dir) / f"MANIFEST_{a.stage}.json").write_text(json.dumps(manifest, indent=1))
    if failed:
        for k, _, msg in failed[:10]:
            print("FAILED", record_name(k), "::", msg)
        raise SystemExit(f"{len(failed)} runs failed")
    print(f"all {len(done)} runs present; manifest written")


if __name__ == "__main__":
    main()

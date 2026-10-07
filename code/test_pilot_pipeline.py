"""
test_pilot_pipeline.py — end-to-end test of run_batch.py and pilot_check.py (v1.1)
==================================================================================
Runs the real engine (run_search.py v1.1) through the pilot stage on a
SYNTHETIC pool (no Materials Project value is read), then the three pilot
checks including the 5 hash-reproducing reruns, and checks the guards of the
confirmatory stage.

  python test_pilot_pipeline.py      (about 3 minutes on 2 cores)
"""

import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

import analyze_runs as ar
import test_analysis as ta

HERE = Path(__file__).resolve().parent


def main():
    root = Path(tempfile.mkdtemp(prefix="p2_pilot_test_"))
    try:
        d, _, ctx, hashes = ta.build(root, n=300, runs=False)
        ids = [r["material_id"] for r in csv.DictReader(open(d / "pool_v1.csv"))]
        X = np.random.default_rng(3).normal(size=(len(ids), 6))
        with open(d / "features_v2_primary.csv", "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["material_id"] + [f"f{j}" for j in range(6)] + ["const"])
            for m, x in zip(ids, X):
                w.writerow([m] + [repr(float(v)) for v in x] + [1.0])
        out = root / "runs_pilot"
        run = lambda *a: subprocess.run([sys.executable, *a], cwd=HERE, capture_output=True, text=True)
        sha = ar.sha256_file(HERE / "analyze_runs.py")

        log_bad, log_ok = root / "log_bad.md", root / "log_ok.md"
        log_bad.write_text("| 35 | analysis code `" + "0" * 64 + "` |")
        log_ok.write_text("| 35 | analysis code `" + sha + "` |")
        p = run("run_batch.py", "--stage", "confirmatory", "--data-dir", str(d), "--out-dir", str(out))
        assert p.returncode != 0 and "is not written in the decision log" in p.stderr
        p = run("run_batch.py", "--stage", "confirmatory", "--data-dir", str(d), "--out-dir", str(out),
                "--decision-log", str(log_bad))
        assert p.returncode != 0 and "is not written in the decision log" in p.stderr
        p = run("run_batch.py", "--stage", "confirmatory", "--data-dir", str(d), "--out-dir", str(out),
                "--decision-log", str(log_ok))
        assert p.returncode != 0 and "--pilot-check must point" in p.stderr
        assert not out.exists() or not list(out.glob("*.json"))
        print("ok  1  confirmatory stage refuses unless the analysis hash is in the decision log and the pilot check passed")

        p = run("run_batch.py", "--stage", "pilot", "--data-dir", str(d), "--out-dir", str(out), "--workers", "2")
        assert p.returncode == 0, p.stdout + p.stderr
        recs = sorted(out.glob("primary_*.json"))
        man = json.loads((out / "MANIFEST_pilot.json").read_text())
        R = [json.loads(r.read_text()) for r in recs]
        assert len(recs) == 80 and man["n_runs"] == 80
        assert {r["run"]["seed"] for r in R} == {10, 11} and {r["run"]["beta"] for r in R} == {0.0}
        assert {r["run"]["prior"] for r in R} == {"P0"} and {r["run"]["objective"] for r in R} == {"Y1", "Y2"}
        assert all(r["registered"] == ar.REGISTRATION and r["engine_version"] == "1.1" for r in R)
        assert all(set(r["env"]["threads"].values()) == {"1"} for r in R)
        print("ok  2  pilot stage: 80 records, beta = 0 only, both objectives, seeds 10-11, DOI and threads recorded")

        p2 = run("run_batch.py", "--stage", "pilot", "--data-dir", str(d), "--out-dir", str(out))
        assert p2.returncode == 0 and json.loads((out / "MANIFEST_pilot.json").read_text())["runs"] == man["runs"]
        p3 = run("run_search.py", *rb_command(d, out)[2:])
        assert p3.returncode != 0 and "not overwritten" in p3.stderr
        print("ok  3  resuming keeps existing records; the engine refuses to overwrite a record")

        p = run("pilot_check.py", "--data-dir", str(d), "--runs-dir", str(out))
        assert p.returncode == 0, p.stdout + p.stderr
        res = json.loads((out / "pilot_check.json").read_text())
        e = res["check3_engine"]
        assert res["n_runs"] == 80 and e["pass"] is True and len(e["rerun_hash_match"]) == 5
        assert set(res) == {"n_runs", "check1_S_dis", "check2_T99", "check3_engine", "fallbacks"}
        assert set(res["check1_S_dis"]) == {"Y1", "Y2"} and len(res["check1_S_dis"]["Y2"]["mean_hits"]) == 4
        print("ok  4  pilot checks run per objective; 5 drawn runs reproduce their run hash")

        p = run("pilot_check.py", "--data-dir", str(d), "--runs-dir", str(out), "--no-rerun")
        assert json.loads((out / "pilot_check.json").read_text())["check3_engine"]["pass"] is None
        p = run("run_batch.py", "--stage", "confirmatory", "--data-dir", str(d), "--out-dir", str(root / "c"),
                "--decision-log", str(log_ok), "--pilot-check", str(out / "pilot_check.json"))
        assert p.returncode != 0 and "--pilot-check must point" in p.stderr
        try:
            ar.fallbacks_from_pilot(out / "pilot_check.json")
            raise AssertionError("analysis accepted a pilot check without the engine reruns")
        except SystemExit:
            pass
        (out / "pilot_check.json").write_text(json.dumps(res))
        t10, n5, h = ar.fallbacks_from_pilot(out / "pilot_check.json")
        assert t10 == tuple(o for o in ("Y1", "Y2") if not res["check1_S_dis"][o]["pass"])
        assert n5 == tuple(o for o in ("Y1", "Y2") if not res["check2_T99"][o]["pass"]) and len(h) == 64
        bad = recs[3].read_text()
        recs[3].write_text("{ not json")
        p = run("run_batch.py", "--stage", "pilot", "--data-dir", str(d), "--out-dir", str(out))
        assert p.returncode != 0 and "cannot be read" in p.stdout + p.stderr and recs[3].read_text() == "{ not json"
        recs[3].write_text(bad)
        print("ok  5  a pilot check without the reruns unlocks nothing; the analysis reads its fallbacks from")
        print("       the pilot check; an unreadable record stops the batch and is left untouched")

        rec = json.loads(recs[7].read_text())
        rec["trajectory"][0], rec["trajectory"][1] = rec["trajectory"][1], rec["trajectory"][0]
        recs[7].write_text(json.dumps(rec))
        p = run("pilot_check.py", "--data-dir", str(d), "--runs-dir", str(out), "--no-rerun")
        assert p.returncode != 0 and "RUN VALIDATION FAILED" in p.stderr
        print("ok  6  a tampered pilot record is rejected")
        print("ALL PILOT PIPELINE TESTS PASSED")
    finally:
        shutil.rmtree(root, ignore_errors=True)


def rb_command(d, out):
    import run_batch as rb
    return rb.command(rb.pilot_runs()[0], d, out, ar.REGISTRATION)


if __name__ == "__main__":
    main()

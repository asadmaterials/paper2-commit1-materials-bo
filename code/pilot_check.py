"""
pilot_check.py — the three registered pilot checks, and nothing else (v1.1)
===========================================================================
Registration: https://doi.org/10.5281/zenodo.23118101 (PREREGISTRATION.md v1.1, Section 9)

The pilot contains beta = 0 runs only (no prior is used in any pilot search),
40 per objective. Checks, each per objective:

  1. S_dis measurable   for every prior-objective cell of the objective, the
                        mean number of S_dis members acquired per run is >= 2
                        (random search gives about 1.0)
  2. T99 measurable     T99 is censored (= 41) in < 50% of runs
  3. Engine             the two objectives share the initial design of each
                        (split, seed); every record validates (engine hash,
                        registration, inputs, design, oracle values, step
                        log); 5 runs drawn with default_rng(20260928) from the
                        sorted record names reproduce their run hash.

Fallbacks apply to the objective whose check failed:
  check 1 -> S_top becomes the top 10% for that objective
  check 2 -> Delta N5 replaces Delta T99 for that objective
  (analyze_runs.py reads both from the pilot_check.json written here)
  check 3 -> fix the engine (logged deviation), rerun the pilot with seeds 12-13

Usage
-----
  python pilot_check.py --data-dir data/derived --runs-dir runs/pilot [--pilot-rerun]
"""

import argparse
import json
import tempfile
from pathlib import Path

import numpy as np

import analyze_runs as ar
import run_batch as rb

MIN_HITS, MAX_CENSORED = 2.0, 0.5


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--runs-dir", required=True)
    ap.add_argument("--pilot-rerun", action="store_true")
    ap.add_argument("--no-rerun", action="store_true",
                    help="skip the 5 reruns; the engine check is then reported as not done")
    a = ap.parse_args(argv)

    keys = rb.pilot_runs(rb.PILOT_RERUN_SEEDS if a.pilot_rerun else rb.PILOT_SEEDS)
    ctx, hashes = ar.load_data(a.data_dir, ("primary",), {"Y1": 0.05, "Y2": 0.05})
    runs = ar.load_runs(a.runs_dir, ctx, hashes, ("primary",), want=keys)

    res = {"n_runs": len(runs), "check1_S_dis": {}, "check2_T99": {}}
    for obj in ("Y1", "Y2"):
        mine = [k for k in keys if k[2] == obj]
        hits = {p: [] for o, p in ar.CELLS["primary"] if o == obj}
        t99 = []
        for k in mine:
            c, r = ctx[("primary", k[1])], runs[k]
            for p in hits:
                hits[p].append(ar.run_outcomes(c, obj, r["init"], r["trajectory"], c.s_dis(obj, p))["dis_hits"])
            t99.append(ar.run_outcomes(c, obj, r["init"], r["trajectory"])["T99"])
        mean_hits = {p: float(np.mean(v)) for p, v in hits.items()}
        cens = float(np.mean(np.array(t99) == 41))
        res["check1_S_dis"][obj] = {"n_runs": len(mine), "mean_hits": mean_hits, "required": MIN_HITS,
                                    "pass": bool(all(v >= MIN_HITS for v in mean_hits.values()))}
        res["check2_T99"][obj] = {"censored_fraction": cens, "required_below": MAX_CENSORED,
                                  "pass": bool(cens < MAX_CENSORED)}

    designs = {}
    for k in keys:
        designs.setdefault((k[1], k[5]), set()).add(tuple(runs[k]["init"]))
    same_init = all(len(v) == 1 for v in designs.values())
    names = sorted(rb.record_name(k) for k in keys)
    by_name = {rb.record_name(k): k for k in keys}
    drawn = [names[i] for i in np.random.default_rng(20260928).choice(len(names), 5, replace=False)]
    rerun = {}
    if not a.no_rerun:
        with tempfile.TemporaryDirectory() as tmp:
            for nm in drawn:
                _, h, _ = rb.run_one(by_name[nm], a.data_dir, tmp, ar.REGISTRATION)
                rerun[nm] = bool(h == runs[by_name[nm]]["run_hash"])
    res["check3_engine"] = {"identical_initial_designs": bool(same_init), "records_validated": True,
                            "reruns_done": not a.no_rerun, "rerun_hash_match": rerun,
                            "pass": None if a.no_rerun else bool(same_init and all(rerun.values()))}

    fb = []
    for obj in ("Y1", "Y2"):
        if not res["check1_S_dis"][obj]["pass"]:
            fb.append(f"{obj}: S_top -> top 10%")
        if not res["check2_T99"][obj]["pass"]:
            fb.append(f"{obj}: efficiency outcome -> N5")
    if res["check3_engine"]["pass"] is False:
        fb.append("engine: fix (logged deviation), rerun the pilot with --pilot-rerun")
    res["fallbacks"] = fb
    (Path(a.runs_dir) / "pilot_check.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

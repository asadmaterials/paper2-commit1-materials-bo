"""
posthoc_reporting.py — descriptive reporting quantities added after the registered analysis
==========================================================================================
NOT part of the registered analysis. Added after entry 40 in response to a review; logged in
the decision log. Reads the frozen data, the confirmatory records and results/run_outcomes.csv.

  t99_censoring.csv   share of runs with no top-1% hit in the 40 acquisitions (T99 = 41),
                      per universe, cell and beta
  init_top1.csv       per universe and objective: number of (split, seed) blocks whose
                      shared 10-point initial design already contains a top-1% material

Usage (from code/):  python posthoc_reporting.py ../data/derived ../runs/confirmatory ../results ../results_posthoc
"""
import json, sys
from pathlib import Path
import pandas as pd
import analyze_runs as ar

data, runs, res, out = map(Path, sys.argv[1:5])
out.mkdir(parents=True, exist_ok=True)

ro = pd.read_csv(res / "run_outcomes.csv")
cen = (ro.assign(no_hit=ro.T99 == 41).groupby(["universe", "cell", "beta"])
       .agg(runs=("no_hit", "size"), no_hit_share=("no_hit", "mean")).reset_index())
cen.to_csv(out / "t99_censoring.csv", index=False, float_format="%.4f")

ctx, _ = ar.load_data(data, ar.UNIVERSES, {"Y1": 0.05, "Y2": 0.05})
rows = []
for u in ar.UNIVERSES:
    for obj in ("Y1", "Y2"):
        n_blocks = n_hit = 0
        per_type = {"random": 0, "chemsys": 0}
        for sp in ar.all_splits():
            c = ctx[(u, sp)]
            for seed in (100, 101, 102):
                f = runs / f"{u}_{sp}_D_{obj}_P0_b0.0_s{seed}.json"
                init = json.loads(f.read_text())["init"]
                hit = any(c.pos[m] in c.top1[obj] for m in init)
                n_blocks += 1; n_hit += hit; per_type[sp.split("_")[0]] += hit
        rows.append({"universe": u, "objective": obj, "blocks": n_blocks, "blocks_init_has_top1": n_hit,
                     "share": n_hit / n_blocks, "random_type": per_type["random"], "chemsys_type": per_type["chemsys"]})
pd.DataFrame(rows).to_csv(out / "init_top1.csv", index=False, float_format="%.4f")
print(pd.DataFrame(rows).to_string(index=False))

# ---- initial-design sensitivity (exploratory): a block whose shared initial design holds a
# top-1% material gets T99 = 0 in every arm, so its paired gain is 0; same split-level t bound.
import numpy as np
from scipy import stats
hitblocks = set()
for u in ar.UNIVERSES:
    for obj in ("Y1", "Y2"):
        for sp in ar.all_splits():
            c = ctx[(u, sp)]
            for seed in (100, 101, 102):
                init = json.loads((runs / f"{u}_{sp}_D_{obj}_P0_b0.0_s{seed}.json").read_text())["init"]
                if any(c.pos[m] in c.top1[obj] for m in init):
                    hitblocks.add((u, obj, sp, seed))

def lower95(v):
    m = v.groupby(level=0).mean(); typ = m.index.str.split("_").str[0]
    vr, vc = m[typ == "random"].var(ddof=1), m[typ == "chemsys"].var(ddof=1)
    se = np.sqrt((vr + vc) / 40); df = (vr + vc) ** 2 / ((vr ** 2 + vc ** 2) / 9)
    return m.mean(), m.mean() - stats.t.ppf(0.95, df) * se

sens = []
for (u, cell), g in ro.groupby(["universe", "cell"]):
    obj = cell.split(":")[0]
    t = g.set_index(["split", "seed", "beta"]).T99
    hit = pd.Series([(u, obj, sp, s) in hitblocks for sp, s, _ in t.index], index=t.index)
    t0 = t.where(~hit, 0)
    for b in sorted(g.beta.unique()):
        if b == 0:
            continue
        reg = t.xs(0.0, level="beta") - t.xs(b, level="beta")
        alt = t0.xs(0.0, level="beta") - t0.xs(b, level="beta")
        (m1, l1), (m2, l2) = lower95(reg), lower95(alt)
        sens.append({"universe": u, "cell": cell, "beta": b, "gain_registered": m1, "lower95_registered": l1,
                     "gain_init_as_time0": m2, "lower95_init_as_time0": l2,
                     "beneficial_registered": l1 > 0, "beneficial_init_as_time0": l2 > 0})
sens = pd.DataFrame(sens)
sens.to_csv(out / "init_design_sensitivity.csv", index=False, float_format="%.4f")
p = sens[sens.universe == "primary"]
print("\nbeneficial decisions changed (primary, beta 0.02-0.75):",
      p[(p.beta < 1) & (p.beneficial_registered != p.beneficial_init_as_time0)][["cell", "beta", "lower95_registered", "lower95_init_as_time0"]].round(3).to_string(index=False))

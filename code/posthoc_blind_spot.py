"""
posthoc_blind_spot.py — interval for the loss of the absolute blind spot (post hoc, exploratory)
================================================================================================
NOT part of the registered analysis. Added after entry 42, prompted by a review of the Results
draft; logged as decision-log entry 43. Specification fixed before any value was computed:

  quantity   L_blind(beta) = R_blind(0) - R_blind(beta), paired on (split, seed), where R_blind is the
             per-run recall of S_blind (top-5% materials below the prior's 80th percentile, minus the
             initial design), as written by the registered analysis in results/run_outcomes.csv
  estimator  the registered split-level one: mean over the 3 seeds of a split, mean of the 20 split
             means, SE = sqrt((v_random + v_chemsys) / 40), Satterthwaite df, two-sided 95% t interval
             (the registered functions analyze_runs.cube, split_t and t_interval)
  also       the level R_blind(beta) with the same interval
  scope      primary universe, all 6 cells, every beta in 0.02 ... 1; everything is reported
  not done   no p-values, no tests, no decision rule; no new runs; registered results unchanged

Usage (from code/):  python posthoc_blind_spot.py ../results ../results_posthoc
"""
import sys
from pathlib import Path

import pandas as pd

import analyze_runs as ar

res, out = Path(sys.argv[1]), Path(sys.argv[2])
ro = pd.read_csv(res / "run_outcomes.csv")
ro = ro[ro.universe == "primary"]
rows = []
for cell, g in ro.groupby("cell"):
    t = {(r.split, r.seed, r.beta): r.R_blind for r in g.itertuples()}
    for b in sorted(g.beta.unique()):
        row = {"universe": "primary", "cell": cell, "beta": b}
        for name, fn in (("R_blind", lambda sp, sd: t[(sp, sd, b)]),
                         ("L_blind", lambda sp, sd: t[(sp, sd, 0.0)] - t[(sp, sd, b)])):
            V = ar.cube(fn)
            est, _, df = ar.split_t(V)
            lo, hi = ar.t_interval(V)
            row.update({name: est, f"{name}_lo": lo, f"{name}_hi": hi, f"{name}_df": df})
        rows.append(row)
out.mkdir(parents=True, exist_ok=True)
df = pd.DataFrame(rows)
df.to_csv(out / "blind_spot_loss.csv", index=False, float_format="%.6f")
print(df[df.beta.isin([0.1, 0.25, 0.5, 0.75])][["cell", "beta", "R_blind", "L_blind", "L_blind_lo", "L_blind_hi"]].round(3).to_string(index=False))

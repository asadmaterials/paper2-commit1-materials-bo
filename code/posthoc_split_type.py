"""
posthoc_split_type.py — gain and coverage loss by split type (post hoc, descriptive)
====================================================================================
NOT part of the registered analysis. Added after entry 44, prompted by a review of the Methods draft
that asked for the main estimates within each split type; logged as decision-log entry 45.
Specification, fixed before any value was computed:

  quantities  gain = T99(0) - T99(beta) and L = R_dis(0) - R_dis(beta), paired on (split, seed), from
              results/run_outcomes.csv
  estimator   within each split type (random, chemsys): the mean of its 10 split means (each the mean of
              the 3 seeds), with a two-sided 95% t interval, SE = sqrt(v_type / 10), 9 degrees of freedom
  scope       primary universe, all 6 cells, every beta in 0.02 ... 1; everything is reported
  not done    no p-values, tests or decisions; no new runs; registered results unchanged

Usage (from code/):  python posthoc_split_type.py ../results ../results_posthoc
"""
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

res, out = Path(sys.argv[1]), Path(sys.argv[2])
ro = pd.read_csv(res / "run_outcomes.csv")
ro = ro[ro.universe == "primary"]
rows = []
for cell, g in ro.groupby("cell"):
    t = g.set_index(["split", "seed", "beta"])
    b0 = t.xs(0.0, level="beta")
    for b in sorted(g.beta.unique()):
        if b == 0:
            continue
        bb = t.xs(b, level="beta")
        for name, v in (("gain", b0.T99 - bb.T99), ("L", b0.R_dis - bb.R_dis)):
            m = v.groupby(level="split").mean()
            for st in ("random", "chemsys"):
                x = m[m.index.str.startswith(st + "_")].to_numpy()
                assert len(x) == 10
                est, se = float(x.mean()), float(x.std(ddof=1) / math.sqrt(10))
                q = float(stats.t.ppf(0.975, 9))
                rows.append({"universe": "primary", "cell": cell, "beta": b, "quantity": name, "split_type": st,
                             "estimate": est, "ci_lo": est - q * se, "ci_hi": est + q * se, "n_splits": 10})
out.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(out / "by_split_type.csv", index=False, float_format="%.6f")
print(len(rows), "rows")

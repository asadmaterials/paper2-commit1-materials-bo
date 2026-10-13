"""
posthoc_h2_adjusted.py — H2 with the beneficial bound also adjusted across weights (post hoc, exploratory)
==========================================================================================================
NOT part of the registered analysis. Added after entry 43, prompted by a review of the Methods draft;
logged as decision-log entry 44. In the registered H2 rule only the safety bound is adjusted across the
six weights (level 1 - 0.05/6); the beneficial bound is an unadjusted one-sided 95% bound, so the joint
"beneficial and safe" call has no controlled error rate. Specification (written after a quick check of
the three registered beneficial-and-safe sets had been seen; recorded in entry 44):

  rule       the registered H2 rule (analyze_runs.h2_rule) unchanged, except that a weight is beneficial
             if the one-sided lower bound of mean gain at level 1 - 0.05/6 (the registered SAFE_LEVEL) is
             above 0. With both bounds at that level, the probability of declaring at least one weight of
             a cell beneficial-and-safe when it is not is at most 0.05 + 0.05 = 0.10 (union bound).
  inputs     results/run_outcomes.csv (gain = T99(0) - T99(beta); loss = R_dis(0) - R_dis(beta); R_dis(0);
             chance recall), the registered efficiency outcome T99 (no pilot fallback was triggered);
             cubes built with analyze_runs.cube; bounds from analyze_runs.t_interval
  margins    20% (registered) and the registered descriptive 10% and 25%
  scope      all three universes, all cells, beta in H2_BETAS; everything is reported
  check      the registered decisions (beneficial at 95%, safe, safe_at_10/25pct) are recomputed from the same
             inputs and must equal results/h2_detail.csv
  not done   no new runs; the registered H2 result is unchanged and remains the registered result

Usage (from code/):  python posthoc_h2_adjusted.py ../results ../results_posthoc
"""
import sys
from pathlib import Path

import pandas as pd

import analyze_runs as ar

res, out = Path(sys.argv[1]), Path(sys.argv[2])
ro = pd.read_csv(res / "run_outcomes.csv")
reg = pd.read_csv(res / "h2_detail.csv")
rows = []
for (u, cell), g in ro.groupby(["universe", "cell"]):
    T = {(r.split, r.seed, r.beta): r for r in g.itertuples()}
    R0 = ar.cube(lambda sp, sd: T[(sp, sd, 0.0)].R_dis)
    for b in ar.H2_BETAS:
        gain = ar.cube(lambda sp, sd: T[(sp, sd, 0.0)].T99 - T[(sp, sd, b)].T99)
        loss = ar.cube(lambda sp, sd: T[(sp, sd, 0.0)].R_dis - T[(sp, sd, b)].R_dis)
        lo95 = ar.t_interval(gain, one_sided=True)[0]
        lo_adj = ar.t_interval(gain, level=ar.SAFE_LEVEL, one_sided=True)[0]
        safe = {f: bool(ar.t_interval(loss - f * R0, level=ar.SAFE_LEVEL, one_sided=True)[1] < 0)
                for f in (0.10, 0.20, 0.25)}
        r = reg[(reg.universe == u) & (reg.cell == cell) & (reg.beta == b)].iloc[0]
        assert (lo95 > 0) == (str(r.beneficial) == "True") and safe[0.20] == (str(r.safe) == "True")
        assert safe[0.10] == (str(r.safe_at_10pct) == "True") and safe[0.25] == (str(r.safe_at_25pct) == "True")
        rows.append({"universe": u, "cell": cell, "beta": b, "gain_mean": float(gain.mean()),
                     "gain_lower95": lo95, "gain_lower_adjusted": lo_adj,
                     "beneficial_registered": lo95 > 0, "beneficial_adjusted": lo_adj > 0,
                     "safe_20": safe[0.20], "safe_10": safe[0.10], "safe_25": safe[0.25],
                     "bs20_registered": lo95 > 0 and safe[0.20], "bs20_adjusted": lo_adj > 0 and safe[0.20],
                     "bs10_adjusted": lo_adj > 0 and safe[0.10], "bs25_adjusted": lo_adj > 0 and safe[0.25]})
out.mkdir(parents=True, exist_ok=True)
df = pd.DataFrame(rows)
assert (pd.read_csv(res / "h2_beneficial_and_safe.csv").applicable.astype(str) == "True").all()   # all cells applicable
df.to_csv(out / "h2_adjusted_beneficial.csv", index=False, float_format="%.6f")
for u, g in df.groupby("universe"):
    print(u)
    for c, h in g.groupby("cell"):
        f = lambda col: "{" + ", ".join(f"{x:g}" for x in h[h[col]].beta) + "}"
        print(f"  {c:<22} registered {f('bs20_registered'):<24} adjusted {f('bs20_adjusted'):<20} "
              f"adj@10 {f('bs10_adjusted'):<12} adj@25 {f('bs25_adjusted')}")

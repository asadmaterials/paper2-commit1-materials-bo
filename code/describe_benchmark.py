"""
describe_benchmark.py — descriptive facts about the frozen benchmark, for the manuscript
=======================================================================================
Not an analysis of search outcomes: reads only the frozen derived data (pool, splits, priors) through the
registered loader (analyze_runs.load_data) and writes counts that the Methods state.

  primary_structures, primary_compositions        structures and distinct reduced compositions (group_id)
  multi_polymorph_compositions, multi_polymorph_structures
                                                  compositions with more than one structure, and their structures
  sdis_tie_splits/<cell>                          splits in which prior scores tied at the S_dis boundary put
                                                  members of S_top on both sides (resolved by the registered
                                                  hash key); sdis_tie_max/<cell>: most such materials outside
  n_splits

Usage (from code/):  python describe_benchmark.py . ../data/derived ../results_posthoc/benchmark_facts.csv
"""
import math
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, sys.argv[1])
import analyze_runs as ar  # noqa: E402

data, out = Path(sys.argv[2]), Path(sys.argv[3])
rows = []
pool = pd.read_csv(data / "pool_v1.csv")
prim = pool[pool.in_primary.astype(str) == "True"]
g = prim.groupby("group_id").size()
rows += [("primary_structures", len(prim)), ("primary_compositions", len(g)),
         ("multi_polymorph_compositions", int((g > 1).sum())), ("multi_polymorph_structures", int(g[g > 1].sum()))]

ctx, _ = ar.load_data(data, ("primary",), {"Y1": 0.05, "Y2": 0.05})
ties = {}
for (u, sp), c in ctx.items():
    for obj, prior in ar.CELLS[u]:
        top, v = c.top[obj], c.prior[prior]
        srt = sorted(top, key=lambda i: (v[i], c.hk[i]))          # the registered S_dis ordering
        k = math.ceil(len(top) / 2)
        assert set(srt[:k]) == c.s_dis(obj, prior)
        outside = sum(1 for i in srt[k:] if v[i] == v[srt[k - 1]])
        n, m = ties.get(f"{obj}:{prior}", (0, 0))
        ties[f"{obj}:{prior}"] = (n + (outside > 0), max(m, outside))
for cell, (n, m) in sorted(ties.items()):
    rows += [(f"sdis_tie_splits/{cell}", n), (f"sdis_tie_max/{cell}", m)]
rows.append(("n_splits", len(ar.all_splits())))
rows.append(("sdis_tie_cells", sum(1 for n, _ in ties.values() if n > 0)))
pd.DataFrame(rows, columns=["key", "value"]).to_csv(out, index=False)
print(pd.DataFrame(rows, columns=["key", "value"]).to_string(index=False))

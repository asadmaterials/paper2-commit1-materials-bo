"""
count_theoretical.py — share of the benchmark carrying MP's `theoretical` flag, from the raw snapshot
====================================================================================================
Descriptive count for the manuscript (Methods, Limitations), which the derived pool file cannot give
because pool_v1.csv does not carry the flag. Reads the raw summary documents (the archived snapshot,
doi:10.5281/zenodo.22978221) and pool_v1.csv; no search outcome is read. Every pool material has exactly
one summary document (a pool filter), and the script asserts this.

Writes <out.csv> (key, value):
  <universe>_structures, <universe>_theoretical   for primary, sens1_no_Pm_Tc, unscreened
  <universe>_theoretical_pct                      percentage, one decimal, half away from zero
  <universe>_ehull_above_0.1                      structures more than 0.1 eV/atom above the hull (pool column)

Usage (from code/):  python count_theoretical.py <snapshot>/summary.jsonl.gz ../data/derived/pool_v1.csv \
                         ../results_posthoc/benchmark_facts_snapshot.csv
"""
import gzip
import json
import sys
from decimal import Decimal, ROUND_HALF_UP

import pandas as pd

summ_path, pool_path, out = sys.argv[1:4]
pool = pd.read_csv(pool_path)
ids = set(pool.material_id)
flag = {}
with gzip.open(summ_path, "rt") as f:
    for line in f:
        d = json.loads(line)
        m = d.get("material_id")
        if m in ids:
            assert m not in flag, f"{m}: more than one summary document"
            flag[m] = bool(d.get("theoretical"))
assert set(flag) == ids, "a pool material has no summary document"

rows = []
for u in ("primary", "sens1_no_Pm_Tc", "unscreened"):
    q = pool[pool[f"in_{u}"].astype(str) == "True"]
    n, t = len(q), sum(flag[m] for m in q.material_id)
    pct = Decimal(100 * t) / Decimal(n)
    rows += [(f"{u}_structures", n), (f"{u}_theoretical", t),
             (f"{u}_theoretical_pct", str(pct.quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))),
             (f"{u}_ehull_above_0.1", int((q.energy_above_hull_eV > 0.1).sum()))]
pd.DataFrame(rows, columns=["key", "value"]).to_csv(out, index=False)
print(pd.DataFrame(rows, columns=["key", "value"]).to_string(index=False))

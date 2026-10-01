"""
make_splits.py — Paper 2 history/discovery splits
=================================================
Assigns every structure in pool_v1.csv to HISTORY (H) or DISCOVERY (D) for
each split type and replicate.

- History: the only data the historical prior P_H may see (with its y).
- Discovery: the search space for the BO campaign; S_top is defined in it.

Split types (grouping unit = the unit that never straddles H and D)
-------------------------------------------------------------------
  random   : composition group (group_id: all polymorphs of a composition)
  chemsys  : chemical system (e.g. Al-Ni); implies composition grouping

Dropped from v1: structural-family hold-out. Family = formula_anonymous +
space group; the Heusler family (ABC2, #225) alone is ~1/3 of the pool and,
after merging polymorph-linked families, forms one block of ~35%. A family
hold-out would then be decided by which side that one block falls on.
Families are still used for the family-level under-discovery analysis.

Design
------
- Splits are defined ONCE on the unscreened universe (4,711) and restricted
  to the primary (3,550) and sensitivity-1 (3,378) universes, so materials
  shared by the universes keep the same assignment in every analysis.
- Groups are shuffled with a seed derived from sha256("paper2-splits-v1:
  <type>:<replicate>") and added to history until the target fraction of
  unscreened structures is reached. Seeding per replicate means adding
  replicates later never changes existing ones.
- The objective is not used, so the same splits serve Y1 (G) and Y2 (G/rho).

Usage:  python make_splits.py --pool pool_v1.csv [--hist-frac 0.30] [--n-rep 3]
"""

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

SPLIT_TYPES = {"random": "group_id", "chemsys": "chemsys"}
UNIVERSES = {"primary": "in_primary", "sens1_no_Pm_Tc": "in_sens1_no_Pm_Tc",
             "unscreened": "in_unscreened"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def seed_for(split_type, rep):
    return int(hashlib.sha256(f"paper2-splits-v1:{split_type}:{rep}".encode()).hexdigest()[:16], 16)


def assign(rows, key, frac, seed):
    size = Counter(r[key] for r in rows)
    groups = sorted(size)
    order = np.random.default_rng(seed).permutation(len(groups))
    target, n, hist = frac * len(rows), 0, set()
    for i in order:
        if n >= target:
            break
        hist.add(groups[i])
        n += size[groups[i]]
    return {r["material_id"]: ("H" if r[key] in hist else "D") for r in rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)
    ap.add_argument("--hist-frac", type=float, default=0.30)
    ap.add_argument("--n-rep", type=int, default=3)
    ap.add_argument("--out-dir", default=".")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.pool)))
    for r in rows:
        r["G"] = float(r["G_vrh_GPa"])
        r["G_rho"] = r["G"] / float(r["density_g_cm3"])
        r["family"] = f'{r["formula_anonymous"]}|{r["spacegroup_number"]}'

    cols, stats = {}, {}
    for st, key in SPLIT_TYPES.items():
        for rep in range(a.n_rep):
            name = f"{st}_r{rep}"
            asg = assign(rows, key, a.hist_frac, seed_for(st, rep))
            cols[name] = asg
            # leakage checks: grouping unit and composition never straddle
            for unit in {key, "group_id"}:
                sides = defaultdict(set)
                for r in rows:
                    sides[r[unit]].add(asg[r["material_id"]])
                bad = [g for g, s in sides.items() if len(s) > 1]
                assert not bad, f"{name}: {unit} straddles H/D: {bad[:3]}"
            s = {}
            for uname, flag in UNIVERSES.items():
                U = [r for r in rows if r[flag] == "True"]
                H = [r for r in U if asg[r["material_id"]] == "H"]
                D = [r for r in U if asg[r["material_id"]] == "D"]
                s[uname] = {
                    "history": len(H), "discovery": len(D),
                    "history_frac": round(len(H) / len(U), 4),
                    "history_top10pct_for_P_H": math.ceil(0.10 * len(H)),
                    "discovery_top5pct": math.ceil(0.05 * len(D)),
                    "discovery_top10pct": math.ceil(0.10 * len(D)),
                    "median_G_history": float(np.median([r["G"] for r in H])),
                    "median_G_discovery": float(np.median([r["G"] for r in D])),
                    "discovery_families": len({r["family"] for r in D}),
                    "discovery_heusler_share": round(
                        sum(r["family"] == "ABC2|225" for r in D) / len(D), 4),
                }
            stats[name] = s

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "splits_v1.csv"
    names = list(cols)
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["material_id"] + names)
        for r in sorted(rows, key=lambda r: r["material_id"]):
            w.writerow([r["material_id"]] + [cols[n][r["material_id"]] for n in names])

    manifest = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "pool_csv_sha256": sha256(a.pool),
        "script_sha256": sha256(__file__),
        "hist_frac_target": a.hist_frac, "n_rep": a.n_rep,
        "split_types": SPLIT_TYPES,
        "seeds": {f"{st}_r{r}": str(seed_for(st, r)) for st in SPLIT_TYPES for r in range(a.n_rep)},
        "defined_on": "unscreened universe; restricted to primary and sens1",
        "stats": stats,
        "splits_csv_sha256": sha256(path),
    }
    (out / "splits_v1.manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"{'split':<12}{'H':>6}{'D':>6}{'H frac':>8}{'P_H ref':>8}{'D top5%':>8}"
          f"{'med G H':>9}{'med G D':>9}{'Heusler D':>10}   (primary universe)")
    for n, s in stats.items():
        p = s["primary"]
        print(f"{n:<12}{p['history']:>6}{p['discovery']:>6}{p['history_frac']:>8.3f}"
              f"{p['history_top10pct_for_P_H']:>8}{p['discovery_top5pct']:>8}"
              f"{p['median_G_history']:>9.1f}{p['median_G_discovery']:>9.1f}"
              f"{p['discovery_heusler_share']:>10.3f}")
    print(f"\nsplits_v1.csv sha256: {manifest['splits_csv_sha256']}")


if __name__ == "__main__":
    main()

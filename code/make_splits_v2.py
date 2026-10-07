"""
make_splits_v2.py — Paper 2 splits for registration v1.1 (10 random + 10 chemical-system)
=========================================================================================
Same algorithm and seed rule as make_splits.py v1, with n_rep = 10. Replicates
r0-r2 of both types are therefore identical to splits_v1.csv (asserted when
--check-v1 is given). r3-r9 are new.

Unlike v1, this script reads no outcome column: the manifest reports sizes
only (v1 also stored median G per side). No discovery-side outcome of the new
splits is computed anywhere in the build.

Output : splits_v2.csv (material_id, random_r0..r9, chemsys_r0..r9), splits_v2.manifest.json
Usage  : python make_splits_v2.py --pool pool_v1.csv --check-v1 splits_v1.csv [--out-dir .]
"""

import argparse
import csv
import datetime as dt
import json
from collections import defaultdict
from pathlib import Path

from make_splits import SPLIT_TYPES, UNIVERSES, assign, seed_for, sha256

N_REP, HIST_FRAC = 10, 0.30


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)
    ap.add_argument("--check-v1", required=True)
    ap.add_argument("--out-dir", default=".")
    a = ap.parse_args()

    keep = ("material_id", "group_id", "chemsys") + tuple(UNIVERSES.values())
    rows = [{k: r[k] for k in keep} for r in csv.DictReader(open(a.pool))]   # no outcome column kept
    cols, stats = {}, {}
    for st, key in SPLIT_TYPES.items():
        for rep in range(N_REP):
            name = f"{st}_r{rep}"
            asg = assign(rows, key, HIST_FRAC, seed_for(st, rep))
            cols[name] = asg
            for unit in {key, "group_id"}:
                sides = defaultdict(set)
                for r in rows:
                    sides[r[unit]].add(asg[r["material_id"]])
                assert all(len(s) == 1 for s in sides.values()), f"{name}: {unit} straddles H/D"
            stats[name] = {u: {"history": sum(asg[r["material_id"]] == "H" for r in rows if r[f] == "True"),
                               "discovery": sum(asg[r["material_id"]] == "D" for r in rows if r[f] == "True")}
                           for u, f in UNIVERSES.items()}

    v1 = {r["material_id"]: r for r in csv.DictReader(open(a.check_v1))}
    old = [c for c in next(iter(v1.values())) if c != "material_id"]
    assert set(v1) == {r["material_id"] for r in rows}
    for c in old:
        assert all(cols[c][m] == v1[m][c] for m in v1), f"{c} differs from splits_v1.csv"

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    names = list(cols)
    with open(out / "splits_v2.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["material_id"] + names)
        for r in sorted(rows, key=lambda r: r["material_id"]):
            w.writerow([r["material_id"]] + [cols[n][r["material_id"]] for n in names])
    manifest = {"created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                "pool_csv_sha256": sha256(a.pool), "splits_v1_csv_sha256": sha256(a.check_v1),
                "script_sha256": sha256(__file__), "make_splits_v1_sha256": sha256(Path(__file__).with_name("make_splits.py")),
                "hist_frac_target": HIST_FRAC, "n_rep": N_REP,
                "identical_to_v1": old, "new": [n for n in names if n not in old],
                "seeds": {n: str(seed_for(*n.rsplit("_r", 1)[0:1], int(n.rsplit("_r", 1)[1]))) for n in names},
                "sizes": stats, "splits_csv_sha256": sha256(out / "splits_v2.csv")}
    (out / "splits_v2.manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"splits_v2.csv: {len(names)} splits; {len(old)} identical to v1; sha256 {manifest['splits_csv_sha256']}")
    for n in names:
        print(f"  {n:11s} primary H/D {stats[n]['primary']['history']}/{stats[n]['primary']['discovery']}")


if __name__ == "__main__":
    main()

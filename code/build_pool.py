"""
build_pool.py — Paper 2 candidate universe from the raw MP snapshot
===================================================================
Reads the raw snapshot written by pull_mp.py (MP database 2026.04.13) and
the frozen elemental reference table, applies the pre-registered filters in
a fixed order, and writes one table of candidates with flags for the
primary universe and the pre-registered sensitivity universes.

Filter order (each step logs how many documents it removes)
-----------------------------------------------------------
 1. elasticity state == "successful"
 2. elasticity document not deprecated
 3. exactly one summary document for the material_id
 4. summary is_metal is True                     (MP PBE-level classification)
 5. nelements >= 2                               (multicomponent scope)
 6. no out-of-scope element (nonmetals, noble gases; SCOPE_EXCLUDED below)
 7. G_VRH finite and > 0
 8. Born stability: all eigenvalues of the symmetrised 6x6 IEEE Voigt
    elastic tensor strictly > 0 (MP's "negative eigenvalue" warning is
    reported as an audit, not used as the criterion)
    -> UNSCREENED universe (sensitivity 2: P0 and P_H only)
 9. every element is eligible in the elemental reference table (rule C)
    -> PRIMARY universe (P0, P_H, P_phys)
    sensitivity 1: additionally exclude materials containing Pm or Tc

No G cutoff, no polymorph collapsing. Polymorphs share a composition group
ID (reduced composition), used later to keep all polymorphs of a
composition on the same side of every split.

Usage
-----
    python build_pool.py --snapshot DIR --reference REF.csv [--out-dir DIR]

DIR must contain elasticity.jsonl.gz, summary.jsonl.gz and manifest.json;
file hashes are verified against the manifest before anything is read.
"""

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np

SCOPE_EXCLUDED = set("H He B C N O F Ne Si P S Cl Ar Se Br Kr I Xe Rn At".split())
SENS1_EXCLUDED = {"Pm", "Tc"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_jsonl_gz(path):
    with gzip.open(path, "rt") as f:
        return [json.loads(l) for l in f]


def born_stable(doc):
    t = (doc.get("elastic_tensor") or {}).get("ieee_format")
    if t is None:
        return False
    a = np.asarray(t, dtype=float)
    if a.shape != (6, 6) or not np.all(np.isfinite(a)):
        return False
    return bool(np.linalg.eigvalsh((a + a.T) / 2).min() > 0)


def group_id(comp_reduced):
    """Canonical reduced-composition key, e.g. 'Al1-Ni1'. Amounts are
    rounded to 6 decimals so float noise cannot split a group."""
    return "-".join(f"{el}{round(float(n), 6):g}" for el, n in sorted(comp_reduced.items()))


def finite_pos(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) and x > 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snapshot", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out-dir", default=".")
    a = ap.parse_args()
    snap = Path(a.snapshot)

    man = json.loads((snap / "manifest.json").read_text())
    for name in ("elasticity.jsonl.gz", "summary.jsonl.gz"):
        want, got = man["files"][name]["sha256"], sha256(snap / name)
        if want != got:
            raise SystemExit(f"{name}: sha256 {got} != manifest {want}")
    print(f"Snapshot verified: MP {man['mp_db_version']}")

    ref = {r["element"]: r for r in csv.DictReader(open(a.reference))}
    eligible = {e for e, r in ref.items() if r["eligible"] == "True"}

    E = read_jsonl_gz(snap / "elasticity.jsonl.gz")
    S = read_jsonl_gz(snap / "summary.jsonl.gz")
    s_count = Counter(d["material_id"] for d in S)
    smap = {d["material_id"]: d for d in S if s_count[d["material_id"]] == 1}

    steps = [
        ("state == successful", lambda e: e.get("state") == "successful"),
        ("not deprecated", lambda e: not e.get("deprecated")),
        ("exactly one summary doc", lambda e: e["material_id"] in smap),
        ("is_metal True", lambda e: smap[e["material_id"]].get("is_metal") is True),
        ("nelements >= 2", lambda e: (e.get("nelements") or 0) >= 2),
        ("no out-of-scope element", lambda e: not set(e["elements"]) & SCOPE_EXCLUDED),
        ("G_VRH finite and > 0", lambda e: finite_pos((e.get("shear_modulus") or {}).get("vrh"))),
        ("Born stable (eigenvalues > 0)", born_stable),
    ]
    funnel = [("all elasticity documents", len(E))]
    cur = E
    for name, f in steps:
        cur = [e for e in cur if f(e)]
        funnel.append((name, len(cur)))
    unscreened = cur

    rows, excl_counter = [], Counter()
    for e in unscreened:
        els = sorted(e["elements"])
        missing = [x for x in els if x not in eligible]
        for x in missing:
            excl_counter[x] += 1
        s = smap[e["material_id"]]
        sym = e.get("symmetry") or {}
        prim = not missing
        rows.append({
            "material_id": e["material_id"],
            "formula_pretty": e.get("formula_pretty"),
            "group_id": group_id(e["composition_reduced"]),
            "chemsys": e.get("chemsys"),
            "elements": ";".join(els),
            "nelements": e["nelements"],
            "nsites": e.get("nsites"),
            "G_vrh_GPa": e["shear_modulus"]["vrh"],
            "K_vrh_GPa": (e.get("bulk_modulus") or {}).get("vrh"),
            "density_g_cm3": e.get("density"),
            "spacegroup_number": sym.get("number"),
            "crystal_system": sym.get("crystal_system"),
            "formula_anonymous": e.get("formula_anonymous"),
            "energy_above_hull_eV": s.get("energy_above_hull"),
            "mp_negative_eigenvalue_warning": any("negative eigenvalue" in w
                                                  for w in (e.get("warnings") or [])),
            "in_primary": prim,
            "in_sens1_no_Pm_Tc": prim and not set(els) & SENS1_EXCLUDED,
            "in_unscreened": True,
            "missing_reference_elements": ";".join(missing),
        })
    rows.sort(key=lambda r: r["material_id"])
    n_prim = sum(r["in_primary"] for r in rows)
    n_s1 = sum(r["in_sens1_no_Pm_Tc"] for r in rows)
    funnel.append(("all elements reference-eligible (PRIMARY)", n_prim))

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    pool_path = out / "pool_v1.csv"
    with open(pool_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    def summarise(sel):
        rs = [r for r in rows if r[sel]]
        g = np.array([r["G_vrh_GPa"] for r in rs])
        groups = Counter(r["group_id"] for r in rs)
        return {"structures": len(rs), "compositions": len(groups),
                "compositions_with_polymorphs": sum(v > 1 for v in groups.values()),
                "elements": len({x for r in rs for x in r["elements"].split(";")}),
                "nelements": dict(sorted(Counter(r["nelements"] for r in rs).items())),
                "G_vrh_GPa": {"min": float(g.min()), "median": float(np.median(g)),
                              "p95": float(np.percentile(g, 95)), "max": float(g.max())},
                "top5pct_size": int(math.ceil(0.05 * len(rs))),
                "top10pct_size": int(math.ceil(0.10 * len(rs)))}

    manifest = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "mp_db_version": man["mp_db_version"],
        "inputs": {"elasticity.jsonl.gz": man["files"]["elasticity.jsonl.gz"]["sha256"],
                   "summary.jsonl.gz": man["files"]["summary.jsonl.gz"]["sha256"],
                   "reference_csv": sha256(a.reference)},
        "script_sha256": sha256(__file__),
        "scope_excluded_elements": sorted(SCOPE_EXCLUDED),
        "sens1_excluded_elements": sorted(SENS1_EXCLUDED),
        "funnel": funnel,
        "structures_removed_by_reference_element": dict(excl_counter.most_common()),
        "audit_mp_neg_eig_warning_in_unscreened": sum(r["mp_negative_eigenvalue_warning"] for r in rows),
        "universes": {"primary": summarise("in_primary"),
                      "sens1_no_Pm_Tc": summarise("in_sens1_no_Pm_Tc"),
                      "unscreened": summarise("in_unscreened")},
        "pool_csv_sha256": sha256(pool_path),
    }
    (out / "pool_v1.manifest.json").write_text(json.dumps(manifest, indent=2))

    print("\nFUNNEL")
    for name, n in funnel:
        print(f"  {name:<46}{n:>6}")
    print(f"  {'sensitivity 1 (also no Pm, Tc)':<46}{n_s1:>6}")
    print("\nStructures removed by reference screen, per element:",
          dict(excl_counter.most_common()))
    for k, v in manifest["universes"].items():
        print(f"\n{k}: {v}")
    print(f"\npool_v1.csv sha256: {manifest['pool_csv_sha256']}")


if __name__ == "__main__":
    main()

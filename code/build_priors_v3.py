"""
build_priors_v3.py — Paper 2 priors for registration v1.1 (20 splits)
=====================================================================
Same prior definitions as build_priors.py v2.1 (imported from it, unchanged):
  P_H_fixed, P_H_matched_Y2   negative distance to the nearest of the top 10% of
                              HISTORY (descriptors standardised on history)
  P_phys_fixed, _matched_Y2   elemental-modulus prior (no MP data)

Differences from v2.1
  - reads splits_v2.csv (20 splits) and the frozen descriptor file
    features_v2_<u>.csv instead of re-featurising (its hash is recorded);
  - writes the alignment table for the HISTORY side only (grouped 5-fold CV for
    P_H). No outcome of any discovery-side material is read into a metric:
    discovery rows are used only as points to be scored by the priors;
  - the retired column S_dis_poolmedian is not written.
For splits r0-r2 the prior scores equal priors_v2_<u>.csv exactly (asserted
with --check-v2).

Outputs: priors_v3_<u>.csv, alignment_v3_<u>.csv, priors_v3_<u>.manifest.json
Usage  : python build_priors_v3.py --pool pool_v1.csv --splits splits_v2.csv
             --reference elemental_reference_springer2005_v1.csv
             --features features_v2_<u>.csv --check-v2 priors_v2_<u>.csv
             [--universe primary] [--out-dir .]
"""

import argparse
import csv
import datetime as dt
import json
from importlib import metadata
from pathlib import Path

import numpy as np

import build_priors as bp

Y2 = "Y2_E_over_rho"


def main():
    ap = argparse.ArgumentParser()
    for k in ("pool", "splits", "reference", "features", "check-v2"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--universe", default="primary", choices=list(bp.UNIVERSE_FLAG))
    ap.add_argument("--out-dir", default=".")
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    u = a.universe

    rows = [r for r in csv.DictReader(open(a.pool)) if r[bp.UNIVERSE_FLAG[u]] == "True"]
    rows.sort(key=lambda r: r["material_id"])
    ids = [r["material_id"] for r in rows]
    ref = {r["element"]: r for r in csv.DictReader(open(a.reference))}
    splits = {r["material_id"]: r for r in csv.DictReader(open(a.splits))}
    split_names = [k for k in next(iter(splits.values())) if k != "material_id"]

    with open(a.features) as f:
        rd = csv.reader(f)
        fnames = next(rd)[1:]
        feat = {r[0]: np.array(r[1:], dtype=float) for r in rd}
    assert set(feat) == set(ids), "descriptor file does not match the universe"
    X = np.array([feat[m] for m in ids])

    G = np.array([float(r["G_vrh_GPa"]) for r in rows])
    K = np.array([float(r["K_vrh_GPa"]) for r in rows])
    rho = np.array([float(r["density_g_cm3"]) for r in rows])
    Y = {"Y1_G": G, Y2: (9 * K * G / (3 * K + G)) / rho}
    fam = np.array([f'{r["formula_anonymous"]}|{r["spacegroup_number"]}' for r in rows])

    have_phys = u != "unscreened"
    if have_phys:
        from pymatgen.core import Composition
        Pfix, Pmat = bp.p_phys([Composition(r["formula_pretty"]) for r in rows], ref)

    prior_rows, align_rows = [], []
    for sp in split_names:
        stype = sp.rsplit("_r", 1)[0]
        side = np.array([splits[m][sp] for m in ids])
        H, D = np.where(side == "H")[0], np.where(side == "D")[0]
        assert len(H) + len(D) == len(ids)
        yH = {k: v[H] for k, v in Y.items()}          # the only outcomes used below are history-side
        units_h = [rows[i][bp.SPLIT_UNIT[stype]] for i in H]
        P_D = {"P_H_fixed": bp.p_hist(X[H], yH["Y1_G"], X[D]),
               "P_H_matched_Y2": bp.p_hist(X[H], yH[Y2], X[D])}
        P_Hs = {"P_H_fixed": bp.p_hist_cv(X[H], yH["Y1_G"], units_h),
                "P_H_matched_Y2": bp.p_hist_cv(X[H], yH[Y2], units_h)}
        if have_phys:
            P_D["P_phys_fixed"], P_Hs["P_phys_fixed"] = Pfix[D], Pfix[H]
            P_D["P_phys_matched_Y2"], P_Hs["P_phys_matched_Y2"] = Pmat[D], Pmat[H]
        for j, i in enumerate(D):
            prior_rows.append({"split": sp, "material_id": ids[i],
                               **{k: repr(float(v[j])) for k, v in P_D.items()}})
        sids = [ids[i] for i in H]
        for yname in Y:
            for pname, p in P_Hs.items():
                if pname.endswith("_matched_Y2") and yname == "Y1_G":
                    continue
                m = bp.metrics(p, yH[yname], sids, sp, fam[H])
                m.pop("S_dis_poolmedian")
                align_rows.append({"universe": u, "side": "history_cv", "split": sp, "split_type": stype,
                                   "objective": yname, "prior": pname, **m})

    # r0-r2 must reproduce the registered v2 priors exactly
    v2 = {(r["split"], r["material_id"]): r for r in csv.DictReader(open(a.check_v2))}
    v3 = {(r["split"], r["material_id"]): r for r in prior_rows}
    old = sorted({k[0] for k in v2})
    assert all(v3[k] == v2[k] for k in v2), "v3 priors differ from priors_v2 on the registered splits"
    assert {k for k in v3 if k[0] in old} == set(v2)

    pf, af = out / f"priors_v3_{u}.csv", out / f"alignment_v3_{u}.csv"
    for path, rws in ((pf, prior_rows), (af, align_rows)):
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rws[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(rws)
    manifest = {"created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), "universe": u,
                "inputs": {k: bp.sha256(getattr(a, k)) for k in ("pool", "splits", "reference", "features")},
                "script_sha256": bp.sha256(__file__), "build_priors_v2_1_sha256": bp.sha256(bp.__file__),
                "versions": {p: metadata.version(p) for p in ("numpy", "scipy", "pymatgen")},
                "n_features": len(fnames), "n_splits": len(split_names),
                "identical_to_priors_v2_on": old, "priors_v2_sha256": bp.sha256(a.check_v2),
                "alignment_sides": ["history_cv"],
                "outputs": {p.name: bp.sha256(p) for p in (pf, af)}}
    (out / f"priors_v3_{u}.manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"{u}: {len(split_names)} splits, {len(prior_rows)} prior rows; identical to v2 on {len(old)} splits")
    for name, h in manifest["outputs"].items():
        print(f"  {h}  {name}")


if __name__ == "__main__":
    main()

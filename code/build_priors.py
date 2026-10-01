"""
build_priors.py — Paper 2 priors (P_H, P_phys) and alignment tables (v2)
========================================================================
Inputs : pool_v1.csv, splits_v1.csv, elemental_reference_springer2005_v1.csv
Outputs: features_v2_<u>.csv   descriptors for every structure (also the GP inputs)
         priors_v2_<u>.csv     raw prior scores, one row per (split, discovery material)
         alignment_v2_<u>.csv  alignment metrics per side x split x prior x objective
         priors_v2_<u>.manifest.json

Changes from v1 (review of 2026-09-27)
  - Y2 is E/rho (conventional specific stiffness), E = 9KG/(3K+G) from MP VRH.
  - nsites and spacegroup_number removed from the descriptors (cell-convention
    dependent / categorical treated as continuous).
  - v2.1: prior scores computed with math.fsum in sorted element order and
    direct-difference distances (no BLAS), then rounded to 12 significant
    digits. Reason: Python >= 3.12 changed float sum(); mathematically tied
    P_phys values (elements with equal tabulated G: Y/Dy 25.5, Er/Ge 29.6)
    were ordered by round-off instead of the hash rule.
  - Alignment is reported on two sides: "history_cv" (a genuine pre-search
    quantity; primary) and "discovery_oracle" (uses discovery ground truth;
    diagnostic only).

Objectives
  Y1 = G_VRH                                   (MP DFT, GPa)
  Y2 = E / rho, E = 9 K_VRH G_VRH / (3 K_VRH + G_VRH), rho = MP density
       (GPa cm^3/g; Spearman with G/rho on the primary pool = 0.998)

Descriptors (same space the GP will use)
  Magpie composition statistics (matminer preset, 132) + structural:
  MP density, volume per atom, crystal-system one-hot. Constant columns (within the universe) are dropped. No elastic
  quantity is used as a descriptor.

Priors (higher = more favoured). Selection later uses rank(P) with a
hash tie-break; raw scores are stored here.
  P_H        historical similarity: -min Euclidean distance, in descriptor
             space standardised on the HISTORY side, to the reference set
             R = top 10% of history by the reference objective.
               fixed   : R ranked by G         (used for Y1 and Y2)
               matched : R ranked by E/rho     (Y2 only)
  P_phys     elemental-physics estimate, atomic-fraction Hill average of
             elemental shear moduli (Springer 2005, screened table v1):
               G_V = sum x_i G_i ; G_R = 1 / sum(x_i / G_i) ; G_H = (G_V+G_R)/2
               fixed   : G_H                (used for Y1 and Y2)
               matched : G_H / rho_mix      (Y2 only), rho_mix volume-additive:
                         rho_mix = sum x_i M_i / sum x_i M_i / rho_i
                         Matched to E/rho up to the factor 2(1+nu), which is
                         not estimated (elemental E is not used: E failed the
                         screen for Sm and Re).
             Uses no MP data. Identical for polymorphs (composition only).
  For Y1, fixed and matched are the same prior.

Leakage rules (asserted)
  - P_H reads y only from HISTORY materials; discovery y is never touched.
  - P_phys reads no MP-derived quantity.
  - Standardisation statistics come from HISTORY descriptors only.

Alignment sides
  history_cv       computed on the HISTORY side of each split. P_phys scores are
                   used directly (they use no MP data). P_H scores are out-of-fold:
                   5-fold CV within history, folds grouped by the split's unit
                   (group_id or chemsys), fold = sha256(unit) mod 5; each fold is
                   scored with reference set and standardisation from the other
                   four folds only.
  discovery_oracle computed on the DISCOVERY side with discovery ground truth.
                   Not available to a real search; diagnostic only.

Alignment metrics (per side x split x prior x objective)
  spearman        Spearman rho(P, y)
  E_5_20          P(y in top5% | P in top20%) / 0.05
  capture_5_20    share of the y-top5% set lying in P's top 20%  (= 0.2 * E_5_20)
  S_top, S_dis    sizes: y-top5% and its members with P below the discovery median
  Neff_top20      effective number of families (exp Shannon entropy) in P's top 20%
  Neff_pool       same for the whole discovery pool (reference)
  Family = formula_anonymous + space group.
Top-k by P uses the pre-registered tie-break: sha256("<split>:<material_id>").

Usage: python build_priors.py --pool pool_v1.csv --splits splits_v1.csv
                              --reference elemental_reference_springer2005_v1.csv
                              [--universe primary] [--out-dir .]
"""

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
from math import fsum
import warnings
from collections import Counter
from importlib import metadata
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr

warnings.filterwarnings("ignore")
UNIVERSE_FLAG = {"primary": "in_primary", "sens1_no_Pm_Tc": "in_sens1_no_Pm_Tc",
                 "unscreened": "in_unscreened"}
CRYSTAL_SYSTEMS = ["Cubic", "Hexagonal", "Monoclinic", "Orthorhombic",
                   "Tetragonal", "Triclinic", "Trigonal"]
TOP_Y, TOP_P, TOP_REF = 0.05, 0.20, 0.10
SIG = 12   # prior scores are rounded to 12 significant digits (see v2.1 note)


def rsig(x):
    """Round to SIG significant digits so that mathematically equal scores are
    bit-identical on every platform and ties go to the hash rule."""
    return np.array([float(f"{v:.{SIG}g}") for v in np.asarray(x, dtype=float)])


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def tiebreak(split, mid):
    return int(hashlib.sha256(f"{split}:{mid}".encode()).hexdigest()[:15], 16)


def top_k(scores, ids, split, frac):
    """indices of the top ceil(frac*n) scores; ties broken by the frozen hash."""
    k = math.ceil(frac * len(scores))
    tb = np.array([tiebreak(split, m) for m in ids])
    order = np.lexsort((tb, -np.asarray(scores)))
    return set(order[:k].tolist())


def featurize(rows):
    from pymatgen.core import Composition
    from matminer.featurizers.composition import ElementProperty
    import pandas as pd
    ep = ElementProperty.from_preset("magpie")
    ep.set_n_jobs(1)
    comps = [Composition(r["formula_pretty"]) for r in rows]
    df = ep.featurize_dataframe(pd.DataFrame({"c": comps}), "c", ignore_errors=False)
    names = list(ep.feature_labels())
    X = df[names].to_numpy(dtype=float)
    struct = []
    for r, c in zip(rows, comps):
        rho = float(r["density_g_cm3"])
        mass_pa = c.weight / c.num_atoms                       # g/mol per atom
        vpa = mass_pa / (rho * 0.602214076)                    # Å^3 per atom
        oh = [1.0 if r["crystal_system"] == cs else 0.0 for cs in CRYSTAL_SYSTEMS]
        struct.append([rho, vpa] + oh)
    X = np.hstack([X, np.array(struct)])
    names += ["mp_density", "volume_per_atom"] + [f"cs_{c}" for c in CRYSTAL_SYSTEMS]
    if not np.all(np.isfinite(X)):
        raise SystemExit("non-finite descriptor values")
    return X, names, comps


def p_phys(comps, ref):
    from pymatgen.core import Element
    G_fix, G_mat = [], []
    for c in comps:
        fr = sorted(c.fractional_composition.get_el_amt_dict().items())
        Gi = {e: float(ref[e]["G_GPa"]) for e, _ in fr}
        rhoi = {e: float(ref[e]["density_g_cm3"]) for e, _ in fr}
        Mi = {e: float(Element(e).atomic_mass) for e, _ in fr}
        gv = fsum(x * Gi[e] for e, x in fr)
        gr = 1.0 / fsum(x / Gi[e] for e, x in fr)
        gh = 0.5 * (gv + gr)
        rho_mix = fsum(x * Mi[e] for e, x in fr) / fsum(x * Mi[e] / rhoi[e] for e, x in fr)
        G_fix.append(gh)
        G_mat.append(gh / rho_mix)
    return rsig(G_fix), rsig(G_mat)


def p_hist(Xh, yh, Xd):
    mu, sd = Xh.mean(0), Xh.std(0)
    sd[sd == 0] = 1.0
    Zh, Zd = (Xh - mu) / sd, (Xd - mu) / sd
    k = math.ceil(TOP_REF * len(yh))
    ref = Zh[np.argsort(-yh, kind="stable")[:k]]
    best = np.full(len(Zd), np.inf)
    for r in ref:                       # direct differences: no BLAS matrix product
        best = np.minimum(best, ((Zd - r) ** 2).sum(1))
    return rsig(-np.sqrt(best))


SPLIT_UNIT = {"random": "group_id", "chemsys": "chemsys"}
N_FOLDS = 5


def p_hist_cv(Xh, yh, units):
    """Out-of-fold P_H scores within history (grouped folds)."""
    fold = np.array([int(hashlib.sha256(u.encode()).hexdigest()[:8], 16) % N_FOLDS
                     for u in units])
    out = np.empty(len(yh))
    for f in range(N_FOLDS):
        te, tr = fold == f, fold != f
        out[te] = p_hist(Xh[tr], yh[tr], Xh[te])
    return out


def metrics(p, y, ids, split, fam):
    ktop = math.ceil(TOP_Y * len(y))
    top_y = set(np.argsort(-y, kind="stable")[:ktop].tolist())
    top_p = top_k(p, ids, split, TOP_P)
    hit = len(top_y & top_p)
    pr = rankdata(p) / len(p)
    return {"n": len(y),
            "spearman": round(float(spearmanr(p, y).statistic), 4),
            "E_5_20": round(hit / len(top_p) / TOP_Y, 4),
            "capture_5_20": round(hit / ktop, 4),
            "S_top": ktop,
            "S_dis_poolmedian": int(sum(pr[i] <= 0.5 for i in top_y)),
            "Neff_top20": round(neff(fam[sorted(top_p)]), 3),
            "Neff_pool": round(neff(fam), 3)}


def neff(labels):
    c = np.array(list(Counter(labels).values()), dtype=float)
    p = c / c.sum()
    return float(np.exp(-(p * np.log(p)).sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)
    ap.add_argument("--splits", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--universe", default="primary", choices=list(UNIVERSE_FLAG))
    ap.add_argument("--out-dir", default=".")
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    rows = [r for r in csv.DictReader(open(a.pool)) if r[UNIVERSE_FLAG[a.universe]] == "True"]
    rows.sort(key=lambda r: r["material_id"])
    ids = [r["material_id"] for r in rows]
    ref = {r["element"]: r for r in csv.DictReader(open(a.reference))}
    if a.universe != "unscreened":
        bad = {e for r in rows for e in r["elements"].split(";") if ref[e]["eligible"] != "True"}
        assert not bad, f"ineligible elements in universe: {bad}"
    splits = {r["material_id"]: r for r in csv.DictReader(open(a.splits))}
    split_names = [k for k in next(iter(splits.values())) if k != "material_id"]

    G = np.array([float(r["G_vrh_GPa"]) for r in rows])
    K = np.array([float(r["K_vrh_GPa"]) for r in rows])
    rho = np.array([float(r["density_g_cm3"]) for r in rows])
    assert np.all(K > 0) and np.all(G > 0) and np.all(rho > 0)
    E = 9 * K * G / (3 * K + G)
    Y = {"Y1_G": G, "Y2_E_over_rho": E / rho}
    fam = np.array([f'{r["formula_anonymous"]}|{r["spacegroup_number"]}' for r in rows])

    print(f"Featurizing {len(rows)} structures ({a.universe})...")
    X, fnames, comps = featurize(rows)
    keep = X.std(0) > 0
    dropped = [n for n, k in zip(fnames, keep) if not k]
    X, fnames = X[:, keep], [n for n, k in zip(fnames, keep) if k]
    with open(out / f"features_v2_{a.universe}.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["material_id"] + fnames)
        for mid, x in zip(ids, X):
            w.writerow([mid] + [repr(float(v)) for v in x])

    have_phys = a.universe != "unscreened"
    if have_phys:
        Pphys_fix, Pphys_mat = p_phys(comps, ref)

    prior_rows, align_rows = [], []
    Y2 = "Y2_E_over_rho"
    for sp in split_names:
        stype = sp.rsplit("_r", 1)[0]
        side = np.array([splits[m][sp] for m in ids])
        H, D = np.where(side == "H")[0], np.where(side == "D")[0]
        assert len(set(H) & set(D)) == 0 and len(H) + len(D) == len(ids)
        units_h = [rows[i][SPLIT_UNIT[stype]] for i in H]
        # discovery-side scores: P_H uses history y only
        P_D = {"P_H_fixed": p_hist(X[H], Y["Y1_G"][H], X[D]),
               "P_H_matched_Y2": p_hist(X[H], Y[Y2][H], X[D])}
        # history-side out-of-fold scores (pre-search alignment)
        P_Hs = {"P_H_fixed": p_hist_cv(X[H], Y["Y1_G"][H], units_h),
                "P_H_matched_Y2": p_hist_cv(X[H], Y[Y2][H], units_h)}
        if have_phys:
            P_D["P_phys_fixed"], P_Hs["P_phys_fixed"] = Pphys_fix[D], Pphys_fix[H]
            P_D["P_phys_matched_Y2"], P_Hs["P_phys_matched_Y2"] = Pphys_mat[D], Pphys_mat[H]
        for j, i in enumerate(D):
            prior_rows.append({"split": sp, "material_id": ids[i],
                               **{k: repr(float(v[j])) for k, v in P_D.items()}})
        for side_name, idx, PP in (("history_cv", H, P_Hs), ("discovery_oracle", D, P_D)):
            sids = [ids[i] for i in idx]
            for yname, yall in Y.items():
                for pname, p in PP.items():
                    if pname.endswith("_matched_Y2") and yname == "Y1_G":
                        continue
                    align_rows.append({
                        "universe": a.universe, "side": side_name, "split": sp,
                        "split_type": stype, "objective": yname, "prior": pname,
                        "variant": "fixed" if (pname.endswith("_fixed") and yname != "Y1_G") else "matched",
                        **metrics(p, yall[idx], sids, sp, fam[idx])})

    with open(out / f"priors_v2_{a.universe}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(prior_rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(prior_rows)
    with open(out / f"alignment_v2_{a.universe}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(align_rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(align_rows)

    vers = {p: metadata.version(p) for p in ["numpy", "scipy", "pymatgen", "matminer", "pandas"]}
    manifest = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "universe": a.universe,
        "inputs": {"pool": sha256(a.pool), "splits": sha256(a.splits), "reference": sha256(a.reference)},
        "script_sha256": sha256(__file__),
        "versions": vers,
        "n_features": len(fnames), "dropped_constant_features": dropped,
        "outputs": {n: sha256(out / n) for n in [f"features_v2_{a.universe}.csv",
                                                  f"priors_v2_{a.universe}.csv",
                                                  f"alignment_v2_{a.universe}.csv"]},
    }
    (out / f"priors_v2_{a.universe}.manifest.json").write_text(json.dumps(manifest, indent=2))

    # summary: mean over replicates, history_cv (pre-search) vs discovery_oracle
    print(f"\n{len(fnames)} descriptors (dropped constant: {dropped})")
    print(f"\n{'objective':<15}{'prior':<19}{'split':<9}{'rho hist':>9}{'rho disc':>9}"
          f"{'E hist':>8}{'E disc':>8}")
    keys = sorted({(r["objective"], r["prior"], r["split_type"]) for r in align_rows})
    for k in keys:
        m = lambda side, f: np.mean([r[f] for r in align_rows if r["side"] == side and
                                      (r["objective"], r["prior"], r["split_type"]) == k])
        print(f"{k[0]:<15}{k[1]:<19}{k[2]:<9}{m('history_cv','spearman'):>9.3f}"
              f"{m('discovery_oracle','spearman'):>9.3f}{m('history_cv','E_5_20'):>8.2f}"
              f"{m('discovery_oracle','E_5_20'):>8.2f}")
    print("\noutputs:", json.dumps(manifest["outputs"], indent=1))


if __name__ == "__main__":
    main()

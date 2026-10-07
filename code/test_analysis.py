"""
test_analysis.py — tests of analyze_runs.py v1.1 on a SYNTHETIC universe
========================================================================
No Materials Project value is read. A synthetic pool (same file formats as the
real inputs) and fabricated run records are written to a temporary directory.
Every reported table is compared with an independent recomputation (plain
Python / pandas, a different code path from analyze_runs.py).

  python test_analysis.py            (about 6 minutes)
"""

import csv
import hashlib
import json
import math
import shutil
import tempfile
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import analyze_runs as ar
import run_search as rs

UNIVERSES = ar.UNIVERSES
SPLITS = ar.all_splits()
FAMS = [("ABC2", 225), ("AB", 221), ("AB3", 139), ("AB2", 227), ("ABC", 216), ("A2B3", 12)]
FLAG = {"primary": "in_primary", "sens1_no_Pm_Tc": "in_sens1_no_Pm_Tc", "unscreened": "in_unscreened"}
ENV = {"python": "3", "numpy": "x", "scipy": "x", "scikit-learn": ar.SKLEARN_VERSION,
       "threads": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}}


def rho_synth(k, j):
    return 0.4 + 0.05 * k + 0.001 * j * j          # not symmetric over splits: mean != median


def build(root, n=640, runs=True):
    """Synthetic data directory + fabricated confirmatory run records and their manifest."""
    rng = np.random.default_rng(1)
    d = root / "data"
    d.mkdir()
    ids = [f"syn-{i:05d}" for i in range(n)]
    G = np.round(rng.gamma(4.0, 20.0, n) / 4.0) * 4.0 + 4.0     # many ties in y (Y1)
    K = np.round(G * rng.uniform(1.2, 3.0, n), 3)
    rho = rng.uniform(2.0, 12.0, n)
    prim = rng.random(n) < 0.8
    sens = prim & (rng.random(n) < 0.95)
    with open(d / "pool_v1.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["material_id", "G_vrh_GPa", "K_vrh_GPa", "density_g_cm3", "spacegroup_number",
                    "formula_anonymous", "in_primary", "in_sens1_no_Pm_Tc", "in_unscreened"])
        rare = set(np.argsort(-(G + 1e6 * prim))[:3].tolist())      # a family with only 3 top materials
        for i, m in enumerate(ids):
            fa, sg = FAMS[int(rng.integers(0, 5))] if G[i] > np.median(G) else FAMS[i % 5]
            if i in rare:
                fa, sg = FAMS[5]
            w.writerow([m, G[i], K[i], repr(float(rho[i])), sg, fa, str(bool(prim[i])),
                        str(bool(sens[i])), "True"])
    side = {sp: np.where(rng.random(n) < 0.7, "D", "H") for sp in SPLITS}
    with open(d / "splits_v2.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["material_id"] + SPLITS)
        for i, m in enumerate(ids):
            w.writerow([m] + [side[sp][i] for sp in SPLITS])
    pri = {}
    for sp in SPLITS:
        for name, noise in (("P_H_fixed", 1.0), ("P_H_matched_Y2", 1.5),
                            ("P_phys_fixed", 0.5), ("P_phys_matched_Y2", 2.0)):
            pri[(sp, name)] = np.round((G - G.mean()) / G.std() + noise * rng.normal(size=n), 1)   # ties
    for u, mask in zip(UNIVERSES, (prim, sens, np.ones(n, bool))):
        with open(d / f"priors_v3_{u}.csv", "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["split", "material_id"] + list(rs.PRIORS[1:]))
            for sp in SPLITS:
                for i, m in enumerate(ids):
                    if mask[i] and side[sp][i] == "D":
                        w.writerow([sp, m] + [pri[(sp, p)][i] for p in rs.PRIORS[1:]])
        (d / f"features_v2_{u}.csv").write_text(f"material_id,f0\n# synthetic {u}\n")
    with open(d / "alignment_v3_primary.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["universe", "side", "split", "objective", "prior", "spearman"])
        for k, (obj, prior) in enumerate(ar.CELLS["primary"]):
            for j, sp in enumerate(SPLITS):
                w.writerow(["primary", "history_cv", sp, ar.ALIGN_OBJ[obj], prior, rho_synth(k, j)])

    ctx, hashes = ar.load_data(d, UNIVERSES, {"Y1": 0.05, "Y2": 0.05})
    rdir = root / "runs"
    rdir.mkdir()
    if runs:
        man = {}
        for key in ar.expected_runs(UNIVERSES):
            f = write_run(rdir, ctx, hashes, key)
            man[f.name] = json.loads(f.read_text())["run_hash"]
        (rdir / "MANIFEST_confirmatory.json").write_text(json.dumps({"runs": man}))
    return d, rdir, ctx, hashes


def write_run(runs, ctx, hashes, key, **override):
    u, sp, obj, prior, beta, seed = key
    c = ctx[(u, sp)]
    init = rs.init_design(c.n, sp, seed)
    r = np.random.default_rng(int(hashlib.sha256(repr(key).encode()).hexdigest()[:12], 16))
    rest = np.array([i for i in range(c.n) if i not in set(init)])
    y = c.y[obj]
    s = (1 - beta) * 1.5 * (y[rest] - y.mean()) / y.std()
    if beta > 0:
        base, sign = (prior.split(":")[1], -1) if prior.startswith("P_wrong:") else (prior, 1)
        p = c.prior[base]
        s = s + min(1.0, 4 * beta) * 3.0 * sign * (p[rest] - p.mean()) / p.std()
    w = np.exp(s - s.max())
    traj_idx = [int(i) for i in r.choice(rest, 40, replace=False, p=w / w.sum())]
    I, T = [c.ids[i] for i in init], [c.ids[i] for i in traj_idx]
    steps = []
    for t, m in enumerate(T, 1):
        nc = c.n - 10 - (t - 1)
        steps.append({"t": t, "material_id": m, "n_cand": nc,
                      "ei_top_ties": None if beta == 1 else (int(0.7 * nc) if t % 10 == 0 else 1),
                      "kernel": None if beta == 1 else [1.0, 3.0 if t % 4 == 0 else (3.5 if t % 4 == 1 else 7.5)]})
    rec = {"engine_version": rs.ENGINE_VERSION,
           "run": {"universe": u, "split": sp, "side": "D", "objective": obj, "prior": prior,
                   "beta": beta, "seed": seed, "pool_size": c.n},
           "registered": ar.REGISTRATION,
           "inputs": {"pool": hashes["pool"], "splits": hashes["splits"],
                      "priors": hashes[f"priors:{u}"], "features": hashes[f"features:{u}"]},
           "script_sha256": ar.ENGINE_SHA256, "env": ENV, "init": I, "trajectory": T,
           "y_init": [float(c.y[obj][i]) for i in init], "y_trajectory": [float(c.y[obj][i]) for i in traj_idx],
           "steps": steps, "run_hash": hashlib.sha256(json.dumps([I, T]).encode()).hexdigest()}
    rec.update(override)
    f = runs / f"{u}_{sp}_D_{obj}_{prior.replace(':', '-')}_b{beta}_s{seed}.json"
    f.write_text(json.dumps(rec))
    return f


# ───────────────────────────── independent re-implementations ────────────────

def full_key(sp, m):
    return hashlib.sha256(f"{sp}:{m}".encode()).hexdigest()        # the whole digest, not rs.hkey


class Brute:
    """Pool-level sets of one (universe, split, objective, prior), computed with dicts and sorts."""

    def __init__(self, pool_rows, split_rows, prior_val, u, sp, obj, frac=Fraction(5, 100)):
        self.D = D = [m for m in pool_rows if pool_rows[m][FLAG[u]] == "True" and split_rows[m][sp] == "D"]

        def yval(m):
            g, k, rho = (float(pool_rows[m][c]) for c in ("G_vrh_GPa", "K_vrh_GPa", "density_g_cm3"))
            return g if obj == "Y1" else 9 * k * g / (3 * k + g) / rho
        self.y = y = {m: yval(m) for m in D}
        ceil = lambda fr: -((-fr.numerator) // fr.denominator)
        ranked = sorted(D, key=lambda m: (-y[m], full_key(sp, m)))
        self.top = ranked[:ceil(frac * len(D))]
        self.top1 = ranked[:ceil(Fraction(1, 100) * len(D))]
        self.sdis = sorted(self.top, key=lambda m: (prior_val[m], full_key(sp, m)))[:ceil(Fraction(len(self.top), 2))]
        vals = sorted(prior_val[m] for m in D)
        self.pct = {}
        for m in self.top:
            below = sum(1 for v in vals if v < prior_val[m])
            equal = sum(1 for v in vals if v == prior_val[m])
            self.pct[m] = (below + (equal + 1) / 2 - 0.5) / len(D)
        self.blind = [m for m in self.top if self.pct[m] < 0.80]
        self.fam = {m: f"{pool_rows[m]['formula_anonymous']}_{pool_rows[m]['spacegroup_number']}" for m in self.top}

    def outcomes(self, rec, tmax=40):
        tr, y, D = rec["trajectory"][:tmax], self.y, self.D
        sdis = [m for m in self.sdis if m not in rec["init"]]
        blind = [m for m in self.blind if m not in rec["init"]]
        best, q = max(y[m] for m in rec["init"]), []
        for m in tr:
            best = max(best, y[m])
            q.append(sum(1 for k in D if y[k] <= best) / len(D))
        return {"T99": min([t for t, m in enumerate(tr, 1) if m in self.top1], default=tmax + 1),
                "T95": min([t for t, m in enumerate(tr, 1) if m in self.top], default=tmax + 1),
                "N5": len([m for m in tr if m in self.top]), "AUC": sum(q) / len(q),
                "S_dis_run": len(sdis), "dis_hits": len([m for m in tr if m in sdis]),
                "R_dis": len([m for m in tr if m in sdis]) / len(sdis),
                "S_blind_run": len(blind), "blind_hits": len([m for m in tr if m in blind])}


def welch(df, col):
    """The stratified split-level test, second implementation (pandas + the general
    Welch-Satterthwaite formula for a linear combination of independent means)."""
    m = df.groupby(["type", "split"])[col].mean().reset_index()
    parts = [m[m["type"] == t][col].to_numpy() for t in ar.SPLIT_TYPES]
    terms = [0.25 * np.var(x, ddof=1) / len(x) for x in parts]
    est = float(sum(0.5 * np.mean(x) for x in parts))
    se = math.sqrt(sum(terms))
    den = sum(t * t / (len(x) - 1) for t, x in zip(terms, parts))
    return est, se, (sum(terms) ** 2 / den if den > 0 else 18.0)


def bounds(df, col):
    est, se, dof = welch(df, col)
    q = lambda lev: stats.t.ppf(lev, dof)
    return {"est": est, "se": se, "df": dof, "lo95": est - q(0.95) * se, "hi95": est + q(0.95) * se,
            "hi_adj": est + q(1 - 0.05 / 6) * se, "lo2": est - q(0.975) * se, "hi2": est + q(0.975) * se,
            "p": float(stats.t.sf(est / se, dof)) if se > 0 else (0.0 if est > 0 else 1.0)}


def paired(run_csv, universe, eff):
    """Paired run table {(cell, beta): frame with L, gain, margin, base} from run_outcomes.csv."""
    df = pd.read_csv(run_csv)
    df = df[df["universe"] == universe].copy()
    df["type"] = df["split"].str.replace(r"_r\d+$", "", regex=True)
    out = {}
    for cell, g in df.groupby("cell"):
        e = eff[cell[:2]]
        base = g[g["beta"] == 0].set_index(["split", "seed"])
        for b in sorted(set(g["beta"]) - {0.0}):
            a = g[g["beta"] == b].set_index(["split", "seed"]).join(base, rsuffix="_0").reset_index()
            a["L"] = a["R_dis_0"] - a["R_dis"]
            a["gain"] = (a["T99_0"] - a["T99"]) if e == "T99" else (a["N5"] - a["N5_0"])
            a["margin"] = a["L"] - 0.20 * a["R_dis_0"]
            a["base"] = a["R_dis_0"] - a["chance_recall_0"]
            out[(cell, float(b))] = a
    return out


def cube_of(frame, col):
    piv = frame.pivot_table(index="split", columns="seed", values=col)
    return np.array([[[piv.loc[ar.split_name(st, r), s] for s in ar.SEEDS] for r in ar.REPS] for st in ar.SPLIT_TYPES])


def expect_fail(fn, text):
    try:
        fn()
    except SystemExit as e:
        assert text in str(e), f"wrong failure: {e}"
        return
    raise AssertionError(f"validation did not fail ({text})")


def close(a, b, tol=1e-8):
    return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(b)))


def main():
    root = Path(tempfile.mkdtemp(prefix="p2_analysis_test_"))
    try:
        d, runs, ctx, hashes = build(root)
        pool_rows = {r["material_id"]: r for r in csv.DictReader(open(d / "pool_v1.csv"))}
        split_rows = {r["material_id"]: r for r in csv.DictReader(open(d / "splits_v2.csv"))}
        prior_tab = {u: {(r["split"], r["material_id"]): r for r in csv.DictReader(open(d / f"priors_v3_{u}.csv"))}
                     for u in UNIVERSES}
        _bc = {}

        def brute(u, sp, obj, prior, frac=Fraction(5, 100)):
            k = (u, sp, obj, prior, frac)
            if k not in _bc:
                pv = {m: float(prior_tab[u][(sp, m)][prior]) for m in ctx[(u, sp)].ids}
                _bc[k] = Brute(pool_rows, split_rows, pv, u, sp, obj, frac)
            return _bc[k]

        loaded = ar.load_runs(runs, ctx, hashes, UNIVERSES)
        assert len(loaded) == 7020 and len(ar.expected_runs(("primary",))) == 3000
        assert len(ar.expected_runs(("sens1_no_Pm_Tc",))) == 2640 and len(ar.expected_runs(("unscreened",))) == 1380
        print("ok  1  7,020 fabricated records load and validate (3,000 + 2,640 + 1,380)")

        n_cmp = 0
        for u, obj, prior in (("primary", "Y1", "P_phys_fixed"), ("unscreened", "Y2", "P_H_matched_Y2")):
            for sp in SPLITS[::3]:
                c, bt = ctx[(u, sp)], brute(u, sp, obj, prior)
                for seed in ar.SEEDS:
                    for b in ar.BETAS:
                        rec = loaded[(u, sp, obj, "P0" if b == 0 else prior, b, seed)]
                        got = ar.run_outcomes(c, obj, rec["init"], rec["trajectory"],
                                              c.s_dis(obj, prior), c.s_blind(obj, prior))
                        for k, v in bt.outcomes(rec).items():
                            assert abs(got[k] - v) < 1e-12, (k, got[k], v)
                        n_cmp += 1
        print(f"ok  2  outcomes match an independent implementation ({n_cmp} runs), incl. S_blind")

        c = ctx[("primary", "random_r0")]
        init = set(rs.init_design(c.n, "random_r0", 100))
        yo = np.lexsort((c.hk, -c.y["Y1"]))
        order = [int(i) for i in yo if int(i) not in init]
        top1 = [i for i in order if i in c.top1["Y1"]]
        non_top = [i for i in order if i not in set(c.top["Y1"])]
        ids_of = lambda idx: [c.ids[i] for i in idx]
        o = ar.run_outcomes(c, "Y1", ids_of(init), ids_of(non_top[:6] + [top1[-1]] + non_top[6:39]))
        assert o["T99"] == 7 and o["T95"] == 7 and o["N5"] == 1
        o = ar.run_outcomes(c, "Y1", ids_of(init), ids_of(non_top[:40]))
        assert o["T99"] == 41 and o["T95"] == 41 and o["N5"] == 0
        seed_ok = next(s for s in range(200, 400) if int(yo[0]) not in rs.init_design(c.n, "random_r0", s))
        i2 = rs.init_design(c.n, "random_r0", seed_ok)
        rest2 = [int(i) for i in yo[::-1] if int(i) not in i2 and int(i) != int(yo[0])]
        assert ar.run_outcomes(c, "Y1", ids_of(i2), ids_of([int(yo[0])] + rest2[:39]))["AUC"] == 1.0
        sd = c.s_dis("Y1", "P_phys_fixed")
        assert len(sd) == math.ceil(len(c.top["Y1"]) / 2) and sd <= set(c.top["Y1"])
        s_run = [i for i in sd if i not in init]
        o = ar.run_outcomes(c, "Y1", ids_of(init), ids_of(s_run[:3] + [i for i in non_top if i not in s_run][:37]), sd)
        assert o["dis_hits"] == 3 and abs(o["R_dis"] - 3 / len(s_run)) < 1e-15
        assert abs(c.chance - 40 / (c.n - 10)) < 1e-15
        print("ok  3  planted trajectories: T99, T95, N5, AUC, S_dis, R_dis and chance recall as constructed")

        rng = np.random.default_rng(5)
        V = rng.normal(0.3, 1.0, size=(2, 10, 3)) + rng.normal(size=(2, 10, 1)) * np.array([1.0, 3.0])[:, None, None]
        fr_ = pd.DataFrame([{"type": ar.SPLIT_TYPES[k], "split": f"{k}_{j}", "v": V[k, j, s]}
                            for k in range(2) for j in range(10) for s in range(3)])
        est, se, dof = ar.split_t(V)
        e2, s2, d2 = welch(fr_, "v")
        sm = V.mean(axis=2)
        assert abs(est - e2) < 1e-12 and abs(se - s2) < 1e-12 and abs(dof - d2) < 1e-9 and 9 < dof < 18
        assert abs(dof - stats.ttest_ind(sm[0], sm[1], equal_var=False).df) < 1e-9
        assert abs(ar.p_greater(V) - stats.t.sf(e2 / s2, d2)) < 1e-12
        lo, hi = ar.t_interval(V)
        assert abs(lo - (e2 - stats.t.ppf(0.975, d2) * s2)) < 1e-12 and abs((lo + hi) / 2 - est) < 1e-12
        l1, h1b = ar.t_interval(V, one_sided=True)
        assert abs(h1b - (e2 + stats.t.ppf(0.95, d2) * s2)) < 1e-12 and lo < l1 < est < h1b < hi
        assert abs(ar.t_interval(V, level=ar.SAFE_LEVEL, one_sided=True)[1]
                   - (e2 + stats.t.ppf(1 - 0.05 / 6, d2) * s2)) < 1e-12 and ar.SAFE_LEVEL == 1 - 0.05 / 6
        V1 = V.copy(); V1[0] += 5.0                     # a pure split-type shift must not change SE or df
        assert abs(ar.split_t(V1)[1] - se) < 1e-12 and abs(ar.split_t(V1)[0] - est - 2.5) < 1e-12
        assert abs(ar.split_t(V1)[2] - dof) < 1e-9
        Ve = np.concatenate([V[:1], V[:1]])             # equal variances -> 18; one type constant -> 9
        Vc = np.concatenate([V[:1], np.zeros((1, 10, 3))])
        assert abs(ar.split_t(Ve)[2] - 18) < 1e-9 and abs(ar.split_t(Vc)[2] - 9) < 1e-9 and ar.DF_MAX == 18
        assert ar.split_t(V[:, 3:, :])[2] <= 12 + 1e-9
        assert ar.p_greater(np.full((2, 10, 3), 0.2)) == 0.0 and ar.p_greater(np.zeros((2, 10, 3))) == 1.0
        idx = ar.boot_indices(300)
        assert idx[0].shape == (300, 2, 10) and idx[1].shape == (300, 2, 10, 3)
        slow = [np.mean([V[k, idx[0][b, k, j], idx[1][b, k, j, i]] for k in range(2) for j in range(10)
                         for i in range(3)]) for b in range(300)]
        assert np.allclose(ar.boot_means(V, idx), slow, atol=1e-13)
        assert np.array_equal(ar.boot_indices(300)[0], idx[0])
        assert ar.boot_p_v1_0(np.zeros(999)) == 1.0 and ar.boot_p_v1_0(np.array([-1.0, 0.0, 1.0, 2.0])) == 3 / 5
        blo, bhi = ar.boot_ci(V, idx)
        assert abs(blo - np.percentile(slow, 2.5)) < 1e-12 and abs(bhi - np.percentile(slow, 97.5)) < 1e-12
        print("ok  4  split-level t-test (Satterthwaite df) equals a second implementation and SciPy's Welch df;")
        print("       bootstrap equals a loop implementation")

        assert np.allclose(ar.holm([0.01, 0.04, 0.03, 0.005]), [0.03, 0.06, 0.06, 0.02])
        assert np.allclose(ar.holm([0.5, 0.5, 0.5]), [1.0, 1.0, 1.0])
        dec = ar.h1_decisions([{"p": v} for v in (0.03, 0.001, 0.2, 0.004, 0.02)])
        assert [x["reject_H0"] for x in dec] == [False, True, False, True, False]
        assert np.allclose([x["p_holm"] for x in dec], [0.06, 0.005, 0.2, 0.016, 0.06])
        dec = ar.h1_decisions([{"p": v} for v in (0.02, 0.2)])           # adjusted 0.04: between alpha/2 and alpha
        assert dec[0]["reject_H0"] is True and abs(dec[0]["p_holm"] - 0.04) < 1e-12
        one = np.ones((2, 10, 3))
        z = 0.01 * rng.normal(size=(2, 10, 1)) * one
        z = z - z.mean()
        _, se_z, df_z = ar.split_t(z)
        q95, qadj = stats.t.ppf(0.95, df_z), stats.t.ppf(ar.SAFE_LEVEL, df_z)
        assert q95 < 2.0 < qadj - 0.2
        R0 = 0.10 * one
        gain = {b: 2.0 * one + 100 * z for b in ar.H2_BETAS}
        gain[0.25] = -1.0 * one + 100 * z                                  # slower: not beneficial
        gain[0.1] = (stats.t.ppf(0.90, df_z) * 100 * se_z) * one + 100 * z  # positive mean, lower bound < 0
        loss = {b: 0.02 * one - (qadj + 0.2) * se_z + z for b in ar.H2_BETAS}    # safe at the adjusted level
        loss[0.5] = 0.02 * one - 2.0 * se_z + z          # safe at 95% only: not confirmatory-safe
        loss[0.75] = 0.02 * one + 2.5 * se_z + z         # demonstrably unsafe
        loss[0.05] = 0.02 * one + z                      # exactly at the margin: inconclusive
        ok, lo_b, rows = ar.h2_rule(R0 - 0.016, gain, loss, R0)
        got = {r["beta"]: (r["beneficial"], r["safe"], r["beneficial_and_safe"], r["state_unadjusted"]) for r in rows}
        assert ok and got[0.02] == (True, True, True, "safe")
        assert got[0.5] == (True, False, False, "safe") and got[0.75] == (True, False, False, "unsafe")
        assert got[0.05] == (True, False, False, "inconclusive") and got[0.25] == (False, True, False, "safe")
        assert got[0.1] == (False, True, False, "safe")
        r02 = [r for r in rows if r["beta"] == 0.02][0]
        assert r02["safe_at_25pct"] and not r02["safe_at_10pct"] and abs(r02["relative_loss"] - loss[0.02].mean() / 0.10) < 1e-12
        ok2, _, rows2 = ar.h2_rule(R0 - 0.2, gain, loss, R0)              # baseline below chance
        assert not ok2 and not any(r["beneficial_and_safe"] for r in rows2) and rows2[0]["safe"]
        ok3, _, _ = ar.h2_rule(0.0 * one + z, gain, loss, R0)             # baseline exactly at chance
        assert not ok3
        print("ok  5  Holm and the H1 decision; H2 rule: gate, benefit, adjusted safety level, three-way state")

        out1, out2 = root / "res1", root / "res2"
        s1 = ar.analyse(d, runs, out1, b=2000)
        s2_ = ar.analyse(d, runs, out2, b=2000)
        assert s1["outputs"] == s2_["outputs"], "analysis is not deterministic"
        rd = lambda f, o=out1: list(csv.DictReader(open(o / f)))
        ga, gb = rd("glmm_fixed_effects.csv"), rd("glmm_fixed_effects.csv", out2)
        assert len(ga) == 6 * 16 and all(abs(float(x["post_mean"]) - float(y["post_mean"])) <= 0.002
                                         for x, y in zip(ga, gb))
        g1 = {r["term"]: float(r["post_mean"]) for r in ga if r["cell"] == "Y1:P_phys_fixed"}
        assert g1["C(beta, Treatment('0'))[T.1]:pctl"] > 1.0 > 0.0 > g1["C(beta, Treatment('0'))[T.1]"]
        fr, h1, h2, h2d = rd("frontier.csv"), rd("h1_coverage_loss.csv"), rd("h2_beneficial_and_safe.csv"), rd("h2_detail.csv")
        assert len(fr) == 15 * 8 and len(h1) == 15 * 2 and len(h2) == 15 and len(h2d) == 15 * 6
        assert len(rd("run_outcomes.csv")) == 15 * 8 * 60 and len(rd("recall_by_prior_percentile.csv")) == 15 * 8 * 5
        assert len(rd("new_splits_only.csv")) == 15 * 7 and len(rd("budget_sensitivity.csv")) == 15 * 7 * 3

        # (a) the run table against the independent implementation, for all 15 cells
        ro = pd.read_csv(out1 / "run_outcomes.csv")
        n_rows = 0
        for u in UNIVERSES:
            for obj, prior in ar.CELLS[u]:
                for sp in (SPLITS[1], SPLITS[14]):
                    bt = brute(u, sp, obj, prior)
                    for seed, b in ((100, 0.0), (101, 0.05), (102, 0.5), (100, 1.0)):
                        rec = loaded[(u, sp, obj, "P0" if b == 0 else prior, b, seed)]
                        row = ro[(ro.universe == u) & (ro.cell == f"{obj}:{prior}") & (ro.split == sp)
                                 & (ro.seed == seed) & (ro.beta == b)]
                        assert len(row) == 1
                        for k, v in bt.outcomes(rec).items():
                            assert close(row.iloc[0][k], v), (u, obj, prior, sp, seed, b, k)
                        assert close(row.iloc[0]["chance_recall"], 40 / (len(bt.D) - 10))
                        n_rows += 1

        # (b) H1, H2 and the frontier from the run table, with the second implementation of the test
        eff0 = {"Y1": "T99", "Y2": "T99"}
        for u in UNIVERSES:
            P = paired(out1 / "run_outcomes.csv", u, eff0)
            cells = sorted({k[0] for k in P})
            if u == "primary":
                conf = [r for r in h1 if r["universe"] == u]
                assert len(conf) == 12 and {r["beta"] for r in conf} == {"0.25", "0.5"} and \
                    all(r["role"] == "confirmatory" for r in conf)
                ps = [bounds(P[(r["cell"], float(r["beta"]))], "L")["p"] for r in conf]
                for r, a in zip(conf, ar.holm(ps)):
                    bd = bounds(P[(r["cell"], float(r["beta"]))], "L")
                    assert close(r["mean_L"], bd["est"]) and close(r["se"], bd["se"]) and close(r["df"], bd["df"])
                    assert close(r["p"], bd["p"]) and close(r["ci_lo"], bd["lo2"]) and close(r["ci_hi"], bd["hi2"])
                    assert close(r["p_holm"], a) and r["reject_H0"] == str(bool(a <= 0.05))
                    bm = ar.boot_means(cube_of(P[(r["cell"], float(r["beta"]))], "L"), ar.boot_indices(2000))
                    assert close(r["boot_ci_lo"], np.percentile(bm, 2.5)) and close(r["boot_ci_hi"], np.percentile(bm, 97.5))
                    assert close(r["boot_p_v1_0"], (1 + np.sum(bm <= 0)) / 2001)
                assert any(r["reject_H0"] == "True" for r in conf)
            else:
                assert all(r["p_holm"] == "" and r["role"].startswith("descriptive") for r in h1 if r["universe"] == u)
            for r in [x for x in h2d if x["universe"] == u]:
                fm = P[(r["cell"], float(r["beta"]))]
                g, mg, ba = bounds(fm, "gain"), bounds(fm, "margin"), bounds(fm, "base")
                assert close(r["gain_lower95"], g["lo95"]) and close(r["margin_upper_adjusted"], mg["hi_adj"])
                assert close(r["margin_upper95"], mg["hi95"]) and close(r["margin_lower95"], mg["lo95"])
                assert close(r["relative_loss"], fm["L"].mean() / fm["R_dis_0"].mean())
                assert (r["beneficial"], r["safe"]) == (str(bool(g["lo95"] > 0)), str(bool(mg["hi_adj"] < 0)))
                assert r["beneficial_and_safe"] == str(bool(ba["lo95"] > 0 and g["lo95"] > 0 and mg["hi_adj"] < 0))
            for r in [x for x in h2 if x["universe"] == u]:
                ba = bounds(P[(r["cell"], 0.02)], "base")
                assert close(r["baseline_minus_chance_lower95"], ba["lo95"]) and r["applicable"] == str(bool(ba["lo95"] > 0))
                want = [b for b in ar.H2_BETAS if bounds(P[(r["cell"], b)], "gain")["lo95"] > 0
                        and bounds(P[(r["cell"], b)], "margin")["hi_adj"] < 0 and ba["lo95"] > 0]
                assert json.loads(r["beneficial_and_safe_betas"]) == want
                assert close(r["baseline_R_dis"], P[(r["cell"], 0.02)]["R_dis_0"].mean())
            for r in [x for x in fr if x["universe"] == u and x["beta"] != "0"]:
                fm = P[(r["cell"], float(r["beta"]))]
                for name in ("gain", "L", "N5", "AUC", "R_dis"):
                    bd = bounds(fm, name)
                    assert close(r[name], bd["est"]) and close(r[f"{name}_lo"], bd["lo2"]) and close(r[f"{name}_hi"], bd["hi2"]), (name, r)
                assert close(r["R_blind_pooled"], fm["blind_hits"].sum() / fm["S_blind_run"].sum())
            r0 = [x for x in fr if x["universe"] == u and x["beta"] == "0"]
            assert len(r0) == len(cells) and all(float(x["gain"]) == 0 and float(x["L"]) == 0 for x in r0)
            bm = ar.boot_means(cube_of(P[(cells[0], 0.5)], "AUC"), ar.boot_indices(2000))
            x = [x for x in fr if x["universe"] == u and x["cell"] == cells[0] and x["beta"] == "0.5"][0]
            assert close(x["AUC_boot_lo"], np.percentile(bm, 2.5)) and close(x["AUC_boot_hi"], np.percentile(bm, 97.5))
            for r in [x for x in rd("new_splits_only.csv") if x["universe"] == u]:      # splits r3-r9 only
                fm = P[(r["cell"], float(r["beta"]))]
                fm = fm[fm["split"].str.extract(r"_r(\d+)$")[0].astype(int) >= 3]
                bl, bg = bounds(fm, "L"), bounds(fm, "gain")
                assert r["splits"] == "14" and close(r["L"], bl["est"]) and close(r["L_lo"], bl["lo2"])
                assert close(r["p_L_greater_0"], bl["p"]) and close(r["gain_hi"], bg["hi2"]) and close(r["df"], bl["df"])
        states = {(x["beneficial"], x["safe"]) for x in h2d}
        assert len(states) >= 3, f"the synthetic data exercise too few H2 states: {states}"

        # (c) H3, manipulation check
        P = paired(out1 / "run_outcomes.csv", "primary", eff0)
        h3 = rd("h3_alignment.csv")
        rho = [np.mean([rho_synth(k, j) for j in range(20)]) for k in range(6)]
        assert len(h3) == 6 and all(close(r["history_rho"], x) for r, x in zip(h3, rho))
        assert all(close(r["gain_b0.5"], P[(r["cell"], 0.5)]["gain"].mean()) for r in h3)
        assert close(s1["h3_kendall_tau"], stats.kendalltau(rho, [P[(r["cell"], 0.5)]["gain"].mean() for r in h3]).statistic)
        c0 = {sp: ctx[("primary", sp)] for sp in SPLITS}

        def wrong_gain(obj, prior, eff, frac=Fraction(5, 100)):
            v = []
            for sp in SPLITS:
                bt = brute("primary", sp, obj, prior, frac)
                for s in ar.SEEDS:
                    o0 = bt.outcomes(loaded[("primary", sp, obj, "P0", 0.0, s)])
                    ow = bt.outcomes(loaded[("primary", sp, obj, f"P_wrong:{prior}", 0.5, s)])
                    v.append(o0["T99"] - ow["T99"] if eff == "T99" else ow["N5"] - o0["N5"])
            return float(np.mean(v))
        mc = rd("manipulation_check.csv")
        assert len(mc) == 6 and close([r for r in mc if r["cell"] == "Y1:P_H_fixed"][0]["gain_P_wrong"],
                                      wrong_gain("Y1", "P_H_fixed", "T99"))

        # (d) recall by prior-percentile bin, family table, engine health, budget, random reference
        edges = [(0.0, 0.5), (0.5, 0.8), (0.8, 0.9), (0.9, 0.95), (0.95, 1.0)]
        for obj, prior, beta in (("Y1", "P_phys_fixed", 0.5), ("Y2", "P_phys_matched_Y2", 0.05)):
            cnt = {e: [0, 0] for e in edges}
            for sp in SPLITS:
                bt = brute("primary", sp, obj, prior)
                for s in ar.SEEDS:
                    rec = loaded[("primary", sp, obj, prior, beta, s)]
                    for m in bt.top:
                        e = [e for e in edges if e[0] <= bt.pct[m] < e[1]][0]
                        cnt[e][0] += m in rec["trajectory"]
                        cnt[e][1] += m not in rec["init"]
            bn = [r for r in rd("recall_by_prior_percentile.csv") if r["universe"] == "primary"
                  and r["cell"] == f"{obj}:{prior}" and float(r["beta"]) == beta]
            assert [(float(r["pctl_lo"]), float(r["pctl_hi"])) for r in bn] == edges
            assert [[int(r["found"]), int(r["available"])] for r in bn] == [cnt[e] for e in edges]
            assert len([v for v in cnt.values() if v[1] > 0]) >= 3
        idx2 = ar.boot_indices(2000)
        for kt, stype in enumerate(ar.SPLIT_TYPES):
          tsplits = SPLITS[10 * kt:10 * kt + 10]
          fam = [r for r in rd("family.csv") if r["universe"] == "primary" and r["cell"] == "Y1:P_H_fixed"
                 and r["split_type"] == stype]
          mats = {}
          for sp in tsplits:
              bt = brute("primary", sp, "Y1", "P_H_fixed")
              for m in bt.top:
                  mats.setdefault(bt.fam[m], set()).add(m)
          big = sorted(f for f, m in mats.items() if len(m) >= 5)
          assert [r["family"] for r in fam] == big + ["Heusler:ABC2_225", "other"]
          assert 0 < len(big) < len(mats), "the synthetic pool must have families on both sides of the threshold"
          for label, test in [("Heusler:ABC2_225", lambda f: f == "ABC2_225"), ("other", lambda f: f != "ABC2_225"),
                              (big[0], lambda f: f == big[0])]:
              N, Dn = {}, {}
              for b, prior in ((0.0, "P0"), (0.5, "P_H_fixed")):
                  N[b], Dn[b] = np.zeros((10, 3)), np.zeros((10, 3))
                  for j, sp in enumerate(tsplits):
                      bt = brute("primary", sp, "Y1", "P_H_fixed")
                      members = [m for m in bt.top if test(bt.fam[m])]
                      for k, s in enumerate(ar.SEEDS):
                          r = loaded[("primary", sp, "Y1", prior, b, s)]
                          N[b][j, k] = len(set(members) & set(r["trajectory"]))
                          Dn[b][j, k] = len(set(members) - set(r["init"]))
              rate = {b: N[b].sum() / Dn[b].sum() for b in N}
              row = [r for r in fam if r["family"] == label][0]
              assert close(row["ratio"], rate[0.5] / rate[0.0]) and close(row["rate_b0"], rate[0.0]) and close(row["rate_b0.5"], rate[0.5])
              assert int(row["distinct_S_top_materials"]) == len(set().union(*[m for f, m in mats.items() if test(f)]))
              bs = []
              for bb in range(2000):
                  tot = {(b, w): 0.0 for b in (0.0, 0.5) for w in "ND"}
                  for j in range(10):
                      js = idx2[0][bb, kt, j]
                      for i in range(3):
                          ks = idx2[1][bb, kt, j, i]
                          for b in (0.0, 0.5):
                              tot[(b, "N")] += N[b][js, ks]
                              tot[(b, "D")] += Dn[b][js, ks]
                  if tot[(0.0, "D")] and tot[(0.5, "D")] and tot[(0.0, "N")]:
                      bs.append((tot[(0.5, "N")] / tot[(0.5, "D")]) / (tot[(0.0, "N")] / tot[(0.0, "D")]))
              assert int(row["valid_resamples"]) == len(bs)
              assert close(row["boot_lo"], np.percentile(bs, 2.5)) and close(row["boot_hi"], np.percentile(bs, 97.5))
        eh = [r for r in rd("engine_health.csv") if r["universe"] == "primary" and r["cell"] == "Y1:P_H_fixed"]
        assert len(eh) == 8 and all(close(r["share_flat_EI"], 0.1) and close(r["share_length_scale_at_lower_bound"], 0.25)
                                    and r["steps"] == "2400" for r in eh if r["beta"] != "1")
        assert [r for r in eh if r["beta"] == "1"][0]["share_flat_EI"] == ""
        bsr = [r for r in rd("budget_sensitivity.csv") if r["universe"] == "primary" and r["cell"] == "Y1:P_H_fixed"
               and r["beta"] == "0.5"]
        assert [r["acquisitions"] for r in bsr] == ["10", "20", "30"]
        for r in bsr:
            t = int(r["acquisitions"])
            dl, dg = [], []
            for sp in SPLITS:
                bt = brute("primary", sp, "Y1", "P_H_fixed")
                for s in ar.SEEDS:
                    o0 = bt.outcomes(loaded[("primary", sp, "Y1", "P0", 0.0, s)], t)
                    o5 = bt.outcomes(loaded[("primary", sp, "Y1", "P_H_fixed", 0.5, s)], t)
                    dl.append(o0["R_dis"] - o5["R_dis"]); dg.append(o0["T99"] - o5["T99"])
            assert close(r["L"], np.mean(dl)) and close(r["gain"], np.mean(dg))
        nc = 0
        for r in [x for x in rd("budget_sensitivity.csv") if x["universe"] == "primary"
                  and x["cell"] == "Y2:P_phys_matched_Y2" and x["acquisitions"] == "10"]:
            dg = []
            for sp in SPLITS:
                bt = brute("primary", sp, "Y2", "P_phys_matched_Y2")
                for s in ar.SEEDS:
                    o0 = bt.outcomes(loaded[("primary", sp, "Y2", "P0", 0.0, s)], 10)
                    ob = bt.outcomes(loaded[("primary", sp, "Y2", "P_phys_matched_Y2", float(r["beta"]), s)], 10)
                    dg.append(o0["T99"] - ob["T99"]); nc += ob["T99"] == 11
            assert close(r["gain"], np.mean(dg))
        assert nc > 20, "the synthetic data must contain runs censored within the truncated budget"
        rr = [r for r in rd("random_search_reference.csv") if r["universe"] == "primary" and r["objective"] == "Y1"][0]
        ns = [c0[sp].n for sp in SPLITS]
        assert close(rr["R_dis"], np.mean([40 / (n - 10) for n in ns]))
        assert close(rr["N5"], np.mean([40 * math.ceil(0.05 * n) / (n - 10) for n in ns]))
        assert close(rr["P_T99_censored"], np.mean([np.prod([(n - 10 - math.ceil(0.01 * n) - j) / (n - 10 - j)
                                                             for j in range(40)]) for n in ns]))
        print(f"ok  6  full analysis: deterministic; the run table ({n_rows} sampled rows over all 15 cells), H1, H2,")
        print("       H3, frontier, never-inspected splits, manipulation check, percentile bins, family table,")
        print("       engine health, budget sensitivity and random reference all match an independent recomputation")

        pc = root / "pilot_check.json"
        pc.write_text(json.dumps({"check1_S_dis": {"Y1": {"pass": False}, "Y2": {"pass": True}},
                                  "check2_T99": {"Y1": {"pass": True}, "Y2": {"pass": False}},
                                  "check3_engine": {"pass": True}}))
        t10, n5, sha = ar.fallbacks_from_pilot(pc)
        assert (t10, n5) == (("Y1",), ("Y2",)) and sha == ar.sha256_file(pc)
        pc2 = root / "pilot_check_norerun.json"
        pc2.write_text(json.dumps({"check1_S_dis": {"Y1": {"pass": True}, "Y2": {"pass": True}},
                                   "check2_T99": {"Y1": {"pass": True}, "Y2": {"pass": True}},
                                   "check3_engine": {"pass": None}}))
        expect_fail(lambda: ar.fallbacks_from_pilot(pc2), "engine check did not pass")
        res3 = root / "res3"
        s3 = ar.analyse(d, runs, res3, universes=("primary",), top10=t10, n5=n5, pilot_check_sha256=sha, b=300, glmm=False)
        st = s3["settings"]
        assert st["top_fraction"] == {"Y1": 0.10, "Y2": 0.05} and st["efficiency_outcome"] == {"Y1": "T99", "Y2": "N5"}
        assert st["pilot_check_sha256"] == sha and math.isnan(s3["h3_kendall_tau"])
        ro3 = pd.read_csv(res3 / "run_outcomes.csv")
        f10 = Fraction(10, 100)
        for obj, prior, frac in (("Y1", "P_phys_fixed", f10), ("Y2", "P_H_fixed", Fraction(5, 100))):
            for sp in (SPLITS[0], SPLITS[17]):
                bt = brute("primary", sp, obj, prior, frac)
                rec = loaded[("primary", sp, obj, prior, 0.25, 101)]
                row = ro3[(ro3.cell == f"{obj}:{prior}") & (ro3.split == sp) & (ro3.seed == 101) & (ro3.beta == 0.25)].iloc[0]
                for k, v in bt.outcomes(rec).items():
                    assert close(row[k], v), (obj, k)
        b10 = brute("primary", SPLITS[0], "Y1", "P_phys_fixed", f10)
        assert len(b10.top1) == math.ceil(0.01 * c0[SPLITS[0]].n) and len(b10.top) == math.ceil(0.10 * c0[SPLITS[0]].n)
        P3 = paired(res3 / "run_outcomes.csv", "primary", {"Y1": "T99", "Y2": "N5"})
        for r in csv.DictReader(open(res3 / "h2_detail.csv")):
            g = bounds(P3[(r["cell"], float(r["beta"]))], "gain")
            assert close(r["gain_lower95"], g["lo95"]) and r["beneficial"] == str(bool(g["lo95"] > 0))
        fr3 = list(csv.DictReader(open(res3 / "frontier.csv")))
        assert {r["efficiency_outcome"] for r in fr3 if r["cell"].startswith("Y2")} == {"N5"}
        assert all(close(r["gain"], P3[(r["cell"], float(r["beta"]))]["gain"].mean()) for r in fr3 if r["beta"] != "0")
        m3 = {r["cell"]: r for r in csv.DictReader(open(res3 / "manipulation_check.csv"))}
        assert close(m3["Y2:P_H_fixed"]["gain_P_wrong"], wrong_gain("Y2", "P_H_fixed", "N5"))
        assert close(m3["Y1:P_phys_fixed"]["gain_P_wrong"], wrong_gain("Y1", "P_phys_fixed", "T99", f10))
        b3 = [r for r in csv.DictReader(open(res3 / "budget_sensitivity.csv"))
              if r["cell"] == "Y2:P_H_fixed" and r["beta"] == "0.5" and r["acquisitions"] == "20"][0]
        dn = []
        for sp in SPLITS:
            bt = brute("primary", sp, "Y2", "P_H_fixed")
            for s in ar.SEEDS:
                dn.append(bt.outcomes(loaded[("primary", sp, "Y2", "P_H_fixed", 0.5, s)], 20)["N5"]
                          - bt.outcomes(loaded[("primary", sp, "Y2", "P0", 0.0, s)], 20)["N5"])
        assert close(b3["gain"], np.mean(dn))
        print("ok  7  fallbacks are read from the pilot check; top 10% for one objective and Delta N5 for the other")
        print("       reach the run table, H2, the frontier, the manipulation check and the budget table")

        f = runs / "primary_random_r0_D_Y1_P_H_fixed_b0.5_s100.json"
        good = f.read_text()
        load = lambda: ar.load_runs(runs, ctx, hashes, UNIVERSES)

        def broken(text, fix=None, **change):
            rec = json.loads(good)
            if fix:
                fix(rec)
            rec.update(change)
            f.write_text(json.dumps(rec))
            expect_fail(load, text)

        def swap(rec):
            for k in ("trajectory", "y_trajectory", "steps"):
                rec[k][0], rec[k][1] = rec[k][1], rec[k][0]

        def repeat(rec):
            for k in ("trajectory", "y_trajectory"):
                rec[k][1] = rec[k][0]
            rec["steps"][1]["material_id"] = rec["trajectory"][0]
        cases = [
            ("engine hash", dict(script_sha256="0" * 64)), ("engine version", dict(engine_version="1.0")),
            ("registration id", dict(registered="10.5281/zenodo.23073511")),
            ("scikit-learn version", dict(env=dict(ENV, **{"scikit-learn": "1.3.0"}))),
            ("thread settings", dict(env=dict(ENV, threads=dict(ENV["threads"], OMP_NUM_THREADS="2")))),
            ("thread settings", dict(env={"scikit-learn": ar.SKLEARN_VERSION})),
            ("oracle values", dict(fix=lambda r: r["y_trajectory"].__setitem__(3, r["y_trajectory"][3] + 1.0))),
            ("initial design", dict(fix=lambda r: r.__setitem__("init", r["init"][::-1]))),
            ("run hash", dict(fix=swap)), ("step log", dict(fix=lambda r: r["steps"].reverse())),
            ("length or repeated material", dict(fix=repeat)),
            ("input hash priors", dict(fix=lambda r: r["inputs"].__setitem__("priors", "0" * 64))),
            ("input hash features", dict(fix=lambda r: r["inputs"].__setitem__("features", None))),
            ("wrong type", dict(fix=lambda r: r["run"].__setitem__("seed", "100"))),
            ("side or pool size", dict(fix=lambda r: r["run"].__setitem__("side", "H"))),
            ("side or pool size", dict(fix=lambda r: r["run"].__setitem__("pool_size", r["run"]["pool_size"] + 1))),
        ]
        for text, kw in cases:
            broken(text, **kw)
        f.write_text(good)
        other = runs / "primary_random_r0_D_Y1_P_H_fixed_b0.25_s100.json"
        og = other.read_text()
        ra, rb = json.loads(good), json.loads(og)                # trajectories swapped between two arms
        for k in ("trajectory", "y_trajectory", "steps", "run_hash"):
            ra[k], rb[k] = rb[k], ra[k]
        f.write_text(json.dumps(ra)); other.write_text(json.dumps(rb))
        expect_fail(load, "batch manifest")
        f.write_text(good); other.write_text(og)
        f.unlink()
        expect_fail(load, "expected runs missing")
        f.write_text(good)
        extra = write_run(runs, ctx, hashes, ("primary", "random_r0", "Y1", "P0", 0.0, 10))
        expect_fail(load, "unexpected runs")
        expect_fail(load, "not in the batch manifest")
        extra.unlink()
        shutil.copy(f, runs / "copy.json")
        expect_fail(load, "duplicate run")
        (runs / "copy.json").unlink()
        man = runs / "MANIFEST_confirmatory.json"
        mg = man.read_text()
        man.unlink()
        expect_fail(load, "no batch manifest")
        man.write_text(mg)
        load()
        print(f"ok  8  validation rejects {len(cases) + 6} kinds of broken record set (engine, registration, environment,")
        print("       values, design, hashes, step log, pool size, repeats, manifest, missing / extra / duplicate runs)")
        print("ALL ANALYSIS TESTS PASSED")
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()

"""
run_search.py — Paper 2 search engine (beta-weighted GP + EI over a fixed pool)
===============================================================================
One run = one (universe, split, side, objective, prior, beta, seed).

Engine (v1.1, registration amendment of 2026-10)
-----------------------------------------------
Pool      : materials of the chosen universe on the chosen split side
            ("D" = discovery for experiments; "H" = history, for testing only).
Inputs    : descriptors from features_v2_<universe>.csv, standardised on the
            searched pool (unsupervised; constant columns dropped).
Surrogate : scikit-learn GaussianProcessRegressor,
            ConstantKernel (amplitude fitted, default bounds) * Matern(nu=2.5,
            one isotropic length-scale, bounds 3.0..1e3, start 3.0),
            alpha = 1e-6, normalize_y = True, n_restarts_optimizer = 2,
            random_state = 0. (Isotropic form: decision-log entry 23.)
            v1.0 had a lower length-scale bound of 1e-2. With two nearly
            identical materials observed (polymorphs) the fit collapsed to
            that bound, EI became equal for every candidate, beta = 0 picked
            in hash order and every beta > 0 picked by the prior alone. The
            bound 3.0 (below the typical nearest-neighbour distance of the
            standardised descriptors) was chosen on history-side runs only.
Acquisition: expected improvement over the best observed y (xi = 0), rounded
            to 9 significant digits before ranking.
Selection : score = (1 - beta) * rank(EI) + beta * rank(P) over the remaining
            candidates (average ranks, 1 = worst); the highest score is
            queried; ties go to the smallest sha256("<split>:<material_id>").
            beta = 0 never reads P; beta = 1 never fits the GP.
            beta in {0, 0.02, 0.05, 0.10, 0.25, 0.5, 0.75, 1} (v1.0: 0, 0.25,
            0.5, 0.75, 1).
Budget    : 10 initial points + 40 acquisitions. The initial design depends
            only on (split, seed): identical across objectives, priors, beta.
Oracle    : direct lookup of the MP value; the engine can obtain y only by
            querying the oracle, and every query is logged.
Priors    : P0 (none; beta = 0 only), P_H_fixed, P_H_matched_Y2, P_phys_fixed,
            P_phys_matched_Y2, and P_wrong:<name> = the reversed prior.

Discovery guard
---------------
Searching a discovery pool requires --registered <Zenodo DOI of the
registration>, which is written into the run record; the analysis accepts
only records carrying the DOI of the registered plan. Without a DOI of the
form 10.5281/zenodo.<number> the CLI refuses (history side is free).
An existing run record is never overwritten.

Usage
-----
  python run_search.py --features F --pool P --splits S --priors PR \
      --split random_r0 --objective Y1 --prior P_phys_fixed --beta 0.5 --seed 100 \
      --registered <id> [--universe primary] [--out-dir runs]
"""

import argparse
import csv
import hashlib
import json
import os
import platform
import re
from importlib import metadata
from pathlib import Path

import numpy as np
from scipy.stats import norm, rankdata

ENGINE_VERSION = "1.1"
N_INIT, N_ACQ = 10, 40
EI_SIG = 9
BETAS = (0.0, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0)
LS_LOWER = 3.0
DOI_RE = re.compile(r"^10\.5281/zenodo\.\d+$")
PRIORS = ("P0", "P_H_fixed", "P_H_matched_Y2", "P_phys_fixed", "P_phys_matched_Y2")
UNIVERSE_FLAG = {"primary": "in_primary", "sens1_no_Pm_Tc": "in_sens1_no_Pm_Tc",
                 "unscreened": "in_unscreened"}


# ───────────────────────────── core (pure functions) ─────────────────────────

def hkey(split, mid):
    return int(hashlib.sha256(f"{split}:{mid}".encode()).hexdigest()[:15], 16)


def rsig(x, sig):
    return np.array([float(f"{v:.{sig}g}") for v in np.asarray(x, dtype=float)])


def init_design(n, split, seed, n_init=N_INIT):
    """Indices of the initial design; depends only on (split, seed, n)."""
    s = int(hashlib.sha256(f"init:{split}:{seed}".encode()).hexdigest()[:16], 16)
    return [int(i) for i in np.random.default_rng(s).choice(n, n_init, replace=False)]


def standardise(X):
    sd = X.std(0)
    keep = sd > 0
    return (X[:, keep] - X[:, keep].mean(0)) / sd[keep]


class Oracle:
    """The only path from the engine to y. Logs every query."""

    def __init__(self, y_by_id):
        self._y = dict(y_by_id)
        self.queries = []

    def __call__(self, mid):
        self.queries.append(mid)
        v = self._y[mid]
        if not np.isfinite(v):
            raise ValueError(f"oracle returned non-finite y for {mid}")
        return float(v)


def fit_gp(Xo, yo):
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel as C, Matern
    k = C(1.0) * Matern(length_scale=LS_LOWER, nu=2.5, length_scale_bounds=(LS_LOWER, 1e3))
    gp = GaussianProcessRegressor(k, alpha=1e-6, normalize_y=True,
                                  n_restarts_optimizer=2, random_state=0)
    return gp.fit(Xo, yo)


def expected_improvement(mu, sd, best):
    sd = np.maximum(sd, 1e-12)
    z = (mu - best) / sd
    return (mu - best) * norm.cdf(z) + sd * norm.pdf(z)


def search(ids, Z, prior, oracle, beta, split, seed, n_init=N_INIT, n_acq=N_ACQ):
    """Run one search. ids: sorted material ids of the pool; Z: standardised
    descriptors (same order); prior: array aligned with ids or None."""
    if beta not in BETAS:
        raise ValueError(f"beta must be one of {BETAS}")
    if beta > 0 and prior is None:
        raise ValueError("beta > 0 needs a prior")
    n = len(ids)
    tb = np.array([hkey(split, m) for m in ids])
    init = init_design(n, split, seed, n_init)
    obs = list(init)
    yo = [oracle(ids[i]) for i in obs]
    steps = []
    for t in range(1, n_acq + 1):
        mask = np.ones(n, bool)
        mask[obs] = False
        cand = np.where(mask)[0]
        r_ei = np.zeros(len(cand))
        r_p = np.zeros(len(cand))
        kern, ties = None, None
        if beta < 1:
            gp = fit_gp(Z[obs], np.array(yo))
            mu, sd = gp.predict(Z[cand], return_std=True)
            ei = rsig(expected_improvement(mu, sd, max(yo)), EI_SIG)
            r_ei = rankdata(ei, method="average")
            ties = int(np.sum(ei == ei.max()))            # size of the top-EI tie group
            kern = [round(float(v), 6) for v in np.exp(gp.kernel_.theta)]
        if beta > 0:
            r_p = rankdata(prior[cand], method="average")
        score = (1.0 - beta) * r_ei + beta * r_p
        j = np.lexsort((tb[cand], -score))[0]
        pick = int(cand[j])
        yo.append(oracle(ids[pick]))
        obs.append(pick)
        steps.append({"t": t, "material_id": ids[pick],
                      "ei_rank": float(r_ei[j]), "prior_rank": float(r_p[j]),
                      "n_cand": int(len(cand)), "kernel": kern, "ei_top_ties": ties})
    traj = [ids[i] for i in obs[n_init:]]
    return {"init": [ids[i] for i in init], "trajectory": traj,
            "y_init": yo[:n_init], "y_trajectory": yo[n_init:], "steps": steps,
            "run_hash": hashlib.sha256(json.dumps([[ids[i] for i in init], traj]).encode()).hexdigest()}


# ───────────────────────────── data loading / CLI ────────────────────────────

def objective_value(r, objective):
    G, K, rho = float(r["G_vrh_GPa"]), float(r["K_vrh_GPa"]), float(r["density_g_cm3"])
    if objective == "Y1":
        return G
    if objective == "Y2":
        return 9 * K * G / (3 * K + G) / rho
    raise ValueError(objective)


def load_pool(features, pool, splits, universe, split, side):
    flag = UNIVERSE_FLAG[universe]
    prow = {r["material_id"]: r for r in csv.DictReader(open(pool)) if r[flag] == "True"}
    srow = {r["material_id"]: r[split] for r in csv.DictReader(open(splits))}
    ids = sorted(m for m in prow if srow[m] == side)
    feat = {}
    with open(features) as f:
        rd = csv.reader(f)
        next(rd)
        for r in rd:
            if r[0] in prow:
                feat[r[0]] = np.array(r[1:], dtype=float)
    X = np.array([feat[m] for m in ids])
    return ids, X, prow


def load_prior(priors, split, ids, name):
    if name == "P0":
        return None
    base, sign = (name.split(":", 1)[1], -1.0) if name.startswith("P_wrong:") else (name, 1.0)
    if base not in PRIORS[1:]:
        raise ValueError(f"unknown prior {name}")
    val = {r["material_id"]: float(r[base]) for r in csv.DictReader(open(priors)) if r["split"] == split}
    missing = [m for m in ids if m not in val]
    if missing:
        raise SystemExit(f"prior {base} missing for {len(missing)} pool materials")
    return sign * np.array([val[m] for m in ids])


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser()
    for a in ("features", "pool", "splits", "priors", "split", "objective", "prior"):
        ap.add_argument(f"--{a}", required=True)
    ap.add_argument("--beta", type=float, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--universe", default="primary", choices=list(UNIVERSE_FLAG))
    ap.add_argument("--side", default="D", choices=["D", "H"])
    ap.add_argument("--registered", default=None)
    ap.add_argument("--out-dir", default="runs")
    a = ap.parse_args(argv)

    if a.side == "D" and not (a.registered and DOI_RE.match(a.registered)):
        raise SystemExit("Refusing to search a DISCOVERY pool without --registered "
                         "<Zenodo DOI of the registration, 10.5281/zenodo.NNN>. "
                         "Register the analysis plan first.")
    if a.prior == "P0" and a.beta != 0:
        raise SystemExit("P0 is only defined at beta = 0")
    if a.side == "H" and a.prior not in ("P0",):
        raise SystemExit("history-side runs support P0 only (prior files cover discovery)")

    out = Path(a.out_dir)
    name = f"{a.universe}_{a.split}_{a.side}_{a.objective}_{a.prior.replace(':', '-')}_b{a.beta}_s{a.seed}.json"
    if (out / name).exists():
        raise SystemExit(f"run record exists, not overwritten: {out / name}")

    ids, X, prow = load_pool(a.features, a.pool, a.splits, a.universe, a.split, a.side)
    prior = load_prior(a.priors, a.split, ids, a.prior)
    oracle = Oracle({m: objective_value(prow[m], a.objective) for m in ids})
    res = search(ids, standardise(X), prior, oracle, a.beta, a.split, a.seed)
    if oracle.queries != res["init"] + res["trajectory"]:
        raise SystemExit("oracle log differs from initial design + acquisitions")

    rec = {"engine_version": ENGINE_VERSION,
           "run": {"universe": a.universe, "split": a.split, "side": a.side,
                   "objective": a.objective, "prior": a.prior, "beta": a.beta,
                   "seed": a.seed, "pool_size": len(ids)},
           "registered": a.registered,
           "inputs": {k: sha256(getattr(a, k)) for k in ("features", "pool", "splits", "priors")},
           "script_sha256": sha256(__file__),
           "env": {"python": platform.python_version(),
                   **{p: metadata.version(p) for p in ("numpy", "scipy", "scikit-learn")},
                   "threads": {v: os.environ.get(v) for v in
                               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}},
           **res}
    out.mkdir(parents=True, exist_ok=True)
    (out / name).write_text(json.dumps(rec, indent=1))
    print(f"{name}  run_hash {res['run_hash'][:16]}")


if __name__ == "__main__":
    main()

"""
engine_design_probe.py — choose the GP input representation using HISTORY data only
===================================================================================
Question: with 140 descriptors and 10-50 observations, which GP configuration
should the Paper 2 engine use?  Candidates (all Matern-5/2, constant scale,
fixed jitter 1e-6, normalised y, 2 optimiser restarts, random_state=0):

  ard140   ARD length-scales on all 140 standardised descriptors
  iso140   one shared length-scale on all 140 standardised descriptors
  pca10    PCA to 10 components, ARD
  pca20    PCA to 20 components, ARD

Standardisation and PCA are fitted on the descriptor matrix of the whole pool
being searched (unsupervised, no y).

Testbed: the HISTORY side of splits random_r0, random_r1, chemsys_r0,
chemsys_r1 (~1,050 materials each). No discovery-side material or outcome is
read. Search: pure EI (beta = 0), 10 random initial points + 40 acquisitions,
seeds 0-4; the initial design is shared by all configurations (paired).

Outcomes per run: N5 = acquisitions in the pool's top 5%; T99 = index of the
first acquisition in the pool's top 1% (41 if none); rho_end = Spearman of the
final posterior mean with y over unobserved pool members.
"""
import csv
import math
import sys
import warnings
from multiprocessing import Pool

import numpy as np
from scipy.stats import norm, spearmanr
from sklearn.decomposition import PCA
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel as C, Matern

warnings.filterwarnings("ignore")
FEAT, POOL, SPLITS = sys.argv[1], sys.argv[2], sys.argv[3]
SPLIT_NAMES = ["random_r0", "random_r1", "chemsys_r0", "chemsys_r1"]
CONFIGS = ["ard140", "iso140", "pca10", "pca20"]
SEEDS = range(5)
N_INIT, N_ACQ = 10, 40

feat = {r[0]: np.array(r[1:], dtype=float) for r in csv.reader(open(FEAT)) if r[0] != "material_id"}
pool = {r["material_id"]: r for r in csv.DictReader(open(POOL)) if r["in_primary"] == "True"}
splits = {r["material_id"]: r for r in csv.DictReader(open(SPLITS))}


def objective(r, which):
    G, K, rho = float(r["G_vrh_GPa"]), float(r["K_vrh_GPa"]), float(r["density_g_cm3"])
    return G if which == "Y1" else 9 * K * G / (3 * K + G) / rho


def represent(X, cfg):
    Z = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1.0)
    if cfg.startswith("pca"):
        Z = PCA(n_components=int(cfg[3:]), random_state=0).fit_transform(Z)
    return Z


def run(job):
    sp, which, cfg, seed = job
    ids = sorted(m for m in pool if splits[m][sp] == "H")
    X = represent(np.array([feat[m] for m in ids]), cfg)
    y = np.array([objective(pool[m], which) for m in ids])
    n = len(y)
    top5 = set(np.argsort(-y)[:math.ceil(0.05 * n)].tolist())
    top1 = set(np.argsort(-y)[:math.ceil(0.01 * n)].tolist())
    rng = np.random.default_rng(1000 * SPLIT_NAMES.index(sp) + seed)
    obs = list(rng.choice(n, N_INIT, replace=False))
    ls = np.ones(X.shape[1]) if cfg != "iso140" else 1.0
    t99, n5 = N_ACQ + 1, 0
    for t in range(1, N_ACQ + 1):
        k = C(1.0) * Matern(length_scale=ls, nu=2.5, length_scale_bounds=(1e-2, 1e3))
        gp = GaussianProcessRegressor(k, alpha=1e-6, normalize_y=True,
                                      n_restarts_optimizer=2, random_state=0)
        gp.fit(X[obs], y[obs])
        cand = np.setdiff1d(np.arange(n), obs)
        mu, sd = gp.predict(X[cand], return_std=True)
        best = y[obs].max()
        z = (mu - best) / np.maximum(sd, 1e-12)
        ei = (mu - best) * norm.cdf(z) + sd * norm.pdf(z)
        pick = int(cand[np.argmax(ei)])
        obs.append(pick)
        if pick in top5:
            n5 += 1
        if pick in top1 and t99 > N_ACQ:
            t99 = t
    cand = np.setdiff1d(np.arange(n), obs)
    rho_end = spearmanr(gp.predict(X[cand]), y[cand]).statistic
    return {"split": sp, "objective": which, "config": cfg, "seed": seed,
            "N5": n5, "T99": t99, "rho_end": round(float(rho_end), 4)}


if __name__ == "__main__":
    jobs = [(sp, w, c, s) for c in ["iso140", "pca10", "pca20", "ard140"] for s in SEEDS for sp in SPLIT_NAMES for w in ("Y1", "Y2")]
    import os
    out = "engine_design_probe_results.csv"
    done = set()
    if os.path.exists(out):
        done = {(r["split"], r["objective"], r["config"], int(r["seed"])) for r in csv.DictReader(open(out))}
    todo = [j for j in jobs if j not in done]
    fields = ["split", "objective", "config", "seed", "N5", "T99", "rho_end"]
    with open(out, "a", newline="") as f, Pool(2) as p:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        if not done:
            w.writeheader()
        for r in p.imap_unordered(run, todo):
            w.writerow(r); f.flush()
            print(f"{r['split']} {r['objective']} {r['config']} s{r['seed']}", flush=True)
    res = list(csv.DictReader(open(out)))
    for r in res:
        r["N5"], r["T99"], r["rho_end"] = int(r["N5"]), int(r["T99"]), float(r["rho_end"])
    print(f"{'objective':<10}{'config':<8}{'N5 mean':>8}{'T99 med':>8}{'T99 cens':>9}{'rho_end':>8}")
    for which in ("Y1", "Y2"):
        for c in CONFIGS:
            g = [r for r in res if r["objective"] == which and r["config"] == c]
            print(f"{which:<10}{c:<8}{np.mean([r['N5'] for r in g]):>8.2f}"
                  f"{np.median([r['T99'] for r in g]):>8.1f}"
                  f"{np.mean([r['T99'] > N_ACQ for r in g]):>9.2f}"
                  f"{np.mean([r['rho_end'] for r in g]):>8.3f}")
    print("random search: N5 = 2.0 expected; P(T99 censored) ~ 0.66 for 11 top-1% members in ~1050")

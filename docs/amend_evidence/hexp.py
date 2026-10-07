"""History-side only experiments for the v1.1 amendment (no discovery pool is touched)."""
import os, sys, json, csv, math, itertools
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"
sys.path.insert(0, "/home/claude")
import numpy as np
from multiprocessing import Pool
import run_search as rs
O = "/mnt/user-data/outputs/"
SPLITS = ["random_r0", "random_r1", "random_r2", "chemsys_r0", "chemsys_r1", "chemsys_r2"]
rs.BETAS = (0.0, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0)

def make_fit(lb):
    def fit_gp(Xo, yo):
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import ConstantKernel as C, Matern
        k = C(1.0) * Matern(length_scale=max(1.0, lb), nu=2.5, length_scale_bounds=(lb, 1e3))
        return GaussianProcessRegressor(k, alpha=1e-6, normalize_y=True, n_restarts_optimizer=2, random_state=0).fit(Xo, yo)
    return fit_gp

_cache = {}
def pool(split):
    if split not in _cache:
        ids, X, prow = rs.load_pool(O + "features_v2_primary.csv", O + "pool_v1.csv", O + "splits_v1.csv", "primary", split, "H")
        from pymatgen.core import Composition, Element
        from math import fsum
        ref = {r["element"]: r for r in csv.DictReader(open(O + "elemental_reference_springer2005_v1.csv"))}
        P = []
        for m in ids:
            fr = sorted(Composition(prow[m]["formula_pretty"]).fractional_composition.get_el_amt_dict().items())
            gv = fsum(x * float(ref[e]["G_GPa"]) for e, x in fr); gr = 1.0 / fsum(x / float(ref[e]["G_GPa"]) for e, x in fr)
            P.append(float(f"{0.5 * (gv + gr):.12g}"))
        _cache[split] = (ids, rs.standardise(X), prow, np.array(P))
    return _cache[split]

def job(a):
    kind, split, obj, seed, lb, beta, pair = a
    ids, Z, prow, P = pool(split)
    y = {m: rs.objective_value(prow[m], obj) for m in ids}
    rs.fit_gp = make_fit(lb)
    if pair is not None:
        base = rs.init_design(len(ids), split, seed)
        init = list(pair) + [i for i in base if i not in pair][:8]
        orig = rs.init_design; rs.init_design = lambda n, s, sd, k=10: init
    res = rs.search(ids, Z, P if beta > 0 else None, rs.Oracle(y), beta, split, seed)
    if pair is not None: rs.init_design = orig
    yv = np.array([y[m] for m in ids]); order = np.argsort(-yv, kind="stable")
    top = set(ids[i] for i in order[:math.ceil(0.05 * len(ids))])
    st = res["steps"]
    flat = sum(abs(s["ei_rank"] - (s["n_cand"] + 1) / 2) < 1e-9 for s in st) if beta < 1 else 0
    return {"kind": kind, "split": split, "obj": obj, "seed": seed, "lb": lb, "beta": beta, "pair": pair is not None,
            "N5": sum(m in top for m in res["trajectory"]), "flat": flat,
            "min_ls": min(s["kernel"][1] for s in st) if beta < 1 else None,
            "ei_pct": [s["ei_rank"] / s["n_cand"] for s in st], "p_pct": [s["prior_rank"] / s["n_cand"] for s in st],
            "ei_argmax": sum(s["ei_rank"] == s["n_cand"] for s in st), "hash": res["run_hash"][:12]}

if __name__ == "__main__":
    jobs = []
    for sp in SPLITS:
        for seed in range(1, 6):
            for obj in ("Y1", "Y2"):
                for lb in (0.01, 1.0, 3.0):
                    jobs.append(("E1", sp, obj, seed, lb, 0.0, None))
    # E2: near-duplicate pairs (descriptor distance < 0.05) with the largest |dy| per split, forced into the design
    for sp in SPLITS:
        ids, Z, prow, P = pool(sp)
        G = np.array([float(prow[m]["G_vrh_GPa"]) for m in ids])
        D = np.sqrt(np.maximum(0, (Z ** 2).sum(1)[:, None] + (Z ** 2).sum(1)[None, :] - 2 * Z @ Z.T)); iu = np.triu_indices(len(ids), 1)
        cand = [(abs(G[i] - G[j]), int(i), int(j)) for i, j in zip(*iu) if D[i, j] < 0.05]
        cand.sort(reverse=True)
        for _, i, j in cand[:2]:
            for lb in (0.01, 3.0):
                for beta in (0.0, 0.05):
                    jobs.append(("E2", sp, "Y1", 1, lb, beta, (i, j)))
    for sp in SPLITS:
        for seed in range(1, 4):
            for beta in (0.02, 0.05, 0.1, 0.25, 0.5):
                jobs.append(("E3", sp, "Y1", seed, 3.0, beta, None))
    print(len(jobs), "jobs", flush=True)
    with Pool(2) as p, open("hexp.jsonl", "w") as f:
        for k, r in enumerate(p.imap_unordered(job, jobs, chunksize=1)):
            f.write(json.dumps(r) + "\n"); f.flush()
            if k % 40 == 0: print(k, flush=True)
    print("DONE", flush=True)

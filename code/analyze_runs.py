"""
analyze_runs.py — Paper 2 registered analysis, v1.1
===================================================
Implements Sections 6-8 of PREREGISTRATION.md v1.1
(registration https://doi.org/10.5281/zenodo.23118101; it amends v1.0,
doi:10.5281/zenodo.23073511, before any discovery-side search was run).

Design: 20 splits (10 random + 10 chemical-system) x 3 seeds (100-102) = 60
blocks; 8 beta values (0, 0.02, 0.05, 0.10, 0.25, 0.5, 0.75, 1).

  Section 6  run-level outcomes (R_dis, R_blind, T99, N5, AUC, T95) and paired
             effects against beta = 0 on the same (split, seed).
  Section 7  the confirmatory test: a t-test on SPLIT MEANS, stratified by
             split type, with Satterthwaite degrees of freedom (9 to 18);
             H1 (12 tests, Holm); H2 (beneficial-and-safe rule); H3 (Kendall
             tau, no test).
  Section 8  frontier at all beta; recall by prior-percentile bin; engine
             health; random-search reference; budget sensitivity; candidate-
             level mixed model; family-level ratios; P_wrong manipulation
             check; sensitivity universes (descriptive). The v1.0 hierarchical
             bootstrap is reported beside every t-based interval.

Definitions fixed in this file (none depends on results)
---------------------------------------------------------
* Orderings. Materials of a discovery pool D are ordered by y descending, ties
  to the smaller hash key sha256("<split>:<material_id>"). S_dis takes the
  lowest prior scores among S_top, ties to the smaller hash key.
* Prior percentile of a material in D = (average rank of its prior score,
  1 = lowest) - 0.5, divided by |D|. S_blind(P) = members of S_top with prior
  percentile < 0.80. Bins for recall by prior percentile: [0, 0.5), [0.5, 0.8),
  [0.8, 0.9), [0.9, 0.95), [0.95, 1].
* Split-level test. For a run-level quantity v: split mean m_ks = mean over the
  3 seeds of split s of type k; estimate = mean of the 20 split means;
  v_k = sample variance (ddof = 1) of the 10 split means of type k;
  SE = sqrt((v_random + v_chemsys) / 40); df = (v_r + v_c)^2 / ((v_r^2 +
  v_c^2) / 9) (Satterthwaite; 18 when the two variances are equal, 9 when one
  type carries all the variance; 18 if both are 0). One-sided p for
  "mean > 0" = P(T_df >= estimate / SE). Two-sided 95% interval = estimate
  +/- t_0.975 SE; one-sided bounds at level q use t_q. If SE = 0, p = 0 when
  the estimate is > 0 and 1 otherwise.
* H2 levels. "Beneficial" uses a one-sided 95% lower bound. "Safe" uses a
  one-sided upper bound at level 1 - 0.05/6, so that the statement "at least
  one of the six beta values is safe" has error at most 0.05. The unadjusted
  95% bounds are written beside it as descriptive, with a three-way state
  (safe / unsafe / inconclusive).
* H3 is computed over the six cells only if both objectives use the same
  efficiency outcome; after a one-objective N5 fallback the six points are
  listed and tau is left empty.
* Splits never inspected. Splits r3-r9 of each type were created for v1.1 and
  no discovery-side outcome of them was seen before registration. The mean
  gain and loss, with split-level intervals and the H1 p-value, are also
  reported on those 14 splits alone (descriptive).
* Chance recall of a run = 40 / (n - 10), the expected recall of 40 uniformly
  random acquisitions from the n - 10 candidates.
* Engine health. A step is "flat" when its top-EI tie group holds more than
  half of the candidates, and "at the lower bound" when the fitted
  length-scale is within 0.1% of its lower bound. Shares are over all steps of an arm.
* Random-search reference (analytic, no runs): recall = chance recall;
  N5 = 40 |S_top| / (n - 10); P(T99 censored) = hypergeometric probability of
  drawing no S_top1 member in 40 draws from n - 10 candidates (initial-design
  overlap ignored). Means over the 20 splits.
* Budget sensitivity. Outcomes recomputed on the first t = 10, 20, 30
  acquisitions: T99 censored at t + 1, N5 and R_dis on the first t.
* H2 margin sensitivity (descriptive): the adjusted safety bound recomputed
  with 10% and 25% of baseline recall in place of 20%.
* Bootstrap (descriptive, as in v1.0): B = 10,000, default_rng(20260928);
  split_draw = rng.integers(0, 10, (B, 2, 10)); seed_draw = rng.integers(0, 3,
  (B, 2, 10, 3)); types in the order (random, chemsys); the same arrays for
  every cell, beta and universe. Percentile intervals (2.5, 97.5).
* Candidate-level model. statsmodels BinomialBayesMixedGLM (logit),
  found ~ C(beta, reference 0) * percentile, independent random intercepts for
  run and material, default priors, variational Bayes (fit_vb, BFGS) with
  numpy.random.seed(20260928). Rows: every S_top member not in the run's
  initial design, for the 8 runs of each (split, seed). Posterior means and
  SDs rounded to 3 decimals. The fit is iterative: repeated fits agree to
  about 3 decimals, not bit for bit, so its table is listed separately from
  the deterministic outputs; mean-field SDs tend to be too small and are not
  confidence intervals. Optimiser warnings are recorded in summary.json.
* Families. family = formula_anonymous + "_" + spacegroup_number. Rate at
  beta = family S_top members acquired / family S_top members not in the
  initial design, summed over the runs of a split type. Ratio = rate(0.5) /
  rate(0). Reported for "Heusler" (ABC2_225), for "other" (all the rest), and
  for every family with at least 5 distinct S_top materials over the type's
  splits. Bootstrap interval within the type; resamples with a zero
  denominator are dropped and counted.

Pre-declared pilot fallbacks, per objective (Section 9), are not chosen by
hand: the analysis reads them from the archived pilot_check.json
(--pilot-check), whose hash is recorded in the output.
  check1_S_dis[obj].pass false  ->  S_top becomes the top 10% for that objective
  check2_T99[obj].pass false    ->  Delta N5 replaces Delta T99 for that objective
An engine fix after the pilot changes run_search.py and therefore the
ENGINE_SHA256 constant below; that is a logged deviation, not a switch.
Every run record must be listed, with its run hash, in a batch manifest
(MANIFEST_*.json written by run_batch.py) in the runs directory.

Usage
-----
  python analyze_runs.py --data-dir data/derived --runs-dir runs \
      --pilot-check runs_pilot/pilot_check.json --out-dir results
"""

import argparse
import csv
import hashlib
import json
import math
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.stats import kendalltau, rankdata

import run_search as rs

REGISTRATION = "10.5281/zenodo.23118101"
ENGINE_SHA256 = "6d900ad9af93c38019e8be2ee7e0a3037e82d1e649678f8f4dededa41c36f13a"
ANALYSIS_VERSION = "1.1"
SKLEARN_VERSION = "1.8.0"

B = 10_000
BOOT_SEED = 20260928
SPLIT_TYPES = ("random", "chemsys")
REPS = tuple(range(10))
SEEDS = (100, 101, 102)
BETAS = rs.BETAS
H1_BETAS = (0.25, 0.5)
H2_BETAS = (0.02, 0.05, 0.1, 0.25, 0.5, 0.75)
ALPHA = 0.05
SAFE_FRACTION = 0.20
SAFE_SENSITIVITY = (0.10, 0.25)
BUDGETS = (10, 20, 30)
BLIND_PCTL = 0.80
PCTL_BINS = (0.0, 0.5, 0.8, 0.9, 0.95, 1.0)
DF_MAX = len(SPLIT_TYPES) * (len(REPS) - 1)
SAFE_LEVEL = 1 - ALPHA / len(H2_BETAS)       # one-sided level of the confirmatory safety bound
NEW_REPS = tuple(range(3, 10))               # splits created for v1.1 (never inspected)

CELLS = {
    "primary": [("Y1", "P_H_fixed"), ("Y1", "P_phys_fixed"),
                ("Y2", "P_H_fixed"), ("Y2", "P_H_matched_Y2"),
                ("Y2", "P_phys_fixed"), ("Y2", "P_phys_matched_Y2")],
    "sens1_no_Pm_Tc": [("Y1", "P_H_fixed"), ("Y1", "P_phys_fixed"),
                       ("Y2", "P_H_fixed"), ("Y2", "P_H_matched_Y2"),
                       ("Y2", "P_phys_fixed"), ("Y2", "P_phys_matched_Y2")],
    "unscreened": [("Y1", "P_H_fixed"), ("Y2", "P_H_fixed"), ("Y2", "P_H_matched_Y2")],
}
UNIVERSES = tuple(CELLS)
WRONG_UNIVERSES = ("primary",)
ALIGN_OBJ = {"Y1": "Y1_G", "Y2": "Y2_E_over_rho"}
HEUSLER = "ABC2_225"
FILES = {"pool": "pool_v1.csv", "splits": "splits_v2.csv", "alignment": "alignment_v3_primary.csv",
         "priors": "priors_v3_{u}.csv", "features": "features_v2_{u}.csv"}


def split_name(stype, rep):
    return f"{stype}_r{rep}"


def all_splits():
    return [split_name(st, rep) for st in SPLIT_TYPES for rep in REPS]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


# ───────────────────────────── expected run set ──────────────────────────────

def expected_runs(universes=UNIVERSES):
    """Every (universe, split, objective, prior, beta, seed) of the confirmatory stage."""
    out = []
    for u in universes:
        for sp in all_splits():
            for seed in SEEDS:
                for obj in ("Y1", "Y2"):
                    out.append((u, sp, obj, "P0", 0.0, seed))
                for obj, prior in CELLS[u]:
                    for b in BETAS[1:]:
                        out.append((u, sp, obj, prior, b, seed))
                    if u in WRONG_UNIVERSES:
                        out.append((u, sp, obj, f"P_wrong:{prior}", 0.5, seed))
    return out


# ───────────────────────────── pool-level sets ───────────────────────────────

class SplitContext:
    """Discovery pool of one (universe, split): ids, y per objective, sets."""

    def __init__(self, universe, split, pool_rows, split_side, prior_rows, top_frac):
        flag = rs.UNIVERSE_FLAG[universe]
        self.universe, self.split = universe, split
        self.ids = sorted(m for m, r in pool_rows.items()
                          if r[flag] == "True" and split_side[m][split] == "D")
        self.n = len(self.ids)
        self.chance = rs.N_ACQ / (self.n - rs.N_INIT)
        self.pos = {m: i for i, m in enumerate(self.ids)}
        self.hk = np.array([rs.hkey(split, m) for m in self.ids])
        self.family = [f"{pool_rows[m]['formula_anonymous']}_{pool_rows[m]['spacegroup_number']}"
                       for m in self.ids]
        self.y, self.top, self.top1 = {}, {}, {}
        for obj in ("Y1", "Y2"):
            y = np.array([rs.objective_value(pool_rows[m], obj) for m in self.ids])
            order = np.lexsort((self.hk, -y))               # y descending, then hash key
            self.y[obj] = y
            self.top[obj] = [int(i) for i in order[:math.ceil(top_frac[obj] * self.n)]]
            self.top1[obj] = set(int(i) for i in order[:math.ceil(0.01 * self.n)])
        self.prior, self.pctl = {}, {}
        for name, col in prior_rows.items():
            v = np.array([col[m] for m in self.ids])
            self.prior[name] = v
            self.pctl[name] = (rankdata(v, method="average") - 0.5) / self.n

    def s_dis(self, obj, prior):
        top = self.top[obj]
        v = self.prior[prior]
        srt = sorted(top, key=lambda i: (v[i], self.hk[i]))  # lowest prior first
        return set(srt[:math.ceil(len(top) / 2)])

    def s_blind(self, obj, prior):
        return {i for i in self.top[obj] if self.pctl[prior][i] < BLIND_PCTL}


def load_data(data_dir, universes, top_frac):
    d = Path(data_dir)
    pool_rows = {r["material_id"]: r for r in csv.DictReader(open(d / FILES["pool"]))}
    split_side = {r["material_id"]: r for r in csv.DictReader(open(d / FILES["splits"]))}
    ctx, hashes = {}, {"pool": sha256_file(d / FILES["pool"]), "splits": sha256_file(d / FILES["splits"])}
    for u in universes:
        pf, ff = d / FILES["priors"].format(u=u), d / FILES["features"].format(u=u)
        hashes[f"priors:{u}"], hashes[f"features:{u}"] = sha256_file(pf), sha256_file(ff)
        by_split = defaultdict(lambda: defaultdict(dict))
        for r in csv.DictReader(open(pf)):
            for name in rs.PRIORS[1:]:
                if r.get(name, "") != "":
                    by_split[r["split"]][name][r["material_id"]] = float(r[name])
        need = sorted({p for _, p in CELLS[u]})
        for sp in all_splits():
            ctx[(u, sp)] = SplitContext(u, sp, pool_rows, split_side,
                                        {p: by_split[sp][p] for p in need}, top_frac)
    return ctx, hashes


# ───────────────────────────── run records ───────────────────────────────────

def load_runs(runs_dir, ctx, hashes, universes, engine_sha=ENGINE_SHA256,
              registration=REGISTRATION, want=None):
    """Load and validate every record; return {(u, split, obj, prior, beta, seed): record}.
    want: the exact set of expected run keys (default: the confirmatory stage)."""
    runs, problems, manifest = {}, [], {}
    for f in sorted(Path(runs_dir).glob("MANIFEST_*.json")):
        manifest.update(json.loads(f.read_text()).get("runs", {}))
    if not manifest:
        raise SystemExit("RUN VALIDATION FAILED\n  no batch manifest (MANIFEST_*.json) in the runs directory")
    for f in sorted(Path(runs_dir).glob("*.json")):
        if f.name.startswith(("MANIFEST_", "pilot_check")):
            continue                                  # written by run_batch.py / pilot_check.py
        rec = json.loads(f.read_text())
        r = rec["run"]
        if r["universe"] not in universes:
            continue
        key = (r["universe"], r["split"], r["objective"], r["prior"], r["beta"], r["seed"])
        c = ctx.get((r["universe"], r["split"]))
        init, traj = rec["init"], rec["trajectory"]

        def bad(msg):
            problems.append(f"{f.name}: {msg}")

        if c is None or not isinstance(r["beta"], float) or not isinstance(r["seed"], int):
            bad("unknown split, or beta / seed of the wrong type")
            continue
        if key in runs:
            bad("duplicate run")
        if rec.get("engine_version") != rs.ENGINE_VERSION:
            bad("engine version")
        if engine_sha is not None and rec.get("script_sha256") != engine_sha:
            bad("engine hash differs from the registered run_search.py")
        if registration is not None and rec.get("registered") != registration:
            bad("registration id")
        env = rec.get("env", {})
        if env.get("scikit-learn") != SKLEARN_VERSION:
            bad("scikit-learn version")
        if set((env.get("threads") or {}).values()) != {"1"}:
            bad("BLAS thread settings are not all 1")
        if r["side"] != "D" or r["pool_size"] != c.n:
            bad("side or pool size")
        inp = rec.get("inputs", {})
        for k in ("pool", "splits"):
            if inp.get(k) != hashes[k]:
                bad(f"input hash {k}")
        for k in ("priors", "features"):
            if inp.get(k) != hashes[f"{k}:{r['universe']}"]:
                bad(f"input hash {k}")
        if manifest.get(f.name) != rec.get("run_hash"):
            bad("record is not in the batch manifest, or its run hash differs from it")
        if len(init) != rs.N_INIT or len(traj) != rs.N_ACQ or len(set(init + traj)) != rs.N_INIT + rs.N_ACQ:
            bad("length or repeated material")
        elif any(m not in c.pos for m in init + traj):
            bad("material outside the discovery pool")
        else:
            if init != [c.ids[i] for i in rs.init_design(c.n, r["split"], r["seed"])]:
                bad("initial design is not the (split, seed) design")
            y = c.y[r["objective"]]
            got = np.array(rec["y_init"] + rec["y_trajectory"])
            if not np.array_equal(got, y[[c.pos[m] for m in init + traj]]):
                bad("oracle values differ from the pool")
            if hashlib.sha256(json.dumps([init, traj]).encode()).hexdigest() != rec["run_hash"]:
                bad("run hash")
            if [s["material_id"] for s in rec.get("steps", [])] != traj:
                bad("step log differs from the trajectory")
        runs[key] = rec
    want = set(expected_runs(universes)) if want is None else set(want)
    have = set(runs)
    if want - have:
        problems.append(f"{len(want - have)} expected runs missing, e.g. {sorted(want - have)[0]}")
    if have - want:
        problems.append(f"{len(have - want)} unexpected runs, e.g. {sorted(have - want)[0]}")
    if problems:
        raise SystemExit("RUN VALIDATION FAILED\n  " + "\n  ".join(problems[:40]))
    return runs


# ───────────────────────────── run-level outcomes ────────────────────────────

def run_outcomes(c, obj, init, traj, s_dis=None, s_blind=None):
    """Outcomes of one run. init / traj: material ids. s_dis, s_blind: sets of pool indices."""
    I = [c.pos[m] for m in init]
    A = [c.pos[m] for m in traj]
    top, top1 = set(c.top[obj]), c.top1[obj]
    t99 = next((t for t, a in enumerate(A, 1) if a in top1), rs.N_ACQ + 1)
    t95 = next((t for t, a in enumerate(A, 1) if a in top), rs.N_ACQ + 1)
    y = c.y[obj]
    best = np.maximum.accumulate(np.concatenate(([y[I].max()], y[A])))[1:]
    auc = float(np.mean(np.searchsorted(np.sort(y), best, side="right") / c.n))
    out = {"T99": t99, "T95": t95, "N5": sum(a in top for a in A), "AUC": auc}
    if s_dis is not None:
        s_run = s_dis - set(I)
        if not s_run:
            raise SystemExit("S_dis,run is empty; the registered outcome is undefined")
        hits = len(s_run & set(A))
        out.update(S_dis_run=len(s_run), dis_hits=hits, R_dis=hits / len(s_run))
    if s_blind is not None:
        b_run = s_blind - set(I)
        hb = len(b_run & set(A))
        out.update(S_blind_run=len(b_run), blind_hits=hb,
                   R_blind=hb / len(b_run) if b_run else float("nan"))
    return out


def cube(fn):
    """Array [type, split, seed] filled by fn(split_name, seed)."""
    return np.array([[[fn(split_name(st, rep), seed) for seed in SEEDS]
                      for rep in REPS] for st in SPLIT_TYPES], dtype=float)


# ───────────────────────────── inference ─────────────────────────────────────

def split_t(V):
    """Stratified split-level summary of a cube [type, split, seed]:
    (estimate, standard error, Satterthwaite degrees of freedom)."""
    m = V.mean(axis=2)                                   # split means [type, split]
    k = m.shape[1]
    v = m.var(axis=1, ddof=1)                            # one variance per split type
    se = math.sqrt(float(v.sum()) / (m.shape[0] ** 2 * k))
    den = float((v ** 2).sum()) / (k - 1)
    df = float(v.sum()) ** 2 / den if den > 0 else float(m.shape[0] * (k - 1))
    return float(m.mean()), se, df


def p_greater(V):
    """One-sided p-value for H0: mean <= 0 against mean > 0."""
    est, se, df = split_t(V)
    if se == 0:
        return 0.0 if est > 0 else 1.0
    return float(stats.t.sf(est / se, df))


def t_interval(V, level=0.95, one_sided=False):
    """(lower, upper) at the given level; one_sided=True gives two one-sided bounds."""
    est, se, df = split_t(V)
    q = stats.t.ppf(level if one_sided else 0.5 + level / 2, df)
    return est - q * se, est + q * se


def boot_indices(b=B, seed=BOOT_SEED):
    rng = np.random.default_rng(seed)
    sp = rng.integers(0, len(REPS), size=(b, len(SPLIT_TYPES), len(REPS)))
    sd = rng.integers(0, len(SEEDS), size=(b, len(SPLIT_TYPES), len(REPS), len(SEEDS)))
    return sp, sd


def boot_means(V, idx):
    sp, sd = idx
    k = np.arange(V.shape[0])[None, :, None, None]
    return V[k, sp[:, :, :, None], sd].reshape(sp.shape[0], -1).mean(1)


def boot_p_v1_0(bm):
    """The v1.0 bootstrap p-value for H0: mean <= 0 (descriptive in v1.1)."""
    bm = np.asarray(bm)
    return float((1 + np.sum(bm <= 0)) / (len(bm) + 1))


def boot_ci(V, idx):
    lo, hi = np.percentile(boot_means(V, idx), [2.5, 97.5])
    return float(lo), float(hi)


def holm(pvals):
    """Holm step-down adjusted p-values (same order as the input)."""
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    adj, run = np.empty(m), 0.0
    for rank, i in enumerate(np.argsort(p, kind="stable")):
        run = max(run, min(1.0, (m - rank) * p[i]))
        adj[i] = run
    return adj


def h1_decisions(rows, alpha=ALPHA):
    """Add the Holm-adjusted p-value and the decision to the rows of the confirmatory family."""
    for r, a in zip(rows, holm([r["p"] for r in rows])):
        r["p_holm"], r["reject_H0"] = float(a), bool(a <= alpha)
    return rows


def h2_rule(R0_minus_chance, gain, loss, R0):
    """H2. Inputs are cubes; gain and loss are {beta: cube}.
    Returns (applicable, baseline lower bound, rows)."""
    base_lo = t_interval(R0_minus_chance, one_sided=True)[0]
    applicable = bool(base_lo > 0)
    rows = []
    for b in H2_BETAS:
        g_lo = t_interval(gain[b], one_sided=True)[0]
        margin = loss[b] - SAFE_FRACTION * R0
        m_adj = t_interval(margin, level=SAFE_LEVEL, one_sided=True)[1]
        m_lo, m_hi = t_interval(margin, one_sided=True)
        ben, safe = bool(g_lo > 0), bool(m_adj < 0)
        alt = {f"safe_at_{int(100 * f)}pct": bool(
            t_interval(loss[b] - f * R0, level=SAFE_LEVEL, one_sided=True)[1] < 0) for f in SAFE_SENSITIVITY}
        rows.append({"beta": b, "gain_mean": float(gain[b].mean()), "gain_lower95": g_lo,
                     "loss_mean": float(loss[b].mean()),
                     "relative_loss": float(loss[b].mean() / R0.mean()) if R0.mean() > 0 else float("nan"),
                     "margin_upper_adjusted": m_adj, "margin_upper95": m_hi, "margin_lower95": m_lo,
                     "state_unadjusted": "safe" if m_hi < 0 else ("unsafe" if m_lo > 0 else "inconclusive"),
                     "beneficial": ben, "safe": safe,
                     "beneficial_and_safe": bool(applicable and ben and safe), **alt})
    return applicable, base_lo, rows


# ───────────────────────────── per-universe analysis ─────────────────────────

def analyse_universe(u, ctx, runs, idx, eff):
    run_rows, frontier, h1, h2, h2d, wrong, bins, health, budget, fresh = [], [], [], [], [], [], [], [], [], []
    gains = {}
    splits = all_splits()

    for obj, prior in CELLS[u]:
        cell = f"{obj}:{prior}"
        sets = {sp: (ctx[(u, sp)].s_dis(obj, prior), ctx[(u, sp)].s_blind(obj, prior)) for sp in splits}

        def outcome(sp, beta, seed, name=None):
            r = runs[(u, sp, obj, name or ("P0" if beta == 0 else prior), beta, seed)]
            return run_outcomes(ctx[(u, sp)], obj, r["init"], r["trajectory"], *sets[sp])

        res = {b: {(sp, seed): outcome(sp, b, seed) for sp in splits for seed in SEEDS} for b in BETAS}
        for b in BETAS:
            for (sp, seed), o in sorted(res[b].items()):
                run_rows.append({"universe": u, "cell": cell, "split": sp, "seed": seed, "beta": b,
                                 "chance_recall": ctx[(u, sp)].chance, **o})
        V = {m: {b: cube(lambda sp, s, b=b, m=m: res[b][(sp, s)][m]) for b in BETAS}
             for m in ("T99", "N5", "AUC", "R_dis")}
        gain = {b: (V["T99"][0.0] - V["T99"][b]) if eff[obj] == "T99" else (V["N5"][b] - V["N5"][0.0])
                for b in BETAS}
        loss = {b: V["R_dis"][0.0] - V["R_dis"][b] for b in BETAS}
        gains[cell] = gain
        blind_n = np.mean([len(sets[sp][1]) for sp in splits])
        for b in BETAS:
            row = {"universe": u, "cell": cell, "beta": b, "efficiency_outcome": eff[obj]}
            for name, arr in (("gain", gain[b]), ("N5", V["N5"][b]), ("AUC", V["AUC"][b]),
                              ("R_dis", V["R_dis"][b]), ("L", loss[b])):
                lo, hi = t_interval(arr)
                blo, bhi = boot_ci(arr, idx)
                row.update({name: float(arr.mean()), f"{name}_lo": lo, f"{name}_hi": hi,
                            f"{name}_boot_lo": blo, f"{name}_boot_hi": bhi})
            num = sum(o["blind_hits"] for o in res[b].values())
            den = sum(o["S_blind_run"] for o in res[b].values())
            row.update({"S_blind_mean_size": float(blind_n), "R_blind_pooled": num / den if den else float("nan")})
            frontier.append(row)
        for b in H1_BETAS:
            lo, hi = t_interval(loss[b])
            est, se, df = split_t(loss[b])
            bm = boot_means(loss[b], idx)
            h1.append({"universe": u, "cell": cell, "beta": b, "mean_L": est, "se": se, "df": df,
                       "ci_lo": lo, "ci_hi": hi, "p": p_greater(loss[b]),
                       "boot_ci_lo": float(np.percentile(bm, 2.5)), "boot_ci_hi": float(np.percentile(bm, 97.5)),
                       "boot_p_v1_0": boot_p_v1_0(bm)})
        for b in BETAS[1:]:
            sub_g, sub_l = gain[b][:, list(NEW_REPS), :], loss[b][:, list(NEW_REPS), :]
            (glo, ghi), (llo, lhi) = t_interval(sub_g), t_interval(sub_l)
            fresh.append({"universe": u, "cell": cell, "beta": b, "splits": 2 * len(NEW_REPS),
                          "efficiency_outcome": eff[obj], "gain": float(sub_g.mean()), "gain_lo": glo,
                          "gain_hi": ghi, "L": float(sub_l.mean()), "L_lo": llo, "L_hi": lhi,
                          "p_L_greater_0": p_greater(sub_l), "df": split_t(sub_l)[2]})
        chance = cube(lambda sp, s: ctx[(u, sp)].chance)
        applicable, base_lo, rows = h2_rule(V["R_dis"][0.0] - chance, gain, loss, V["R_dis"][0.0])
        h2.append({"universe": u, "cell": cell, "efficiency_outcome": eff[obj],
                   "baseline_R_dis": float(V["R_dis"][0.0].mean()), "chance_recall": float(chance.mean()),
                   "baseline_minus_chance_lower95": base_lo, "applicable": applicable,
                   "beneficial_and_safe_betas": [r["beta"] for r in rows if r["beneficial_and_safe"]]})
        h2d += [{"universe": u, "cell": cell, **r} for r in rows]
        # recall of S_top by prior-percentile bin
        for b in BETAS:
            for lo_b, hi_b in zip(PCTL_BINS[:-1], PCTL_BINS[1:]):
                num = den = 0
                for sp in splits:
                    c = ctx[(u, sp)]
                    p = c.pctl[prior]
                    members = [i for i in c.top[obj] if lo_b <= p[i] < hi_b]      # p is always < 1
                    for seed in SEEDS:
                        r = runs[(u, sp, obj, "P0" if b == 0 else prior, b, seed)]
                        I = {c.pos[m] for m in r["init"]}
                        A = {c.pos[m] for m in r["trajectory"]}
                        num += sum(i in A for i in members)
                        den += sum(i not in I for i in members)
                bins.append({"universe": u, "cell": cell, "beta": b, "pctl_lo": lo_b, "pctl_hi": hi_b,
                             "available": den, "found": num, "recall": num / den if den else float("nan")})
        # engine health per arm
        for b in BETAS:
            flat = bound = steps = 0
            if b < 1:
                for sp in splits:
                    for seed in SEEDS:
                        for st in runs[(u, sp, obj, "P0" if b == 0 else prior, b, seed)]["steps"]:
                            steps += 1
                            flat += st["ei_top_ties"] > 0.5 * st["n_cand"]
                            bound += st["kernel"][1] <= rs.LS_LOWER * 1.001
            health.append({"universe": u, "cell": cell, "beta": b, "steps": steps,
                           "share_flat_EI": flat / steps if steps else float("nan"),
                           "share_length_scale_at_lower_bound": bound / steps if steps else float("nan")})
        # budget sensitivity
        for tmax in BUDGETS:
            def trunc(sp, beta, seed):
                r = runs[(u, sp, obj, "P0" if beta == 0 else prior, beta, seed)]
                o = run_outcomes(ctx[(u, sp)], obj, r["init"], r["trajectory"][:tmax], sets[sp][0])
                return min(o["T99"], tmax + 1), o["N5"], o["R_dis"]
            tr = {b: {(sp, seed): trunc(sp, b, seed) for sp in splits for seed in SEEDS} for b in BETAS}
            for b in BETAS[1:]:
                Vg = cube(lambda sp, s: (tr[0.0][(sp, s)][0] - tr[b][(sp, s)][0]) if eff[obj] == "T99"
                          else (tr[b][(sp, s)][1] - tr[0.0][(sp, s)][1]))
                Vl = cube(lambda sp, s: tr[0.0][(sp, s)][2] - tr[b][(sp, s)][2])
                (glo, ghi), (llo, lhi) = t_interval(Vg), t_interval(Vl)
                budget.append({"universe": u, "cell": cell, "beta": b, "acquisitions": tmax,
                               "efficiency_outcome": eff[obj], "gain": float(Vg.mean()), "gain_lo": glo,
                               "gain_hi": ghi, "L": float(Vl.mean()), "L_lo": llo, "L_hi": lhi})
        if u in WRONG_UNIVERSES:
            def g(sp, seed):
                o0, ow = res[0.0][(sp, seed)], outcome(sp, 0.5, seed, f"P_wrong:{prior}")
                return (o0["T99"] - ow["T99"]) if eff[obj] == "T99" else (ow["N5"] - o0["N5"])
            W = cube(g)
            lo, hi = t_interval(W)
            wrong.append({"universe": u, "cell": cell, "beta": 0.5, "efficiency_outcome": eff[obj],
                          "gain_P_wrong": float(W.mean()), "ci_lo": lo, "ci_hi": hi})
    return run_rows, frontier, h1, h2, h2d, wrong, bins, health, budget, fresh, gains


def random_reference(u, ctx):
    rows = []
    for obj in ("Y1", "Y2"):
        cs = [ctx[(u, sp)] for sp in all_splits()]
        rows.append({"universe": u, "objective": obj,
                     "R_dis": float(np.mean([c.chance for c in cs])),
                     "N5": float(np.mean([rs.N_ACQ * len(c.top[obj]) / (c.n - rs.N_INIT) for c in cs])),
                     "P_T99_censored": float(np.mean([stats.hypergeom.pmf(0, c.n - rs.N_INIT, len(c.top1[obj]), rs.N_ACQ)
                                                      for c in cs]))})
    return rows


def h3_table(alignment_csv, gains, eff):
    rho = defaultdict(dict)
    for r in csv.DictReader(open(alignment_csv)):
        if r["universe"] == "primary" and r["side"] == "history_cv":
            rho[(r["objective"], r["prior"])][r["split"]] = float(r["spearman"])
    pts = []
    for obj, prior in CELLS["primary"]:
        vals = rho[(ALIGN_OBJ[obj], prior)]
        if sorted(vals) != sorted(all_splits()):
            raise SystemExit(f"alignment table: need one history_cv row per split for {obj}:{prior}")
        pts.append({"cell": f"{obj}:{prior}", "history_rho": float(np.mean(list(vals.values()))),
                    "gain_b0.5": float(gains[f"{obj}:{prior}"][0.5].mean())})
    if len(set(eff.values())) > 1:                       # the six gains are not in one unit
        return pts, float("nan")
    tau = kendalltau([p["history_rho"] for p in pts], [p["gain_b0.5"] for p in pts]).statistic
    return pts, float(tau)


# ───────────────────────────── family level (exploratory) ────────────────────

def family_table(u, ctx, runs, idx):
    out = []
    sp_i, sd_i = idx
    for obj, prior in CELLS[u]:
        for k, st in enumerate(SPLIT_TYPES):
            mats = defaultdict(set)
            for rep in REPS:
                c = ctx[(u, split_name(st, rep))]
                for i in c.top[obj]:
                    mats[c.family[i]].add(c.ids[i])
            fams = sorted(f for f, m in mats.items() if len(m) >= 5)
            groups = [(f, {f}) for f in fams] + [("Heusler:" + HEUSLER, {HEUSLER}), ("other", None)]
            for label, members in groups:
                def nd(beta):
                    N = np.zeros((len(REPS), len(SEEDS)))
                    D = np.zeros_like(N)
                    for j, rep in enumerate(REPS):
                        sp = split_name(st, rep)
                        c = ctx[(u, sp)]
                        fam = [i for i in c.top[obj]
                               if (c.family[i] in members if members is not None else c.family[i] != HEUSLER)]
                        for s, seed in enumerate(SEEDS):
                            r = runs[(u, sp, obj, "P0" if beta == 0 else prior, beta, seed)]
                            I = {c.pos[m] for m in r["init"]}
                            A = {c.pos[m] for m in r["trajectory"]}
                            N[j, s] = sum(i in A for i in fam)
                            D[j, s] = sum(i not in I for i in fam)
                    return N, D
                (n0, d0), (n5, d5) = nd(0.0), nd(0.5)

                def ratio(a0, b0, a5, b5):
                    if b0 == 0 or b5 == 0 or a0 == 0:
                        return float("nan")
                    return (a5 / b5) / (a0 / b0)
                bs = np.array([ratio(n0[j, s].sum(), d0[j, s].sum(), n5[j, s].sum(), d5[j, s].sum())
                               for j, s in zip(sp_i[:, k, :, None], sd_i[:, k])])
                ok = bs[np.isfinite(bs)]
                lo, hi = (float(v) for v in np.percentile(ok, [2.5, 97.5])) if len(ok) else (float("nan"),) * 2
                size = len(set().union(*[m for f, m in mats.items()
                                         if (f in members if members is not None else f != HEUSLER)] or [set()]))
                out.append({"universe": u, "cell": f"{obj}:{prior}", "split_type": st, "family": label,
                            "distinct_S_top_materials": size,
                            "rate_b0": float(n0.sum() / d0.sum()) if d0.sum() else float("nan"),
                            "rate_b0.5": float(n5.sum() / d5.sum()) if d5.sum() else float("nan"),
                            "ratio": ratio(n0.sum(), d0.sum(), n5.sum(), d5.sum()),
                            "boot_lo": lo, "boot_hi": hi, "valid_resamples": int(len(ok))})
    return out


# ───────────────────────────── candidate-level model (secondary) ─────────────

def glmm_rows(u, ctx, runs, obj, prior):
    rows = []
    for sp in all_splits():
        c = ctx[(u, sp)]
        for seed in SEEDS:
            for b in BETAS:
                r = runs[(u, sp, obj, "P0" if b == 0 else prior, b, seed)]
                I = {c.pos[m] for m in r["init"]}
                A = {c.pos[m] for m in r["trajectory"]}
                for i in c.top[obj]:
                    if i not in I:
                        rows.append({"found": int(i in A), "beta": f"{b:g}", "pctl": float(c.pctl[prior][i]),
                                     "run": f"{sp}|{seed}|{b:g}", "material": c.ids[i]})
    return rows


def fit_glmm(rows):
    import pandas as pd
    from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM
    df = pd.DataFrame(rows)
    np.random.seed(BOOT_SEED)
    model = BinomialBayesMixedGLM.from_formula(
        "found ~ C(beta, Treatment('0')) * pctl",
        {"run": "0 + C(run)", "material": "0 + C(material)"}, df)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = model.fit_vb()
    k = len(model.exog_names)
    return [{"term": t, "post_mean": round(float(m), 3), "post_sd": round(float(s), 3)}
            for t, m, s in zip(model.exog_names, fit.fe_mean[:k], fit.fe_sd[:k])], \
        {"n_rows": len(df), "n_runs": int(df["run"].nunique()), "n_materials": int(df["material"].nunique()),
         "vcp_mean": [round(float(v), 3) for v in fit.vcp_mean],
         "optimizer_warnings": sorted({str(w.message) for w in caught})}


# ───────────────────────────── output ────────────────────────────────────────

def fmt(v):
    if isinstance(v, (float, np.floating)):
        return "" if math.isnan(v) else f"{float(v):.10g}"
    return json.dumps(v) if isinstance(v, (list, dict)) else str(v)


def write_csv(path, rows):
    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([fmt(r.get(k, "")) for k in cols])


def fallbacks_from_pilot(path):
    """(top10 objectives, N5 objectives, file hash) from the archived pilot_check.json."""
    pc = json.loads(Path(path).read_text())
    if pc["check3_engine"]["pass"] is not True:
        raise SystemExit("pilot_check.json: the engine check did not pass")
    top10 = tuple(o for o in ("Y1", "Y2") if pc["check1_S_dis"][o]["pass"] is not True)
    n5 = tuple(o for o in ("Y1", "Y2") if pc["check2_T99"][o]["pass"] is not True)
    return top10, n5, sha256_file(path)


def analyse(data_dir, runs_dir, out_dir, universes=UNIVERSES, top10=(), n5=(), pilot_check_sha256=None,
            engine_sha=ENGINE_SHA256, registration=REGISTRATION, b=B, glmm=True):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    top_frac = {o: (0.10 if o in top10 else 0.05) for o in ("Y1", "Y2")}
    eff = {o: ("N5" if o in n5 else "T99") for o in ("Y1", "Y2")}
    ctx, hashes = load_data(data_dir, universes, top_frac)
    runs = load_runs(runs_dir, ctx, hashes, universes, engine_sha, registration)
    idx = boot_indices(b)
    T = defaultdict(list)
    gains, glmm_meta = {}, {}
    for u in universes:
        rr, fr, a1, a2, a2d, wr, bn, hl, bd, fs, g = analyse_universe(u, ctx, runs, idx, eff)
        for name, rows in (("run_outcomes", rr), ("frontier", fr), ("h1", a1), ("h2", a2), ("h2_detail", a2d),
                           ("manipulation_check", wr), ("recall_by_prior_percentile", bn),
                           ("engine_health", hl), ("budget_sensitivity", bd), ("new_splits_only", fs),
                           ("random_search_reference", random_reference(u, ctx)),
                           ("family", family_table(u, ctx, runs, idx))):
            T[name] += rows
        gains[u] = g
    h1_decisions([r for r in T["h1"] if r["universe"] == "primary"])
    for name in ("h1", "h2", "h2_detail"):
        for r in T[name]:
            r["role"] = "confirmatory" if r["universe"] == "primary" else "descriptive (sensitivity)"
    h3, tau = h3_table(Path(data_dir) / FILES["alignment"], gains["primary"], eff) \
        if "primary" in universes else ([], float("nan"))
    T["h3_alignment"] = h3
    if glmm and "primary" in universes:
        for obj, prior in CELLS["primary"]:
            terms, meta = fit_glmm(glmm_rows("primary", ctx, runs, obj, prior))
            T["glmm_fixed_effects"] += [{"cell": f"{obj}:{prior}", **t} for t in terms]
            glmm_meta[f"{obj}:{prior}"] = meta
    names = {"h1": "h1_coverage_loss", "h2": "h2_beneficial_and_safe"}
    for name in ("run_outcomes", "frontier", "h1", "h2", "h2_detail", "h3_alignment", "manipulation_check",
                 "recall_by_prior_percentile", "engine_health", "budget_sensitivity", "new_splits_only",
                 "random_search_reference", "family", "glmm_fixed_effects"):
        write_csv(out / f"{names.get(name, name)}.csv", T[name])
    det = sorted(p for p in out.glob("*.csv") if p.name != "glmm_fixed_effects.csv")
    summary = {"analysis_version": ANALYSIS_VERSION, "registration": registration,
               "analysis_sha256": sha256_file(__file__), "engine_sha256": engine_sha, "n_runs": len(runs),
               "settings": {"top_fraction": top_frac, "efficiency_outcome": eff, "safe_fraction": SAFE_FRACTION,
                            "safe_level": SAFE_LEVEL, "blind_percentile": BLIND_PCTL, "alpha": ALPHA,
                            "df_max": DF_MAX, "B": b, "boot_seed": BOOT_SEED,
                            "fallback_top10": sorted(top10), "fallback_N5": sorted(n5),
                            "pilot_check_sha256": pilot_check_sha256},
               "input_hashes": hashes, "h3_kendall_tau": tau, "glmm": glmm_meta,
               "outputs": {p.name: sha256_file(p) for p in det},
               "outputs_iterative_fit": {"glmm_fixed_effects.csv": sha256_file(out / "glmm_fixed_effects.csv")}}
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--runs-dir", required=True)
    ap.add_argument("--pilot-check", required=True,
                    help="pilot_check.json of the last pilot; the per-objective fallbacks are read from it")
    ap.add_argument("--out-dir", default="results")
    a = ap.parse_args(argv)
    top10, n5, sha = fallbacks_from_pilot(a.pilot_check)
    s = analyse(a.data_dir, a.runs_dir, a.out_dir, top10=top10, n5=n5, pilot_check_sha256=sha)
    print(f"analysed {s['n_runs']} runs; fallbacks: top 10% {list(top10)}, N5 {list(n5)}; outputs in {a.out_dir}")
    for name, h in s["outputs"].items():
        print(f"  {h[:16]}  {name}")


if __name__ == "__main__":
    main()

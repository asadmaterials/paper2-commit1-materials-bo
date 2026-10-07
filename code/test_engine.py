"""
test_engine.py — tests for run_search.py v1.1. Uses HISTORY pools and synthetic
data only; no discovery-side material is searched.

Usage: python test_engine.py FEATURES POOL SPLITS REFERENCE PRIORS
Prints one line per test. The reference run hashes (one 40-step search per
beta value) are ASSERTED, so a change to the surrogate, the acquisition or
the mixing weights fails the suite.
"""
import sys
import time
import warnings

import numpy as np

import run_search as rs

warnings.filterwarnings("ignore")
FEAT, POOL, SPLITS, REF, PRIORS = sys.argv[1:6]
SPLIT = "random_r0"
RESULTS, HASHES = [], {}


def test(fn):
    t = time.time()
    try:
        fn()
        RESULTS.append((fn.__name__, "PASS", time.time() - t))
    except Exception as e:                      # any error counts as a failure
        RESULTS.append((fn.__name__, f"FAIL: {type(e).__name__}: {e}", time.time() - t))
    print(f"{RESULTS[-1][1]:<6} {fn.__name__}  ({RESULTS[-1][2]:.1f}s)", flush=True)
    return fn


# fixture: history pool of random_r0, Y1, with P_phys computed for history materials
ids, X, prow = rs.load_pool(FEAT, POOL, SPLITS, "primary", SPLIT, "H")
Z = rs.standardise(X)
y = {m: rs.objective_value(prow[m], "Y1") for m in ids}


def pphys():
    import csv
    import importlib.util
    spec = importlib.util.spec_from_file_location("bp", "build_priors.py")
    bp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bp)
    from pymatgen.core import Composition
    ref = {r["element"]: r for r in csv.DictReader(open(REF))}
    fix, _ = bp.p_phys([Composition(prow[m]["formula_pretty"]) for m in ids], ref)
    return fix


P = pphys()
R = np.random.default_rng(7).normal(size=len(ids))      # synthetic unrelated prior
NA = 12                                                  # shorter budget for most tests


def run(prior, beta, seed=0, n_acq=NA, yy=None):
    o = rs.Oracle(yy if yy is not None else y)
    r = rs.search(ids, Z, prior, o, beta, SPLIT, seed, n_acq=n_acq)
    return r, o


@test
def test_init_depends_only_on_split_and_seed():
    a = rs.init_design(len(ids), SPLIT, 3)
    assert a == rs.init_design(len(ids), SPLIT, 3)
    assert a != rs.init_design(len(ids), SPLIT, 4)
    assert a != rs.init_design(len(ids), "chemsys_r0", 3)
    r1, _ = run(P, 1.0, seed=3, n_acq=1)
    r2, _ = run(R, 0.5, seed=3, n_acq=1)
    r3, _ = run(None, 0.0, seed=3, n_acq=1)
    assert r1["init"] == r2["init"] == r3["init"] == [ids[i] for i in a]


@test
def test_beta0_ignores_prior():
    a, _ = run(None, 0.0)
    b, _ = run(P, 0.0)
    c, _ = run(-R, 0.0)
    assert a["trajectory"] == b["trajectory"] == c["trajectory"]


@test
def test_beta1_is_prior_order_with_hash_ties():
    r, _ = run(P, 1.0, n_acq=40)
    init = set(r["init"])
    rest = [(-P[i], rs.hkey(SPLIT, m), m) for i, m in enumerate(ids) if m not in init]
    expect = [m for _, _, m in sorted(rest)[:40]]
    assert r["trajectory"] == expect, "beta=1 must follow prior order, ties by hash"


@test
def test_tie_break_all_equal_prior():
    r, _ = run(np.zeros(len(ids)), 1.0)
    init = set(r["init"])
    expect = [m for m in sorted((m for m in ids if m not in init), key=lambda m: rs.hkey(SPLIT, m))][:NA]
    assert r["trajectory"] == expect


@test
def test_wrong_prior_picks_bottom():
    r, _ = run(-P, 1.0)
    init = set(r["init"])
    rest = sorted(((P[i], rs.hkey(SPLIT, m), m) for i, m in enumerate(ids) if m not in init))
    assert r["trajectory"] == [m for _, _, m in rest[:NA]]


@test
def test_engine_reads_only_queried_y():
    r, o = run(P, 0.5)
    q = set(o.queries)
    assert o.queries == r["init"] + r["trajectory"], "oracle log must equal init + trajectory"
    blind = {m: (y[m] if m in q else float("nan")) for m in ids}
    r2, o2 = run(P, 0.5, yy=blind)                       # would raise if any other y were read
    assert r2["trajectory"] == r["trajectory"]


@test
def test_candidates_valid():
    r, _ = run(P, 0.25, n_acq=40)
    allq = r["init"] + r["trajectory"]
    assert len(set(allq)) == 50, "no repeats"
    assert set(allq) <= set(ids), "only pool materials"
    assert all(s["n_cand"] == len(ids) - 10 - (s["t"] - 1) for s in r["steps"])


@test
def test_prior_mixing_changes_path():
    a, _ = run(P, 0.0)
    b, _ = run(P, 0.5)
    c, _ = run(P, 1.0)
    t = [a["trajectory"], b["trajectory"], c["trajectory"]]
    assert t[0] != t[1] and t[1] != t[2] and t[0] != t[2]


@test
def test_guard_refuses_discovery_without_registration():
    # Nonexistent inputs: if the guard ever fails, the call crashes on loading
    # instead of searching a discovery pool.
    base = ["--features", "/nonexistent/f.csv", "--pool", "/nonexistent/p.csv",
            "--splits", "/nonexistent/s.csv", "--priors", "/nonexistent/pr.csv",
            "--split", SPLIT, "--objective", "Y1", "--prior", "P_phys_fixed",
            "--beta", "0.5", "--seed", "0", "--out-dir", "/nonexistent/out"]
    for extra in ([], ["--registered", "yes"], ["--registered", "zenodo.123"]):
        try:
            rs.main(base + extra)
        except SystemExit as e:
            assert "Refusing" in str(e), f"{extra}: {e}"
            continue
        raise AssertionError(f"discovery run was not refused with {extra}")


@test
def test_rejects_invalid_inputs():
    assert rs.BETAS == (0.0, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0)
    for bad in [(None, 0.5), (P, 0.3)]:
        try:
            run(*bad, n_acq=1)
        except ValueError:
            continue
        raise AssertionError(f"accepted invalid {bad[1]}")


REFERENCE = {   # history pool random_r0, Y1, elemental-modulus prior, seed 1, 40 steps
    0.0: "79a48ce1eb7b349f4e823f628e57edc0523f6441be09165ab8db25a7c6a062ed",
    0.02: "d212dce81e7fc194a03736290786cc371c8794be18f50204cbcd2cc9cd4fe37b",
    0.05: "1677aa4616b2998201005cf2117cd401d0203bbe0a483aeb47ae7d0a037c9c40",
    0.1: "27d9d2756e1333d817e6466568b2fb0dde3efaa5a8f1a5c5b69f2d842ef4e607",
    0.25: "3d038019f5fbf255d2c2ba79c4d0e3d5d2501716d347631eb61cca6b6b1e67da",
    0.5: "28330f33123d57aa9d4ab450f2295e80089448e5479bd3e61fb3ee4884291107",
    0.75: "40a378e11125d0b189729620de62e46025089489e6cd913b7896ef3d9010a657",
    1.0: "0cbe8e2100adaf1f8eab65d6a679e022bc204435bbeef720595a31a1e696d386",
}


@test
def test_reference_hashes_every_beta():
    for beta in rs.BETAS:
        a, _ = run(P, beta, seed=1, n_acq=40)
        HASHES[f"H_random_r0_Y1_Pphys_b{beta}_s1"] = a["run_hash"]
        assert a["run_hash"] == REFERENCE[beta], f"beta={beta}: {a['run_hash']}"
    b, _ = run(P, 0.05, seed=1, n_acq=40)
    assert b["run_hash"] == REFERENCE[0.05], "not deterministic"
    assert len(set(REFERENCE.values())) == len(rs.BETAS), "two beta values gave the same run"


@test
def test_expected_improvement_against_quadrature():
    from scipy import integrate, stats
    for mu, sd, best in [(0.0, 1.0, 0.5), (1.2, 0.3, 1.0), (-2.0, 0.7, 0.0), (3.0, 2.0, 3.0), (0.4, 1e-3, 0.5)]:
        ref, _ = integrate.quad(lambda v: max(v - best, 0.0) * stats.norm.pdf(v, mu, sd),
                                mu - 12 * sd, mu + 12 * sd, points=[best], limit=200)
        got = float(rs.expected_improvement(np.array([mu]), np.array([sd]), best)[0])
        assert abs(got - ref) <= 1e-9 + 1e-7 * abs(ref), (mu, sd, best, got, ref)
    hi = rs.expected_improvement(np.array([1.0, 2.0]), np.array([1.0, 1.0]), 0.0)
    assert hi[1] > hi[0], "EI must increase with the mean (maximisation)"


@test
def test_selection_rule_matches_logged_ranks():
    """Recompute step 1 of a beta = 0.25 run from an independent GP fit and the
    written rule: score = (1 - beta) * rank(EI) + beta * rank(P), ties by hash."""
    from scipy.stats import rankdata
    beta = 0.25
    r, _ = run(P, beta, seed=2, n_acq=1)
    init = [ids.index(m) for m in r["init"]]
    cand = [i for i in range(len(ids)) if i not in set(init)]
    yo = np.array([y[ids[i]] for i in init])
    gp = rs.fit_gp(Z[init], yo)
    mu, sd = gp.predict(Z[cand], return_std=True)
    ei = rs.rsig(rs.expected_improvement(mu, sd, yo.max()), rs.EI_SIG)
    score = (1 - beta) * rankdata(ei) + beta * rankdata(P[cand])
    best = max(range(len(cand)), key=lambda j: (score[j], -rs.hkey(SPLIT, ids[cand[j]])))
    assert r["trajectory"] == [ids[cand[best]]]
    assert r["steps"][0]["ei_rank"] == rankdata(ei)[best] and r["steps"][0]["prior_rank"] == rankdata(P[cand])[best]


@test
def test_no_collapse_with_near_duplicate_pair():
    """v1.0 defect: a near-identical pair with different y in the design drove the
    length-scale to its bound and made EI flat. v1.1 must keep EI informative."""
    d2 = ((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    yv = np.array([y[m] for m in ids])
    close = np.argwhere(d2 < 0.05 ** 2)
    i, j = max(((a, b) for a, b in close if a < b), key=lambda ab: abs(yv[ab[0]] - yv[ab[1]]))
    base = rs.init_design(len(ids), SPLIT, 1)
    init = [int(i), int(j)] + [k for k in base if k not in (i, j)][:8]
    orig = rs.init_design
    rs.init_design = lambda n, split, seed, n_init=rs.N_INIT: init
    try:
        r = rs.search(ids, Z, None, rs.Oracle(y), 0.0, SPLIT, 1, n_acq=15)
    finally:
        rs.init_design = orig
    assert all(s["kernel"][1] >= rs.LS_LOWER - 1e-9 for s in r["steps"])
    assert all(s["ei_top_ties"] < 0.5 * s["n_cand"] for s in r["steps"]), "EI is flat"


def twin_pool(eps):
    """Synthetic pool where the two best-EI candidates are twins (identical, or
    differing by eps in one descriptor). Everything else copies an initial point."""
    split, n = "synthetic", 40
    init = rs.init_design(n, split, 0)
    rng = np.random.default_rng(11)
    Zs = np.zeros((n, 4))
    Zs[init] = rng.normal(size=(len(init), 4))
    others = [i for i in range(n) if i not in init]
    a, b = others[0], others[1]
    Zs[a] = np.array([6.0, 6.0, 6.0, 6.0])
    Zs[b] = Zs[a] + np.array([eps, 0, 0, 0])
    for k, i in enumerate(others[2:]):
        Zs[i] = Zs[init[k % len(init)]]
    sids = [f"m{i:03d}" for i in range(n)]
    ys = {m: float(np.sin(Zs[i]).sum()) for i, m in enumerate(sids)}
    winner = min((sids[a], sids[b]), key=lambda m: rs.hkey(split, m))
    loser = sids[b] if winner == sids[a] else sids[a]
    return sids, Zs, ys, split, winner, loser


@test
def test_ei_tie_goes_to_hash_not_prior():
    sids, Zs, ys, split, win, lose = twin_pool(0.0)
    prior = np.array([1.0 if m == lose else 0.0 for m in sids])      # prior favours the hash loser
    r = rs.search(sids, Zs, prior, rs.Oracle(ys), 0.0, split, 0, n_acq=1)
    assert r["trajectory"] == [win], f"beta=0 tie decided by prior: got {r['trajectory']}"
    st = r["steps"][0]
    assert st["ei_top_ties"] == 2 and st["ei_rank"] == st["n_cand"] - 0.5, "tied EI must get the average rank"


@test
def test_standardise():
    Xs = np.column_stack([np.arange(6.0), np.full(6, 3.0), np.array([1.0, 1, 1, 2, 2, 5]) * 1e3])
    Zs = rs.standardise(Xs)
    assert Zs.shape == (6, 2), "constant column must be dropped"
    assert np.allclose(Zs.mean(0), 0) and np.allclose(Zs.std(0), 1)


@test
def test_ei_rounding_makes_near_ties_exact():
    for eps in (1e-13, -1e-13):
        sids, Zs, ys, split, win, _ = twin_pool(eps)
        r = rs.search(sids, Zs, None, rs.Oracle(ys), 0.0, split, 0, n_acq=1)
        assert r["trajectory"] == [win], f"eps={eps}: near-tie not resolved by hash"


if __name__ == "__main__":
    n_fail = sum(r[1] != "PASS" for r in RESULTS)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} passed")
    for k, v in HASHES.items():
        print(f"RUNHASH {k} {v}")
    sys.exit(1 if n_fail else 0)

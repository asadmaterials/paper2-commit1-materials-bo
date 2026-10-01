"""
test_engine.py — tests for run_search.py. Uses HISTORY pools and synthetic data
only; no discovery-side material is searched.

Usage: python test_engine.py FEATURES POOL SPLITS REFERENCE PRIORS
Prints one line per test and, at the end, the run hashes used for the
cross-environment determinism check.
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
    assert a["trajectory"] != b["trajectory"] != c["trajectory"]


@test
def test_guard_refuses_discovery_without_registration():
    try:
        # Nonexistent inputs: if the guard ever fails, the call crashes on loading
        # instead of searching a discovery pool.
        rs.main(["--features", "/nonexistent/f.csv", "--pool", "/nonexistent/p.csv",
                 "--splits", "/nonexistent/s.csv", "--priors", "/nonexistent/pr.csv",
                 "--split", SPLIT, "--objective", "Y1", "--prior", "P_phys_fixed",
                 "--beta", "0.5", "--seed", "0", "--out-dir", "/nonexistent/out"])
    except SystemExit as e:
        assert "Refusing" in str(e)
        return
    raise AssertionError("discovery run was not refused")


@test
def test_rejects_invalid_inputs():
    for bad in [(None, 0.5), (P, 0.3)]:
        try:
            run(*bad, n_acq=1)
        except ValueError:
            continue
        raise AssertionError(f"accepted invalid {bad[1]}")


@test
def test_deterministic_full_runs():
    for beta in (0.0, 0.5):
        a, _ = run(P, beta, seed=1, n_acq=40)
        b, _ = run(P, beta, seed=1, n_acq=40)
        assert a["run_hash"] == b["run_hash"]
        HASHES[f"H_random_r0_Y1_Pphys_b{beta}_s1"] = a["run_hash"]


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

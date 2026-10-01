# Search-engine cross-machine check: Colab result (2026-09-27)

**Result: PASS.** `test_engine.py` passed 13/13 tests, and both full-search run hashes matched the reference values in `ENGINE_CHECK_COLAB.md` exactly.

Every value below is copied from the Colab console output of this run.

## Purpose

To check that `run_search.py` produces identical search paths on different hardware. The check runs two full 40-step GP+EI searches on a HISTORY pool (split `random_r0`, objective `Y1`, prior `P_phys`, seed 1, β = 0.0 and 0.5) and compares their run hashes with reference values produced beforehand. It also runs 11 behavioural tests of the engine (initial design, tie-breaking, prior mixing, oracle access, input validation, and the discovery-run guard). No discovery-side material is searched.

## Environment

| Component | Version |
|---|---|
| Platform | Google Colab (fresh VM; the runtime had been reset earlier in the day) |
| Python | 3.13.15 |
| numpy | 2.4.4 |
| scipy | 1.17.1 |
| scikit-learn | 1.8.0 (upgraded from Colab's preinstalled 1.6.1) |
| Threading | `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1` |

Packages were installed from `requirements-build.txt` (see the notes). pip reported dependency conflicts with Colab's own packages (`google-colab` pins pandas 2.2.3 and requests 2.32.4; `numba` 0.61.2 requires numpy < 2.3). None of these packages is used by the engine, and all 13 tests passed.

## Files

### Scripts

| File | SHA-256 |
|---|---|
| `run_search.py` | `ec4fe544b478ef8e7a236b860eff1aa93644c3e3acc620cc8491caaf1fdbca16` |
| `test_engine.py` | `c0c1c1b75b1267cf12c50a128bbdf0c48e50e2bb67733300ee45b04ff0ec0cfe` |
| `build_priors.py` (v2.1) | `d2fa186eb7488c1091d45a56f479cf7ff5d9f7085c6a58972c692f88154bf8b5` |
| `requirements-build.txt` | `d7606a8aba2909f0add4c68c2c681e9f196c418aa1f7c5c558ce33f968dfabc8` |

### Input data

Restored from the `rebuild_out/` folder saved to Google Drive after the verified rebuild. All hashes match the committed versions.

| File | SHA-256 |
|---|---|
| `features_v2_primary.csv` | `8ff6c835c1a3a087ede1cd91a9c047857612f6a0ef79b608697df1f6b2d9b30f` |
| `pool_v1.csv` | `fb5873dabdc0cf462d6678fe47be30ab2f21e8fb659604cb208cecd3162f9540` |
| `splits_v1.csv` | `d4b473a7900c60363c7c93b0aef7633a2f330b533ed93c4eb64954a4e5f8d92a` |
| `elemental_reference_springer2005_v1.csv` | `fb9363a5b741596de99e902a2cd6a79c216145afe85180b0365c9641ebc70898` |
| `priors_v2_primary.csv` (v2.1) | `de86901c3b806cc6ad236e7849f429928ffab7a4f9443cbf7f5cb3238a7b6b0a` |

## Command

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python test_engine.py \
    rebuild_out/features_v2_primary.csv rebuild_out/pool_v1.csv rebuild_out/splits_v1.csv \
    rebuild_out/elemental_reference_springer2005_v1.csv rebuild_out/priors_v2_primary.csv
```

## Test results

| Test | Result | Time (s) |
|---|---|---|
| `test_init_depends_only_on_split_and_seed` | PASS | 0.5 |
| `test_beta0_ignores_prior` | PASS | 1.6 |
| `test_beta1_is_prior_order_with_hash_ties` | PASS | 0.0 |
| `test_tie_break_all_equal_prior` | PASS | 0.0 |
| `test_wrong_prior_picks_bottom` | PASS | 0.0 |
| `test_engine_reads_only_queried_y` | PASS | 1.1 |
| `test_candidates_valid` | PASS | 2.3 |
| `test_prior_mixing_changes_path` | PASS | 0.9 |
| `test_guard_refuses_discovery_without_registration` | PASS | 0.0 |
| `test_rejects_invalid_inputs` | PASS | 0.0 |
| `test_deterministic_full_runs` | PASS | 7.3 |
| `test_ei_tie_goes_to_hash_not_prior` | PASS | 0.0 |
| `test_ei_rounding_makes_near_ties_exact` | PASS | 0.1 |
| **Total** | **13/13** | **13.8** |

## Run hashes

| Run | Reference (`ENGINE_CHECK_COLAB.md`) | Colab | Match |
|---|---|---|---|
| `H_random_r0_Y1_Pphys_b0.0_s1` | `79a48ce1eb7b349f4e823f628e57edc0523f6441be09165ab8db25a7c6a062ed` | `79a48ce1eb7b349f4e823f628e57edc0523f6441be09165ab8db25a7c6a062ed` | ✓ |
| `H_random_r0_Y1_Pphys_b0.5_s1` | `28330f33123d57aa9d4ab450f2295e80089448e5479bd3e61fb3ee4884291107` | `28330f33123d57aa9d4ab450f2295e80089448e5479bd3e61fb3ee4884291107` | ✓ |

## Notes

1. **`requirements-build.txt` is newer than the committed version.** The file used here (`d7606a8a…`, 328 bytes) is the committed file (`10a4b43e…`, 308 bytes) with one added line, `scikit-learn==1.8.0`. All other pins are identical. The repository copy needs updating.
2. **Provenance of `run_search.py`.** Its hash was recorded only in Colab, at the time of the run. `test_engine.py` was also hashed independently outside Colab and matches (`c0c1c1b7…`).
3. **Scope.** This result shows that the two reference search paths reproduce byte-for-byte on this Colab machine, with the environment above and single-threaded BLAS. The reference environment for the expected run hashes is not recorded in `ENGINE_CHECK_COLAB.md`, so this document makes no claim about which differences between the two environments were tested. Reproduction with multi-threaded BLAS was not tested.
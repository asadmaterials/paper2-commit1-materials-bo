# Registered analysis plan — Paper 2

**The efficiency–coverage trade-off of scientific priors in Bayesian search:
a registered benchmark on metallic multi-element MP elasticity data**

- Author: Muhammad Asad, independent researcher, Pakistan
- Version: 1.0, 2026-09-28
- Registration record: https://doi.org/10.5281/zenodo.23073511
- Data record (raw snapshot): https://doi.org/10.5281/zenodo.22978221
- Code: https://github.com/asadmaterials/paper2-commit1-materials-bo (the
  registered commit is linked from the registration record)

---

## 0. What this registration does and does not claim

The study has two phases.

**Phase A — exploratory development (2026-09-24 to 2026-09-28).** The data
were acquired, the benchmark was curated, and the design, metrics and
hypotheses were developed. Many choices in Phase A were made **after
inspecting the data**, including outcome data (MP elastic values). All 28
decisions are listed, dated and flagged in `DECISION_LOG.md`: seven were made
after seeing outcome data of the whole pool or the discovery side, three after
seeing history-side outcomes only (engine design and testing), and two after
seeing aggregate statistics only. Phase A is not presented as preregistered.

**Phase B — the registered experiment (from this registration on).** No
discovery-side search has been run for results before this registration. One
discovery-side run was executed accidentally during software testing; it was
deleted unread and its seed is retired (Section 9). From this point, the
benchmark, the search engine, the outcomes, the hypotheses and the tests below
are fixed. The engine refuses to search a discovery pool without this
registration's identifier, which it writes into every run record.

Suggested wording for the paper: *"The benchmark and analysis framework were
developed through an explicitly documented exploratory phase. The discovery
experiment was prospectively registered and locked before any discovery-side
results were generated."*

## 1. Scope

All claims concern **the metallic multi-element Materials Project elasticity
benchmark studied here**: computed (PBE DFT) elastic moduli of materials
classified as metallic by MP, containing at least two elements, all from a
metallic-element list, and restricted to elements with screened reference data
(Section 2). About 45% of the primary universe carries MP's `theoretical`
flag. Results are not claimed to generalise to materials discovery at large.

## 2. Frozen inputs

All files are listed with SHA-256 in Appendix A; `rebuild_and_verify.sh`
regenerates every derived file from the raw snapshot and the handbook PDF and
checks the hashes (verified identical on Python 3.11 and 3.13 and on Colab).

| Input | Content |
|---|---|
| Raw data | MP database 2026.04.13, all 13,283 elasticity documents and matching summary documents, unfiltered |
| Elemental reference | Springer Handbook of Condensed Matter and Materials Data (Martienssen & Warlimont, 2005), mechanical-property tables 2.1-6B(b)–2.1-26B(b), compiled mainly from Landolt–Börnstein. Scalar G is kept for an element if at least one computable internal-consistency check passes: implied Poisson ratio E/2G − 1 in [0, 0.5] and within 0.10 of the printed value, or \|ρv_t² − G\|/G ≤ 0.15. No value is repaired or imputed. All 355 extracted values of the 71 candidate elements were verified against the printed tables (no corrections) |
| Universes | **Primary: 3,550 structures, 60 elements.** Sensitivity 1: 3,378 (primary without Pm, Tc). Sensitivity 2: 4,711 (no reference screen; P0 and P_H only) |
| Splits | 6 splits: random ×3 (grouped by composition) and chemical-system ×3 (grouped by chemical system); history ≈ 30%, discovery ≈ 70%. All polymorphs of a composition are always on the same side |

### Information boundary

- Every prior is built from **history-side** materials and elemental data only.
  In each universe, history = that universe's materials on the history side;
  in the primary universe, P_H's history is a subset of the 3,550 benchmark
  structures (materials removed by the reference screen are not used).
- The elemental-modulus prior uses no MP data.
- The search sees a y value only by querying the oracle; the oracle log must
  equal the initial design plus the acquisitions (asserted in every run).
- "Oracle alignment" (computed on the discovery side with ground truth) is a
  diagnostic only and is never presented as knowledge available before search.

## 3. Objectives

- **Y1** = G_VRH (GPa), MP value.
- **Y2** = E/ρ, with E = 9 K_VRH G_VRH / (3 K_VRH + G_VRH) and ρ the MP density.

## 4. Priors

Higher score = more favoured. Scores are rounded to 12 significant digits.

| Name in code | Name in text | Definition |
|---|---|---|
| `P_H_fixed` | historical prior | minus the Euclidean distance (descriptors standardised on history) to the nearest of the top 10% of history materials ranked by **G** |
| `P_H_matched_Y2` | historical prior, matched | as above, top 10% ranked by **E/ρ** (Y2 only) |
| `P_phys_fixed` | elemental-modulus prior | atomic-fraction Hill average of elemental shear moduli: ½(Σxᵢ Gᵢ + 1/Σ(xᵢ/Gᵢ)) |
| `P_phys_matched_Y2` | elemental-modulus prior, matched | the above divided by a volume-additive mixture density Σxᵢ Mᵢ / Σ(xᵢ Mᵢ/ρᵢ) (Y2 only) |
| `P_wrong:<prior>` | reversed prior | the negative of a prior (manipulation check only) |

The elemental-modulus prior is a **heuristic, structure-blind** prior
(composition only; identical for polymorphs), not a physical prediction of
compound moduli.

**Prior–objective cells (6):** Y1 × {P_H_fixed, P_phys_fixed};
Y2 × {P_H_fixed, P_H_matched_Y2, P_phys_fixed, P_phys_matched_Y2}.

## 5. Search engine

`run_search.py` v1.0 (hash in Appendix A); tests in `test_engine.py`.

- Pool: the discovery side of one split in one universe.
- Descriptors: 132 Magpie composition statistics + MP density, volume per
  atom, crystal-system one-hot; constant columns dropped; standardised on the
  searched pool.
- Surrogate: scikit-learn 1.8.0 `GaussianProcessRegressor`,
  ConstantKernel(1.0) × Matérn(ν = 2.5, one isotropic length-scale, bounds
  1e-2–1e3), alpha = 1e-6, `normalize_y=True`, `n_restarts_optimizer=2`,
  `random_state=0`.
- Acquisition: expected improvement over the best observed y (ξ = 0), rounded
  to 9 significant digits.
- Selection: score = (1 − β)·rank(EI) + β·rank(P) over remaining candidates
  (average ranks); highest score queried; ties to the smallest
  sha256("<split>:<material_id>"). β = 0 never reads P; β = 1 never fits the
  GP.
- Budget: 10 initial points + 40 acquisitions. The initial design depends only
  on (split, seed) and is identical across objectives, priors and β.
- All runs are launched with `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`
  (single-threaded BLAS, the configuration verified across machines).
- β ∈ {0, 0.25, 0.5, 0.75, 1}. β = 1 is labelled **prior-only**. The β = 0 run
  of a (split, seed, objective) is shared by all cells.

## 6. Outcomes

For split s, objective Y and the discovery pool D (n = 2,479–2,525 in the
primary universe), materials are ordered by y descending, ties broken by the
hash key above.

- **S_top**: the top ⌈0.05 n⌉ (124–127 materials).
- **S_top1**: the top ⌈0.01 n⌉ (25–26).
- **S_dis(P)**: the ⌈|S_top| / 2⌉ members of S_top with the **lowest** prior
  score under the cell's prior P (ties by hash key); 62–64 materials.

For one run with initial design I and acquisitions a₁…a₄₀:

- **S_dis,run** = S_dis(P) minus the members of I.
- **R_dis** = |{a₁…a₄₀} ∩ S_dis,run| / |S_dis,run|.
- **T99** = the smallest t with aₜ ∈ S_top1; 41 if none.
- **N5** = |{a₁…a₄₀} ∩ S_top|.
- **AUC** = (1/40) Σₜ qₜ, where qₜ = the fraction of D with y ≤ the best y
  among I ∪ {a₁…aₜ}.
- **T95** = the smallest t with aₜ ∈ S_top; 41 if none (descriptive only).

**Paired effects.** Runs are paired on (split, seed); β = 0 is the reference.

- Efficiency gain: ΔT99(β) = T99(β = 0) − T99(β).
- Coverage loss: L(β) = R_dis(β = 0) − R_dis(β), with S_dis from the cell's prior.

**Estimator.** The mean over the 60 confirmatory runs (6 splits × 10 seeds),
every split weighted equally.

**Hierarchical bootstrap.** B = 10,000 resamples, NumPy `default_rng(20260928)`.
In each resample: for each split type, draw 3 splits with replacement from its
3; for each drawn split, draw 10 seeds with replacement; the statistic is the
mean over the 60 drawn runs. The same resample indices are used for every
cell and β.

## 7. Hypotheses and tests

### H1 — coverage loss at moderate prior strength (confirmatory)

For each of the 6 cells and each β ∈ {0.25, 0.5} (12 tests):
H₀: mean L(β) ≤ 0; H₁: mean L(β) > 0.

p = (1 + #{b : bootstrap mean L*_b ≤ 0}) / (B + 1). Holm correction across the
12 tests, familywise α = 0.05. Reported for each test: the estimate, the
two-sided 95% bootstrap interval, and the Holm-adjusted decision.
β = 0.75 and β = 1 are reported descriptively only, because there the prior
dominates the selection and the outcome is largely implied by alignment.

### H2 — does a safe prior strength exist? (confirmatory decision rule)

For each cell, let G₁ = mean ΔT99(1). If G₁ ≤ 0, H2 is **not applicable** in
that cell (prior-only search gives no speed-up to retain). Otherwise, for each
β ∈ {0.25, 0.5, 0.75}:

- Efficiency retention r(β) = mean ΔT99(β) / G₁ (point estimate).
- Coverage bound U(β) = the 95th percentile of the bootstrap distribution of
  mean L(β) (one-sided upper 95% bound).
- **β is safe if r(β) ≥ 0.80 and U(β) < δ, with δ = 0.016** (about one S_dis
  material per run).

Reported per cell: the set of safe β (possibly empty). No multiplicity
adjustment, because this is a pre-specified decision rule, not a family of
significance tests. Both outcomes, "a safe β exists" and "no safe β was
demonstrated", are reported as results. Because recall per run moves in steps
of about 1/62 ≈ δ, H2 can fail to demonstrate safety even when the true loss
is near zero; this strictness was accepted before any result was seen.

### H3 — the frontier and pre-search alignment (exploratory)

For each cell c: xᶜ = the history-side Spearman ρ between prior and objective
(`alignment_v2_primary.csv`, `side = history_cv`, mean over the 6 splits);
yᶜ = mean ΔT99(0.5). Reported: the 6 points and Kendall's τ. No test (n = 6).
A change of sign in a prior's net value across objectives is not predicted,
because no prior is anti-aligned (lowest history-side ρ 0.45 for a single
split, 0.48 as a cell mean).

## 8. Secondary, exploratory and sensitivity analyses

- **Frontier (primary descriptive result):** for every cell, mean ΔT99, N5,
  AUC and L at all five β values, with bootstrap intervals.
- **Candidate-level model (secondary):** logistic mixed model of discovery of
  each S_top member on β (categorical) × the member's prior percentile in D,
  with random intercepts for run and material. The software and settings are
  fixed in the analysis code, which is committed and hashed before the
  confirmatory runs (Section 10).
- **Family level (exploratory):** family = anonymous formula + space group;
  for families with ≥ 5 S_top members pooled within a split type, the ratio of
  discovery rates at β = 0.5 vs β = 0, with bootstrap intervals; Heusler (ABC2,
  #225) and other families reported separately.
- **Manipulation check:** P_wrong at β = 0.5 for each cell; mean ΔT99 < 0
  expected; descriptive.
- **Sensitivity 1** (no Pm/Tc) and **sensitivity 2** (unscreened; P_H cells
  only): the H1 and H2 estimands recomputed with the same δ = 0.016 and
  reported descriptively beside the primary results.

## 9. Runs

**Pilot (feasibility only).** Universe primary; Y1; β = 0, and P_H_fixed and
P_phys_fixed at β ∈ {0.5, 1}; all 6 splits; **seeds 10–14**; 150 runs.

The pilot computes only these checks, at β = 0 over its 30 runs, and no
comparison between arms:

| Check | Pass criterion | Pre-declared fallback if it fails |
|---|---|---|
| S_dis measurable | mean \|A ∩ S_dis,run\| ≥ 1, for each of P_H_fixed and P_phys_fixed | S_top becomes the top 10% (S_dis 124–127, δ = 0.008, one material per run) |
| T99 measurable | T99 censored (= 41) in < 50% of runs | primary efficiency outcome becomes N5; ΔN5(β) = N5(β) − N5(0) replaces ΔT99 everywhere, including H2 retention |
| Engine | identical initial designs across arms per (split, seed); oracle log = I + acquisitions; 5 runs drawn with NumPy `default_rng(20260928)` reproduce their run hash when rerun | fix the code (logged), rerun the pilot with seeds 15–19 |

Pilot runs are not part of any confirmatory analysis. Pilot outcomes may not
change hypotheses, thresholds or metrics beyond the fallbacks above.

**Confirmatory.** **Seeds 100–109** on all 6 splits (60 paired runs per arm).
Runs per (split, seed): primary 32 (Y1: 1 + 2 cells × 4 β + 2 P_wrong;
Y2: 1 + 4 × 4 + 4 P_wrong); sensitivity 1: 26; sensitivity 2: 14. In total
1,920 + 1,560 + 840 = 4,320 runs.

**Retired seed.** Seed 0 is never used. One discovery run with seed 0
(random_r0, Y1, P_phys_fixed, β = 0.5) was executed during software testing
before registration; it was deleted unread (`EXCLUDED_RUNS.json`,
decision-log entry 25).

## 10. Analysis code and deviations

- The analysis code (`analyze_runs.py`) implementing Sections 6–8 exactly is
  written and tested on pilot-format records **before** the confirmatory runs,
  then committed; its hash is added to the decision log before the first
  confirmatory run.
- Any departure from this plan is logged in `DECISION_LOG.md` with its date and
  reason, and reported in the paper next to the registered analysis. Registered
  analyses are always reported as registered, whatever their result.

## 11. Statement on AI assistance

The study design, the code (data extraction, pool, splits, priors, search
engine and tests) and the text of this plan were developed with an AI
assistant (Claude, Anthropic). Separate AI-assisted reviews were used to
critique the design before registration. The author reviewed and approved
every design decision recorded in `DECISION_LOG.md`, ran the cross-machine
rebuilds and engine checks on Colab, and checked the extracted elemental reference values against
the printed handbook tables (with a separate AI-assisted read). The author is
responsible for the content of this registration.

---

## Appendix A — SHA-256 of all frozen inputs

Verify from the repository root with `sha256sum -c` (raw snapshot and
handbook files at the paths shown).

```
4db3c2a0b600da330eeafd18525021971c585e753e051b94cd73e2f5c64763fc  data/raw/mp_2026.04.13/elasticity.jsonl.gz
b485561ade56cd7c15c6e819ac707276114945cb9d6a00ecce3dcb6b70ce91e6  data/raw/mp_2026.04.13/summary.jsonl.gz
dbee15d30fcd8426f3de0acf1b12801923cf46d04b2bb2e37110108ecf208153  handbook/springer-handbook-of-condensed-matter-and-materials-data_compress.pdf
3a31c2e65d463e6bb20c4b8552f710f8b1e7a3a993b4754e4c35e8ad64908e7e  code/pull_mp.py
f7146cece6936b437f73a54ccd078624b50549094c0baa7e18f9b8d7c8517f30  code/extract_reference.py
9d0f4f31ec46d1bd3a6489fc818936e1c387dabdec5126b9b8a9e74cbbf243e6  code/build_pool.py
8d30b0b10f37ff5f2e2a60916622cef167ccb7059d403506a50a2d95c5065926  code/make_splits.py
d2fa186eb7488c1091d45a56f479cf7ff5d9f7085c6a58972c692f88154bf8b5  code/build_priors.py
ec4fe544b478ef8e7a236b860eff1aa93644c3e3acc620cc8491caaf1fdbca16  code/run_search.py
c0c1c1b75b1267cf12c50a128bbdf0c48e50e2bb67733300ee45b04ff0ec0cfe  code/test_engine.py
00a2bc777432a4ae3853ff6d41d92df6b0325bfb5fd6deddc6366cd4b507f0f0  code/engine_design_probe.py
d7606a8aba2909f0add4c68c2c681e9f196c418aa1f7c5c558ce33f968dfabc8  code/requirements-build.txt
f49cec604475ea0ccab9063f2b030cdce4c9cc069a32e7e97bc3f6f05b757565  code/rebuild_and_verify.sh
65454f097428a155d7ecf8b0b0100d8c1f287476c2aa4fb224ca6df129291123  code/expected_hashes.txt
fb9363a5b741596de99e902a2cd6a79c216145afe85180b0365c9641ebc70898  data/derived/elemental_reference_springer2005_v1.csv
fb5873dabdc0cf462d6678fe47be30ab2f21e8fb659604cb208cecd3162f9540  data/derived/pool_v1.csv
d4b473a7900c60363c7c93b0aef7633a2f330b533ed93c4eb64954a4e5f8d92a  data/derived/splits_v1.csv
8ff6c835c1a3a087ede1cd91a9c047857612f6a0ef79b608697df1f6b2d9b30f  data/derived/features_v2_primary.csv
de86901c3b806cc6ad236e7849f429928ffab7a4f9443cbf7f5cb3238a7b6b0a  data/derived/priors_v2_primary.csv
644598ad65d161db0398854092db38f4097be68439c62995f2524c1fc35ff2f5  data/derived/alignment_v2_primary.csv
33526d4b743bc1514746e26cea485df5849fbddf06a5b020458987ed6c63e256  data/derived/features_v2_sens1_no_Pm_Tc.csv
6f9dc7928ba6f2011abc27f84da81ea5bc0b911b93ffa3a8ec89d1bbeec604a7  data/derived/priors_v2_sens1_no_Pm_Tc.csv
7e93481ca8b5f5dceace6cc417875e30a3cb7e0a9f64e1de803928f5c2d28de8  data/derived/alignment_v2_sens1_no_Pm_Tc.csv
4039be545a0ac4c8665fd625c69b788680f387269c9749ca40625dd742a1dd26  data/derived/features_v2_unscreened.csv
4f33ec37a24e3786a284ad88a354ee4884bd98a394e23c9b909c3bdfb0defd0e  data/derived/priors_v2_unscreened.csv
1413095c8763e7f73afa7460154333379ddd77450024da746c587e8c3af4442e  data/derived/alignment_v2_unscreened.csv
7ffeceecd2c85c749156f4ea33b04499f3cb5e3d394113cf4517ccda7ebbfeef  data/derived/engine_design_probe_results.csv
24721262f9e0501d8074af11566368d72d6fc63fad9819bc1ba591bd182273a3  docs/reference_verification_sheet.csv
048c7cde729add2323fa97896a859511531ee44572b224e089f71b30d3431c36  docs/DECISION_LOG.md
5caa1e2425ef405bbb7c48383f26eca09d407b86ab0154d6a318b5a9fb505469  docs/EXCLUDED_RUNS.json
```

The raw snapshot and the handbook are not in the repository: the snapshot is
on Zenodo (data record above); the handbook is cited, not redistributed.

## Appendix B — Reference run hashes (engine reproducibility)

Two 40-step searches on a history pool (random_r0, Y1, elemental-modulus
prior, seed 1), reproduced on Python 3.11 and 3.13 and on Colab hardware:

```
β = 0.0  79a48ce1eb7b349f4e823f628e57edc0523f6441be09165ab8db25a7c6a062ed
β = 0.5  28330f33123d57aa9d4ab450f2295e80089448e5479bd3e61fb3ee4884291107
```

# Registered analysis plan — Paper 2 (version 1.1)

**The efficiency–coverage trade-off of scientific priors in Bayesian search:
a registered benchmark on metallic multi-element MP elasticity data**

- Author: Muhammad Asad, independent researcher, Pakistan
- Version: 1.1, 2026-10-06. Amends version 1.0 (dated 2026-09-28, published
  2026-10-02) **before any discovery-side search was run for results** (one
  accidental run, deleted unread, is described in Section 9)
- Registration record, this version: https://doi.org/10.5281/zenodo.23118101
- Registration record, version 1.0: https://doi.org/10.5281/zenodo.23073511
- Data record (raw snapshot): https://doi.org/10.5281/zenodo.22978221
- Code: https://github.com/asadmaterials/paper2-commit1-materials-bo (the
  registered commit is linked from the registration record)

Every change against version 1.0 is listed, with its reason, in Appendix C.

---

## 0. What this registration does and does not claim

The study has two phases.

**Phase A — exploratory development (2026-09-24 to 2026-10-06).** The data
were acquired, the benchmark was curated, and the design, metrics and
hypotheses were developed, registered as version 1.0, reviewed, and amended.
Many choices in Phase A were made **after inspecting the data**, including
outcome data (MP elastic values). All 36 decisions are listed, dated and
flagged in `DECISION_LOG.md`. Phase A is not presented as preregistered.
What had been seen:

- Aggregate outcome statistics of the whole pool (range, median, 95th
  percentile of G) were seen when the pool was curated (entries 5, 6).
- Outcome values of the whole pool and of the discovery sides of splits r0–r2
  were inspected when the benchmark, the coverage set S_dis, the objectives
  and the hypotheses were developed (entries 12, 13, 15, 16, 17, 19, 21; also
  18, whose flag is corrected in entry 31). The prior-only arm (β = 1) is a
  deterministic function of those data and was computed for two splits
  (entry 16).
- The search engine was chosen on history-side outcomes of four splits
  (entry 23), which cover 75% of the primary pool and 65–77% of every
  discovery pool, because all splits partition the same pool.
- The engine tests (entries 24, 26) and the review used the history sides of
  all six splits r0–r2, which together cover 87% of the primary pool, 82% of
  each r0–r2 discovery pool and 86–88% of each r3–r9 discovery pool. Between
  version 1.0 and this version, the review ran about 1,500 searches on those
  **history** pools, some with the elemental-modulus prior at β = 0.02–0.75.
  They gave engine-health statistics, the rank geometry of picks, top-5% hit
  counts and, for one reviewer, history-side efficiency contrasts between
  arms. No coverage outcome was computed (entries 31–33, 35).
- For the 14 splits added in this version (r3–r9), no discovery-side outcome
  or alignment has been computed.

**Phase B — the registered experiment (from this registration on).** Apart
from one accidental run, no search driven by the Gaussian-process model has
been run on a discovery pool. That run (β = 0.5) was executed during software
testing before version 1.0; it was deleted unread and its seed is retired
(Section 9).
From this point, the benchmark, the search engine, the outcomes, the
hypotheses and the tests below are fixed. The analysis accepts only run
records that carry this version's DOI.

Wording for the paper: *"The benchmark, metrics and hypotheses were developed
with access to the outcome labels of all pools, in an explicitly documented
exploratory phase. The analysis plan for the model-driven search arms was
registered, and amended once after an independent review, before any such
search was run on a discovery pool for results."*

## 1. Scope and limitations

All claims concern **the metallic multi-element Materials Project elasticity
benchmark studied here**: computed (PBE DFT) elastic moduli of materials
classified as metallic by MP, containing at least two elements, none of H, B,
C, N, O, F, Si, P, S, Cl, Se, Br, I, At or a noble gas, and restricted to
elements with screened reference data (Section 2). Ge and Sb are in scope
(330 and 214 structures of the primary universe). About 45% of the primary
universe carries MP's `theoretical` flag, and 304 structures lie more than
0.1 eV/atom above the convex hull; no stability filter is applied. Results
are not claimed to generalise to materials discovery at large.

Limitations that this plan cannot remove and the paper will state:

1. The exposure to outcome data described in Section 0.
2. The universe is filtered by the quality of the elemental prior's own input
   data: the reference screen removes 1,161 of 4,711 structures (24.6%).
3. All splits re-partition one pool; any two discovery pools share about 70%
   of their materials. Inference is conditional on this pool: it describes
   variation over partitions and initial designs, not over benchmarks.
4. The chemical-system splits are only mildly out of distribution: 24–31% of
   their discovery materials have a proper sub-system in history.
5. Descriptors include the DFT-relaxed density and volume per atom of each
   candidate, so the setting is the screening of already-relaxed structures.
6. The historical prior uses about 1,050 labelled history materials, while
   the no-prior search starts from 10 points; for that prior, "prior against
   no prior" is partly "labelled data against none".
7. The prior enters through fixed-weight rank mixing. β is a rank weight; it
   is not comparable to the decaying multiplicative weights of πBO-type
   methods, and no decaying arm is run.
8. One surrogate (an isotropic Gaussian process), one budget (10 + 40), batch
   size 1, a noiseless lookup oracle, one property family.
9. S_dis is selected on low prior score, so some coverage loss at high β is
   built into its definition; H1 asks whether the loss is already present at
   moderate β, and the frontier at small β and the absolute set S_blind carry
   the less constructed evidence.
10. Six of the 20 splits (r0–r2 of each type) were inspected in Phase A. The
    main results are therefore also reported on the 14 splits created for this
    version alone (Section 8).

## 2. Frozen inputs

All files are listed with SHA-256 in Appendix A. `rebuild_and_verify.sh`
regenerates the 12 derived data files of version 1.0 from the raw snapshot
and the handbook PDF; `rebuild_v1_1.sh` regenerates the 7 files added in this
version from those. Hashes are identical on Python 3.11 and 3.13 and on Colab
(`docs/COLAB_CHECK_v1.1_2026-10-06.txt`, decision-log entry 36).

| Input | Content |
|---|---|
| Raw data | MP database 2026.04.13, all 13,283 elasticity documents and matching summary documents, unfiltered |
| Elemental reference | Springer Handbook of Condensed Matter and Materials Data (Martienssen & Warlimont, 2005), mechanical-property tables 2.1-6B(b)–2.1-26B(b). Scalar G is kept for an element if at least one computable internal-consistency check passes: implied Poisson ratio E/2G − 1 in [0, 0.5] and within 0.10 of the printed value, or \|ρv_t² − G\|/G ≤ 0.15. No value is repaired or imputed. The screen checks consistency of a handbook row, not its accuracy. All 355 extracted values of the 71 candidate elements were checked against the printed tables (AI-assisted read plus the author's manual look; no corrections). The handbook PDF is not redistributed; the values extracted from its tables are included with citation |
| Universes | **Primary: 3,550 structures, 60 elements.** Sensitivity 1: 3,378 (primary without Pm, Tc). Sensitivity 2: 4,711 (no reference screen; historical prior only). Unchanged from version 1.0 |
| Splits | **20 splits: random ×10 (grouped by composition) and chemical-system ×10 (grouped by chemical system)**; history ≈ 30%, discovery ≈ 70%; all polymorphs of a composition are always on the same side. Replicates r0–r2 of each type are byte-identical to version 1.0; r3–r9 are generated by the same seeded algorithm (`splits_v2.csv`) |
| Descriptors | 131 Magpie composition statistics, the density of the structure in the MP elasticity document, volume per atom, and a crystal-system one-hot (140 columns; 141 in sensitivity 2). Unchanged from version 1.0 |

Known defects of the frozen inputs, left unchanged and disclosed (decision-log
entry 31): the stability filter tests MP's integer-rounded elastic tensor
with zero tolerance, so 2 primary structures (mp-aaaabhfw, mp-aaacdoei) are
admitted and 19 otherwise-stable structures are rejected on round-off; the
silicon row of the reference table is wrong (silicon is in no pool).

### Information boundary

- Every prior is built from **history-side** outcomes and elemental data
  only. In each universe, history = that universe's materials on the history
  side of the split.
- The elemental-modulus prior uses no MP outcome data.
- The descriptors of a discovery material (including its MP density and
  volume) are used to score it; its elastic outcome is not.
- The search sees a y value only by querying the oracle; the oracle log must
  equal the initial design plus the acquisitions (checked in every run and
  again by the analysis).

## 3. Objectives

- **Y1** = G_VRH (GPa), MP value.
- **Y2** = E/ρ, with E = 9 K_VRH G_VRH / (3 K_VRH + G_VRH) and ρ the density
  of the structure in the MP elasticity document.

## 4. Priors

Higher score = more favoured. Scores are rounded to 12 significant digits.
Definitions are unchanged from version 1.0.

| Name in code | Name in text | Definition |
|---|---|---|
| `P_H_fixed` | historical prior | minus the Euclidean distance (descriptors standardised on history) to the nearest of the top 10% of history materials ranked by **G** |
| `P_H_matched_Y2` | historical prior, matched | as above, top 10% ranked by **E/ρ** (Y2 only) |
| `P_phys_fixed` | elemental-modulus prior | atomic-fraction Hill average of elemental shear moduli: ½(Σxᵢ Gᵢ + 1/Σ(xᵢ/Gᵢ)) |
| `P_phys_matched_Y2` | elemental-modulus prior, matched | the above divided by a volume-additive mixture density Σxᵢ Mᵢ / Σ(xᵢ Mᵢ/ρᵢ) from handbook densities (Y2 only) |
| `P_wrong:<prior>` | reversed prior | the negative of a prior (manipulation check only) |

The elemental-modulus prior is a **heuristic, structure-blind** prior
(composition only; identical for polymorphs), not a physical prediction of
compound moduli.

**Prior–objective cells (6):** Y1 × {P_H_fixed, P_phys_fixed};
Y2 × {P_H_fixed, P_H_matched_Y2, P_phys_fixed, P_phys_matched_Y2}.

## 5. Search engine

`run_search.py` **v1.1** (hash in Appendix A); tests in `test_engine.py`.

- Pool: the discovery side of one split in one universe.
- Inputs: the descriptors of Section 2; constant columns of the searched pool
  are dropped (139 or 140 inputs remain in the primary universe; 140 or 141
  in sensitivity 2); standardised on the searched pool.
- Surrogate: scikit-learn 1.8.0 `GaussianProcessRegressor`, ConstantKernel
  (amplitude fitted) × Matérn(ν = 2.5, one isotropic length-scale, **bounds
  3.0–1e3, optimiser start 3.0**), alpha = 1e-6, `normalize_y=True`,
  `n_restarts_optimizer=2`, `random_state=0`. Version 1.0 had a lower bound of
  0.01; with two nearly identical materials observed the fit collapsed to it
  and the search became prior-only (Appendix C, change 1). On history pools
  the new bound is reached in at least one step of about half of the runs and
  costs 0.2 ± 0.1 top-5% hits per run at β = 0.
- Acquisition: expected improvement over the best observed y (ξ = 0), rounded
  to 9 significant digits.
- Selection: score = (1 − β)·rank(EI) + β·rank(P) over remaining candidates
  (average ranks); highest score queried; ties to the smallest
  sha256("<split>:<material_id>"). β = 0 never reads P; β = 1 never fits the
  GP.
- Budget: 10 initial points + 40 acquisitions. The initial design depends only
  on (split, seed) and is identical across objectives, priors and β.
- All runs are launched with `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and
  `MKL_NUM_THREADS` set to 1; the settings are written into every run record.
- **β ∈ {0, 0.02, 0.05, 0.10, 0.25, 0.5, 0.75, 1}.** β = 1 is labelled
  **prior-only**. The β = 0 run of a (split, seed, objective) is shared by
  all cells.
- The command line refuses to search a discovery pool without a Zenodo DOI in
  `--registered`, and never overwrites a run record; the batch driver stops
  on a record it cannot read and deletes nothing.

## 6. Outcomes

For split s, objective Y and the discovery pool D (n = 2,447–2,525 in the
primary universe), materials are ordered by y descending, ties broken by the
hash key above.

- **S_top**: the top ⌈0.05 n⌉ (123–127 materials).
- **S_top1**: the top ⌈0.01 n⌉ (25–26).
- **S_dis(P)**, "the lower-prior half of the top 5%": the ⌈|S_top| / 2⌉
  members of S_top with the **lowest** prior score under the cell's prior P
  (ties by hash key); 62–64 materials. S_dis is defined relative to S_top, not
  to the pool: for a well-aligned prior its members can still lie in the upper
  part of the pool's prior ranking.
- **S_blind(P)** (secondary): the members of S_top whose prior percentile in D
  is below 0.80, where the percentile of a material is (average rank of its
  prior score, 1 = lowest, − 0.5) / n. It may be small or empty.

For one run with initial design I and acquisitions a₁…a₄₀:

- **S_dis,run** = S_dis(P) minus the members of I.
- **R_dis** = |{a₁…a₄₀} ∩ S_dis,run| / |S_dis,run|.
- **Chance recall** = 40 / (n − 10), the expected recall of 40 uniformly
  random acquisitions (0.0159–0.0164 in the primary universe).
- **T99** = the smallest t with aₜ ∈ S_top1; 41 if none.
- **N5** = |{a₁…a₄₀} ∩ S_top|.
- **AUC** = (1/40) Σₜ qₜ, where qₜ = the fraction of D with y ≤ the best y
  among I ∪ {a₁…aₜ}.
- **T95** = the smallest t with aₜ ∈ S_top; 41 if none (descriptive only).

**Paired effects.** Runs are paired on (split, seed); β = 0 is the reference.

- Efficiency gain: ΔT99(β) = T99(β = 0) − T99(β).
- Coverage loss: L(β) = R_dis(β = 0) − R_dis(β), with S_dis from the cell's prior.

**Estimator.** The mean over the 60 confirmatory runs (20 splits × 3 seeds),
every split weighted equally.

**Confirmatory test (split-level t-test).** For a run-level quantity v, the
split mean is the mean over the 3 seeds of a split. The estimate is the mean
of the 20 split means. With v_r and v_c the sample variances of the 10 split
means of the random and the chemical-system type, the standard error is
√((v_r + v_c) / 40). The reference distribution is Student's t with
Satterthwaite degrees of freedom, df = (v_r + v_c)² / ((v_r² + v_c²) / 9),
which lies between 9 and 18. Split type is treated as a fixed stratum. (With
a fixed 18 degrees of freedom the test rejects a true null too often when the
two types differ in between-split variance: 0.0077 at a nominal 0.0042 in
simulation; the Satterthwaite form holds 0.0039–0.0042.)

**Hierarchical bootstrap (descriptive).** As in version 1.0, adapted to the
new structure: B = 10,000, NumPy `default_rng(20260928)`; per split type draw
10 splits with replacement, then 3 seeds with replacement per drawn split.
Its percentile intervals are reported beside the t-based intervals. With few
splits it understates between-split variance, which is why it is no longer
the confirmatory test (Appendix C, change 3).

## 7. Hypotheses and tests

### H1 — coverage loss at moderate prior strength (confirmatory)

For each of the 6 cells and each β ∈ {0.25, 0.5} (12 tests):
H₀: mean L(β) ≤ 0; H₁: mean L(β) > 0.

p = P(T_df ≥ estimate / standard error), one-sided. Holm correction across the
12 tests, familywise α = 0.05. Reported for each test: the estimate, the
two-sided 95% t interval, the bootstrap interval, and the Holm-adjusted
decision. H1's content is unchanged from version 1.0; only the test is.
At β = 1 the loss equals the baseline recall by construction (barring large
tie groups in the prior, an S_dis member has at least 61 members of S_top
above it in prior score, and the prior-only search queries the 40 highest),
so β = 1 and β = 0.75 are descriptive.

### H2 — is there a prior strength that is beneficial and safe? (confirmatory decision rule)

Per cell:

- **Applicable** only if the no-prior search finds S_dis members above
  chance: the one-sided lower 95% bound of mean [R_dis(0) − chance recall] is
  above 0. Otherwise the cell is reported as "baseline at chance; H2 not
  applicable".
- For each β ∈ {0.02, 0.05, 0.10, 0.25, 0.5, 0.75}:
  - **beneficial** if the one-sided lower 95% bound of mean ΔT99(β) is above 0;
  - **safe** if the one-sided upper bound, at level 1 − 0.05/6, of mean
    [L(β) − 0.20 · R_dis(0)] is below 0, that is, the search demonstrably
    loses less than 20% of the baseline's recall of S_dis.
- Reported: the set of β that are both beneficial and safe (possibly empty),
  with every bound.

All bounds use the split-level t distribution. The safety bound is taken at
level 1 − 0.05/6 so that the statement "at least one of the six β values is
safe" is wrong with probability at most 0.05 when every β truly loses 20% or
more (with unadjusted 95% bounds that probability is 0.10–0.26 in simulation,
depending on how correlated the six arms are). The
unadjusted 95% bounds are reported beside it as descriptive, with a three-way
state per β: safe (upper bound below 0), unsafe (lower bound above 0) or
inconclusive. Both outcomes, "such a β exists" and "none was demonstrated",
are reported as results. The 20% margin is a judgement fixed before any
result; the rule recomputed with 10% and 25% is reported as descriptive.

### H3 — the frontier and pre-search alignment (exploratory)

For each cell c: xᶜ = the history-side Spearman ρ between prior and objective
(`alignment_v3_primary.csv`, mean over the 20 splits); yᶜ = mean ΔT99(0.5).
Reported: the 6 points and Kendall's τ. No test (n = 6). If a pilot fallback
leaves the two objectives with different efficiency outcomes, the six points
are listed and τ is not computed. No prior is anti-aligned: the lowest cell
mean is 0.53 and the lowest single split 0.45.

## 8. Secondary, exploratory and sensitivity analyses

- **Frontier (primary descriptive result):** for every cell, mean ΔT99, N5,
  AUC, R_dis and L at all eight β values, with t and bootstrap intervals; the
  chance level of R_dis on every plot; pooled recall of S_blind with its size.
- **Recall by prior percentile:** recall of S_top members in the prior-
  percentile bins [0, 0.5), [0.5, 0.8), [0.8, 0.9), [0.9, 0.95), [0.95, 1],
  at every β.
- **Engine health (diagnostic):** per arm, the share of steps whose top-EI
  tie group exceeds half of the candidates, and the share of steps with the
  length-scale at its lower bound.
- **Splits never inspected:** the mean gain and loss at every β, with
  split-level intervals and the H1 p-value, on the 14 splits r3–r9 alone
  (up to 12 degrees of freedom); descriptive.
- **Candidate-level model (secondary):** logistic mixed model of discovery of
  each S_top member on β (categorical) × the member's prior percentile, with
  random intercepts for run and material; software and settings fixed in the
  analysis code. The fit is iterative and reproducible to about 3 decimals;
  its mean-field standard deviations are not confidence intervals.
- **Family level (exploratory):** family = anonymous formula + space group;
  ratio of discovery rates at β = 0.5 vs β = 0 for Heusler (ABC2, #225), for
  all other families together, and for every family with at least 5 distinct
  S_top materials within a split type; bootstrap intervals.
- **Manipulation check:** P_wrong at β = 0.5 for each cell; mean ΔT99 < 0
  expected; descriptive.
- **Random-search reference:** analytic expectations (R_dis = chance recall;
  N5 = 40 |S_top| / (n − 10) ≈ 2.0; P(T99 censored) ≈ 0.66, the probability
  of no S_top1 member in 40 uniform draws), no runs.
- **Budget sensitivity:** the efficiency gain and the coverage loss recomputed
  on the first 10, 20 and 30 acquisitions (T99 censored at t + 1); descriptive.
- **Sensitivity 1** (no Pm/Tc) and **sensitivity 2** (unscreened; P_H cells
  only): the H1 and H2 quantities recomputed and reported descriptively beside
  the primary results. In sensitivity 2, S_dis has 83 members.

## 9. Runs

**Pilot (feasibility only).** Primary universe; **β = 0 only**, so no prior
is used in any pilot search; Y1 and Y2; all 20 splits; **seeds 10–11**;
80 runs.

The pilot computes only these checks, per objective (40 runs each):

| Check | Pass criterion | Pre-declared fallback if it fails, for that objective |
|---|---|---|
| S_dis measurable | for every cell of the objective, mean \|A ∩ S_dis,run\| ≥ 2 (random search gives about 1.0) | S_top becomes the top 10% (S_dis then has about 124 members, random search gives about 2 hits, and N5 counts top-10% hits) |
| T99 measurable | T99 censored (= 41) in < 50% of runs | ΔN5(β) = N5(β) − N5(0) replaces ΔT99 everywhere, including H2 |
| Engine | the two objectives share each (split, seed) initial design; every record validates; 5 runs drawn with NumPy `default_rng(20260928)` reproduce their run hash when rerun | fix the code (a logged deviation that changes the engine and analysis hashes), rerun the pilot with seeds 12–13 |

Pilot runs are not part of any confirmatory analysis. Pilot outcomes may not
change hypotheses, thresholds or metrics beyond the fallbacks above. The
pilot records and `pilot_check.json` are archived with their hashes, and the
analysis reads the fallbacks from that file; they are not chosen by hand.

**Confirmatory.** **Seeds 100–102** on all 20 splits (60 paired runs per
arm). Runs per (split, seed): primary 50 (Y1: 1 + 2 cells × 7 β + 2 P_wrong;
Y2: 1 + 4 × 7 + 4 P_wrong); sensitivity 1: 44; sensitivity 2: 23. In total
3,000 + 2,640 + 1,380 = 7,020 runs.

**Retired seed.** Seed 0 is never used. One discovery run with seed 0
(random_r0, Y1, P_phys_fixed, β = 0.5) was executed during software testing
before version 1.0; it was deleted unread (`EXCLUDED_RUNS.json`,
decision-log entry 25). No run planned in version 1.0 (pilot seeds 10–14,
confirmatory seeds 100–109) was ever executed.

## 10. Analysis code and deviations

- `analyze_runs.py` v1.1 implements Sections 6–8; its hash is in Appendix A
  and in the decision log (entry 35). The batch driver refuses to start the
  confirmatory runs unless the hash of the analysis file it finds is written
  in the decision log and the pilot's engine check has passed; the analysis
  refuses run records that are not listed in a batch manifest.
- Details that this text leaves to the code (tie-breaks, the bootstrap
  indices, the mixed-model settings, the family-level rates) are fixed in the
  header of `analyze_runs.py`; none depends on results.
- Any departure from this plan is logged in `DECISION_LOG.md` with its date and
  reason, and reported in the paper next to the registered analysis. Registered
  analyses are always reported as registered, whatever their result.

## 11. Statement on AI assistance

The study design, the code (data extraction, pool, splits, priors, search
engine, analysis and tests) and the text of this plan were developed with an
AI assistant (Claude, Anthropic). Version 1.0 was audited before any run by
seven AI-assisted reviews (plan, analysis code, data pipeline, priors, engine,
whole-study referee, documents); their findings led to this amendment, and an
eighth review audited the amendment itself. The
author reviewed and approved every design decision recorded in
`DECISION_LOG.md`, ran the cross-machine rebuilds and engine checks of
versions 1.0 and 1.1 on Colab, and
checked the extracted elemental reference values against the printed
handbook tables (with a separate AI-assisted read). The author is responsible
for the content of this registration.

---

## Appendix A — SHA-256 of all frozen inputs

Verify from the repository root with `sha256sum -c` (raw snapshot and
handbook files at the paths shown).

```
4db3c2a0b600da330eeafd18525021971c585e753e051b94cd73e2f5c64763fc  data/raw/mp_2026.04.13/elasticity.jsonl.gz
b485561ade56cd7c15c6e819ac707276114945cb9d6a00ecce3dcb6b70ce91e6  data/raw/mp_2026.04.13/summary.jsonl.gz
dbee15d30fcd8426f3de0acf1b12801923cf46d04b2bb2e37110108ecf208153  handbook/springer-handbook-of-condensed-matter-and-materials-data_compress.pdf
8592e2b369eb8458c0eb78570d715f2ccf733276d6757a60e4f65dfa13d8be70  code/analyze_runs.py
9d0f4f31ec46d1bd3a6489fc818936e1c387dabdec5126b9b8a9e74cbbf243e6  code/build_pool.py
d2fa186eb7488c1091d45a56f479cf7ff5d9f7085c6a58972c692f88154bf8b5  code/build_priors.py
073225e454a5c211c5b69e033f5daed13432c22c5e73c82783a8cda199352db1  code/build_priors_v3.py
00a2bc777432a4ae3853ff6d41d92df6b0325bfb5fd6deddc6366cd4b507f0f0  code/engine_design_probe.py
65454f097428a155d7ecf8b0b0100d8c1f287476c2aa4fb224ca6df129291123  code/expected_hashes.txt
88d6b326da1932a464ef958073a9b434effe4dbf94047753c76c25d66e2c890f  code/expected_hashes_v1_1.txt
f7146cece6936b437f73a54ccd078624b50549094c0baa7e18f9b8d7c8517f30  code/extract_reference.py
8d30b0b10f37ff5f2e2a60916622cef167ccb7059d403506a50a2d95c5065926  code/make_splits.py
f397ab295a6b2783cf2069542a5759a2f037c838bc642126b3bcdef51b22cb49  code/make_splits_v2.py
4ccd25f9a73845fc4f2525241c0d20ab7977bb4011b4ccb7baf23ea4d8daa65a  code/pilot_check.py
3a31c2e65d463e6bb20c4b8552f710f8b1e7a3a993b4754e4c35e8ad64908e7e  code/pull_mp.py
f49cec604475ea0ccab9063f2b030cdce4c9cc069a32e7e97bc3f6f05b757565  code/rebuild_and_verify.sh
6c7f0bddd77147a45e857283fcd82a45b42e12f91ecf278f7300cdb689656171  code/rebuild_v1_1.sh
862a40d5b6300e404cc60c74b764050b921e536af8eef1d68dfca0dffed19141  code/requirements-analysis.txt
d7606a8aba2909f0add4c68c2c681e9f196c418aa1f7c5c558ce33f968dfabc8  code/requirements-build.txt
97a17fe4f7daa68da09db31b603218286006e337d9997cdbe0d19142b4be59fd  code/run_batch.py
6d900ad9af93c38019e8be2ee7e0a3037e82d1e649678f8f4dededa41c36f13a  code/run_search.py
67620ab5d6e27ed3b9fa66f59eb96a680f23f61b46319fdc84fca104d019373d  code/test_analysis.py
a58cd38b3c745ce8d2f23a8a3db4bb496a5f68a0d94344a44e6227b8a162be31  code/test_engine.py
07f6b3221744a0397fc6f348862d78c904cca0dd9db678e569cfda08322f1ece  code/test_pilot_pipeline.py
644598ad65d161db0398854092db38f4097be68439c62995f2524c1fc35ff2f5  data/derived/alignment_v2_primary.csv
7e93481ca8b5f5dceace6cc417875e30a3cb7e0a9f64e1de803928f5c2d28de8  data/derived/alignment_v2_sens1_no_Pm_Tc.csv
1413095c8763e7f73afa7460154333379ddd77450024da746c587e8c3af4442e  data/derived/alignment_v2_unscreened.csv
398f82345a39cbd20caff1fe99bd02c699dfd49cafd8e59bb2a041fad2e87302  data/derived/alignment_v3_primary.csv
c6db38dc0b718134b39c48732f16dcdd9688022020780078e0432065eb2c601c  data/derived/alignment_v3_sens1_no_Pm_Tc.csv
cc75e1d6dbedd64c48555294b482421880bfd5c3189b5aeb0b660a5b8076b795  data/derived/alignment_v3_unscreened.csv
fb9363a5b741596de99e902a2cd6a79c216145afe85180b0365c9641ebc70898  data/derived/elemental_reference_springer2005_v1.csv
0d2aa4206286f44ed528b8bdea1518fa37b388bec1c246d9e9150a4dbd598b6b  data/derived/elemental_reference_springer2005_v1.manifest.json
7ffeceecd2c85c749156f4ea33b04499f3cb5e3d394113cf4517ccda7ebbfeef  data/derived/engine_design_probe_results.csv
8ff6c835c1a3a087ede1cd91a9c047857612f6a0ef79b608697df1f6b2d9b30f  data/derived/features_v2_primary.csv
33526d4b743bc1514746e26cea485df5849fbddf06a5b020458987ed6c63e256  data/derived/features_v2_sens1_no_Pm_Tc.csv
4039be545a0ac4c8665fd625c69b788680f387269c9749ca40625dd742a1dd26  data/derived/features_v2_unscreened.csv
fb5873dabdc0cf462d6678fe47be30ab2f21e8fb659604cb208cecd3162f9540  data/derived/pool_v1.csv
ef052f539c38409d2632dcfdc799e4b75d670d5cff890050d5df784772a017c3  data/derived/pool_v1.manifest.json
de86901c3b806cc6ad236e7849f429928ffab7a4f9443cbf7f5cb3238a7b6b0a  data/derived/priors_v2_primary.csv
ce08cd5faf7fc29674bd2c4476d1896eec5b9a1f94ad038e3a2c5bae23bce9ba  data/derived/priors_v2_primary.manifest.json
6f9dc7928ba6f2011abc27f84da81ea5bc0b911b93ffa3a8ec89d1bbeec604a7  data/derived/priors_v2_sens1_no_Pm_Tc.csv
16d587f8af31992d0cb9cef1e5e853815cca79603d860d73f42cb294bca6b36a  data/derived/priors_v2_sens1_no_Pm_Tc.manifest.json
4f33ec37a24e3786a284ad88a354ee4884bd98a394e23c9b909c3bdfb0defd0e  data/derived/priors_v2_unscreened.csv
1a13301b25c9908390aca2537e27784e9f331a04cf0b37b99a2de0554e4fd3df  data/derived/priors_v2_unscreened.manifest.json
03ef265b21cc530f270d06e4f4f4e7456f2f1771b7701e032f02dce62cb95043  data/derived/priors_v3_primary.csv
3a81ee7813a0f98a5a6d0b33d8397ae2dfabcaf9b7b96cbb689a1c1e63da6195  data/derived/priors_v3_primary.manifest.json
6e23b19c6258ffbaeb6a7141942db5462b5ecdef0d109ff762fcfb277e2d150b  data/derived/priors_v3_sens1_no_Pm_Tc.csv
052f19c7df3add07ba98a3a1a7cad28b57ddaaf6cdbf96bde0b8738c35897736  data/derived/priors_v3_sens1_no_Pm_Tc.manifest.json
8e2e8398129128234fa9ba6a2c3ecdc70a65a4ccc09203b0d31eacdbc783fb5b  data/derived/priors_v3_unscreened.csv
006e8fe2fd26b3bda1709c4c756c813340ba3180354bea81b330523451ed489c  data/derived/priors_v3_unscreened.manifest.json
d4b473a7900c60363c7c93b0aef7633a2f330b533ed93c4eb64954a4e5f8d92a  data/derived/splits_v1.csv
d857fe0c615bb8305e095d7f6de26c2374a115f2f99c2b66e4ec9f1ef18e9bb7  data/derived/splits_v1.manifest.json
6b6142ddeaeab2a9be14a58e88c34f4e506cdefd739101c3ac1208d945b673e1  data/derived/splits_v2.csv
62cdc533403aebb0f5eeec1c4fc0db47386c61a0236f97c0303c6a0889be2708  data/derived/splits_v2.manifest.json
24721262f9e0501d8074af11566368d72d6fc63fad9819bc1ba591bd182273a3  docs/reference_verification_sheet.csv
5caa1e2425ef405bbb7c48383f26eca09d407b86ab0154d6a318b5a9fb505469  docs/EXCLUDED_RUNS.json
dc948f8bae81526fe6aa750c0e49f5aa7bc1568ae4c030f761950115be30b4e5  docs/DECISION_LOG.md
4e913ffa6ff63c104072e72256249a0db5d5dc0542bea303c98a1c6379eec8ca  docs/registered_v1.0/PREREGISTRATION.md
048c7cde729add2323fa97896a859511531ee44572b224e089f71b30d3431c36  docs/registered_v1.0/DECISION_LOG.md
```

The raw snapshot and the handbook are not in the repository: the snapshot is
on Zenodo (data record above). The handbook hash identifies the copy used;
the committed reference table, not the PDF, is the verifiable input.
`docs/registered_v1.0/` holds byte-identical copies of the version 1.0 plan
and decision log.

## Appendix B — Reference run hashes (engine v1.1)

One 40-step search per β on a history pool (random_r0, Y1, elemental-modulus
prior, seed 1). They are asserted by `test_engine.py`, and are identical on
Python 3.11 and 3.13 and on Colab. The β = 0 and β = 0.5 hashes equal those of engine v1.0,
because that run never fits a length-scale below 3.

```
β = 0.0   79a48ce1eb7b349f4e823f628e57edc0523f6441be09165ab8db25a7c6a062ed
β = 0.02  d212dce81e7fc194a03736290786cc371c8794be18f50204cbcd2cc9cd4fe37b
β = 0.05  1677aa4616b2998201005cf2117cd401d0203bbe0a483aeb47ae7d0a037c9c40
β = 0.1   27d9d2756e1333d817e6466568b2fb0dde3efaa5a8f1a5c5b69f2d842ef4e607
β = 0.25  3d038019f5fbf255d2c2ba79c4d0e3d5d2501716d347631eb61cca6b6b1e67da
β = 0.5   28330f33123d57aa9d4ab450f2295e80089448e5479bd3e61fb3ee4884291107
β = 0.75  40a378e11125d0b189729620de62e46025089489e6cd913b7896ef3d9010a657
β = 1.0   0cbe8e2100adaf1f8eab65d6a679e022bc204435bbeef720595a31a1e696d386
```

## Appendix C — Changes from version 1.0, and errata

Every change was decided before any model-driven search on a discovery pool,
on the evidence named. Evidence files: `docs/amend_evidence/`,
`docs/INDEPENDENT_REVIEW_2026-10-03.md`, decision-log entries 31–36.

| # | What | Version 1.0 | Version 1.1 | Reason and evidence |
|---|---|---|---|---|
| 1 | GP length-scale lower bound | 0.01 | 3.0 | With a near-identical pair of materials observed, the fit collapsed to the bound, EI was equal for all candidates, β = 0 picked in hash order and β > 0 by the prior alone. History pools, β = 0: 40 of 2,400 steps flat in normal runs at 0.01, none at 1, 2 or 3; with such a pair forced into the design, top-5% hits were 1.8 (0.01), 10.9 (1), 15.6 (2), 16.1 (3) against 13–14 without the pair. One version 1.0 confirmatory design (primary, chemsys_r1, seed 104) contained such a pair |
| 2 | β grid | 0, 0.25, 0.5, 0.75, 1 | adds 0.02, 0.05, 0.10 (secondary arms) | Rank mixing takes hold early: on history pools the β = 0.02 search already leaves the EI maximiser in 58% of steps. The frontier needs points below 0.25 |
| 3 | Replicates and confirmatory test | 3 + 3 splits × 10 seeds; hierarchical bootstrap p-value | 10 + 10 splits × 3 seeds; split-level t-test with Satterthwaite degrees of freedom (9–18); bootstrap descriptive | Simulation: with 3 splits per type the bootstrap rejects a true null at 0.020 (nominal 0.0042) when 20% of the variance is between splits; the new design holds 0.004 with power 0.61 at a standardised effect of 0.5, against 0.12 for a split-level test on 6 splits. Satterthwaite rather than a fixed 18: see Section 6 |
| 4 | H1 | 12 tests, bootstrap p | 12 tests, same cells and β, split-level t | Content unchanged; test as in change 3 |
| 5 | H2 | retention ≥ 0.80 of the β = 1 speed-up and upper bound of L < δ = 0.016 | applicability gate (baseline above chance); beneficial (lower 95% bound of ΔT99 > 0) and safe (upper bound of L − 0.20·R_dis(0) < 0 at level 1 − 0.05/6) | δ = 0.016 equals the recall of random search; the β = 1 anchor rests on about six independent values; retention had no uncertainty; six β values are examined, so the safety level is adjusted |
| 6 | Coverage sets | S_dis | S_dis unchanged (described as the lower-prior half of the top 5%); S_blind and recall by prior-percentile bin added as secondary | S_dis is relative to S_top; for well-aligned priors its members are not low in the pool's prior ranking |
| 7 | Pilot | 150 runs incl. two priors at β = 0.5 and 1; Y1; S_dis pass mark 1 hit; rerun seeds 15–19 | 80 runs, β = 0 only; Y1 and Y2; pass mark 2 hits; fallbacks per objective, read by the analysis from `pilot_check.json`; rerun seeds 12–13 | The checks need β = 0 only; 1 hit per run is what random search gives; 4 of 6 cells are Y2 |
| 8 | Runs | 4,320 (seeds 100–109) | 7,020 (seeds 100–102) | Changes 2 and 3 |
| 9 | Engine record and guard | — | thread settings and top-EI tie count recorded; DOI-format guard; no overwrite; explicit oracle check | Review findings |
| 10 | Analysis | — | stricter record validation (scikit-learn version, thread settings, step log, batch manifest, mandatory descriptor hash); frontier reports R_dis and t intervals | Review findings |
| 11 | Other | optimiser start 1.0; two thread variables set; family table for families with ≥ 5 S_top memberships | optimiser start 3.0 (it must lie inside the bounds); `MKL_NUM_THREADS` also set to 1; families with ≥ 5 distinct S_top materials; new descriptive analyses: S_blind, recall by prior percentile, engine health, random-search reference, budget sensitivity, 10% and 25% margins, splits never inspected | Consequences of changes 1–3 and review findings |

Unchanged: the raw snapshot, the pool and its three universes, splits r0–r2,
descriptors, the priors' definitions and their scores on r0–r2, the
objectives, the budget, rank mixing, S_top, S_top1, S_dis, T99, H1's content,
H3, the retired seed 0.

**Errata to the text of version 1.0** (no analysis was affected):

- Section 5 said "132 Magpie composition statistics" and implied 140 GP
  inputs. The descriptor files hold 131 Magpie columns (one constant column
  was dropped at build time), and 139 or 140 inputs remain after constant
  columns of a discovery pool are dropped.
- Section 7 gave the lowest history-side ρ as "0.48 as a cell mean". The
  lowest mean over the six splits was 0.52; 0.48 was the mean over the three
  chemical-system splits of that cell.
- Section 7 said H2 could fail because recall moves in steps of about 1/62.
  The mean over 60 runs has far finer steps; the standard error is what
  limits it.
- Section 5 wrote "ConstantKernel(1.0)"; the amplitude is fitted.
- Section 2 said the rebuild script regenerates "every derived file"; it
  covers the 12 derived data files, not the engine-probe results or manifests.
- Sections 0 and 2 said the handbook is "cited, not redistributed"; values
  extracted from its tables are included in the reference table.
- Section 1 said "all from a metallic-element list"; the list is the
  exclusion list now printed in Section 1.
- The registered copy of `DECISION_LOG.md` carried an unfilled drafting note
  about the AI-assistance statement.

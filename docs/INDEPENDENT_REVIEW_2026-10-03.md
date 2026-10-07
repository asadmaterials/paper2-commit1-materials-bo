> Note added 2026-10-03: this report led to registration v1.1. Its Section 2.2 ("H1 is close to true by
> construction") was later corrected by measurements on history pools; see AMENDMENT_PLAN_v1.1.md,
> Section 1. H1 is unchanged in v1.1.

# Paper 2 — independent whole-project review (2026-10-03)

Five fresh reviewers each audited one part of the project without seeing the
earlier findings: (1) data pipeline and splits, (2) priors and information
boundary, (3) search engine, (4) the study as a journal referee, (5)
documents, repository and records. Their main claims were then checked
directly; the "Checked" column says which.

- **Verified** = reproduced by the assistant after the review.
- **Two reviewers** = found independently by two reviewers, not rerun.
- **Reviewer** = one reviewer's result, not rerun.

No discovery-side search was run and no discovery-side outcome was computed
by anyone. The reviewers did run searches on **history** pools, including
arms with the elemental-modulus prior; see Section 6.

---

## 1. Verdict

The engineering is sound: every hash matches, the funnel and splits were
reproduced by independent code, the priors were recomputed exactly, and no
discovery-side outcome reaches any prior. **The confirmatory layer is not
sound as registered.** Four problems, each sufficient on its own to weaken
the paper, sit in the engine, the main hypothesis, the safety rule and the
statistics. All four can still be fixed, because no discovery-side search
result exists.

## 2. Problems that need an amendment (v1.1) before any run

### 2.1 The engine collapses, and one registered block already has the trigger

| Fact | Checked |
|---|---|
| When two near-identical materials with different values are both observed, the GP length-scale falls to its lower bound (0.01), expected improvement is equal for all candidates, β = 0 picks in hash order and every β > 0 picks by the prior alone | Verified (2 of 60 history runs; one collapsed for 39 of 40 steps) |
| The confirmatory initial design for primary / chemsys_r1 / seed 104 contains two Li2ZnGe polymorphs (mp-aaabghir, mp-aaaaasjj) at descriptor distance 0.025. It is the only one of the 360 registered designs with a pair closer than 0.27 | Verified |
| History pairs at that distance collapse the GP in about 70% of fits; a collapsed block compares prior-only search against hash-order search (history test: ΔT99(0.5) of 31 against a normal 9.7) | Reviewer |
| A second, transient route needs no near-duplicate: with few points the likelihood prefers "no structure" (20% of Y1 history runs at step 1) | Reviewer |
| Prior arms do **not** co-select polymorphs more often than β = 0 (0.40–0.55 pairs per run against 0.80), so the feared arm asymmetry was not found; the bias comes through the shared initial design and spontaneous episodes | Reviewer |
| Remedy tested on history data only: length-scale lower bound 3.0 instead of 0.01 removed every flat state in 276 runs, left 45 of 60 healthy Y1 trajectories bit-identical and did not change β = 0 quality (N5 −0.07 ± 0.18). A noise term or larger `alpha` did not remove it | Reviewer |

Direction of the bias: toward larger ΔT99 and larger coverage loss, in favour
of H1.

### 2.2 H1 is close to true by construction, and the β grid starts too high

| Fact | Checked |
|---|---|
| An S_dis member has at least 62 top-set members above it in prior score, so it lies outside the prior's top 2.5%; β = 1 queries only the top 40 (1.6%). L(1) = R_dis(0) exactly | Two reviewers (arithmetic) |
| With rank mixing, β = 0.25 already behaves like a near-hard prior filter: on history runs the β = 0.25 pick sits at the median 96.7th prior percentile and abandons the EI maximiser in 88% of steps | Reviewer (measured, history) |
| Consequence: at β ≥ 0.25, L(β) ≈ R_dis(0), so H1 mostly tests whether the no-prior arm finds any S_dis member. The trade-off itself lives at β of roughly 0.01–0.10, which the design does not run | Reviewer (synthetic rank simulation) |

### 2.3 The H2 margin and the pilot criterion equal random search

| Fact | Checked |
|---|---|
| Random search recalls 40 / 2,490 = 0.016 of S_dis, which is δ; its expected S_dis hits per run are 1.0, which is the pilot pass mark. Random search would pass the pilot check about half the time | Two reviewers (arithmetic) |
| "Safe" is declared vacuously if the baseline is near chance, and is nearly impossible otherwise (mean loss must be under about half a material per run). Retention is a ratio of point estimates with no interval; G₁ rests on about six independent values | Two reviewers |

### 2.4 The bootstrap rejects a true null too often

| Share of variance between splits | Registered bootstrap, P(p ≤ 0.05/12), nominal 0.0042 | 95% interval coverage |
|---|---|---|
| 0 | 0.001–0.002 | 0.97 |
| 0.2 | 0.015–0.019 | 0.90 |
| 0.5 | about 0.03–0.04 | 0.85–0.87 |

Two reviewers ran this simulation independently and agree. A t-test on split
means stratified by split type (4 degrees of freedom) held its nominal level
at every setting. The cost is power: the critical t at 0.05/12 is 4.85.
More replicate splits (for example 10 per type with 3 seeds, the same number
of runs) would restore power but changes frozen inputs.

### 2.5 Smaller items that an amendment should carry

- Text errors in the published plan (Verified): 131 Magpie columns, not 132;
  139–141 inputs after constant-column removal, not always 140; lowest cell
  mean of history ρ is 0.52, not 0.48; "ConstantKernel(1.0)" is fitted;
  "regenerates every derived file" omits the probe results and manifests.
- The pilot runs β = 0.5 and β = 1 on the confirmatory pools; only β = 0 is
  needed for its three checks. The pilot checks Y1 only, although 4 of 6
  cells are Y2, and the plan does not say whether a fallback applies per cell.
- The engine's guard accepts any non-empty string, not only this
  registration's identifier (the analysis does check the exact DOI).

## 3. Cannot be fixed; must be disclosed

| # | Limitation | Checked |
|---|---|---|
| 1 | The benchmark, metrics and hypotheses were developed with outcome labels of all pools in view. The β = 1 arm is a deterministic function of data already inspected (decision 16 computed it for 2 splits). "Registered before any discovery-side result" is accurate only for the GP-driven arms | Two reviewers |
| 2 | The engine probe used history sides of 4 splits, which cover 75.4% of the pool and 65–75% of every discovery pool | Two reviewers |
| 3 | The universe is filtered by the quality of the elemental prior's own input data (1,161 structures, 24.6%, removed) | Two reviewers |
| 4 | The 6 splits re-partition one pool; any two discovery pools share about 70% of materials. Inference is conditional on this pool | Two reviewers |
| 5 | The chemical-system split is only mildly out-of-distribution: 29–31% of its discovery materials have a sub-system in history, and 2,501 of 3,299 systems hold one structure | Reviewer |
| 6 | Descriptors include the DFT-relaxed density and volume of each candidate, so the setting is "screening relaxed structures". Volume per atom alone reaches \|ρ\| of 0.62–0.72 with the objectives on history | Reviewer |
| 7 | P_H uses about 1,050 labelled history materials while the β = 0 GP starts from 10 points: "prior against no prior" is partly "labels against none" | Reviewer |
| 8 | A plain random forest on the descriptors reaches history-side ρ of 0.88, above both priors (0.66 and 0.79) | Reviewer |
| 9 | No thermodynamic-stability filter: 304 primary structures lie more than 0.1 eV/atom above the hull | Reviewer |
| 10 | "MP density" is the elasticity document's density; it differs by more than 1% from the summary density for 3,444 of 4,709 structures | Reviewer |
| 11 | Scope wording: B, Si and Se are excluded but Ge and Sb are in (330 and 214 primary structures) | Two reviewers |
| 12 | The reference screen checks internal consistency of one handbook row, not accuracy. Six eligible elements pass one check and fail the other; some kept values look low against the handbook's own elastic constants (Ge) | Reviewer |
| 13 | Fixed-weight rank mixing is not what the literature uses (πBO and successors decay a multiplicative prior); β is not comparable to those weights | Two reviewers |

## 4. Defects in frozen files (fixable only by amendment; small effect)

| # | Defect | Checked |
|---|---|---|
| 1 | The Born test runs on MP's integer-rounded tensor with zero tolerance. 19 structures are rejected and 2 admitted (LiIn mp-aaaabhfw, Li3Mg mp-aaacdoei, eigenvalues about 1e-15) on round-off; on the unrounded tensor the pool would be 4,730, not 4,711 | Verified for the 2 admitted cases and the integer rounding |
| 2 | The reference table's Si row is wrong (graphite's modification and density; a column collision in the carbon table). Si is not in any pool, so nothing downstream changes, but the file is hashed and "61 eligible" counts it | Verified |
| 3 | `rebuild_and_verify.sh` aborts intermittently at its `pdftotext -v \| head -1` line (6 of 30 trials) | Reviewer |
| 4 | The registration hashes a handbook file named `…_compress.pdf`. A reader with a publisher copy will get a different hash and possibly a different text layout, so step 1 of the rebuild is not reproducible by others; the committed CSV should be the verifiable input | Reviewer |
| 5 | "Exactly one summary document" really means "has a summary document"; 30 of the 45 structures dropped there contain Yb | Reviewer |
| 6 | Stale comments in hashed scripts ("pre-registered", "Y2 (G/rho)", the retired S_dis definition, "v2") and an unfilled AI-statement note in the registered decision log | Two reviewers |

## 5. Free to fix now (uncommitted or metadata)

**Analysis code and tests**
- H2 retention of exactly 0.80 can evaluate to 0.7999999999999999 and be
  marked not safe (Verified).
- The tests do not check the H1/H2/H3 wiring (11 of 14 planted bugs passed);
  the code itself matched an independent recomputation.
- `test_engine.py` never asserts the reference hashes and has no test or
  reference hash for β = 0.25 and 0.75: mutants that swap those two arms pass
  all 13 tests and reproduce both Appendix B hashes (Reviewer).
- Validation gaps: scikit-learn version, mandatory descriptor file, manifest
  cross-check; `pilot_check.py --no-rerun` reports a pass.
- Decision-log entry 30 is wrong on three points: "frozen" (the file is
  uncommitted), "no file changes under a fallback" (false for an engine fix),
  and the mixed-model precision. The Python 3.13 claim does hold for the
  deterministic outputs of the current file (Verified today).

**Repository and records**
- No LICENSE; no Materials Project attribution (CC BY 4.0; Jain 2013, de Jong
  2015) although `pool_v1.csv` redistributes MP values.
- README at the registered commit lists the wrong contents; `SHA256SUMS.txt`
  does not cover the plan; the superseded hypotheses draft is in the repo
  unmarked and contradicts the plan (pilot seeds 0–4, H1 at β = 0.5 only).
- Commit 4 as planned would make `sha256sum -c` of Appendix A fail on the
  decision log without a visible explanation; keep a byte-identical v1.0 copy
  or tag the registered commit.
- "The handbook is cited, not redistributed" should say that extracted values
  are included with citation.
- Zenodo metadata: affiliation (PIEAS on the snapshot record, "independent
  researcher" in the plan); check the snapshot README for unfilled
  placeholders; note that the record's files are flat.
- A stale `pull_mp.py` (the version that did not run) sits in the outputs
  folder next to the registered one.

## 6. What the review itself exposed (must be logged)

To test the engine, the reviewers and the assistant ran about 1,000 searches on **history** pools,
including β = 0.25–0.75 arms with the elemental-modulus prior computed for
history materials. These produced history-side efficiency comparisons between
arms (for example ΔT99(0.5) of about 9.7 and ΔN5 of about +5.3 on Y1). No
coverage outcome was computed. History materials are 65–75% of other splits'
discovery pools, so this is information about likely results and belongs in
the decision log with an "outcome data seen: history, by arm" flag.

## 7. What was checked and found correct

- Every SHA-256 in Appendix A, the manifests and the checksum lists matches;
  the rebuild from a fresh clone gives ALL MATCH.
- An independent funnel script reproduces every count and the three universe
  id sets; splits regenerate from their seeds with 0 mismatches; no
  composition or chemical system straddles a split.
- All 14,978 stored elemental-prior scores equal an exact-arithmetic
  recomputation; P_H recomputed bit-identically; no discovery-side outcome
  reaches a prior, a descriptor or a standardisation.
- Expected-improvement formula, score arithmetic, tie-break, initial design,
  oracle log, P_wrong: correct. β = 0 is a credible baseline (N5 of 14
  against 2 for random search on history pools).
- Run counts (150; 1,920 + 1,560 + 840 = 4,320) and set sizes match the plan.
- `analyze_runs.py` implements Sections 6–8 as written in all four switch
  settings.
- No API key or secret in any committed file.

## 8. Literature for positioning (titles confirmed by search; author lists
not all verified)

No study was found that pre-empts this one. Closest: πBO (Hvarfner et al.,
ICLR 2022, arXiv:2204.11051); BOPrO (Souza et al., 2021, arXiv:2006.14608);
ColaBO (Hvarfner et al., ICLR 2024, arXiv:2311.14645); PriorBand (Mallik et
al., NeurIPS 2023); DynaBO (arXiv:2511.02570); dynamic mean decay (ICLR
2026); Xu et al., NeurIPS 2024 (arXiv:2410.10452); Liang et al., npj Comput.
Mater. 2021; Rohr et al., Chem. Sci. 2020; Terayama et al. (BLOX), Chem.
Sci. 2020; Borg et al., Digital Discovery 2023. Surviving contribution: on a
fixed, reproducible DFT elasticity pool, measure how much a well-aligned
prior narrows **which** top materials a GP search finds, not only how fast it
finds one — provided coverage is measured on absolute blind-spot sets and at
small β.

## 9. Options

| | What changes | What it buys | Cost |
|---|---|---|---|
| **A. Keep v1.0** | Nothing registered. Log errata; pre-specify diagnostics and supplementary analyses | Fastest | H1 is near-tautological, H2 near-vacuous, the test is anti-conservative, and at least one block is contaminated. A referee can dismiss the confirmatory part |
| **B. Minimal amendment v1.1** | Engine length-scale bound; split-level stratified t-test as the confirmatory test (bootstrap reported beside it); pilot at β = 0 only; errata. Hypotheses, sets, δ, β grid, seeds unchanged | Valid engine and valid error rate | About 2–3 days. H1 and H2 stay weak as constructs |
| **C. Amendment v1.1 with a repaired confirmatory layer** | B, plus: small-β arms (for example 0.02, 0.05, 0.10); coverage also on absolute blind-spot sets; H2 gated on baseline recall above chance, with a relative margin and an interval for retention; pilot criterion above chance | A trade-off result that is not implied by construction | About a week; more runs; the largest change to what was registered, all of it made before any discovery-side search |

All three keep the pool, splits, priors and the documented exploratory phase.
None can undo Section 3.

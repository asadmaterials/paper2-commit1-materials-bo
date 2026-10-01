# Paper 2 — Hypotheses, metrics and analysis plan (DRAFT v0.2)

**Status:** draft, not yet registered. To be registered **before any optimization
run**. Every design choice made after inspecting data is listed in
`DECISION_LOG.md`; this plan does not claim to precede data inspection.

**Working title:** *The efficiency–coverage trade-off of scientific priors in
Bayesian materials search.*

---

## 1. What is and is not being tested

The pre-search alignment table (`alignment_v2_primary.csv`) already fixes two
outcomes, so they are **not** hypotheses:

- A prior that ranks well-aligned materials first accelerates search when it
  is followed exclusively (β = 1). P_phys alone, taking its top 50, recovers
  ~36% of the top-5% set versus ~2% for random selection.
- At β = 1, top-set materials that the prior ranks low are not found. With
  S_dis defined relative to the prior (Sect. 3), their recall at β = 1 is ~0
  by construction.

These endpoints are reported descriptively as the **prior-only** condition.
The confirmatory questions concern the **interior** of the frontier
(0 < β < 1), where the GP's expected improvement and the prior compete and
the outcome is not implied by the alignment table.

## 2. Setting (frozen inputs)

| Item | Value | Source |
|---|---|---|
| Data | MP database 2026.04.13, raw snapshot | Zenodo DOI [snapshot] |
| Universe (primary) | 3,550 structures, 60 elements | `pool_v1.csv` (`in_primary`) |
| Sensitivity universes | no Pm/Tc (3,378); unscreened (4,711; P0 and P_H only) | `pool_v1.csv` |
| Splits | random ×3, chemical-system ×3; history ≈ 30%, composition-grouped | `splits_v1.csv` |
| Objectives | Y1 = G_VRH; Y2 = E/ρ, E = 9KG/(3K+G) (MP VRH, MP density) | `build_priors.py` v2 |
| Priors | P0 (none); P_H fixed/matched; P_phys fixed/matched; P_wrong = reversed rank (control, full grid only) | `priors_v2_primary.csv` |
| Ground truth | PBE DFT elastic tensors (MP). Results concern computed, not measured, moduli. | |

"Materials" in the text means computed MP entries; ~45% of the primary
universe carries MP's `theoretical` flag.

## 3. Search and outcome definitions

**Engine** (`run_search.py` v1.0, frozen). Scikit-learn Gaussian process,
ConstantKernel × Matérn-5/2 with **one isotropic length-scale** over the 140
descriptors standardised on the searched pool (bounds 1e-2–1e3), jitter 1e-6,
normalised y, 2 optimiser restarts, random_state 0 (choice: decision-log entry
23, history-only probe). Expected improvement over the best observed y
(ξ = 0), rounded to 9 significant digits before ranking. Selection score
S_β = (1 − β)·rank(EI) + β·rank(P) over the remaining candidates (average
ranks), greedy, one query per step; β = 0 never reads P and β = 1 never fits
the GP. The β = 0 run of a (split, seed, objective) is shared by all priors.
The engine obtains y only through a logged oracle, and refuses to search a
discovery pool without a registration identifier, which it records.
Prior scores are rounded to 12 significant digits, so mathematically equal
scores are exactly equal on every platform. Ties in rank(P) (polymorphs under
P_phys; different compositions whose elements share a tabulated G, e.g. Y/Dy,
Er/Ge) are broken by sha256("<split>:<material_id>"). β ∈ {0, 0.25, 0.5, 0.75, 1}; β = 1 is
labelled *prior-only*.

**Budget.** 10 initial points + 40 acquisitions. The initial design is drawn
uniformly from the discovery pool with a seed determined by (split, seed
index) and is **identical across all priors and β**. Initial points are
excluded from every outcome below.

**Sets (per split, per objective, computed on the discovery pool).**
- S_top: top 5% by y (124–127 materials).
- S_top1: top 1% by y (25–26 materials).
- S_dis(P): members of S_top whose prior score is at or below the median
  prior score within S_top (hash tie-break); about half of S_top. Members
  in the initial design are removed from S_dis for that run.

**Run-level outcomes.**
- **T99**: index (1–40) of the first acquisition that is in S_top1;
  censored at 41 if none.
- **N5**: number of acquisitions in S_top.
- **R_dis**: fraction of S_dis acquired.
- Paired effects against β = 0 on the same split, seed and initial design:
  ΔT99(β) = T99(0) − T99(β) (positive = faster),
  ΔR_dis(β) = R_dis(β) − R_dis(0) (negative = coverage loss).

## 4. Confirmatory hypotheses

Cells: Y1 × {P_H, P_phys} and Y2 × {P_H fixed, P_H matched, P_phys fixed,
P_phys matched} = **6 prior–objective cells**.

**H1 — coverage cost at moderate prior strength.**
In each cell, ΔR_dis(0.5) < 0.
*Why it can fail:* for P_phys on Y1, S_dis members still sit between the
51st and 96th prior percentile of the discovery pool, so a 50/50 mix with EI
can reach them; the GP may also find them on its own.
*Test:* one-sided, hierarchical bootstrap (resample splits within split
type, then seeds within split) of the mean paired difference; Holm
correction across the 6 cells; α = 0.05.

**H2 — does a safe prior strength exist?** (two-outcome question; no
directional prediction)
In each cell, a β ∈ {0.25, 0.5, 0.75} is **safe** if both
(a) median ΔT99(β) ≥ 0.8 × median ΔT99(1), and
(b) ΔR_dis(β) is non-inferior to 0 with margin δ: the lower 95% bootstrap
bound of the mean ΔR_dis(β) exceeds −δ.
δ = 1 / |S_dis| (one S_dis material per run, ≈ 0.016). **[to confirm]**
*Reported outcome:* for each cell, the set of safe β (possibly empty). Both
"a safe β exists" and "no safe β exists" are informative results.

**H3 — the frontier depends on pre-search alignment.**
Across the 6 cells, the ordering of median ΔT99(0.5) follows the ordering of
**history-side** Spearman ρ (`side = history_cv`). Reported as Kendall τ
with the six points shown; n = 6, so **descriptive, no significance claim.**
*Retired from v0.1:* the prediction that a prior's net value changes sign
when the objective changes. The alignment table has no anti-aligned prior
(lowest history-side ρ: 0.45 for a single split, 0.48 as a cell mean), so a
reversal is not expected. If one occurs it
is reported as unexpected.

## 5. Secondary and exploratory analyses

- **Candidate-level model (secondary).** Logistic GLMM,
  found ~ β × prior percentile + (1 | run) + (1 | candidate), over S_top
  members; the β × percentile interaction relative to β = 0.
- **Secondary outcomes:** N5, best-so-far AUC, T99 at all β.
- **Family level (exploratory).** Family = anonymous formula + space group.
  For families with ≥ 5 S_top members (pooled within split type), rate ratio
  of discovery at β = 0.5 vs β = 0 with run-level bootstrap CIs. Heusler
  (ABC2, #225) and non-Heusler reported separately for H1 as well.
- **P_wrong (full grid only).** Manipulation check: ΔT99(0.5) < 0 expected.
- **Sensitivity analyses:** (i) primary analyses in the no-Pm/Tc universe;
  (ii) P0 and P_H on the unscreened universe.

## 6. Pilot and GO/KILL (feasibility only)

**Pilot:** Y1; P0, P_H fixed, P_phys fixed; β ∈ {0, 0.5, 1}; all 6 splits;
seeds 0–4. **Pilot seeds are not reused** in the confirmatory runs
(confirmatory seeds start at 100).

The pilot checks that the outcomes are measurable. It does **not** test
effects, and pilot effect sizes may not be used to change hypotheses,
margins or metrics.

| Check | Criterion | If it fails (pre-declared) |
|---|---|---|
| Engine tests | identical initial designs across arms; no discovery y reaches any prior; reproducible run hashes | fix code, rerun pilot |
| S_dis measurability | at β = 0, mean ≥ 1 S_dis acquisition per run | S_top becomes top 10% |
| T99 measurability | at β = 0, censoring < 50% of runs | primary efficiency outcome becomes N5 |

## 7. Inference notes

- Replicate splits share materials; runs are **not independent**. All CIs
  come from the hierarchical bootstrap, never from tests that treat runs as
  i.i.d.
- The frontier and all claims apply to the **screened universe** defined by
  the elemental reference screen.
- Wording: "compiled elemental shear moduli (Springer Handbook 2005, mainly
  from Landolt–Börnstein), screened with internal-consistency checks";
  "registered before any optimization run"; "computed materials"; "PBE DFT
  ground truth".

## 8. Open items before registration

1. Confirm δ (H2 margin) and the 0.8 efficiency fraction.
2. Number of confirmatory seeds per split (proposal: 10 → 60 runs per arm per
   split type).
3. Engine details (Week 3): GP hyperparameter priors, handling of 140
   descriptors with 10–50 points, EI ties at numerical zero.

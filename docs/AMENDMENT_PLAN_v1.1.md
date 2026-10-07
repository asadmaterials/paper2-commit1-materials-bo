> Note added 2026-10-03: the proposal below was adopted with three refinements after a further audit
> (decision-log entry 35): Satterthwaite degrees of freedom instead of a fixed 18; the H2 safety bound at
> level 1 - 0.05/6; pilot fallbacks read from pilot_check.json. PREREGISTRATION.md v1.1 is authoritative.

# Paper 2 — proposed amendment to the registration (v1.0 → v1.1)

**Status:** proposal for the author's decision, 2026-10-03. Nothing has been
changed yet. No discovery-side search has been run.

**What informed this proposal.** The independent review of 2026-10-03;
synthetic simulations; and searches on **history** pools only (about 470 new
runs for this document, on top of the reviewers' runs). The history-side runs
gave engine-health and selection-geometry measurements, and top-5% hit counts
at β = 0 and β = 0.05. They are listed in Section 7 so they can be logged.

---

## 1. One correction to the review

The review said H1 is "close to true by construction" because rank mixing
makes the lower-prior half of the top set almost unreachable at β ≥ 0.25.
That rested on a synthetic rank simulation. Measured on history pools (Y1,
elemental-modulus prior, 18 runs per β, 720 steps each) the picture is more
gradual:

| β | Pick is the EI maximiser | Prior percentile of the pick: median | 10th percentile | Picks below the 97.5th prior percentile |
|---|---|---|---|---|
| 0.02 | 42% | 0.916 | 0.583 | 65% |
| 0.05 | 27% | 0.945 | 0.662 | 62% |
| 0.10 | 18% | 0.951 | 0.736 | 61% |
| 0.25 | 10% | 0.970 | 0.798 | 53% |
| 0.50 | 6% | 0.976 | 0.852 | 49% |

About half of the picks at β = 0.25 and 0.5 still fall in the prior range
where S_dis members sit, so H1 at those β is **not** vacuous and can stay as
registered. What the table does confirm is that the prior takes hold very
early: at β = 0.02 the search already abandons the EI maximiser in 58% of
steps. The frontier therefore needs arms below 0.25, as secondary arms.

## 2. Proposed changes

### 2.1 Engine (v1.0 → v1.1)

| Change | Before | After | Evidence |
|---|---|---|---|
| Length-scale lower bound | 0.01 | **3.0** | see below |
| β grid | 0, 0.25, 0.5, 0.75, 1 | adds **0.02, 0.05, 0.10** | Section 1 |
| Run record | — | adds BLAS thread settings and the size of the top-EI tie group per step | reviewer request |
| Guard | any non-empty `--registered` string | must equal the v1.1 registration DOI | reviewer finding |
| Oracle check | `assert` | explicit error; refuses to overwrite an existing record | reviewer finding |

Length-scale bound, history pools, β = 0:

| Lower bound | Normal runs (60): flat-EI steps | Change in top-5% hits vs v1.0 | Identical trajectories | Near-duplicate pair forced into the design (12 blocks): top-5% hits |
|---|---|---|---|---|
| 0.01 (v1.0) | 40 of 2,400, in 2 runs | — | — | **1.8** (10 of 12 runs collapsed) |
| 1.0 | 0 | +0.03 ± 0.03 | 58 of 60 | 10.9 |
| 2.0 | 0 | +0.20 ± 0.33 | 53 of 60 | 15.6 |
| 3.0 | 0 | −0.22 ± 0.14 | 49 of 60 | 16.1 |

The same 12 blocks without the forced pair give 12.8–14.0 hits. A bound of
1.0 stops the collapse but leaves a degraded model; 2.0 and 3.0 both restore
it. **3.0 is proposed** because discovery pools are about 2.4 times denser
than the history pools tested here, so the larger margin is the safer choice,
and because the reviewer's 276 runs at 3.0 showed no flat state. Its cost on
normal runs is within noise.

### 2.2 Replicate structure and the confirmatory test

| | v1.0 | Proposed |
|---|---|---|
| Splits | 3 random + 3 chemical-system | **10 + 10** (r0–r2 unchanged, r3–r9 added with the same seeded code) |
| Seeds per split | 10 (100–109) | **3 (100–102)** |
| Blocks | 60 | 60 |
| Confirmatory test | hierarchical bootstrap p-value | **t-test on split means, stratified by split type (18 degrees of freedom)**; the v1.0 bootstrap is reported beside it |

Simulation (2,000–4,000 data sets per row; total variance 1; "share" is the
share of variance between splits):

| Design and test | Share | False positives at 0.05/12 (nominal 0.0042) | Power at effect 0.5, level 0.05/12 |
|---|---|---|---|
| 3 × 10, bootstrap (v1.0) | 0.2 | **0.020** | 0.46 |
| 3 × 10, bootstrap (v1.0) | 0.5 | **0.035** | 0.37 |
| 3 × 10, split-level t (4 df) | 0.2 | 0.003 | 0.12 |
| 3 × 10, split-level t (4 df) | 0.5 | 0.002 | 0.05 |
| **10 × 3, split-level t (18 df)** | 0.2 | 0.004 | **0.61** |
| **10 × 3, split-level t (18 df)** | 0.5 | 0.004 | **0.43** |

Keeping 6 splits and only changing the test would be valid but would have
almost no power. Ten splits per type give a valid test with more power than
the invalid one had, for the same number of blocks.

Consequences: new files `splits_v2.csv` and priors for the 14 new splits. For
the new splits **no discovery-side alignment is computed**; only the
history-side (cross-validated) alignment is.

### 2.3 Hypotheses

| | v1.0 | Proposed |
|---|---|---|
| **H1** | L(β) > 0 at β = 0.25 and 0.5; 6 cells × 2 = 12 tests, Holm | **Unchanged in content.** Tested with the split-level t-test |
| **Secondary arms** | — | β = 0.02, 0.05, 0.10 for every cell; frontier reported at all 8 β values |
| **H2** | retention ≥ 0.80 of the β = 1 speed-up, and upper bound of L(β) < δ = 0.016 | **Replaced**, see below |
| **H3** | exploratory | unchanged |

Why H2 must change: δ = 0.016 is what random search recalls (40 of about
2,490), the retention anchor G₁ rests on about six independent values, and
retention had no uncertainty attached.

Proposed H2 ("is there a β that speeds the search up without a material loss
of coverage?"), per cell:

- **Applicable** only if the baseline finds S_dis members above chance: the
  one-sided lower 95% bound of mean [R_dis(0) − chance] is above 0, where
  chance = 40 / (n − 10) for that pool.
- For each β in {0.02, 0.05, 0.10, 0.25, 0.5, 0.75}:
  - **beneficial** if the lower 95% bound of mean ΔT99(β) is above 0;
  - **safe** if the upper 95% bound of mean [L(β) − 0.20 · R_dis(0)] is below
    0, that is, the search loses less than 20% of the baseline's recall.
- Reported: the set of β that are both beneficial and safe (possibly empty).
  All bounds from the split-level t distribution. A pre-specified decision
  rule; no multiplicity adjustment; stated as such.

The 20% margin is a judgement call (Section 3, decision 3).

### 2.4 Coverage sets

- S_dis keeps its definition and is renamed in the text "the lower-prior half
  of the top 5%". It stays the confirmatory coverage set.
- Added, secondary: **S_blind(P)** = members of S_top whose prior score is
  below the 80th percentile of the discovery pool. Reported per cell with its
  size; recall at every β. It can be small or empty for well-aligned priors,
  and that is itself a reported fact.
- Added, secondary: recall of S_top by prior-percentile bin at every β.

### 2.5 Pilot

| | v1.0 | Proposed |
|---|---|---|
| Arms | β = 0, and two priors at β = 0.5 and 1 | **β = 0 only** |
| Objectives | Y1 | Y1 and Y2 |
| Runs | 150 (seeds 10–14, 6 splits) | **80** (seeds 10–11, 20 splits, 2 objectives) |
| S_dis check | mean hits ≥ 1 (equal to random search) | mean hits **≥ 2** for every cell's S_dis (twice the random expectation), checked per objective |
| T99 check | censored in < 50% of runs | unchanged, per objective |
| Engine check | 5 reruns reproduce their hash | unchanged |

Fallbacks stay as in v1.0 and apply per objective.

### 2.6 Runs

Per (split, seed): primary 50, sensitivity 1: 44, sensitivity 2: 23.
With 60 blocks: 3,000 + 2,640 + 1,380 = **7,020 runs** (v1.0: 4,320).
At about 7 s per run on two cores this is roughly 7 hours.

### 2.7 Corrections to the text

131 Magpie columns; 139–141 inputs after constant-column removal; lowest
6-split cell mean of the history-side ρ is 0.52; the kernel amplitude is
fitted; the rebuild script covers 12 derived files; "MP density" is the
elasticity document's density; the element scope is stated as an explicit
list; "the handbook PDF is not redistributed; values extracted from its
tables are included with citation".

## 3. Decisions that are yours

1. **Replicates:** 10 + 10 splits with 3 seeds (proposed), or keep the 6
   splits and accept a valid but weak test.
2. **H2:** replace it as in 2.3 (proposed), or keep the v1.0 rule and report
   the new one as supplementary.
3. **The safety margin** in the new H2: 20% of baseline recall (proposed),
   10%, or 25%.
4. **Known defects in frozen inputs** (the Born test on rounded tensors: 2
   structures admitted and 19 rejected on round-off; the wrong silicon row in
   the reference table): leave the pool and reference table exactly as they
   are and disclose (proposed), or rebuild them.

## 4. What does not change

The raw snapshot, the pool and its three universes, splits r0–r2, descriptors,
both priors and their matched variants, the objectives, the 10 + 40 budget,
rank mixing, S_top / S_top1 / S_dis, T99 as the efficiency outcome, the
registered β arms, H1's content, H3, the retired seed 0, and the Phase A /
Phase B framing.

## 5. What an amendment cannot fix

Listed in Section 3 of the review: the exploratory phase saw outcome labels;
the β = 1 arm is computable from data already inspected; the engine was
chosen on materials that are 65–75% of every discovery pool; the universe is
filtered by the prior's own data quality; the splits re-partition one pool.
v1.1 states these plainly. Proposed wording: *"The benchmark, metrics and
hypotheses were developed with access to outcome labels of all pools. The
analysis plan for the GP-driven search arms was registered, and amended once,
before any such search was run on a discovery pool."*

## 6. Work after your decisions (about a week)

1. Engine v1.1; engine tests extended (reference hashes for every β, asserted).
2. `make_splits.py` v2 and priors for the new splits; rebuild check on two
   Python versions and on Colab.
3. `analyze_runs.py` v1.1 (split-level test, new H2, new sets, fixes from the
   review) with tests of the full H1/H2/H3 wiring; pilot driver.
4. `PREREGISTRATION.md` v1.1 with a table of every change against v1.0;
   decision-log entries; README, licence, Materials Project attribution.
5. You: run the rebuild and engine check on Colab; publish v1.1 as a new
   version of the same Zenodo record; commit.
6. Pilot, then the confirmatory runs.

## 7. History-side runs made for this proposal (to be logged)

| Set | Runs | What was measured |
|---|---|---|
| Engine health, β = 0, bounds 0.01 / 1 / 2 / 3 | 240 | flat-EI steps, fitted length-scale, top-5% hits |
| Forced near-duplicate designs, β = 0 and 0.05 | 96, plus 48 reference | the same; prior percentile of picks |
| Selection geometry, bound 3, β = 0.02–0.5, elemental prior, Y1 | 90 | EI rank and prior rank of each pick (no outcome summarised) |
| Statistical simulation | synthetic | error rate and power |

Scripts and raw logs: `amend_evidence/` (`hexp.py`, `hexp2.py`, `statsim.py`,
`hexp.jsonl`, `hexp2.jsonl`).

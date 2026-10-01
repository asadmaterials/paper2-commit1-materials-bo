# Paper 2 — Decision log

Every design decision, in order, with the data that had been **seen** when it
was made. Its purpose is to let a reader judge which choices could have been
influenced by the data. Entries are append-only; corrections are new entries.

"Outcome data" = MP elastic values (G, K, derived E/ρ) of materials, which are
what the search is scored on.

[Author: add a statement on AI assistance per the target journal's policy.
Design discussions and code were developed with an AI assistant and reviewed
by separate AI-assisted reviews.]

| # | Date | Decision | Data seen at the time | Outcome data seen? |
|---|---|---|---|---|
| 1 | 2026-09-24 | Candidate scope: metallic (MP `is_metal`) + metallic-element list; N ≥ 2; borides/carbides/metalloid-anion compounds excluded; all polymorphs kept; composition-grouped splits | none (before the pull) | no |
| 2 | 2026-09-24 | Fixed priors primary for the objective-change comparison; matched priors added | none | no |
| 3 | 2026-09-24 | `young_modulus` removed from requested fields | `--inspect` output (3 documents) | no |
| 4 | 2026-09-24 | Raw snapshot acquired, MP 2026.04.13; deposited on Zenodo, unfiltered | — | — |
| 5 | ≤ 2026-09-26 | Born stability = all eigenvalues of the symmetrised IEEE Voigt tensor > 0 (MP warning used as audit only) | provisional funnel counts; 161 disagreements between the eigenvalue check and the MP warning | aggregate only (G range, median, p95) |
| 6 | ≤ 2026-09-26 | Elemental reference source: Springer Handbook 2005 (Paper 1 `ELEM_G` table rejected: no provenance) | per-element structure counts in the provisional pool | aggregate only |
| 7 | ≤ 2026-09-26 | Consistency checks and thresholds (ν within 0.10; ρv_t² within 15%) **proposed in the same step in which the check values for every element were computed** | all handbook values; per-element structure counts | no per-material outcomes |
| 8 | 2026-09-26 | Rule C adopted (keep G if ≥ 1 computable check passes); no repairs; Pm/Tc kept, with C-without-Pm/Tc as a sensitivity analysis; Ca kept (348 K) | pool size under options A–D | no |
| 9 | 2026-09-26 | Removed the undeclared 200 m/s v_t cutoff; stated the ν_implied ∈ [0, 0.5] clause (changes no element) | as above | no |
| 10 | 2026-09-26 | **Implementation corrected to the written sound rule** (v0 code used G/ρv_t² ∈ [0.85, 1.15]). Na moved from excluded to eligible; pool 3,473 → 3,550 | provisional counts under v0 code | no |
| 11 | 2026-09-26 | One shared universe for all priors; P0/P_H on the unscreened pool as sensitivity analysis | pool counts | no |
| 12 | 2026-09-26 | Structural-family hold-out dropped | family sizes (Heusler ABC2 #225 = 33% of pool); **top-5% family composition by G and by G/ρ on the whole primary pool** | **yes** |
| 13 | 2026-09-26 | History fraction 0.30; 3 replicates per split type; splits defined on the unscreened universe | per-split size and median G by side | yes (medians) |
| 14 | 2026-09-26 | Descriptors v1 (Magpie + density, volume/atom, space group, nsites, crystal system); P_H = nearest-neighbour distance to top 10% of history | none beyond pool | no |
| 15 | 2026-09-26 | Alignment table v1 computed on the discovery side; S_dis (pool-median split) found empty for P_phys; own-median S_dis **proposed** | discovery-side y for all splits | **yes** |
| 16 | 2026-09-27 | Own-median S_dis shown to give zero recall at β = 1 by construction; H1/H2 as prior-only tests **retired**; hypotheses reframed to the frontier interior (draft v0.2) | prior-only recall of S_top on discovery y (2 splits) | **yes** |
| 17 | 2026-09-27 | Y2 changed from G/ρ to E/ρ (convention) | Spearman(E/ρ, G/ρ) = 0.998; top-5% overlap 173/178 on the primary pool | **yes** |
| 18 | 2026-09-27 | `nsites` and `spacegroup_number` removed from descriptors | none (methodological) | no |
| 19 | 2026-09-27 | Alignment reported on the history side (grouped 5-fold CV for P_H) as the pre-search quantity; discovery side kept as oracle diagnostic | history vs discovery ρ for P_phys (2 splits) | yes |
| 20 | 2026-09-27 | Primary efficiency outcome T95 → T99 (T95 degenerate: 40% of random initial designs already contain a top-5% material) | pool size only (combinatorial calculation) | no |
| 21 | 2026-09-27 | H3 sign-reversal prediction retired (no anti-aligned prior in the alignment table) | alignment table | yes |
| 22 | 2026-09-27 | `build_priors.py` v2.1: prior scores computed with `math.fsum` in sorted element order and direct-difference distances (no BLAS), rounded to 12 significant digits. Cause: cross-environment rebuild (Colab, Python 3.13) mismatched `priors_v2_primary.csv`. Python ≥ 3.12 changed float `sum()`, and BLAS results vary by machine; 24–36 pairs of different compositions per split have mathematically equal P_phys (Y/Dy share G = 25.5 GPa, Er/Ge share 29.6 GPa), and v2.0 ordered some of them by round-off instead of the hash tie-break (e.g. TiGe2 vs ErTiGe; 4 rank changes). Alignment table unchanged. v2.1 verified identical on Python 3.11 and 3.13 and across BLAS thread counts | rebuild hashes; prior values of the two environments | no (ranks of prior scores only) |

## Threshold sensitivity (descriptive; the frozen rule is unchanged)

| ν tolerance | ρv_t² tolerance | Primary universe | Pool elements excluded |
|---|---|---|---|
| 0.05 | 10% | 3,286 | Ac As Ce Cr Ga Hg **Na** Nb Np Pa Pu **Ta Te Tm** |
| **0.10** | **15%** (frozen) | **3,550** | Ac As Ce Cr Ga Hg Nb Np Pa Pu Te |
| 0.15 | 20% | 3,609 | Ac As Cr Ga Hg Nb Np Pa Te |

## Known limitations recorded at decision time

- Extraction verification: density and G of the 71 pool elements were compared
  visually against the rendered handbook pages; **E, ν and v_t were visually
  verified only for flagged elements**. A full independent check of these
  three columns is pending (`reference_verification_sheet.csv`).
- Artifacts to date were generated in the assistant's environment (Python 3.11.15,
  numpy 2.4.4, scipy 1.17.1, pandas 2.3.3, matminer 0.10.1, pymatgen 2026.9.24,
  poppler pdftotext 24.02.0). Only `build_priors.py` uses pandas/matminer.
  Colab (Python 3.13.15) reproduced 6/7 artifacts; the seventh led to entry 22.
  v2.1 is identical on Python 3.11 and 3.13 in the assistant's environment, and
  the Colab rerun (Python 3.13.15, 2026-09-27) matched all 7 artifacts.

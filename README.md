# Paper 2 — the efficiency–coverage trade-off of scientific priors in Bayesian search

A registered benchmark on metallic multi-element Materials Project (MP)
elasticity data. This repository holds the code, the derived data and the
protocol documents. **Apart from one accidental run that was deleted unread,
no model-driven search was run on a discovery pool before registration v1.1.**

| | |
|---|---|
| Registered plan, current | `docs/PREREGISTRATION.md` **v1.1**, doi:10.5281/zenodo.23118101 |
| Registered plan, first version | v1.0, doi:10.5281/zenodo.23073511; byte-identical copies of its plan and decision log are in `docs/registered_v1.0/`; git commit `21bc34dee85938d143c52f283f1b31f0474b8531` |
| Raw data | MP database 2026.04.13 snapshot, doi:10.5281/zenodo.22978221 |
| What changed between the versions, and why | `docs/PREREGISTRATION.md`, Appendix C; `docs/DECISION_LOG.md`, entries 31–36; `docs/INDEPENDENT_REVIEW_2026-10-03.md` |

## Layout

```
code/           pipeline, search engine, analysis, drivers, tests
data/derived/   elemental reference, pool, splits (v1, v2), descriptors,
                priors (v2, v3), alignment tables, manifests
docs/           registered plan, decision log, review, evidence for the amendment
```

All scripts are run from `code/`.

| Script | Role |
|---|---|
| `pull_mp.py` | acquired the raw snapshot (needs an MP API key; see the data record for its environment) |
| `extract_reference.py`, `build_pool.py`, `make_splits.py`, `build_priors.py` | v1.0 inputs: reference table, pool, 6 splits, descriptors, priors |
| `make_splits_v2.py`, `build_priors_v3.py` | v1.1 inputs: 20 splits and their priors (the first 6 splits are unchanged) |
| `run_search.py` | search engine v1.1 |
| `run_batch.py`, `pilot_check.py` | launch the pilot and confirmatory runs; the three pilot checks |
| `analyze_runs.py` | the registered analysis v1.1 |
| `test_engine.py`, `test_analysis.py`, `test_pilot_pipeline.py` | tests (history pools and synthetic data only) |
| `rebuild_and_verify.sh`, `rebuild_v1_1.sh` | rebuild the derived data and compare hashes |
| `engine_design_probe.py` | history-side comparison of GP forms (decision 23) |

## Verify

```
# hashes of every frozen file (Appendix A of the plan); from the repository root
sha256sum -c SHA256SUMS.txt

# rebuild the v1.1 inputs from the v1.0 inputs (about 5 minutes); from code/
pip install -r requirements-build.txt
bash rebuild_v1_1.sh ../data/derived

# tests; from code/
D=../data/derived
python test_engine.py $D/features_v2_primary.csv $D/pool_v1.csv $D/splits_v1.csv \
       $D/elemental_reference_springer2005_v1.csv $D/priors_v2_primary.csv
pip install -r requirements-analysis.txt
python test_analysis.py
python test_pilot_pipeline.py
```

`rebuild_and_verify.sh` rebuilds the v1.0 inputs from the raw snapshot and the
handbook PDF. The handbook is not distributed here; a reader without the same
PDF file should treat `data/derived/elemental_reference_springer2005_v1.csv`
as the verifiable input. Set `PYTHONUTF8=1` where the default encoding is not
UTF-8. `SHA256SUMS.txt` lists the current files; Appendix A of the v1.0 plan
verifies at commit `21bc34d`, where `docs/DECISION_LOG.md` had 28 entries;
Appendix A of the v1.1 plan verifies at commit `4ba3afa`, where it had 36.
The confirmatory stage is started with `run_batch.py --stage confirmatory
--decision-log ../docs/DECISION_LOG.md --pilot-check ../runs/pilot/pilot_check.json`.

## Run (after registration v1.1 is published)

```
python run_batch.py --stage pilot --data-dir ../data/derived --out-dir ../runs/pilot --workers 2
python pilot_check.py --data-dir ../data/derived --runs-dir ../runs/pilot
```

## Notes on the frozen v1.0 scripts

Their file hashes are registered, so their comments were not edited. Read
them with these corrections: "pre-registered" in their docstrings means
"fixed before the search" (the exploratory phase is not presented as
preregistered); `make_splits.py` mentions "Y2 (G/rho)" but Y2 is E/ρ;
`build_priors.py` describes a retired definition of S_dis and its title says
v2 (it is v2.1); `docs/PREREG_HYPOTHESES_DRAFT.md` is superseded.

## Licences and attribution

- **Code:** MIT (see `LICENSE`).
- **Data derived from the Materials Project** (`data/derived/pool_v1.csv` and
  files built from it): Creative Commons Attribution 4.0, as for the MP data.
  Please cite A. Jain et al., *APL Materials* 1, 011002 (2013), and M. de Jong
  et al., *Scientific Data* 2, 150009 (2015).
- **Elemental reference values** were extracted from W. Martienssen and
  H. Warlimont (eds.), *Springer Handbook of Condensed Matter and Materials
  Data* (Springer, 2005), mechanical-property tables 2.1-6B(b)–2.1-26B(b).
  The handbook itself is not redistributed; the extracted numerical values
  are included with this citation.
- Design, code and documents were developed with an AI assistant; see the
  plan, Section 11.

# Paper 2 — registration bundle v1

Frozen code, derived data and protocol documents for the registered analysis
plan (PREREGISTRATION.md, https://doi.org/10.5281/zenodo.23073511).

- `code/`  pipeline, search engine, tests, rebuild script, expected hashes
- `data/derived/`  elemental reference, pool, splits, descriptors, priors and
  alignment tables (primary and both sensitivity universes), with manifests
- `docs/`  decision log, excluded-runs record, reference verification sheet,
  engine cross-machine check

Every file listed in Appendix A of PREREGISTRATION.md is included and matches
its hash (`sha256sum -c` from this directory). `SHA256SUMS.txt` covers all files.

Not included: the raw MP 2026.04.13 snapshot (separate Zenodo record) and the
Springer Handbook (cited, not redistributed). `rebuild_and_verify.sh`
regenerates every derived file from those two inputs.

No discovery-side optimization run had been inspected when this bundle was
frozen (see docs/DECISION_LOG.md, entry 25, and docs/EXCLUDED_RUNS.json).

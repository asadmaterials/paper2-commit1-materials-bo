#!/usr/bin/env bash
# rebuild_v1_1.sh — rebuild the inputs added in registration v1.1 (splits_v2,
# priors_v3 and history-side alignment for the three universes) from the
# registered v1.0 inputs, and compare SHA-256 hashes with expected_hashes_v1_1.txt.
#
# Usage (from the code/ directory):
#   bash rebuild_v1_1.sh ../data/derived
#
# DATA_DIR must hold pool_v1.csv, splits_v1.csv, elemental_reference_springer2005_v1.csv,
# features_v2_<u>.csv and priors_v2_<u>.csv (all in the repository; they are
# themselves rebuilt from the raw snapshot by rebuild_and_verify.sh).
# No outcome of a discovery-side material enters any statistic in this build.
# Needs numpy, scipy and pymatgen (requirements-build.txt). Set PYTHONUTF8=1 on
# systems whose default encoding is not UTF-8.
set -euo pipefail

DATA=${1:?data dir}
OUT=rebuild_v1_1_out
rm -rf "$OUT"; mkdir -p "$OUT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUTF8=1

echo "== environment"
python3 --version
python3 - <<'PY'
import importlib.metadata as m
for p in ["numpy", "scipy", "pymatgen"]:
    print(p, m.version(p))
PY

echo "== 1/2 splits (r0-r2 must equal splits_v1.csv)"
python3 make_splits_v2.py --pool "$DATA/pool_v1.csv" --check-v1 "$DATA/splits_v1.csv" --out-dir "$OUT" > "$OUT/make_splits_v2.log"
echo "== 2/2 priors and history-side alignment (scores on r0-r2 must equal priors_v2)"
for U in primary sens1_no_Pm_Tc unscreened; do
    python3 build_priors_v3.py --pool "$DATA/pool_v1.csv" --splits "$OUT/splits_v2.csv" \
        --reference "$DATA/elemental_reference_springer2005_v1.csv" --features "$DATA/features_v2_$U.csv" \
        --check-v2 "$DATA/priors_v2_$U.csv" --universe "$U" --out-dir "$OUT" > "$OUT/build_priors_v3_$U.log"
done

echo "== compare"
fail=0
while read -r want name; do
    [[ -z "$want" || "$want" == \#* ]] && continue
    got=$(sha256sum "$OUT/$name" | cut -d' ' -f1)
    if [[ "$got" == "$want" ]]; then echo "MATCH     $name"; else echo "MISMATCH  $name  got $got"; fail=1; fi
done < expected_hashes_v1_1.txt
if [[ $fail -eq 0 ]]; then echo "ALL MATCH"; else echo "SOME MISMATCH: see above"; exit 1; fi

#!/usr/bin/env bash
# rebuild_and_verify.sh — rebuild every Paper 2 artifact from its inputs and
# compare SHA-256 hashes against expected_hashes.txt.
#
# Usage (from the directory holding the scripts):
#   bash rebuild_and_verify.sh SNAPSHOT_DIR HANDBOOK_PDF
#
# SNAPSHOT_DIR must contain elasticity.jsonl.gz, summary.jsonl.gz and
# manifest.json from the Zenodo deposit. HANDBOOK_PDF is the Springer Handbook
# file whose SHA-256 is recorded below.
#
# Colab:
#   !apt-get -qq install poppler-utils
#   !pip install -q -r requirements-build.txt
#   !bash rebuild_and_verify.sh /content/snapshot /content/handbook.pdf
set -euo pipefail

SNAP=${1:?snapshot dir}
PDF=${2:?handbook pdf}
OUT=rebuild_out
rm -rf "$OUT"; mkdir -p "$OUT"

echo "== environment"
python3 --version
pdftotext -v 2>&1 | head -1
python3 - <<'EOF'
import importlib.metadata as m
for p in ["numpy","scipy","pandas","matminer","pymatgen","pymatgen-core"]:
    print(f"{p} {m.version(p)}")
EOF

echo "== input hashes"
sha256sum "$PDF" "$SNAP/elasticity.jsonl.gz" "$SNAP/summary.jsonl.gz"

echo "== 1/4 elemental reference"
python3 extract_reference.py "$PDF" --out-dir "$OUT"
echo "== 2/4 candidate pool"
python3 build_pool.py --snapshot "$SNAP" --reference "$OUT/elemental_reference_springer2005_v1.csv" --out-dir "$OUT" > "$OUT/build_pool.log"
echo "== 3/4 splits"
python3 make_splits.py --pool "$OUT/pool_v1.csv" --out-dir "$OUT" > "$OUT/make_splits.log"
echo "== 4/4 priors and alignment (primary and both sensitivity universes)"
for U in primary sens1_no_Pm_Tc unscreened; do
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 build_priors.py --pool "$OUT/pool_v1.csv" --splits "$OUT/splits_v1.csv" \
        --reference "$OUT/elemental_reference_springer2005_v1.csv" --universe "$U" --out-dir "$OUT" > "$OUT/build_priors_$U.log" 2>/dev/null
done

echo "== compare"
fail=0
while read -r want name; do
    [[ -z "$want" || "$want" == \#* ]] && continue
    f="$OUT/$name"; [[ "$name" == handbook.pdf ]] && f="$PDF"
    got=$(sha256sum "$f" | cut -d' ' -f1)
    if [[ "$got" == "$want" ]]; then echo "MATCH     $name"; else echo "MISMATCH  $name  got $got"; fail=1; fi
done < expected_hashes.txt
if [[ $fail -eq 0 ]]; then echo "ALL MATCH"; else echo "SOME MISMATCH: see above"; exit 1; fi

# Cross-machine check of the search engine (Colab)

Run in the same Colab folder used for `rebuild_and_verify.sh` (packages already
installed; add scikit-learn first if the session is new).

```
!pip install -q scikit-learn==1.8.0
!OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python test_engine.py \
    rebuild_out/features_v2_primary.csv rebuild_out/pool_v1.csv rebuild_out/splits_v1.csv \
    rebuild_out/elemental_reference_springer2005_v1.csv rebuild_out/priors_v2_primary.csv
```

Expected: `13/13 passed` and exactly these two lines:

```
RUNHASH H_random_r0_Y1_Pphys_b0.0_s1 79a48ce1eb7b349f4e823f628e57edc0523f6441be09165ab8db25a7c6a062ed
RUNHASH H_random_r0_Y1_Pphys_b0.5_s1 28330f33123d57aa9d4ab450f2295e80089448e5479bd3e61fb3ee4884291107
```

Reference environment (where the expected hashes were produced): the
assistant's Linux VM, numpy 2.4.4, scipy 1.17.1, scikit-learn 1.8.0; identical
hashes on Python 3.11.15 and 3.13.7, with 1 and 2 BLAS threads.

These are full 40-step GP+EI searches on a HISTORY pool. Identical hashes on
Colab hardware mean search paths reproduce across machines.
Files needed next to the test: run_search.py, test_engine.py, build_priors.py.

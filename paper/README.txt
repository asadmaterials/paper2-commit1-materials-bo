Paper 2 -- Overleaf project (Digital Discovery / RSC article template)

Main document: main.tex (two-column RSC layout). Overleaf compiles it by default.
Supplementary Information: supplement.tex. The latexmkrc file builds it automatically before main.tex, so that the
cross-references between the two documents resolve. To view the Supplement's PDF, set supplement.tex as the main
document (Menu > Main document) and recompile.

Do not edit by hand (generated from the committed results by the repository scripts):
  numbers.tex                      code/make_numbers.py
  supp_tables_registered.tex,
  supp_tables_posthoc.tex          code/make_supplement_tables.py
  figures/*.pdf                    code/make_figures.py
Text files: abstract.tex, intro.tex, methods.tex, results.tex, discussion.tex, em_*.tex (end matter).
Template files: head_foot/, rsc.bst. Bibliography: references.bib.

This folder holds the manuscript sources as submitted to Digital Discovery (repository tag dd-submission-v1, decision-log
entry 48). The ESI DOI in the main.tex footnote is assigned by the journal.
claims.tex is the claims ledger (what each result sentence may claim, with its source file); it is not part of the submission.

Regenerating the generated files (from code/):
  python make_numbers.py ../results ../results_posthoc ../paper/numbers.tex
  python make_supplement_tables.py ../results ../results_posthoc ../paper/supp_tables.tex
  python make_figures.py ../results ../paper/figures

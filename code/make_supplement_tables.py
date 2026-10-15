"""
make_supplement_tables.py — the Supplement's tables, generated from the result tables
=====================================================================================
Writes the Supplement's tables: one LaTeX table per registered or post hoc output that the main text does not
show in full. No number is typed by hand; every cell is formatted from results/ or results_posthoc/.

Writes <stem>_registered.tex and <stem>_posthoc.tex next to the given path.
Usage (from code/):  python make_supplement_tables.py ../results ../results_posthoc ../paper/supp_tables.tex
"""
import json
import math
from decimal import Decimal, ROUND_HALF_UP
import sys
from pathlib import Path

import pandas as pd

RES, POST, OUT = (Path(a) for a in sys.argv[1:4])
CELLS = ["Y1:P_H_fixed", "Y1:P_phys_fixed", "Y2:P_H_fixed", "Y2:P_H_matched_Y2", "Y2:P_phys_fixed", "Y2:P_phys_matched_Y2"]
NAME = {"Y1:P_H_fixed": "Y1 hist.", "Y1:P_phys_fixed": "Y1 elem.", "Y2:P_H_fixed": "Y2 hist.",
        "Y2:P_H_matched_Y2": "Y2 hist.\\ m.", "Y2:P_phys_fixed": "Y2 elem.", "Y2:P_phys_matched_Y2": "Y2 elem.\\ m."}
UNI = {"primary": "primary", "sens1_no_Pm_Tc": "without Pm, Tc", "unscreened": "unscreened"}


def f(v, d=3):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "--"
    q = Decimal(repr(float(v))).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP)   # round half up
    s = f"{q:.{d}f}".lstrip("-")
    return r"$-$" + s if float(v) < 0 else s                  # a negative that rounds to zero keeps its sign


def p(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "--"
    return r"$<$0.001" if v < 0.001 else f"{v:.2g}"


def ci(lo, hi, d=3):
    return r"{\scriptsize[" + f(lo, d) + ", " + f(hi, d) + "]}"


def b(v):
    return f"{float(v):g}"


def yes(v):
    return "yes" if str(v) == "True" else "no"


def table(label, caption, colspec, header, rows, long=False, size=r"\footnotesize"):
    body = "\n".join(" & ".join(r) + r" \\" for r in rows)
    if long:
        return (f"{{{size}\\setlength{{\\tabcolsep}}{{3pt}}\n\\begin{{longtable}}{{{colspec}}}\n"
                f"\\caption{{{caption}}}\\label{{{label}}}\\\\\n\\toprule\n{header} \\\\\n\\midrule\n\\endfirsthead\n"
                f"\\toprule\n{header} \\\\\n\\midrule\n\\endhead\n\\bottomrule\n\\endfoot\n{body}\n\\end{{longtable}}}}\n")
    return (f"\\begin{{table}}[htbp]\n\\centering{size}\n\\caption{{{caption}}}\\label{{{label}}}\n"
            f"\\setlength{{\\tabcolsep}}{{3pt}}\n\\begin{{tabular}}{{{colspec}}}\n\\toprule\n{header} \\\\\n\\midrule\n"
            f"{body}\n\\bottomrule\n\\end{{tabular}}\n\\end{{table}}\n")


parts = []
fr = pd.read_csv(RES / "frontier.csv")
h1 = pd.read_csv(RES / "h1_coverage_loss.csv")
h2 = pd.read_csv(RES / "h2_detail.csv")
h2b = pd.read_csv(RES / "h2_beneficial_and_safe.csv")
ns = pd.read_csv(RES / "new_splits_only.csv")
bud = pd.read_csv(RES / "budget_sensitivity.csv")
eh = pd.read_csv(RES / "engine_health.csv")
rr = pd.read_csv(RES / "random_search_reference.csv")
fam = pd.read_csv(RES / "family.csv")
gl = pd.read_csv(RES / "glmm_fixed_effects.csv")
summ = json.loads((RES / "summary.json").read_text())
cen = pd.read_csv(POST / "t99_censoring.csv")
it1 = pd.read_csv(POST / "init_top1.csv")
isv = pd.read_csv(POST / "init_design_sensitivity.csv")
bl = pd.read_csv(POST / "blind_spot_loss.csv")
ha = pd.read_csv(POST / "h2_adjusted_beneficial.csv")
st = pd.read_csv(POST / "by_split_type.csv")
P = lambda df: df[df.universe == "primary"]

# ---------------------------------------------------------------- S1 budget (registered, descriptive)
rows = []
for c in CELLS:
    for beta in (0.25, 0.5):
        r = [NAME[c], b(beta)]
        for t in (10, 20, 30):
            x = bud[(bud.universe == "primary") & ((bud.cell == c) & (bud.beta == beta) & (bud.acquisitions == t))].iloc[0]
            r += [f(x.gain, 1) + " " + ci(x.gain_lo, x.gain_hi, 1), f(x.L) + " " + ci(x.L_lo, x.L_hi)]
        rows.append(r)
parts.append(table("tab:s-budget", r"Budget sensitivity (registered, descriptive; primary universe). Gain (reduction in restricted mean "
    r"sequential $T_{99}$, restricted at $t+1$) and coverage loss $L$ after the first $t$ acquisitions, with two-sided 95\% "
    r"split-level $t$ intervals.", "@{}llcccccc@{}",
    r"Cell & $\beta$ & gain, $t{=}10$ & $L$, $t{=}10$ & gain, $t{=}20$ & $L$, $t{=}20$ & gain, $t{=}30$ & $L$, $t{=}30$", rows))

# ---------------------------------------------------------------- S2 censoring (post hoc)
betas = sorted(cen.beta.unique())
rows = [[NAME[c]] + [f(100 * cen[(cen.universe == "primary") & ((cen.cell == c) & (cen.beta == x))].iloc[0].no_hit_share, 1) for x in betas] for c in CELLS]
parts.append(table("tab:s-censor", r"Share of runs (\%) without a top-1\% material among the 40 acquisitions ($T_{99}=41$), by prior "
    r"weight (post hoc, decision-log entry 42; primary universe; 60 runs per arm; the $\beta=0$ runs are shared by the cells of an "
    r"objective). The sensitivity universes are in \texttt{results\_posthoc/t99\_censoring.csv}.", "@{}l" + "c" * len(betas) + "@{}", "Cell & " + " & ".join(f"$\\beta={b(x)}$" for x in betas), rows))

# ---------------------------------------------------------------- S3 initial design (post hoc)
rows = [[UNI[r.universe], r.objective, str(int(r.blocks_init_has_top1)), str(int(r.blocks)), str(int(r.random_type)),
         str(int(r.chemsys_type))] for r in it1.itertuples()]
pv = P(isv)
flips = pv[(pv.beta < 1) & (pv.beneficial_registered != pv.beneficial_init_as_time0)]
note = (r" Random and chemical-system: blocks with a hit, by split type. In the primary universe, scoring those blocks as "
        r"discovered at time zero (every arm of the block given $T_{99}=0$) changes the gains at $\beta\le0.75$ by at most "
        + f((pv.gain_init_as_time0 - pv.gain_registered)[pv.beta < 1].abs().max(), 2) + " acquisitions and changes "
        + str(len(flips)) + " beneficial decision" + ("" if len(flips) == 1 else "s") + " from not beneficial to beneficial"
        * bool(len(flips)) * bool((~flips.beneficial_registered.astype(bool)).all())
        + ("" if not len(flips) else " (" + "; ".join(f"{NAME[x.cell]}, $\\beta={b(x.beta)}$" for x in flips.itertuples()) + ")")
        + "; " + ("no beneficial and safe set changes, because "
                  + ("no changed weight is" if len(flips) != 1 else "the changed weight is not") + " safe."
                  if not flips.merge(h2, on=["universe", "cell", "beta"]).safe.astype(str).eq("True").any()
                  else r"\textbf{a beneficial and safe set changes}."))
parts.append(table("tab:s-init", r"Split--seed blocks whose shared 10-point initial design already contains a top-1\% material (post "
    r"hoc, decision-log entry 42)." + note, "@{}llcccc@{}", r"Universe & Objective & blocks with a hit & blocks & random & chemical-system",
    rows))

# ---------------------------------------------------------------- S4 blind spot (post hoc)
rows = []
for c in CELLS:
    g = bl[bl.cell == c].set_index("beta")
    rows.append([NAME[c], f(g.loc[0.0].R_blind)] + [f(g.loc[x].L_blind) + " " + ci(g.loc[x].L_blind_lo, g.loc[x].L_blind_hi)
                                                     for x in (0.1, 0.25, 0.5, 0.75)])
parts.append(table("tab:s-blind", r"Loss of recall of the absolute blind spot $S_{\mathrm{blind}}$ (top-5\% materials below the prior's "
    r"80th percentile), $L_{\mathrm{blind}}=R_{\mathrm{blind}}(0)-R_{\mathrm{blind}}(\beta)$, per-run recall with the registered "
    r"split-level estimator and two-sided 95\% $t$ intervals (post hoc, decision-log entry 43; primary universe). "
    r"$R_{\mathrm{blind}}(0)$ is the mean of per-run recall and differs slightly from the pooled recall in Supplementary "
    r"Table~\ref{tab:s-frontier}, which weights runs by the size of their set. All weights are in "
    r"\texttt{results\_posthoc/blind\_spot\_loss.csv}.", "@{}lccccc@{}",
    r"Cell & $R_{\mathrm{blind}}(0)$ & $\beta=0.1$ & $\beta=0.25$ & $\beta=0.5$ & $\beta=0.75$", rows))

# ---------------------------------------------------------------- S5 H1 in full (registered)
rows = []
for c in CELLS:
    for beta in (0.25, 0.5):
        x = h1[(h1.universe == "primary") & ((h1.cell == c) & (h1.beta == beta))].iloc[0]
        rows.append([NAME[c], b(beta), f(x.mean_L), ci(x.ci_lo, x.ci_hi), ci(x.boot_ci_lo, x.boot_ci_hi), f(x.df, 1), p(x.p),
                     p(x.p_holm), yes(x.reject_H0), p(x.boot_p_v1_0)])
parts.append(table("tab:s-h1", r"H1 in full (registered, confirmatory; primary universe): mean coverage loss $L$, two-sided 95\% split-level $t$ "
    r"interval, hierarchical-bootstrap 95\% interval (descriptive), Satterthwaite degrees of freedom, one-sided $p$, Holm-adjusted $p$ "
    r"over the 12 tests and decision, and the version-1.0 bootstrap $p$ (descriptive).", "@{}llcccccccc@{}",
    r"Cell & $\beta$ & $L$ & $t$ interval & bootstrap & df & $p$ & $p_{\mathrm{Holm}}$ & reject & boot.\ $p$", rows))

# ---------------------------------------------------------------- S6 H2 in full (registered + post hoc column)
rows = []
for c in CELLS:
    for x in h2[(h2.universe == "primary") & (h2.cell == c)].sort_values("beta").itertuples():
        a = ha[(ha.universe == "primary") & ((ha.cell == c) & (ha.beta == x.beta))].iloc[0]
        rows.append([NAME[c], b(x.beta), f(x.gain_mean, 1), f(x.gain_lower95, 2), f(a.gain_lower_adjusted, 2), f(x.relative_loss),
                     f(x.margin_upper_adjusted), x.state_unadjusted, yes(x.beneficial), yes(x.safe), yes(x.safe_at_10pct),
                     yes(x.safe_at_25pct)])
parts.append(table("tab:s-h2", r"H2 in full (registered confirmatory decision rule; primary universe). Gain and its one-sided lower 95\% bound (registered) and "
    r"lower bound at level $1-0.05/6$ (post hoc, entry 44); relative loss $L/R_{\mathrm{dis}}(0)$; upper bound at level $1-0.05/6$ of "
    r"$L-0.2R_{\mathrm{dis}}(0)$; the unadjusted three-way state; the registered beneficial and safe calls; safety at 10\% and 25\% "
    r"margins (descriptive).", "@{}llcccccccccc@{}",
    r"Cell & $\beta$ & gain & lower 95 & lower adj. & rel.\ loss & margin up. & state & B & S & S 10\% & S 25\%", rows, long=True))

# ---------------------------------------------------------------- S7 frontier (registered)
rows = []
for c in CELLS:
    for x in fr[(fr.universe == "primary") & (fr.cell == c)].sort_values("beta").itertuples():
        rows.append([NAME[c], b(x.beta), f(x.gain, 1), ci(x.gain_lo, x.gain_hi, 1), ci(x.gain_boot_lo, x.gain_boot_hi, 1),
                     f(x.R_dis), f(x.L) + " " + ci(x.L_lo, x.L_hi), ci(x.L_boot_lo, x.L_boot_hi), f(x.N5, 1), f(x.AUC),
                     f(x.R_blind_pooled)])
parts.append(table("tab:s-frontier", r"Frontier in full (registered primary descriptive result, primary universe): gain with $t$ and "
    r"bootstrap 95\% intervals, recall of $S_{\mathrm{dis}}$, coverage loss $L$ with $t$ and bootstrap intervals, $N_5$ (top-5\% "
    r"materials acquired), area under the best-so-far quantile curve, and pooled recall of $S_{\mathrm{blind}}$.",
    "@{}llccccccccc@{}", r"Cell & $\beta$ & gain & $t$ int. & boot.\ int. & $R_{\mathrm{dis}}$ & $L$ & boot.\ int. & $N_5$ & AUC & "
    r"$R_{\mathrm{blind}}$", rows, long=True))

# ---------------------------------------------------------------- S8 r3--r9 at all weights (registered, descriptive)
rows = []
for c in CELLS:
    for x in ns[(ns.universe == "primary") & (ns.cell == c)].sort_values("beta").itertuples():
        rows.append([NAME[c], b(x.beta), f(x.gain, 1) + " " + ci(x.gain_lo, x.gain_hi, 1), f(x.L) + " " + ci(x.L_lo, x.L_hi),
                     p(x.p_L_greater_0), f(x.df, 1)])
parts.append(table("tab:s-newsplits", r"Previously unanalysed splits r3--r9 only (registered, descriptive; 14 splits): gain and coverage "
    r"loss with two-sided 95\% $t$ intervals, and the unadjusted one-sided $p$ of $L>0$ with its degrees of freedom.",
    "@{}llcccc@{}", r"Cell & $\beta$ & gain & $L$ & $p(L>0)$ & df", rows, long=True))

# ---------------------------------------------------------------- S9 by split type (post hoc)
rows = []
for c in CELLS:
    for beta in (0.25, 0.5):
        r = [NAME[c], b(beta)]
        for q in ("gain", "L"):
            for t in ("random", "chemsys"):
                x = st[(st.cell == c) & (st.beta == beta) & (st.quantity == q) & (st.split_type == t)].iloc[0]
                d = 1 if q == "gain" else 3
                r.append(f(x.estimate, d) + " " + ci(x.ci_lo, x.ci_hi, d))
        rows.append(r)
parts.append(table("tab:s-splittype", r"Gain and coverage loss within each split type (post hoc, decision-log entry 45; primary "
    r"universe; mean of 10 split means with a two-sided 95\% $t$ interval on 9 degrees of freedom). All weights are in "
    r"\texttt{results\_posthoc/by\_split\_type.csv}.", "@{}llcccc@{}",
    r"Cell & $\beta$ & gain, random & gain, chem.-system & $L$, random & $L$, chem.-system", rows))

# ---------------------------------------------------------------- S10 sensitivity universes (registered, descriptive)
rows = []
for u in ("sens1_no_Pm_Tc", "unscreened"):
    for c in CELLS:
        g = h1[(h1.universe == u) & (h1.cell == c)]
        if g.empty:
            continue
        s = h2b[(h2b.universe == u) & (h2b.cell == c)].iloc[0]
        sel = [b(x) for x in h2[(h2.universe == u) & (h2.cell == c) & (h2.beneficial_and_safe.astype(str) == "True")].beta]
        a = ha[(ha.universe == u) & (ha.cell == c) & ha.bs20_adjusted]
        x25, x50 = g[g.beta == 0.25].iloc[0], g[g.beta == 0.5].iloc[0]
        rows.append([UNI[u], NAME[c], f(s.baseline_R_dis), f(x25.mean_L) + " " + ci(x25.ci_lo, x25.ci_hi),
                     f(x50.mean_L) + " " + ci(x50.ci_lo, x50.ci_hi), r"\{" + ", ".join(sel) + r"\}" if sel else "none",
                     r"\{" + ", ".join(b(v) for v in a.beta) + r"\}" if len(a) else "none"])
parts.append(table("tab:s-universes", r"Sensitivity universes (registered, descriptive): baseline recall of $S_{\mathrm{dis}}$, coverage "
    r"loss at $\beta=0.25$ and $0.5$ with two-sided 95\% $t$ intervals, and the beneficial and safe weights under the registered rule "
    r"and with the beneficial bound also adjusted (post hoc, entry 44). Sensitivity universe 2 (unscreened) has the historical "
    r"priors only.", "@{}llccccc@{}",
    r"Universe & Cell & $R_{\mathrm{dis}}(0)$ & $L$, $\beta=0.25$ & $L$, $\beta=0.5$ & B\,+\,S & B\,+\,S, adj.", rows, long=True, size=r"\scriptsize"))

# ---------------------------------------------------------------- S11 engine health and random-search reference (registered)
rows = []
for c in CELLS:
    g = eh[(eh.universe == "primary") & (eh.cell == c)].dropna().sort_values("beta")
    rows.append([NAME[c]] + [f(100 * g[g.beta == x].iloc[0].share_length_scale_at_lower_bound, 1) if (g.beta == x).any() else "--"
                             for x in betas] + [f(100 * g.share_flat_EI.max(), 1)])
parts.append(table("tab:s-health", r"Engine health (registered diagnostic, primary universe): share of GP steps (\%) with the "
    r"length-scale within 0.1\% of its lower bound, by prior weight ($\beta=1$ fits no GP), and the largest share of steps in any arm "
    r"whose top-EI tie group holds more than half of the candidates.", "@{}l" + "c" * (len(betas) + 1) + "@{}",
    "Cell & " + " & ".join(f"$\\beta={b(x)}$" for x in betas) + " & flat EI", rows))
rows = [[UNI[r.universe], r.objective, f(r.R_dis, 4), f(r.N5, 2), f(r.P_T99_censored)] for r in rr.itertuples()]
parts.append(table("tab:s-random", r"Analytic random-search reference (registered; no runs): expected recall of any target set with "
    r"the initial design removed, expected $N_5$, and probability that no top-1\% material is drawn in 40 acquisitions (mean over "
    r"the 20 splits).", "@{}llccc@{}", r"Universe & Objective & recall & $N_5$ & $P(T_{99}=41)$", rows))

# ---------------------------------------------------------------- S12 family level (registered, exploratory)
rows = []
for c in CELLS:
    for t in ("random", "chemsys"):
        for fam_name in ("Heusler:ABC2_225", "other"):
            x = fam[(fam.universe == "primary") & ((fam.cell == c) & (fam.split_type == t) & (fam.family == fam_name))].iloc[0]
            rows.append([NAME[c], "random" if t == "random" else "chem.-system",
                         "Heusler (ABC$_2$, 225)" if fam_name.startswith("Heusler") else "all other families",
                         str(int(x.distinct_S_top_materials)), f(x.rate_b0), f(x["rate_b0.5"]), f(x.ratio, 2),
                         ci(x.boot_lo, x.boot_hi, 2)])
parts.append(table("tab:s-family", r"Family level (registered, exploratory; primary universe): discovery rate of top-5\% materials at "
    r"$\beta=0$ and $0.5$ and their ratio, with a bootstrap 95\% interval, for Heusler compounds and all other families together. "
    r"Top-5\% materials: distinct members over the 10 splits of the type. Families with at least 5 distinct top-5\% materials "
    r"within a split type are in \texttt{results/family.csv}.",
    "@{}lllccccc@{}", r"Cell & Split type & Family & top-5\% materials & rate, $\beta=0$ & rate, $\beta=0.5$ & ratio & bootstrap", rows,
    long=True))

# ---------------------------------------------------------------- S13 mixed model (registered, secondary)
rows = []
for c in CELLS:
    g = gl[gl.cell == c]
    r = [NAME[c]]
    for x in ("0.02", "0.05", "0.1", "0.25", "0.5", "0.75", "1"):
        t = g[g.term == f"C(beta, Treatment('0'))[T.{x}]:pctl"].iloc[0]
        r.append(f(t.post_mean, 2))
    meta = summ["glmm"][c]
    r.append("yes" if "VB fitting did not converge" in meta["optimizer_warnings"] else "no")
    rows.append(r)
parts.append(table("tab:s-glmm", r"Candidate-level mixed model (registered, secondary; primary universe): posterior mean of the "
    r"$\beta\times$percentile interaction (logit scale, relative to $\beta=0$) and whether the variational fit reported "
    r"non-convergence. Shown for direction only; see Supplementary Note~\ref{note:glmm}.", "@{}lcccccccc@{}",
    r"Cell & 0.02 & 0.05 & 0.1 & 0.25 & 0.5 & 0.75 & 1 & not converged", rows))

# two files, registered outputs first, then post hoc ones; tables are numbered in this order
by_label = {t.split("\\label{")[1].split("}")[0]: t for t in parts}
ORDER = {"registered": ["tab:s-frontier", "tab:s-h1", "tab:s-h2", "tab:s-universes", "tab:s-newsplits", "tab:s-budget",
                        "tab:s-health", "tab:s-random", "tab:s-family", "tab:s-glmm"],
         "posthoc": ["tab:s-censor", "tab:s-init", "tab:s-blind", "tab:s-splittype"]}
assert sorted(sum(ORDER.values(), [])) == sorted(by_label)
for part, labels in ORDER.items():
    path = OUT.parent / f"{OUT.stem}_{part}.tex"
    path.write_text("% generated by code/make_supplement_tables.py -- do not edit by hand\n"
                    + "\n".join(by_label[x] for x in labels), encoding="utf-8")
    print(len(labels), "tables written to", path)

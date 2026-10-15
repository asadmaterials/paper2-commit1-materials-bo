"""
make_numbers.py — every number quoted in the paper, as LaTeX macros, read from the result tables
=================================================================================================
Writes paper/numbers.tex. The paper never types a result by hand: it writes \R{<key>}, and the key
names the table, the row and the column, e.g.

    \R{frontier/primary/Y1:P_H_fixed/0.1/gain}      mean ΔT99 at β = 0.1
    \R{h1/primary/Y2:P_phys_fixed/0.5/p_holm}       Holm-adjusted p of an H1 test
    \R{h2set/primary/Y1:P_phys_fixed}               beneficial-and-safe β set
    \R{posthoc/censor/primary/Y2:P_H_fixed/1.0}     share of runs with T99 = 41

A missing key prints a bold "??key" and a LaTeX warning, so a typo cannot pass silently.
Registered outputs come from results/ (analysis v1.1); post-hoc ones from results_posthoc/
(decision-log entry 42) and are kept under the separate prefix "posthoc/".

Usage (from code/):  python make_numbers.py ../results ../results_posthoc ../paper/numbers.tex
"""
import json
import math
from decimal import Decimal, ROUND_HALF_UP
import sys
from pathlib import Path

import pandas as pd

RES, POST, OUT = (Path(a) for a in sys.argv[1:4])

PCOLS = ("p", "p_holm", "p_L_greater_0", "boot_p_v1_0")
DEC = {"gain": 1, "N5": 1, "df": 1, "AUC": 3, "tau": 2, "pct": 1, "S_blind_mean_size": 1}       # decimals by column prefix; default 3


def fmt(col, v):
    if isinstance(v, (bool,)) or str(v) in ("True", "False"):
        return "yes" if str(v) == "True" else "no"
    if isinstance(v, str):
        return None
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "--"
    if col in PCOLS:                                                   # p values: 2 significant digits
        if v < 0.001:
            return r"\ensuremath{<}0.001"
        s = f"{v:.2g}"
        return s if "e" not in s else f"{v:.3f}"
    d = next((n for k, n in DEC.items() if col.startswith(k)), 3)
    if col in ("runs", "blocks", "blocks_init_has_top1", "random_type", "chemsys_type", "splits", "steps",
               "acquisitions", "available", "found", "distinct_S_top_materials", "valid_resamples"):
        d = 0
    s = half_up(v, d)                                 # round half up on the decimal value; a negative that rounds
    return r"\ensuremath{-}" + s[1:] if s.startswith("-") else s   # to zero keeps its sign (-0.000)


def half_up(v, d):
    q = Decimal(repr(float(v))).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP)
    s = f"{q:.{d}f}"
    s = s.lstrip("-")
    return ("-" + s) if float(v) < 0 else s          # an exact zero is unsigned, as in make_supplement_tables.py


def bkey(b):
    return f"{float(b):g}" if float(b) != 1.0 else "1"


lines = []


def put(key, val):
    if val is not None:
        lines.append(rf"\expandafter\def\csname R:{key}\endcsname{{{val}}}")


def table(name, df, keys):
    for _, r in df.iterrows():
        k = "/".join(bkey(r[c]) if c == "beta" else str(r[c]) for c in keys)
        for col in df.columns:
            if col in keys:
                continue
            put(f"{name}/{k}/{col}", fmt(col, r[col]))
            if col in PCOLS and not pd.isna(r[col]):          # "= 0.062" / "< 0.001", for use after $p$
                v = fmt(col, r[col])
                put(f"{name}/{k}/{col}_rel", r"\ensuremath{{}<{}}" + v[len("\\ensuremath{<}"):]
                    if v.startswith("\\ensuremath{<}") else r"\ensuremath{{}={}}" + v)
            if col == "relative_loss" and not pd.isna(r[col]):     # the same value as a percentage
                put(f"{name}/{k}/{col}_pct", fmt("pct", 100 * r[col]))


for name, keys in (("frontier", ["universe", "cell", "beta"]), ("h1", ["universe", "cell", "beta"]),
                   ("h2detail", ["universe", "cell", "beta"]), ("manip", ["universe", "cell", "beta"]),
                   ("newsplits", ["universe", "cell", "beta"]), ("health", ["universe", "cell", "beta"]),
                   ("randref", ["universe", "objective"]), ("h3", ["cell"])):
    f = {"frontier": "frontier", "h1": "h1_coverage_loss", "h2detail": "h2_detail", "manip": "manipulation_check",
         "newsplits": "new_splits_only", "health": "engine_health", "randref": "random_search_reference",
         "h3": "h3_alignment"}[name]
    table(name, pd.read_csv(RES / f"{f}.csv"), keys)
rb = pd.read_csv(RES / "recall_by_prior_percentile.csv")
table("rbins", rb.drop(columns=["pctl_hi"]), ["universe", "cell", "beta", "pctl_lo"])
for (u, c, bb), g in rb.groupby(["universe", "cell", "beta"]):
    put(f"rbelow95/{u}/{c}/{bkey(bb)}/found", str(int(g[g.pctl_lo < 0.95].found.sum())))
b = pd.read_csv(RES / "budget_sensitivity.csv")
table("budget", b, ["universe", "cell", "beta", "acquisitions"])

# derived: sets and counts (selection only, no new statistics)
h2 = pd.read_csv(RES / "h2_detail.csv")
for (u, c), g in h2.groupby(["universe", "cell"]):
    for col, name in (("beneficial_and_safe", "h2set"), ("safe_at_10pct", "safe10"), ("safe_at_25pct", "safe25"),
                      ("safe", "safe20")):
        sel = [f"{x:g}" for x in g[g[col].astype(str) == "True"].beta]
        put(f"{name}/{u}/{c}", "\\{" + ", ".join(sel) + "\\}" if sel else "none")
    for col, name in (("safe_at_10pct", "bs10"), ("safe_at_25pct", "bs25")):
        sel = [f"{x:g}" for x in g[(g[col].astype(str) == "True") & (g.beneficial.astype(str) == "True")].beta]
        put(f"{name}/{u}/{c}", "\\{" + ", ".join(sel) + "\\}" if sel else "none")
for (u, bb), g in h2.groupby(["universe", "beta"]):
    put(f"benefcount/{u}/{bkey(bb)}", str(int((g.beneficial.astype(str) == "True").sum())))
    put(f"safecount/{u}/{bkey(bb)}", str(int((g.safe.astype(str) == "True").sum())))
for u, g in h2.groupby("universe"):          # largest number of cells in which one weight is beneficial and safe
    put(f"h2maxcellsperbeta/{u}", str(int(g[g.beneficial_and_safe.astype(str) == "True"].groupby("beta").cell.nunique().max())) if (g.beneficial_and_safe.astype(str) == "True").any() else "0")
for u, g in h2.groupby("universe"):          # largest weight that is beneficial and safe in any cell
    q = g[g.beneficial_and_safe.astype(str) == "True"]
    put(f"h2maxbeta/{u}", f"{q.beta.max():g}" if len(q) else "none")
h1 = pd.read_csv(RES / "h1_coverage_loss.csv")
for u, g in h1.groupby("universe"):
    put(f"h1count/{u}", str(int((g.reject_H0.astype(str) == "True").sum())) if u == "primary" else None)
    put(f"h1cells/{u}", str(g[g.reject_H0.astype(str) == "True"].cell.nunique()) if u == "primary" else None)
    put(f"h1n/{u}", str(len(g)))
h2b = pd.read_csv(RES / "h2_beneficial_and_safe.csv")
for u, g in h2b.groupby("universe"):
    put(f"h2nonempty/{u}", str(int((g.beneficial_and_safe_betas != "[]").sum())))
    put(f"h2applicable/{u}", str(int((g.applicable.astype(str) == "True").sum())))
    put(f"h2cells/{u}", str(len(g)))
    for _, r in g.iterrows():
        put(f"h2base/{u}/{r.cell}/baseline_R_dis", fmt("R", r.baseline_R_dis))
        put(f"h2base/{u}/{r.cell}/chance", fmt("R", r.chance_recall))
eh = pd.read_csv(RES / "engine_health.csv").dropna()
for u, g in eh.groupby("universe"):
    put(f"healthrange/{u}/lb_min", fmt("pct", 100 * g.share_length_scale_at_lower_bound.min()))
    put(f"healthrange/{u}/lb_max", fmt("pct", 100 * g.share_length_scale_at_lower_bound.max()))
    put(f"healthrange/{u}/flat_max", fmt("pct", 100 * g.share_flat_EI.max()))
s = json.loads((RES / "summary.json").read_text())
put("summary/n_runs", f"{s['n_runs']:,}".replace(",", "{,}"))
put("summary/tau", fmt("tau", s["h3_kendall_tau"]))
# H3 is registered without a test; an exact Kendall p is given for orientation only (manuscript says so)
from scipy.stats import kendalltau
h3 = pd.read_csv(RES / "h3_alignment.csv")
kt = kendalltau(h3.history_rho, h3["gain_b0.5"], method="exact")
assert abs(kt.statistic - s["h3_kendall_tau"]) < 1e-9
put("summary/tau_p_exact", half_up(kt.pvalue, 2))
# N5 (registered secondary): cells whose mean N5 at beta exceeds that at beta = 0, descriptive count only
frn = pd.read_csv(RES / "frontier.csv")
for (u, bb), g in frn[frn.beta > 0].groupby(["universe", "beta"]):
    base = frn[(frn.universe == u) & (frn.beta == 0)].set_index("cell").N5
    up = g.set_index("cell").N5 > base.reindex(g.cell).values
    put(f"n5up/{u}/{bkey(bb)}", str(int(up.sum())))
    h1r = pd.read_csv(RES / "h1_coverage_loss.csv")
    h1r = h1r[(h1r.universe == u) & (h1r.beta == bb)].set_index("cell").reject_H0.astype(str) == "True"
    if len(h1r):                                     # of the cells with higher N5, those whose loss H1 detected
        put(f"n5up_h1/{u}/{bkey(bb)}", str(int(sum(1 for c in up[up].index if h1r.get(c, False)))))
NONCONV = "VB fitting did not converge"
put("summary/glmm_nonconverged", str(sum(1 for m in s["glmm"].values() if NONCONV in m["optimizer_warnings"])))
put("summary/analysis_sha", s["analysis_sha256"][:8])
for c, meta in s["glmm"].items():
    for k in ("n_rows", "n_runs", "n_materials"):
        put(f"summary/glmm/{c}/{k}", f"{meta[k]:,}".replace(",", "{,}"))
rows_ = [m["n_rows"] for m in s["glmm"].values()]
put("summary/glmm/n_rows_min", f"{min(rows_):,}".replace(",", "{,}"))
put("summary/glmm/n_rows_max", f"{max(rows_):,}".replace(",", "{,}"))
put("summary/glmm/n_runs_all", str(int(max(m["n_runs"] for m in s["glmm"].values())))
     if len({m["n_runs"] for m in s["glmm"].values()}) == 1 else None)
put("summary/glmm/n_materials_all", str(int(max(m["n_materials"] for m in s["glmm"].values())))
     if len({m["n_materials"] for m in s["glmm"].values()}) == 1 else None)
CELLNAME = {"Y1:P_H_fixed": r"\cellA{}", "Y1:P_phys_fixed": r"\cellB{}", "Y2:P_H_fixed": r"\cellC{}",
            "Y2:P_H_matched_Y2": r"\cellD{}", "Y2:P_phys_fixed": r"\cellE{}", "Y2:P_phys_matched_Y2": r"\cellF{}"}
nc = [CELLNAME[c] for c, m in s["glmm"].items() if NONCONV in m["optimizer_warnings"]]
put("summary/glmm_nonconverged_cells", " and ".join([", ".join(nc[:-1]), nc[-1]]) if len(nc) > 1 else (nc[0] if nc else "none"))

# post-hoc (decision-log entry 42) — separate prefix
c = pd.read_csv(POST / "t99_censoring.csv")
for _, r in c.iterrows():
    put(f"posthoc/censor/{r.universe}/{r.cell}/{bkey(r.beta)}", fmt("pct", 100 * r.no_hit_share))
i = pd.read_csv(POST / "init_top1.csv")
for _, r in i.iterrows():
    put(f"posthoc/inittop1/{r.universe}/{r.objective}", str(int(r.blocks_init_has_top1)))
    put(f"posthoc/initblocks/{r.universe}/{r.objective}", str(int(r.blocks)))
sv = pd.read_csv(POST / "init_design_sensitivity.csv")
for u, g in sv.groupby("universe"):
    put(f"posthoc/initsens/maxshift/{u}", f"{(g.gain_init_as_time0 - g.gain_registered)[g.beta < 1].abs().max():.2f}")
put("posthoc/initsens/maxshift", f"{(sv.gain_init_as_time0 - sv.gain_registered)[sv.beta < 1].abs().max():.2f}")
flip = sv[(sv.beta < 1) & (sv.beneficial_registered != sv.beneficial_init_as_time0)]
put("posthoc/initsens/nflips/primary", str(int((flip.universe == "primary").sum())))

bl = pd.read_csv(POST / "blind_spot_loss.csv")
table("posthoc/blind", bl, ["universe", "cell", "beta"])
for (u, bb), g in bl.groupby(["universe", "beta"]):
    put(f"posthoc/blindcount/{u}/{bkey(bb)}", str(int((g.L_blind_lo > 0).sum())))

ha = pd.read_csv(POST / "h2_adjusted_beneficial.csv")
table("posthoc/h2adj", ha, ["universe", "cell", "beta"])
for (u, c), g in ha.groupby(["universe", "cell"]):
    for col in ("bs20_adjusted", "bs10_adjusted", "bs25_adjusted"):
        sel = [f"{x:g}" for x in g[g[col].astype(str) == "True"].beta]
        put(f"posthoc/{col}/{u}/{c}", "\\{" + ", ".join(sel) + "\\}" if sel else "none")
for u, g in ha.groupby("universe"):
    put(f"posthoc/h2adjnonempty/{u}", str(int(g[g.bs20_adjusted.astype(str) == "True"].cell.nunique())))

# post-hoc (decision-log entry 45): estimates within each split type, descriptive
st = pd.read_csv(POST / "by_split_type.csv")
for _, r in st.iterrows():
    k = f"posthoc/splittype/{r.universe}/{r.cell}/{bkey(r.beta)}/{r.quantity}/{r.split_type}"
    for col in ("estimate", "ci_lo", "ci_hi"):
        put(f"{k}/{col}", fmt("gain" if r.quantity == "gain" else "L", r[col]))
fr = pd.read_csv(RES / "frontier.csv")
m = st.merge(fr[["universe", "cell", "beta", "gain", "L"]], on=["universe", "cell", "beta"])
m["registered"] = [r.gain if r.quantity == "gain" else r.L for r in m.itertuples()]
for u, g in m.groupby("universe"):
    mod = g[g.beta.isin([0.25, 0.5])]
    put(f"posthoc/splittype/signagree/{u}", str(int((mod.estimate * mod.registered > 0).sum())))
    put(f"posthoc/splittype/n/{u}", str(len(mod)))
    put(f"posthoc/splittype/includezero/{u}", str(int(((mod.ci_lo <= 0) & (mod.ci_hi >= 0)).sum())))
    opp = g[((g.registered > 0) & (g.ci_hi < 0)) | ((g.registered < 0) & (g.ci_lo > 0))]
    put(f"posthoc/splittype/opposite/{u}", str(len(opp)))
    put(f"posthoc/splittype/n_all/{u}", str(len(g)))

# descriptive facts about the frozen benchmark (code/describe_benchmark.py; data only, no search outcomes)
for r in pd.read_csv(POST / "benchmark_facts.csv").itertuples():
    put(f"bench/{r.key}", f"{int(r.value):,}".replace(",", "{,}"))

head = [r"% generated by code/make_numbers.py from results/ and results_posthoc/ -- do not edit by hand",
        r"\makeatletter",
        r"\newcommand{\R}[1]{\ifcsname R:#1\endcsname\csname R:#1\endcsname\else"
        r"\textbf{??\detokenize{#1}}\PackageWarning{numbers}{missing result key \detokenize{#1}}\fi}",
        r"\makeatother"]
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text("\n".join(head + lines) + "\n", encoding="utf-8")
print(len(lines), "macros written to", OUT)

"""
make_figures.py — Paper 2 figures, drawn only from the committed results/ tables
================================================================================
Reads the registered outputs of analyze_runs.py v1.1 (results/*.csv) and draws:

  fig1_frontier       efficiency-coverage frontier per objective (mean gain in T99
                      vs mean recall of S_dis across beta, 95% split-level t intervals)
  fig2_h1_forest      H1: mean coverage loss L at beta = 0.25 and 0.5, 95% t intervals,
                      Holm decisions
  fig3_h2_grid        H2: beneficial / safe per cell and beta (registered rule)
  fig4_recall_bins    recall of S_top members by prior-percentile bin, pooled over runs

No number is computed here beyond selecting rows and the 0.8 x baseline reference.
Usage (from code/):  python make_figures.py ../results ../figures
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from decimal import Decimal, ROUND_HALF_UP


def r1(v):  # one decimal, half away from zero, as in numbers.tex
    return f"{Decimal(repr(float(v))).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP):.1f}".replace('-', '\u2212')

RES = Path(sys.argv[1] if len(sys.argv) > 1 else "../results")
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "../figures")
OUT.mkdir(parents=True, exist_ok=True)

# palette (validated: categorical blue/orange; ordinal blue ramp 250/400/550/700)
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#ffffff"
FAMILY = {"P_H": "#2a78d6", "P_phys": "#eb6834"}
RAMP = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
CELLS = ["Y1:P_H_fixed", "Y1:P_phys_fixed", "Y2:P_H_fixed", "Y2:P_H_matched_Y2",
         "Y2:P_phys_fixed", "Y2:P_phys_matched_Y2"]
LABEL = {"Y1:P_H_fixed": "historical", "Y1:P_phys_fixed": "elemental-modulus",
         "Y2:P_H_fixed": "historical", "Y2:P_H_matched_Y2": "historical, matched",
         "Y2:P_phys_fixed": "elemental-modulus", "Y2:P_phys_matched_Y2": "elemental-modulus, matched"}
OBJ = {"Y1": "Y1: shear modulus $G_{\\mathrm{VRH}}$", "Y2": "Y2: specific modulus $E/\\rho$"}


def style(c):
    fam = "P_phys" if "P_phys" in c else "P_H"
    matched = "matched" in c
    return dict(color=FAMILY[fam], ls="--" if matched else "-", marker="s" if matched else "o")


plt.rcParams.update({"font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
                     "savefig.facecolor": SURF, "axes.facecolor": SURF})


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def grid(ax, axis="y"):
    ax.grid(axis=axis, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


fr = pd.read_csv(RES / "frontier.csv")
fr = fr[fr.universe == "primary"]
ref = pd.read_csv(RES / "random_search_reference.csv")
chance = float(ref[ref.universe == "primary"].R_dis.mean())

# ---------------------------------------------------------------- Figure 1
XMIN, XMAX = {"Y1": -5, "Y2": -6}, {"Y1": 21, "Y2": 13}
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.4), sharey=True)
for ax, obj in zip(axes, ("Y1", "Y2")):
    grid(ax, "both")
    ax.axhline(chance, color=MUTED, lw=0.8, ls=":")
    ax.set_xlim(XMIN[obj], XMAX[obj])
    ax.text(0.99, chance + 0.002, "random search", transform=ax.get_yaxis_transform(), ha="right", va="bottom",
            fontsize=6.5, color=INK2)
    lo = XMIN[obj]
    offscale = []
    for c in [c for c in CELLS if c.startswith(obj)]:
        d = fr[fr.cell == c].sort_values("beta").copy()
        off = d.gain < lo
        offscale += [f"{LABEL[c]}: {r1(r.gain)}" for r in d[off].itertuples()]
        s = style(c)
        on = d[~off]                                  # off-scale points are not drawn at a false position
        ax.errorbar(on.gain, on.R_dis, xerr=[on.gain - on.gain_lo, on.gain_hi - on.gain],
                    yerr=[on.R_dis - on.R_dis_lo, on.R_dis_hi - on.R_dis], fmt="none", ecolor=s["color"],
                    elinewidth=0.6, alpha=0.45, capsize=0)
        ax.plot(on.gain, on.R_dis, color=s["color"], ls=s["ls"], lw=1.4, marker=s["marker"], ms=4.5,
                mec=SURF, mew=0.8, label=LABEL[c])
        if off.any():                                 # faint dotted lead to an open arrow at the axis edge
            prev = on.iloc[-1]
            for r in d[off].itertuples():
                ax.plot([prev.gain, lo], [prev.R_dis, r.R_dis], color=s["color"], ls=":", lw=0.9, alpha=0.6)
                k = len(offscale)                     # k-th off-scale point of this panel: k-th arrow left of the axis
                ax.plot(lo - 0.03 * (XMAX[obj] - XMIN[obj]) * k, r.R_dis, marker="<",
                        ms=6, mfc="none", mec=s["color"], mew=1.0, ls="", clip_on=False)
        last = d.iloc[-1]
        if not off.iloc[-1]:
            ax.annotate("β = 1", (last.gain, last.R_dis), xytext=(3, -2), textcoords="offset points",
                        fontsize=6, color=INK2)
    if offscale:
        ax.text(0.02, 0.98, "◁ off scale at β = 1 (open arrow):\n" + "\n".join(offscale), transform=ax.transAxes,
                fontsize=6, color=INK2, va="top")
    ax.set_title(OBJ[obj], loc="left")
    ax.set_xlabel("reduction in restricted mean $T_{99}$ (acquisitions)")
axes[0].set_ylabel("mean recall of $S_{\\mathrm{dis}}$")
axes[0].annotate("β = 0", (0, fr[(fr.cell == "Y1:P_H_fixed") & (fr.beta == 0)].R_dis.iloc[0]),
                 xytext=(4, 4), textcoords="offset points", fontsize=6, color=INK2)
fig.tight_layout(rect=(0, 0.07, 1, 1))
h, l = axes[1].get_legend_handles_labels()
fig.legend(h, l, frameon=False, loc="lower center", ncol=4, handlelength=2.6, bbox_to_anchor=(0.5, 0))
save(fig, "fig1_frontier")

# ---------------------------------------------------------------- Figure 2
h1 = pd.read_csv(RES / "h1_coverage_loss.csv")
h1 = h1[h1.universe == "primary"].copy()
h1["order"] = h1.cell.map(CELLS.index) * 2 + (h1.beta == 0.5)
h1 = h1.sort_values("order", ascending=False).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(4.6, 3.6))
ax.axvline(0, color=MUTED, lw=0.8)
ax.grid(axis="x", color=GRID, lw=0.6)
ax.set_axisbelow(True)
for i, r in h1.iterrows():
    s = style(r.cell)
    ax.plot([r.ci_lo, r.ci_hi], [i, i], color=s["color"], lw=1.4, solid_capstyle="round")
    rej = str(r.reject_H0) == "True"
    ax.plot(r.mean_L, i, marker=s["marker"], ms=6, color=s["color"], mfc=s["color"] if rej else SURF,
            mec=s["color"], mew=1.2)
    ax.text(1.02, i, f"{r.p_holm:.2g}" if r.p_holm >= 1e-3 else "<0.001", transform=ax.get_yaxis_transform(),
            va="center", fontsize=6.5, color=INK if rej else INK2, fontweight="bold" if rej else "normal")
ax.text(1.02, len(h1) - 0.3, "Holm p", transform=ax.get_yaxis_transform(), fontsize=6.5, color=INK2, va="bottom")
ax.set_yticks(range(len(h1)))
ax.set_yticklabels([f"{r.cell.split(':')[0]} {LABEL[r.cell]}, β = {r.beta:g}" for r in h1.itertuples()])
ax.set_xlabel("mean coverage loss L(β) = R_dis(0) − R_dis(β)")
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=INK2, mfc=INK2, label="H₀ rejected (Holm, α = 0.05)"),
                   Line2D([], [], marker="o", ls="", color=INK2, mfc=SURF, label="not rejected")],
          frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, fontsize=6.5)
fig.tight_layout()
save(fig, "fig2_h1_forest")

# ---------------------------------------------------------------- Figure 3
h2 = pd.read_csv(RES / "h2_detail.csv")
h2 = h2[h2.universe == "primary"]
betas = [0.02, 0.05, 0.1, 0.25, 0.5, 0.75]
fill = {2: "#2a78d6", 1: "#b7d3f6", 0: "#f0efec"}
fig, ax = plt.subplots(figsize=(4.8, 2.6))
for i, c in enumerate(CELLS):
    for j, b in enumerate(betas):
        r = h2[(h2.cell == c) & (h2.beta == b)].iloc[0]
        ben, safe = str(r.beneficial) == "True", str(r.safe) == "True"
        n = ben + safe
        ax.add_patch(plt.Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, color=fill[n], lw=0))
        txt = "B + S" if n == 2 else ("B" if ben else ("S" if safe else "–"))
        ax.text(j + 0.5, i + 0.5, txt, ha="center", va="center", fontsize=7,
                color="#ffffff" if n == 2 else (INK if n == 1 else MUTED), fontweight="bold" if n == 2 else "normal")
ax.set_xlim(0, len(betas)); ax.set_ylim(len(CELLS), 0)
ax.set_xticks([j + 0.5 for j in range(len(betas))]); ax.set_xticklabels([f"{b:g}" for b in betas])
ax.set_yticks([i + 0.5 for i in range(len(CELLS))])
ax.set_yticklabels([f"{c.split(':')[0]} {LABEL[c]}" for c in CELLS])
ax.set_xlabel("prior weight β")
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
fig.text(0.02, -0.02, "B = beneficial: lower 95% bound of mean ΔT99 > 0.   S = safe: loses < 20% of baseline "
         "recall of $S_{\\mathrm{dis}}$ (upper bound at level 1 − 0.05/6).", fontsize=6, color=INK2, va="top")
fig.tight_layout()
save(fig, "fig3_h2_grid")

# ---------------------------------------------------------------- Figure 4
rb = pd.read_csv(RES / "recall_by_prior_percentile.csv")
rb = rb[rb.universe == "primary"]
show = [0.0, 0.1, 0.5, 1.0]
bins = sorted(rb.pctl_lo.unique())
blab = ["< 0.5", "0.5–0.8", "0.8–0.9", "0.9–0.95", "≥ 0.95"]
fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.0), sharey=True)
for ax, c in zip(axes.flat, CELLS):
    grid(ax)
    d = rb[rb.cell == c]
    for b, col in zip(show, RAMP):
        e = d[d.beta == b].set_index("pctl_lo").loc[bins]
        ax.plot(range(len(bins)), e.recall, color=col, lw=1.4, marker="o", ms=3.5, mec=SURF, mew=0.6,
                label=f"β = {b:g}")
    ax.set_title(f"{c.split(':')[0]} {LABEL[c]}", loc="left")
    ax.set_xticks(range(len(bins))); ax.set_xticklabels(blab, rotation=30, ha="right")
for ax in axes[:, 0]:
    ax.set_ylabel("recall of top-5% materials")
for ax in axes[1, :]:
    ax.set_xlabel("prior percentile of the material")
fig.tight_layout(rect=(0, 0, 1, 0.95))
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.0))
save(fig, "fig4_recall_bins")
print("figures written to", OUT)

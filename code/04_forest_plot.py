# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

"""04_forest_plot.py — generate forest plot for MUC5 VNTR haplogroup paper.

Inputs: none (numbers hard-coded from 02_primary_analysis.py and
        03_supplementary_analyses.py outputs).
Output: muc5_forest_plot.png

All ORs are from logistic regression models adjusted for age, ever-smoking,
PC1–5, and recruitment centre (N=798). Joint MUC5AC model uses H1 as reference.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Data ─────────────────────────────────────────────────────────────────────
# (label, OR, CI_lo, CI_hi, p-value, group)
rows = [
    # rs35705950 reference
    ("rs35705950 T allele",      3.96, 2.93, 5.34, 3.01e-19, "rs35705950"),

    # MUC5AC main effects
    ("MUC5AC H3 (unadj.)",       1.69, 1.31, 2.17, 4.2e-5,   "MUC5AC main"),
    ("MUC5AC H3 (+rs35705950)",  1.17, 0.89, 1.54, 0.27,     "MUC5AC main"),

    # MUC5AC joint model — H1 as reference
    ("Joint: H2 vs H1",          0.62, 0.47, 0.82, 7.6e-4,   "MUC5AC joint"),
    ("Joint: H3 vs H1",          1.50, 1.15, 1.94, 2.4e-3,   "MUC5AC joint"),
    ("Joint+rs: H2 vs H1",       0.71, 0.53, 0.95, 0.022,    "MUC5AC joint"),
    ("Joint+rs: H3 vs H1",       1.08, 0.81, 1.44, 0.59,     "MUC5AC joint"),

    # MUC5B
    ("MUC5B H1 (unadj.)",        0.72, 0.52, 0.99, 0.042,    "MUC5B"),
    ("MUC5B H1 (+rs35705950)",   0.81, 0.58, 1.14, 0.23,     "MUC5B"),
]

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 6.5))

colors = {
    "rs35705950":   "#222222",
    "MUC5AC main":  "#1a6faf",
    "MUC5AC joint": "#4ca3dd",
    "MUC5B":        "#c0392b",
}
sig_marker = "D"
ns_marker  = "o"

y_positions = list(range(len(rows)))[::-1]

for i, (label, OR, lo, hi, p, group) in enumerate(rows):
    y   = y_positions[i]
    col = colors[group]
    marker = sig_marker if p < 0.05 else ns_marker
    msize  = 7           if p < 0.05 else 6
    mfc    = col         if p < 0.05 else "white"
    ax.errorbar(OR, y, xerr=[[OR - lo], [hi - OR]],
                fmt=marker, color=col, markerfacecolor=mfc,
                markeredgecolor=col, markersize=msize,
                elinewidth=1.2, capsize=3, capthick=1.2)
    pstr = f"p={p:.2e}" if p < 0.001 else f"p={p:.3f}"
    ax.text(hi + 0.05, y, pstr, va='center', ha='left', fontsize=7.5, color='#444')

ax.set_yticks(y_positions)
ax.set_yticklabels([r[0] for r in rows], fontsize=8.5)
ax.axvline(1.0, color='grey', linewidth=0.8, linestyle='--')
ax.set_xscale('log')
ax.set_xlabel("Odds ratio (95% CI)", fontsize=9)
ax.set_xlim(0.28, 10)
ax.set_xticks([0.3, 0.5, 1.0, 2.0, 4.0])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

for y_div in [6.5, 5.5, 2.5, 1.5]:
    ax.axhline(y_div, color='#cccccc', linewidth=0.6)

patches = [mpatches.Patch(color=v, label=k) for k, v in colors.items()]
patches += [
    plt.Line2D([0], [0], marker=sig_marker, color='grey',
               markerfacecolor='grey', markersize=7, label='p<0.05 (filled ◆)'),
    plt.Line2D([0], [0], marker=ns_marker, color='grey',
               markerfacecolor='white', markeredgecolor='grey',
               markersize=6, label='p≥0.05 (open ●)'),
]
ax.legend(handles=patches, fontsize=7.5, loc='lower right', framealpha=0.9)

ax.set_title(
    "MUC5AC and MUC5B VNTR haplogroup associations with IPF\n"
    "(N=798; adjusted for age, smoking, PC1–5, centre; joint model H1 as reference)",
    fontsize=8.5, pad=8)

plt.tight_layout()
plt.savefig('./muc5_forest_plot.png', dpi=150,
            bbox_inches='tight')
print("Saved muc5_forest_plot.png")

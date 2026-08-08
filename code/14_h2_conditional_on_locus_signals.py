# ---
# jupyter:
#   jupytext:
#     formats: py:percent,ipynb
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 14. Does the MUC5AC H2 residual survive conditioning on ALL independent 11p15.5 signals?
#
# The joint model conditions the H2 residual on rs35705950 alone. But 11p15.5 is a
# multi-signal locus: Kousathanas et al. (2025) GCTA-COJO identified several
# conditionally-independent IPF signals here beyond the promoter (in MUC6, MUC2,
# and MUC5B). "Independent of rs35705950" is therefore not the same as "independent
# of the locus" — the H2 tag could be indexing one of those other signals.
#
# This tests it directly: the H2 tag conditioned on rs35705950 plus the eight
# independent 11p15.5 signals from Kousathanas Table S5. If H2 attenuates, it was
# tagging a neighbour; if it survives, it is distinct from every known signal here.
#
# Requires `kous_independent_signals.raw` from REPRODUCE.md step 3c.

# %%
import warnings
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")
ROOT = "/home/drcjar/ipfjes-omics"

raw = pd.read_csv(f"{ROOT}/kous_independent_signals.raw", sep=r"\s+")
raw["pid"] = raw["IID"].astype(str).str.split("_").str[0]
SNP = {c: "S_" + c.split(":")[1] for c in raw.columns if c.startswith("chr11:")}
raw = raw.rename(columns=SNP)
MUC5B = ["S_1208601", "S_1214934", "S_1219775", "S_1224068"]
MUC2_6 = ["S_1013992", "S_1057948", "S_1078249", "S_1081750"]   # MUC6 + MUC2 x3

clin = pd.read_csv(f"{ROOT}/flat_data_genotype_subset.csv"); clin["pid"] = clin.participant_id.astype(str)
hap = pd.read_csv(f"{ROOT}/muc5_haplogroups.csv"); hap["pid"] = hap.IID.str.split("_").str[0]
pcs = pd.read_csv(f"{ROOT}/code/ipfjes_pca.eigenvec", sep=r"\s+"); pcs["pid"] = pcs.IID.astype(str)
# Primary H2 tag rs879201833 (chr11:1,182,279, r2=0.888) — matches the manuscript's primary analysis.
h2raw = pd.read_csv(f"{ROOT}/muc5ac_h2_tag.raw", sep=r"\s+"); h2raw["pid"] = h2raw.IID.astype(str).str.split("_").str[0]
H2COL = [c for c in h2raw.columns if c.startswith("chr11:1182279")][0]

d = (clin
     .merge(hap[["pid", "H3_count", "rs35705950_T_count"]], on="pid")
     .merge(h2raw[["pid", H2COL]], on="pid")
     .merge(raw[["pid"] + MUC5B + MUC2_6], on="pid")
     .merge(pcs[["pid", "PC1", "PC2", "PC3", "PC4", "PC5"]], on="pid")
     .dropna(subset=["age", "ever_smoked", "peto_exposed", "centre"]))
d = d[d.centre != 16].copy()
f0 = d.loc[d.case == 0, H2COL].mean() / 2
d["H2d"] = d[H2COL] if f0 < 0.5 else 2 - d[H2COL]   # orient to the H2-tracking allele
COVS = "age + ever_smoked + PC1+PC2+PC3+PC4+PC5 + C(centre)"


def fit(f):
    try:
        return smf.logit(f, data=d).fit(disp=0)
    except Exception:
        return smf.logit(f, data=d).fit(method="bfgs", disp=0)


def orci(m, t):
    b, se, p = m.params[t], m.bse[t], m.pvalues[t]
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), p


# %%
print(f"N={len(d)} (cases={int(d.case.sum())}, controls={int((d.case == 0).sum())})\n")
print("MUC5AC H2 (direct tag rs28519516) under progressive conditioning on Kousathanas 2025")
print("independent 11p15.5 signals:\n")
for lab, adj in [("+ rs35705950 only",               "rs35705950_T_count"),
                 ("+ rs35705950 + MUC5B signals",     "rs35705950_T_count+" + "+".join(MUC5B)),
                 ("+ rs35705950 + MUC2/MUC6 signals",  "rs35705950_T_count+" + "+".join(MUC2_6)),
                 ("+ rs35705950 + ALL 8 signals",      "rs35705950_T_count+" + "+".join(MUC5B + MUC2_6))]:
    o, l, h, p = orci(fit(f"case~H2d+H3_count+{adj}+{COVS}"), "H2d")
    print(f"  {lab:34s} H2 OR {o:.2f} ({l:.2f}-{h:.2f})  p={p:.3f}")

print("\nLD (r) of the H2 tag with each independent signal (all |r|<0.4 -> H2 is orthogonal):")
for s in MUC5B + MUC2_6:
    print(f"  {s}: r={np.corrcoef(d.H2d, d[s])[0, 1]:+.3f}")
print("DONE")

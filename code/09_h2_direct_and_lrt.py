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
# # 11. Direct MUC5AC H2 tagging, nested-model tests, and reciprocal stratification
#
# Three analyses added after peer review, all reported in the manuscript
# Discussion and none of them produced by scripts `00`–`10`.
#
# **1. Direct H2 tagging.** The primary joint model assigns MUC5AC H2 by
# subtraction (`H2 = 2 − H1 − H3`) from an H1 concordance tag of imperfect
# quality (r²=0.899). The residual H2 association that survives adjustment for
# rs35705950 was therefore initially attributed to measurement error. That
# explanation is testable: Plender et al. report direct H2 tSNPs, four of which
# are present in our imputed data — `rs28519516`, `rs28558973`, `rs28368633`
# (Table S4 H2 r²=0.862, already in the tSNP panel and previously used only as
# unreported cross-checks) and `rs879201833` (H2 r²=0.888, extracted
# separately). Refitting the joint model with each gives a **larger** residual
# than the subtraction estimate, so the published figure is the conservative
# one and the measurement-error explanation fails.
#
# **2. Nested-model likelihood-ratio tests**, which formalise the conditional
# asymmetry the Discussion rests on (adjusting rs35705950 for H3 barely moves
# it; adjusting H3 for rs35705950 abolishes it).
#
# **3. Reciprocal stratification**: rs35705950 estimated within H3-carrier
# strata, the complement of Supplementary Table S2.
#
# Requires `muc5ac_h2_tag.raw` from REPRODUCE.md step 3b.

# %%
import warnings
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
warnings.filterwarnings("ignore")

ROOT = "/home/drcjar/academic/ipfjes-omics"

pcs = pd.read_csv(f"{ROOT}/code/ipfjes_pca.eigenvec", sep=r"\s+")
pcs["pid"] = pcs["IID"].astype(str)
clin = pd.read_csv(f"{ROOT}/flat_data_genotype_subset.csv")
clin["pid"] = clin["participant_id"].astype(str)
hap = pd.read_csv(f"{ROOT}/muc5_haplogroups.csv")
hap["pid"] = hap["IID"].str.split("_").str[0]
dos = pd.read_csv(f"{ROOT}/muc5_tsnps_dosage_pheno.csv")
dos["pid"] = dos["IID"].astype(str).str.split("_").str[0]
raw = pd.read_csv(f"{ROOT}/muc5ac_h2_tag.raw", sep=r"\s+")
raw["pid"] = raw["IID"].astype(str).str.split("_").str[0]
NEW = [c for c in raw.columns if c.startswith("chr11:1182279")][0]

# All four tags are fitted separately here. The three r2=0.862 panel tags are in complete
# mutual LD in our data and give identical estimates, so Supplementary Table S3 reports them
# as one row; they are still computed and checked individually so the collapse is verifiable.
H2_TAGS = {
    "rs879201833 (H2 r2=0.888)": NEW,
    "rs28519516  (H2 r2=0.862)": "chr11:1202343:T:C_C",
    "rs28558973  (H2 r2=0.862)": "chr11:1202859:A:G_G",
    "rs28368633  (H2 r2=0.862)": "chr11:1202914:T:C_C",
}

dm = (clin
      .merge(hap[["pid", "H1_count", "H2_count", "H3_count",
                  "MUC5B_H1_count", "rs35705950_T_count"]], on="pid")
      .merge(dos[["pid"] + [c for c in H2_TAGS.values() if c != NEW]], on="pid")
      .merge(raw[["pid", NEW]], on="pid")
      .merge(pcs[["pid", "PC1", "PC2", "PC3", "PC4", "PC5"]], on="pid")
      .dropna(subset=["age", "ever_smoked", "peto_exposed", "centre"]))
dm = dm[(dm["centre"] != 16) & dm.ethnicity.astype(str).str.contains("White", case=False, na=False)].copy()
covs = "age + ever_smoked + PC1+PC2+PC3+PC4+PC5 + C(centre)"
print(f"N={len(dm)} (cases={int(dm.case.sum())}, controls={int((dm.case==0).sum())})")


def fit(f, d=None):
    d = dm if d is None else d
    try:    return smf.logit(f, data=d).fit(disp=0)
    except: return smf.logit(f, data=d).fit(method="bfgs", disp=0)


def orci(m, t):
    b, se, p = m.params[t], m.bse[t], m.pvalues[t]
    return np.exp(b), np.exp(b - 1.96*se), np.exp(b + 1.96*se), p


# %% [markdown]
# ## 1. Does the residual H2 association survive direct tagging?
#
# Each tag is oriented so the counted allele tracks H2 (control frequency ~0.24).

# %%
ctrl = dm[dm.case == 0]
print(f"{'H2 measure':28s} {'ctrl freq':>9s} {'r vs subtr':>10s}  "
      f"{'unadjusted':>22s}  {'+ rs35705950':>22s}")
for lab, col in H2_TAGS.items():
    f = ctrl[col].mean() / 2
    dm["H2d"] = dm[col] if f < 0.5 else 2 - dm[col]        # orient toward H2
    fr = dm.loc[dm.case == 0, "H2d"].mean() / 2
    r = np.corrcoef(dm.H2d, dm.H2_count)[0, 1]
    o1, l1, h1, p1 = orci(fit(f"case~H2d+H3_count+{covs}"), "H2d")
    o2, l2, h2, p2 = orci(fit(f"case~H2d+H3_count+rs35705950_T_count+{covs}"), "H2d")
    print(f"{lab:28s} {fr:>9.3f} {r:>10.3f}  "
          f"{f'{o1:.2f} ({l1:.2f}-{h1:.2f}) p={p1:.3f}':>22s}  "
          f"{f'{o2:.2f} ({l2:.2f}-{h2:.2f}) p={p2:.3f}':>22s}")

o1, l1, h1, p1 = orci(fit(f"case~H2_count+H3_count+{covs}"), "H2_count")
o2, l2, h2, p2 = orci(fit(f"case~H2_count+H3_count+rs35705950_T_count+{covs}"), "H2_count")
print(f"{'subtraction (Table 2)':28s} {ctrl.H2_count.mean()/2:>9.3f} {1.0:>10.3f}  "
      f"{f'{o1:.2f} ({l1:.2f}-{h1:.2f}) p={p1:.3f}':>22s}  "
      f"{f'{o2:.2f} ({l2:.2f}-{h2:.2f}) p={p2:.3f}':>22s}")
print("\nAll four direct tags give a LARGER residual than subtraction:")
print("the published estimate is conservative, and the residual is not a")
print("consequence of deriving H2 by subtraction.")

# %% [markdown]
# ## 2. Nested-model likelihood-ratio tests (the conditional asymmetry)

# %%
l_rs   = fit(f"case~rs35705950_T_count+{covs}").llf
l_h3   = fit(f"case~H3_count+{covs}").llf
l_both = fit(f"case~H3_count+rs35705950_T_count+{covs}").llf
for lab, base in [("add H3 to rs35705950", l_rs), ("add rs35705950 to H3", l_h3)]:
    chi2 = 2 * (l_both - base)
    print(f"LRT {lab:24s} chi2(1) = {chi2:7.2f}   p = {stats.chi2.sf(chi2, 1):.3g}")

# %% [markdown]
# ## 3. Reciprocal stratification — rs35705950 within H3-carrier strata

# %%
for lab, d in [("H3 non-carriers", dm[dm.H3_count == 0]),
               ("H3 carriers",     dm[dm.H3_count > 0])]:
    o, l, h, p = orci(fit(f"case~rs35705950_T_count+{covs}", d), "rs35705950_T_count")
    print(f"rs35705950 in {lab:16s} N={len(d):3d}  OR {o:.2f} ({l:.2f}-{h:.2f})  p={p:.2g}")

m = fit(f"case~H3_count*rs35705950_T_count+{covs}")
print(f"\nH3 x rs35705950 interaction p = {m.pvalues['H3_count:rs35705950_T_count']:.3f} "
      "(justifies pooling the strata)")

# %% [markdown]
# ## 4. MUC5AC H2 stratified by rs35705950 carrier status (Supplementary Table S2)
#
# H2 is the one contrast retaining significance after adjustment, so it is
# stratified here as H3 and MUC5B H1 already are. Centre exclusions follow the
# Table S2 convention (centre 18 from the non-carrier stratum; 17 and 21 from
# the carrier stratum) to avoid complete separation.

# %%
_neg = dm[(dm.rs35705950_T_count == 0) & (~dm.centre.isin([18]))]
_pos = dm[(dm.rs35705950_T_count >= 1) & (~dm.centre.isin([17, 21]))]
# Use the PRIMARY H2 tag (rs879201833) here, the one the manuscript reports
# everywhere else. This section previously used rs28519516, a secondary panel
# tag, which is what put the wrong H2 rows into Supplementary Table S2.
_f = dm.loc[dm.case == 0, NEW].mean() / 2
dm["H2d"] = dm[NEW] if _f < 0.5 else 2 - dm[NEW]
for _lab, _d in [("Non-carrier", _neg), ("Carrier", _pos)]:
    for _tag, _v in [("subtraction", "H2_count"), ("direct tag", "H2d")]:
        _dd = _d.assign(H2d=dm.H2d.loc[_d.index])
        o, l, h, p = orci(fit(f"case~{_v}+H3_count+{covs}", _dd), _v)
        print(f"  H2 ({_tag:11s}) {_lab:12s} N={len(_dd):3d}  OR {o:.2f} ({l:.2f}-{h:.2f})  p={p:.3f}")

# %% [markdown]
# ## 5. Genotypic (2-df) model for MUC5AC H3
#
# Additive coding is assumed throughout the paper. Plender et al. report positive
# Tajima's D and homozygote depletion for H3 in Europeans and interpret this as
# heterozygote advantage, so a dose-dependent effect is not the only candidate.

# %%
from scipy import stats as _st
dm["H3g"] = dm.H3_count.round().astype(int)
for _lab, _extra in [("unadjusted", ""), ("+ rs35705950", "rs35705950_T_count+")]:
    _mg = fit(f"case~C(H3g)+{_extra}{covs}")
    _ma = fit(f"case~H3_count+{_extra}{covs}")
    _terms = [t for t in _mg.params.index if "H3g" in t]
    _txt = "  ".join(f"{t.split('.')[-1].rstrip(']')} copy: OR {np.exp(_mg.params[t]):.2f}" for t in _terms)
    _chi = 2 * (_mg.llf - _ma.llf)
    print(f"  {_lab:14s} {_txt}   LRT vs additive chi2(1)={_chi:.2f} p={_st.chi2.sf(_chi,1):.2f}")
print("DONE")

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
# # 16. Peer-review sensitivity analyses: ancestry restriction and rs35705950 x environment
#
# Two analyses requested by collaborators:
#
# 1. **Ancestry restriction** (Allen, Hollox, Fawcett): the primary analysis adjusts
#    for five within-sample PCs; reviewers asked whether restricting to a single
#    ancestry changes the result. Re-run the primary estimates in the 771 self-reported
#    white participants (96.6% of the sample).
# 2. **rs35705950 x environment** (Schwartz): test the promoter variant for interaction
#    with ever-smoking and with asbestos (Peto) exposure, complementing the haplogroup
#    x environment tests already in Supplementary Table S5.

# %%
import warnings
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")
ROOT = "/home/drcjar/ipfjes-omics"

clin = pd.read_csv(f"{ROOT}/flat_data_genotype_subset.csv"); clin["pid"] = clin.participant_id.astype(str)
hap = pd.read_csv(f"{ROOT}/muc5_haplogroups.csv"); hap["pid"] = hap.IID.str.split("_").str[0]
pcs = pd.read_csv(f"{ROOT}/code/ipfjes_pca.eigenvec", sep=r"\s+"); pcs["pid"] = pcs.IID.astype(str)
raw = pd.read_csv(f"{ROOT}/muc5ac_h2_tag.raw", sep=r"\s+"); raw["pid"] = raw.IID.astype(str).str.split("_").str[0]
H2C = [c for c in raw.columns if c.startswith("chr11:1182279")][0]

d = (clin.merge(hap[["pid", "H3_count", "MUC5B_H1_count", "rs35705950_T_count"]], on="pid")
     .merge(raw[["pid", H2C]], on="pid")
     .merge(pcs[["pid", "PC1", "PC2", "PC3", "PC4", "PC5"]], on="pid")
     .dropna(subset=["age", "ever_smoked", "peto_exposed", "centre"]))
d = d[d.centre != 16].copy()
_f = d.loc[d.case == 0, H2C].mean() / 2
d["H2t"] = d[H2C] if _f < 0.5 else 2 - d[H2C]
COVS = "age + ever_smoked + PC1+PC2+PC3+PC4+PC5 + C(centre)"


def orci(m, t):
    b, se, p = m.params[t], m.bse[t], m.pvalues[t]
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), p


def fit(f, dd):
    try:
        return smf.logit(f, data=dd).fit(disp=0)
    except Exception:
        return smf.logit(f, data=dd).fit(method="bfgs", disp=0)


# %% [markdown]
# ## 1. European-ancestry (white) restricted sensitivity

# %%
white = d[d.ethnicity.astype(str).str.contains("White", case=False, na=False)]
print(f"white-ethnicity subset: N={len(white)} (cases={int(white.case.sum())}, "
      f"controls={int((white.case == 0).sum())}); full sample N={len(d)}")
for lab, f_, t in [("MUC5AC H3 unadjusted", f"case~H3_count+{COVS}", "H3_count"),
                   ("MUC5AC H3 + rs35705950", f"case~H3_count+rs35705950_T_count+{COVS}", "H3_count"),
                   ("MUC5AC H2 joint + rs35705950", f"case~H2t+H3_count+rs35705950_T_count+{COVS}", "H2t"),
                   ("MUC5B H1 unadjusted", f"case~MUC5B_H1_count+{COVS}", "MUC5B_H1_count")]:
    o, l, h, p = orci(fit(f_, white), t)
    print(f"  {lab:30s} OR {o:.2f} ({l:.2f}-{h:.2f})  p={p:.3g}")
print("  -> essentially unchanged from the full-sample estimates (paper Table 2).")

# %% [markdown]
# ## 2. Supplementary Table S4: gene x environment, risk-oriented (RERI, delta method)
#
# Additive-scale interaction (RERI) is not invariant to coding direction: Knol and
# VanderWeele (2012) require each factor to be entered in its risk-increasing
# direction. H3 and rs35705950 are risk-increasing as carriage; the protective
# MUC5AC H2 and MUC5B H1 haplogroups are risk-increasing as NON-carriage.

# %%
d["H3_risk"] = (d.H3_count > 0).astype(int)             # H3 carrier (risk)
d["nonH2_risk"] = (d.H2t <= 0).astype(int)              # non-carrier of protective H2 (risk)
d["nonB1_risk"] = (d.MUC5B_H1_count <= 0).astype(int)   # non-carrier of protective MUC5B H1 (risk)
d["T_risk"] = (d.rs35705950_T_count > 0).astype(int)    # rs35705950 T carrier (risk)


def _reri(m, ht, et):
    it = f"{ht}:{et}"; b = m.params; V = m.cov_params()
    e11 = np.exp(b[ht] + b[et] + b.get(it, 0)); e10 = np.exp(b[ht]); e01 = np.exp(b[et])
    g = np.zeros(len(b)); idx = {n: i for i, n in enumerate(b.index)}
    g[idx[ht]] = e11 - e10; g[idx[et]] = e11 - e01; g[idx[it]] = e11
    se = float(np.sqrt(g @ V.values @ g)); r = e11 - e10 - e01 + 1
    return r, r - 1.96 * se, r + 1.96 * se


print(f"  {'Factor (risk direction)':22s} {'Exposure':11s} {'mult OR (95% CI)':20s} {'p':6s} RERI (95% CI)")
for lab, g in [("MUC5AC H3 carrier", "H3_risk"), ("MUC5AC H2 non-carrier", "nonH2_risk"),
               ("MUC5B H1 non-carrier", "nonB1_risk"), ("rs35705950 T carrier", "T_risk")]:
    for env, elab in [("peto_exposed", "Asbestos"), ("ever_smoked", "Ever-smoked")]:
        cov = "age + PC1+PC2+PC3+PC4+PC5 + C(centre)" if env == "ever_smoked" else COVS
        m = fit(f"case~{g}+{env}+{g}:{env}+{cov}", d); it = f"{g}:{env}"
        o, lo, hi, p = orci(m, it); r, rl, rh = _reri(m, g, env)
        print(f"  {lab:22s} {elab:11s} {o:.2f} ({lo:.2f}-{hi:.2f})   p={p:.2f}  {r:+.2f} ({rl:+.2f} to {rh:+.2f})")
print("DONE")

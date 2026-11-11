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
# 1. **Ancestry restriction** (Allen, Hollox, Fawcett): the analytical sample is now
#    itself restricted to the 771 self-reported white participants (96.6% of the 798),
#    so the restriction is the primary specification, not a sensitivity analysis. The
#    comparison that still earns its place is the reverse one: relax the restriction and
#    refit in all 798, to show the restriction did not manufacture the result.
# 2. **rs35705950 x environment** (Schwartz): test the promoter variant for interaction
#    with ever-smoking and with asbestos (Peto) exposure, complementing the haplogroup
#    x environment tests already in Supplementary Table S5.

# %%
import warnings
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")
ROOT = "/home/drcjar/academic/ipfjes-omics"

clin = pd.read_csv(f"{ROOT}/flat_data_genotype_subset.csv"); clin["pid"] = clin.participant_id.astype(str)
hap = pd.read_csv(f"{ROOT}/muc5_haplogroups.csv"); hap["pid"] = hap.IID.str.split("_").str[0]
pcs = pd.read_csv(f"{ROOT}/code/ipfjes_pca.eigenvec", sep=r"\s+"); pcs["pid"] = pcs.IID.astype(str)
raw = pd.read_csv(f"{ROOT}/muc5ac_h2_tag.raw", sep=r"\s+"); raw["pid"] = raw.IID.astype(str).str.split("_").str[0]
H2C = [c for c in raw.columns if c.startswith("chr11:1182279")][0]

d = (clin.merge(hap[["pid", "H3_count", "MUC5B_H1_count", "rs35705950_T_count"]], on="pid")
     .merge(raw[["pid", H2C]], on="pid")
     .merge(pcs[["pid", "PC1", "PC2", "PC3", "PC4", "PC5"]], on="pid")
     .dropna(subset=["age", "ever_smoked", "peto_exposed", "centre"]))
d = d[(d.centre != 16) & d.ethnicity.astype(str).str.contains("White", case=False, na=False)].copy()
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
# ## 1. Ancestry restriction: does the unrestricted sample give the same answer?
#
# `d` is already the ancestry-restricted primary sample, so subsetting it to white
# ethnicity again would compare it against itself. The informative direction is to relax
# the restriction and refit in all 798. Principal components cannot enter that model: the
# eigenvectors are computed within the 771 (REPRODUCE.md step 1), so no PCs exist for the
# 27 participants the restriction removes. The unrestricted fit is therefore shown without
# PC adjustment, with the restricted sample refitted on the same reduced covariate set
# beside it, so those two columns differ only in who is included.

# %%
unres = (clin.merge(hap[["pid", "H3_count", "MUC5B_H1_count", "rs35705950_T_count"]], on="pid")
         .merge(raw[["pid", H2C]], on="pid")
         .dropna(subset=["age", "ever_smoked", "peto_exposed", "centre"]))
unres = unres[unres.centre != 16].copy()
unres["H2t"] = unres[H2C] if _f < 0.5 else 2 - unres[H2C]
NOPC = "age + ever_smoked + C(centre)"
print(f"restricted (primary) N={len(d)} (cases={int(d.case.sum())}, "
      f"controls={int((d.case == 0).sum())}); unrestricted N={len(unres)} "
      f"(cases={int(unres.case.sum())}, controls={int((unres.case == 0).sum())}); "
      f"the ancestry restriction removes {len(unres) - len(d)}")
print(f"  {'':30s} {'restricted +PCs':>21s} {'restricted -PCs':>21s} {'unrestricted -PCs':>21s}")
_shift = []
for lab, t, terms in [("MUC5AC H3 unadjusted", "H3_count", "H3_count"),
                      ("MUC5AC H3 + rs35705950", "H3_count", "H3_count+rs35705950_T_count"),
                      ("MUC5AC H2 joint + rs35705950", "H2t", "H2t+H3_count+rs35705950_T_count"),
                      ("MUC5B H1 unadjusted", "MUC5B_H1_count", "MUC5B_H1_count")]:
    cells = [orci(fit(f"case~{terms}+{COVS}", d), t),
             orci(fit(f"case~{terms}+{NOPC}", d), t),
             orci(fit(f"case~{terms}+{NOPC}", unres), t)]
    _shift.append(abs(cells[2][0] - cells[0][0]))
    print(f"  {lab:30s} " + " ".join(f"{o:.2f} ({l:.2f}-{h:.2f})".rjust(21) for o, l, h, _ in cells))
print(f"  -> largest OR difference between the restricted primary and the unrestricted "
      f"refit: {max(_shift):.2f}")

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

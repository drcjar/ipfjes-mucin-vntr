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

"""05_h2ref_sensitivity.py — MUC5AC joint-model reference-category sensitivity.

Reviewer concern: the joint MUC5AC model uses H1 as the reference, and the H1
concordance tSNP is imperfect (r²=0.899) so some H2 may be misclassified as H1.
Re-fit the joint model with H2 (and H3) as the reference to check robustness.

Note: H1+H2+H3 = 2 per person, so these are reparameterisations of one model
(the fit is identical; only the contrasts change). The H3 count is tagged at
r²=1.0, so the H3-vs-non-H3 main effect (model 1a) is the estimate that does
NOT depend on the H1/H2 split at all. Runs on the regression sample N=798.
"""
import pandas as pd, numpy as np, warnings
import statsmodels.formula.api as smf
warnings.filterwarnings("ignore")
R = "."

pcs = pd.read_csv(f"{R}/code/ipfjes_pca.eigenvec", sep=r"\s+"); pcs["pid"] = pcs.IID.astype(str)
clin = pd.read_csv(f"{R}/flat_data_genotype_subset.csv"); clin["pid"] = clin.participant_id.astype(str)
hap = pd.read_csv(f"{R}/muc5_haplogroups.csv"); hap["pid"] = hap.IID.str.split("_").str[0]
df = (clin.merge(hap[["pid","H1_count","H2_count","H3_count","rs35705950_T_count"]], on="pid")
      .merge(pcs[["pid","PC1","PC2","PC3","PC4","PC5"]], on="pid")
      .dropna(subset=["age","ever_smoked","peto_exposed","centre"]))
dm = df[df.centre != 16].copy()
covs = "age + ever_smoked + PC1+PC2+PC3+PC4+PC5 + C(centre)"
fit = lambda f: smf.logit(f, data=dm).fit(disp=0)
def orci(m, t):
    b, se, p = m.params[t], m.bse[t], m.pvalues[t]
    return np.exp(b), np.exp(b-1.96*se), np.exp(b+1.96*se), p

s = (dm.H1_count + dm.H2_count + dm.H3_count)
print(f"N={len(dm)}; H1+H2+H3 per person min/max/mean = {s.min():.2f}/{s.max():.2f}/{s.mean():.3f}\n")

def show(title, formula, terms):
    m = fit(formula)
    print(title)
    for t, lab in terms:
        o, l, h, p = orci(m, t)
        print(f"    {lab:14s} OR {o:.2f} ({l:.2f}–{h:.2f})  p={p:.2e}")
    print()

show("Joint, H1 reference (published):", f"case~H2_count+H3_count+{covs}",
     [("H2_count","H2 vs H1"), ("H3_count","H3 vs H1")])
show("Joint, H2 reference (sensitivity):", f"case~H1_count+H3_count+{covs}",
     [("H1_count","H1 vs H2"), ("H3_count","H3 vs H2")])
show("Joint, H3 reference (sensitivity):", f"case~H1_count+H2_count+{covs}",
     [("H1_count","H1 vs H3"), ("H2_count","H2 vs H3")])
show("Joint, H2 reference + rs35705950:", f"case~H1_count+H3_count+rs35705950_T_count+{covs}",
     [("H1_count","H1 vs H2"), ("H3_count","H3 vs H2"), ("rs35705950_T_count","rs35705950 T")])
show("H3 main effect (model 1a; r²=1.0, independent of H1/H2 split):",
     f"case~H3_count+{covs}", [("H3_count","H3 vs non-H3")])

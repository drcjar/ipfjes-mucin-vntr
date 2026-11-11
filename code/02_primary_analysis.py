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

import pandas as pd, numpy as np, warnings
from scipy import stats
import statsmodels.formula.api as smf
warnings.filterwarnings('ignore')

pcs = pd.read_csv('/home/drcjar/academic/ipfjes-omics/code/ipfjes_pca.eigenvec', sep=r'\s+')
pcs['pid'] = pcs['IID'].astype(str)
clin = pd.read_csv('/home/drcjar/academic/ipfjes-omics/flat_data_genotype_subset.csv')
clin['pid'] = clin['participant_id'].astype(str)
hap  = pd.read_csv('/home/drcjar/academic/ipfjes-omics/muc5_haplogroups.csv')
hap['pid'] = hap['IID'].str.split('_').str[0]

df = (clin
      .merge(hap[['pid','H1_count','H2_count','H3_count','MUC5B_H1_count','MUC5B_H2_count','rs35705950_T_count']], on='pid')
      .merge(pcs[['pid','PC1','PC2','PC3','PC4','PC5']], on='pid')
      .dropna(subset=['age','ever_smoked','peto_exposed','centre']))
# Primary MUC5AC H2 is the direct tag rs879201833 (H2_count is the diplotype complement).
_h2 = pd.read_csv('/home/drcjar/academic/ipfjes-omics/muc5ac_h2_tag.raw', sep=r'\s+')
_h2['pid'] = _h2['IID'].astype(str).str.split('_').str[0]
_h2c = [c for c in _h2.columns if c.startswith('chr11:1182279')][0]
df = df.merge(_h2[['pid', _h2c]], on='pid')
_f = df.loc[df.case == 0, _h2c].mean() / 2
df['H2t'] = df[_h2c] if _f < 0.5 else 2 - df[_h2c]

# Centre 16 has no controls (31/31 cases) — exclude from regression models only
dm = df[(df['centre'] != 16) & df.ethnicity.astype(str).str.contains("White", case=False, na=False)].copy()

cases   = dm[dm['case']==1]        # Table 1 is the N=798 analytical sample (centre 16 excluded)
controls= dm[dm['case']==0]
print(f"Merged N={len(df)} -> analytical N={len(dm)} (centre 16 excluded); cases={len(cases)}, controls={len(controls)}")

# TABLE 1 (N=798)
age_p = stats.mannwhitneyu(cases['age'], controls['age'], alternative='two-sided').pvalue
es_c,es_ct = int(cases['ever_smoked'].sum()), int(controls['ever_smoked'].sum())
_,es_p,_,_ = stats.chi2_contingency([[es_c,len(cases)-es_c],[es_ct,len(controls)-es_ct]])
sc  = cases[cases['ever_smoked']==1]['packyrs'].dropna()
sct = controls[controls['ever_smoked']==1]['packyrs'].dropna()
py_p = stats.mannwhitneyu(sc, sct, alternative='two-sided').pvalue
ab_c,ab_ct = int(cases['peto_exposed'].sum()), int(controls['peto_exposed'].sum())
_,ab_p,_,_ = stats.chi2_contingency([[ab_c,len(cases)-ab_c],[ab_ct,len(controls)-ab_ct]])

print(f"Age: {cases['age'].median():.1f}({cases['age'].quantile(.25):.1f}-{cases['age'].quantile(.75):.1f}) vs {controls['age'].median():.1f}({controls['age'].quantile(.25):.1f}-{controls['age'].quantile(.75):.1f}) p={age_p:.3f}")
print(f"Ever-smoked: {es_c}({100*es_c/len(cases):.1f}%) vs {es_ct}({100*es_ct/len(controls):.1f}%) p={es_p:.3f}")
print(f"Pack-years: {sc.median():.1f}({sc.quantile(.25):.1f}-{sc.quantile(.75):.1f}) vs {sct.median():.1f}({sct.quantile(.25):.1f}-{sct.quantile(.75):.1f}) p={py_p:.3f}")
print(f"Asbestos: {ab_c}({100*ab_c/len(cases):.1f}%) vs {ab_ct}({100*ab_ct/len(controls):.1f}%) p={ab_p:.3f}")
print(f"Peto exposed all: {ab_c+ab_ct}/{len(dm)} = {(ab_c+ab_ct)/len(dm):.1%}")

for col,lab in [('H1_count','MUC5AC H1'),('H2t','MUC5AC H2'),('H3_count','MUC5AC H3'),('MUC5B_H1_count','MUC5B H1')]:
    cn  = int((cases[col]>0).sum());   ctn = int((controls[col]>0).sum())
    _,p,_,_ = stats.chi2_contingency([[cn,len(cases)-cn],[ctn,len(controls)-ctn]])
    print(f"{lab}: {cn}({100*cn/len(cases):.1f}%) vs {ctn}({100*ctn/len(controls):.1f}%) p={p:.3e}")

rs_c  = int((cases['rs35705950_T_count']>0).sum())
rs_ct = int((controls['rs35705950_T_count']>0).sum())
tt_c  = int((cases['rs35705950_T_count']==2).sum())
tt_ct = int((controls['rs35705950_T_count']==2).sum())
_,rs_p,_,_ = stats.chi2_contingency([[rs_c,len(cases)-rs_c],[rs_ct,len(controls)-rs_ct]])
_,tt_p,_,_ = stats.chi2_contingency([[tt_c,len(cases)-tt_c],[tt_ct,len(controls)-tt_ct]])
print(f"T-carrier: {rs_c}({100*rs_c/len(cases):.1f}%) vs {rs_ct}({100*rs_ct/len(controls):.1f}%) p={rs_p:.2e}")
print(f"TT: {tt_c}({100*tt_c/len(cases):.1f}%) vs {tt_ct}({100*tt_ct/len(controls):.1f}%) p={tt_p:.2e}")

# LD (N=798)
r_h3,p_h3 = stats.pearsonr(dm['H3_count'], dm['rs35705950_T_count'])
h3c = (dm['H3_count']>0)
h3T = int((dm.loc[h3c,'rs35705950_T_count']>0).sum())
r_h1  = np.corrcoef(dm['H1_count'], dm['rs35705950_T_count'])[0,1]
r_mh1 = np.corrcoef(dm['MUC5B_H1_count'], dm['rs35705950_T_count'])[0,1]
print(f"r(H3,rs35705950)={r_h3:.3f} p={p_h3:.2e}; H3 with T: {h3T}/{h3c.sum()}={100*h3T/h3c.sum():.1f}%")
print(f"r(H1,rs35705950)={r_h1:.3f}; r(MUC5B_H1,rs35705950)={r_mh1:.3f}")

# MODELS (regression N=798, excluding centre 16)
covs = 'age + ever_smoked + PC1+PC2+PC3+PC4+PC5 + C(centre)'
def fit(f, d=dm):
    return smf.logit(f, data=d).fit(disp=0)
def fmt(m, t):
    b,se,p = m.params[t], m.bse[t], m.pvalues[t]
    return np.exp(b), np.exp(b-1.96*se), np.exp(b+1.96*se), p

m1a = fit(f'case~H3_count+{covs}')
o,l,h,p = fmt(m1a,'H3_count')
print(f"1a H3: OR {o:.2f}({l:.2f}-{h:.2f}) p={p:.2e}")

m1b = fit(f'case~H3_count+rs35705950_T_count+{covs}')
o,l,h,p = fmt(m1b,'H3_count'); print(f"1b H3+rs35705950 | H3: OR {o:.2f}({l:.2f}-{h:.2f}) p={p:.2e}")
o,l,h,p = fmt(m1b,'rs35705950_T_count'); print(f"1b H3+rs35705950 | rs35705950: OR {o:.2f}({l:.2f}-{h:.2f}) p={p:.2e}")

m1c = fit(f'case~H2t+H3_count+{covs}')
o,l,h,p = fmt(m1c,'H2t'); print(f"1c joint (H1 ref) | H2: OR {o:.2f}({l:.2f}-{h:.2f}) p={p:.2e}")
o,l,h,p = fmt(m1c,'H3_count'); print(f"1c joint (H1 ref) | H3: OR {o:.2f}({l:.2f}-{h:.2f}) p={p:.2e}")

m4a = fit(f'case~MUC5B_H1_count+{covs}')
o,l,h,p = fmt(m4a,'MUC5B_H1_count')
print(f"4a MUC5B H1: OR {o:.2f}({l:.2f}-{h:.2f}) p={p:.4f}")

m4b = fit(f'case~MUC5B_H1_count+rs35705950_T_count+{covs}')
o,l,h,p = fmt(m4b,'MUC5B_H1_count'); print(f"4b MUC5B H1+rs | H1: OR {o:.2f}({l:.2f}-{h:.2f}) p={p:.2e}")
o,l,h,p = fmt(m4b,'rs35705950_T_count'); print(f"4b MUC5B H1+rs | rs35705950: OR {o:.2f}({l:.2f}-{h:.2f}) p={p:.2e}")

mrs = fit(f'case~rs35705950_T_count+{covs}')
o,l,h,p = fmt(mrs,'rs35705950_T_count')
print(f"rs35705950 alone: OR {o:.2f}({l:.2f}-{h:.2f}) p={p:.2e}")

# Gene x environment interactions / RERI are reported in Supplementary Table S4, computed in
# 14_ancestry_and_promoter_gxe.py (risk-oriented, carrier-coded per Knol & VanderWeele) and
# cross-checked in 03_supplementary_analyses.py. Not duplicated here: the earlier additive
# 0/1/2-count RERI parameterisation was superseded by the carrier-level one.
print("DONE")

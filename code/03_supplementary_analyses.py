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

"""03_supplementary_analyses.py — supplementary analyses for MUC5 VNTR paper.

Runs in sequence after 01_haplogroup_assignment.py and 02_primary_analysis.py.
Requires: muc5_haplogroups.csv, flat_data_genotype_subset.csv, ipfjes_pca.eigenvec

Tasks:
  1. rs35705950-stratified haplogroup associations (Supplementary Table S2)
  2. Joint MUC5AC model + rs35705950 adjustment; attenuation quantification
  3. Bootstrap RERI 95% CIs (500 resamples, seed=42)
  4. Hardy-Weinberg equilibrium in controls
  5. Post-hoc interaction power simulations (2,000 replicates)
"""
import pandas as pd, numpy as np, warnings
from scipy import stats
import statsmodels.formula.api as smf
warnings.filterwarnings('ignore')

# Set QUICK=True for a fast (~few-second) preview run with reduced resampling.
# The PUBLISHED Supplementary Table S5 RERI CIs and power figures need QUICK=False
# (full counts below); a QUICK=True run prints PREVIEW-ONLY numbers.
QUICK  = True
N_BOOT = 20 if QUICK else 500   # Task 3 — bootstrap RERI resamples
N_SIM  = 50 if QUICK else 2000  # Task 5 — power-simulation replicates
if QUICK:
    print("*** QUICK PREVIEW MODE (N_BOOT=20, N_SIM=50): numbers below are NOT the "
          "published S5/power values — set QUICK=False for those. ***\n")

# ── Load data (same merge as 02_primary_analysis.py) ─────────────────────
pcs  = pd.read_csv('/home/drcjar/academic/ipfjes-omics/code/ipfjes_pca.eigenvec', sep=r'\s+')
pcs['pid'] = pcs['IID'].astype(str)
clin = pd.read_csv('/home/drcjar/academic/ipfjes-omics/flat_data_genotype_subset.csv')
clin['pid'] = clin['participant_id'].astype(str)
hap  = pd.read_csv('/home/drcjar/academic/ipfjes-omics/muc5_haplogroups.csv')
hap['pid'] = hap['IID'].str.split('_').str[0]

df = (clin
      .merge(hap[['pid','H1_count','H2_count','H3_count','MUC5B_H1_count','MUC5B_H2_count','rs35705950_T_count']], on='pid')
      .merge(pcs[['pid','PC1','PC2','PC3','PC4','PC5']], on='pid')
      .dropna(subset=['age','ever_smoked','peto_exposed','centre']))
# Primary MUC5AC H2 is the direct tag rs879201833 (H2_count is the diplotype complement,
# kept for the descriptive diplotype and the tag-sensitivity table).
_h2 = pd.read_csv('/home/drcjar/academic/ipfjes-omics/muc5ac_h2_tag.raw', sep=r'\s+')
_h2['pid'] = _h2['IID'].astype(str).str.split('_').str[0]
_h2c = [c for c in _h2.columns if c.startswith('chr11:1182279')][0]
df = df.merge(_h2[['pid', _h2c]], on='pid')
_f = df.loc[df.case == 0, _h2c].mean() / 2
df['H2t'] = df[_h2c] if _f < 0.5 else 2 - df[_h2c]
dm = df[(df['centre'] != 16) & df.ethnicity.astype(str).str.contains("White", case=False, na=False)].copy()   # regression dataset (N=771, European-ancestry restricted)

covs = 'age + ever_smoked + PC1+PC2+PC3+PC4+PC5 + C(centre)'

def fit(f, d=dm):
    try:    return smf.logit(f, data=d).fit(disp=0)
    except: return smf.logit(f, data=d).fit(method='bfgs', disp=0)

def fmt(m, t):
    b,se,p = m.params[t], m.bse[t], m.pvalues[t]
    return np.exp(b), np.exp(b-1.96*se), np.exp(b+1.96*se), p

print('='*70)
print('TASK 1 — rs35705950-Stratified Analysis')
print('='*70)

# Centre 18 all-cases in non-carrier stratum; centres 17+21 all-cases in carrier stratum
dm_neg = dm[(dm['rs35705950_T_count'] == 0) & (~dm['centre'].isin([18]))].copy()
dm_pos = dm[(dm['rs35705950_T_count'] >= 1) & (~dm['centre'].isin([17, 21]))].copy()
print(f'Non-carriers (T_count=0, excl c18): N={len(dm_neg)}, cases={int(dm_neg.case.sum())}, controls={int((dm_neg.case==0).sum())}')
print(f'Carriers (T_count>=1, excl c17,21): N={len(dm_pos)}, cases={int(dm_pos.case.sum())}, controls={int((dm_pos.case==0).sum())}')

results_s1 = []
for label, d in [('Non-carrier', dm_neg), ('Carrier', dm_pos), ('Overall', dm)]:
    for hc, hl in [('H3_count','MUC5AC H3'), ('MUC5B_H1_count','MUC5B H1')]:
        try:
            m = fit(f'case~{hc}+{covs}', d)
            o,l,h,p = fmt(m, hc)
            print(f"{hl} | {label} | N={len(d)} | OR {o:.2f}({l:.2f}-{h:.2f}) p={p:.3f}")
            results_s1.append({'Haplogroup': hl, 'Stratum': label, 'N': len(d),
                                'OR': f'{o:.2f}', 'CI': f'{l:.2f}–{h:.2f}', 'p': f'{p:.3f}'})
        except Exception as e:
            print(f"{hl} | {label} | FAILED: {e}")
            results_s1.append({'Haplogroup': hl, 'Stratum': label, 'N': len(d),
                                'OR': 'NC', 'CI': '—', 'p': '—'})

print()
print('='*70)
print('TASK 2 — rs35705950-Adjusted Joint MUC5AC Model')
print('='*70)

m_joint_adj = fit(f'case~H2t+H3_count+rs35705950_T_count+{covs}')
o2,l2,h2,p2 = fmt(m_joint_adj, 'H2t')
o3,l3,h3,p3 = fmt(m_joint_adj, 'H3_count')
or_,l_,h_,p_ = fmt(m_joint_adj, 'rs35705950_T_count')
print(f"Joint+rs35705950 (H1 ref) | H2: OR {o2:.2f}({l2:.2f}-{h2:.2f}) p={p2:.2e}")
print(f"Joint+rs35705950 (H1 ref) | H3: OR {o3:.2f}({l3:.2f}-{h3:.2f}) p={p3:.2e}")
print(f"Joint+rs35705950 (H1 ref) | rs35705950: OR {or_:.2f}({l_:.2f}-{h_:.2f}) p={p_:.2e}")

# Attenuation on log-OR scale
import math
m_joint = fit(f'case~H2t+H3_count+{covs}')
o2u,_,_,_ = fmt(m_joint,'H2t'); o3u,_,_,_ = fmt(m_joint,'H3_count')
atten_h2 = 100*(1 - math.log(o2)/math.log(o2u))
atten_h3 = 100*(1 - math.log(o3)/math.log(o3u))
print(f"H2 attenuation on log-OR: {atten_h2:.0f}%  ({math.log(o2u):.3f} -> {math.log(o2):.3f})")
print(f"H3 attenuation on log-OR: {atten_h3:.0f}%  ({math.log(o3u):.3f} -> {math.log(o3):.3f})")

# Also compute H3 attenuation (main model 1a vs 1b)
m1a = fit(f'case~H3_count+{covs}')
m1b = fit(f'case~H3_count+rs35705950_T_count+{covs}')
o1a,_,_,_ = fmt(m1a,'H3_count'); o1b,_,_,_ = fmt(m1b,'H3_count')
atten_main = 100*(1 - math.log(o1b)/math.log(o1a))
m4a = fit(f'case~MUC5B_H1_count+{covs}')
m4b = fit(f'case~MUC5B_H1_count+rs35705950_T_count+{covs}')
o4a,_,_,_ = fmt(m4a,'MUC5B_H1_count'); o4b,_,_,_ = fmt(m4b,'MUC5B_H1_count')
atten_mh1 = 100*(1 - math.log(o4b)/math.log(o4a))
print(f"H3 main effect attenuation (1a->1b): {atten_main:.0f}%")
print(f"MUC5B H1 attenuation (4a->4b): {atten_mh1:.0f}%")

print()
print('='*70)
print('TASK 3 — Bootstrap CIs for RERI')
print('='*70)

def reri_point(m, ht, et):
    it = f'{ht}:{et}'
    b = m.params
    lor11 = b[ht]+b[et]+b.get(it,0)
    return np.exp(lor11)-np.exp(b[ht])-np.exp(b[et])+1

# RERI is defined for BINARY exposures. Earlier versions passed the additive
# 0/1/2 haplogroup count, so b[ht] was a one-copy-versus-zero log-OR and the
# homozygote stratum was ignored — not the carrier-level RERI a reader assumes.
# Both terms are now carrier indicators (Knol & VanderWeele 2012).
dm['H3_carrier']      = (dm.H3_count > 0).astype(int)
dm['MUC5B_H1_carrier'] = (dm.MUC5B_H1_count > 0).astype(int)

results_boot = []
for hc,ec,hl,el in [('H3_carrier','peto_exposed','H3','Asbestos'),
                     ('H3_carrier','ever_smoked','H3','Smoking'),
                     ('MUC5B_H1_carrier','peto_exposed','MUC5B_H1','Asbestos'),
                     ('MUC5B_H1_carrier','ever_smoked','MUC5B_H1','Smoking')]:
    it = f'{hc}:{ec}'
    mi = fit(f'case~{hc}+{ec}+{hc}:{ec}+{covs}')
    rp = reri_point(mi, hc, ec)

    boot = []
    rng = np.random.default_rng(42)
    for _ in range(N_BOOT):
        idx = rng.integers(0, len(dm), len(dm))
        bd  = dm.iloc[idx].copy()
        try:
            mb = smf.logit(f'case~{hc}+{ec}+{hc}:{ec}+{covs}', data=bd).fit(disp=0)
            boot.append(reri_point(mb, hc, ec))
        except: pass

    boot = np.array(boot)
    lo,hi = np.percentile(boot, [2.5,97.5])
    results_boot.append((hl,el,rp,lo,hi))
    print(f"{hl}×{el}: RERI {rp:+.3f}  Bootstrap 95% CI [{lo:.3f}, {hi:.3f}]  (n_boot={len(boot)})")

print()
print('='*70)
print('TASK 3b — Supplementary Table S5 as published (delta-method CIs)')
print('='*70)
# The published S5 uses delta-method CIs, not the percentile bootstrap above.
# Bootstrap replicates fail silently for ~10% of resamples (the separated fits),
# so the surviving percentile interval is biased inward; the delta method is
# reported instead. Both haplogroup and exposure are carrier indicators.
def reri_delta_ci(m, ht, et):
    it = f'{ht}:{et}'; b = m.params; V = m.cov_params()
    e11 = np.exp(b[ht]+b[et]+b.get(it,0)); e10 = np.exp(b[ht]); e01 = np.exp(b[et])
    g = np.zeros(len(b)); idx = {n:i for i,n in enumerate(b.index)}
    g[idx[ht]] = e11-e10; g[idx[et]] = e11-e01; g[idx[it]] = e11
    se = float(np.sqrt(g @ V.values @ g)); r = e11-e10-e01+1
    return r, r-1.96*se, r+1.96*se

print(f"{'Haplogroup':12s} {'Exposure':12s} {'multiplicative OR (95% CI)':28s} {'p':>6s}  RERI (95% CI)")
for hc, ec, hl, el in [('H3_carrier','peto_exposed','MUC5AC H3','Asbestos'),
                       ('H3_carrier','ever_smoked','MUC5AC H3','Ever-smoked'),
                       ('MUC5B_H1_carrier','peto_exposed','MUC5B H1','Asbestos'),
                       ('MUC5B_H1_carrier','ever_smoked','MUC5B H1','Ever-smoked')]:
    cv = covs if ec != 'ever_smoked' else 'age + PC1+PC2+PC3+PC4+PC5 + C(centre)'
    m = fit(f'case~{hc}+{ec}+{hc}:{ec}+{cv}')
    b, se, pv = m.params[f'{hc}:{ec}'], m.bse[f'{hc}:{ec}'], m.pvalues[f'{hc}:{ec}']
    r, lo, hi = reri_delta_ci(m, hc, ec)
    print(f"{hl:12s} {el:12s} "
          f"{f'{np.exp(b):.2f} ({np.exp(b-1.96*se):.2f}-{np.exp(b+1.96*se):.2f})':28s} "
          f"{pv:>6.2f}  {r:+.2f} ({lo:+.2f} to {hi:+.2f})")

print()
print('='*70)
print('TASK 4 — Hardy-Weinberg Equilibrium in Controls')
print('='*70)

controls = dm[dm['case']==0]
from scipy.stats import chi2 as chi2_dist

def hwe_test(col, ctrldf):
    n = len(ctrldf)
    cnt = ctrldf[col]
    hom   = int((cnt==2).sum())
    het   = int((cnt==1).sum())
    wt    = int((cnt==0).sum())
    p_hat = cnt.mean() / 2   # haplotype frequency
    q_hat = 1 - p_hat
    exp_hom = n * p_hat**2
    exp_het = n * 2*p_hat*q_hat
    exp_wt  = n * q_hat**2
    obs = np.array([wt, het, hom], dtype=float)
    exp = np.array([exp_wt, exp_het, exp_hom])
    stat = np.sum((obs-exp)**2 / exp)
    p = 1 - chi2_dist.cdf(stat, df=1)
    return stat, p, obs.astype(int), exp.round(1), p_hat

for col, label in [('H3_count','MUC5AC H3'), ('H1_count','MUC5AC H1'), ('MUC5B_H1_count','MUC5B H1')]:
    stat, p, obs, exp, freq = hwe_test(col, controls)
    print(f"{label}: haplotype freq={freq:.3f}  obs[WT,het,hom]={obs}  exp={exp}  chi2={stat:.3f} p={p:.3f}")

print()
print('='*70)
print('TASK 5 — Interaction Power Calculations')
print('='*70)

def interaction_power_sim(n_cases, n_controls, p_hap, p_exp,
                          or_hap, or_exp, or_int, alpha=0.05, n_sim=2000, seed=42):
    rng = np.random.default_rng(seed)
    rejected = 0
    for _ in range(n_sim):
        n = n_cases + n_controls
        hap = rng.binomial(2, p_hap, n)   # additive 0/1/2 dosage under HWE (p_hap = allele freq), matching the primary coding
        exp = rng.binomial(1, p_exp, n)
        lo  = (np.log(or_hap)*hap + np.log(or_exp)*exp + np.log(or_int)*hap*exp)
        prob= 1/(1+np.exp(-lo))
        y   = rng.binomial(1, prob, n)
        d   = pd.DataFrame({'y':y,'hap':hap,'exp':exp})
        try:
            m = smf.logit('y~hap+exp+hap:exp', data=d).fit(disp=0)
            if m.pvalues.get('hap:exp',1.0) < alpha: rejected+=1
        except: pass
    return rejected/n_sim

# Parameters from N=798 analysis
n_cases, n_controls = 397, 401
p_h3 = dm['H3_count'].mean()/2   # haplotype freq
p_h1 = dm['MUC5B_H1_count'].mean()/2
p_asb= dm['peto_exposed'].mean()
p_smk= dm['ever_smoked'].mean()

print(f'Parameters: H3 hap freq={p_h3:.3f}, MUC5B H1 hap freq={p_h1:.3f}')
print(f'            asbestos {p_asb:.3f}, ever-smoked {p_smk:.3f}')
print(f'Running power simulations (n_sim={N_SIM} each)...')

scenarios = [
    ('MUC5AC H3 × Asbestos (OR_int=1.5)', p_h3, p_asb, 1.69, 1.3, 1.5),
    ('MUC5AC H3 × Asbestos (OR_int=2.0)', p_h3, p_asb, 1.69, 1.3, 2.0),
    ('MUC5B H1 × Asbestos (OR_int=0.5)',  p_h1, p_asb, 0.72, 1.3, 0.5),
    ('MUC5B H1 × Asbestos (OR_int=0.3)',  p_h1, p_asb, 0.72, 1.3, 0.3),
    ('MUC5AC H3 × Smoking (OR_int=1.5)',  p_h3, p_smk, 1.69, 1.5, 1.5),
    ('MUC5B H1 × Smoking (OR_int=0.5)',   p_h1, p_smk, 0.72, 1.5, 0.5),
]

power_results = []
for label, ph, pe, oh, oe, oi in scenarios:
    pw = interaction_power_sim(n_cases, n_controls, ph, pe, oh, oe, oi, n_sim=N_SIM)
    print(f"  {label}: power={pw:.1%}")
    power_results.append((label, pw))

print()
print('DONE — all analyses complete')

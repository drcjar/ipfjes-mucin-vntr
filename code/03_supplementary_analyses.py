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
  1. rs35705950-stratified haplogroup associations (Supplementary Table S1)
  2. Joint MUC5AC model + rs35705950 adjustment; attenuation quantification
  3. TaqMan sensitivity: H3 re-modelled with TaqMan-derived rs35705950 genotype
  4. Bootstrap RERI 95% CIs (500 resamples, seed=42)
  5. Hardy-Weinberg equilibrium in controls
  6. Post-hoc interaction power simulations (2,000 replicates)
"""
import pandas as pd, numpy as np, warnings
from scipy import stats
import statsmodels.formula.api as smf
warnings.filterwarnings('ignore')

# Set QUICK=True for a fast (~few-second) preview run with reduced resampling.
# The PUBLISHED Supplementary Table S5 RERI CIs and power figures need QUICK=False
# (full counts below); a QUICK=True run prints PREVIEW-ONLY numbers.
QUICK  = True
N_BOOT = 20 if QUICK else 500   # Task 4 — bootstrap RERI resamples
N_SIM  = 50 if QUICK else 2000  # Task 6 — power-simulation replicates
if QUICK:
    print("*** QUICK PREVIEW MODE (N_BOOT=50, N_SIM=200): numbers below are NOT the "
          "published S5/power values — set QUICK=False for those. ***\n")

# ── Load data (same merge as 02_primary_analysis.py) ─────────────────────
pcs  = pd.read_csv('./code/ipfjes_pca.eigenvec', sep=r'\s+')
pcs['pid'] = pcs['IID'].astype(str)
clin = pd.read_csv('./flat_data_genotype_subset.csv')
clin['pid'] = clin['participant_id'].astype(str)
hap  = pd.read_csv('./muc5_haplogroups.csv')
hap['pid'] = hap['IID'].str.split('_').str[0]

df = (clin
      .merge(hap[['pid','H1_count','H2_count','H3_count','MUC5B_H1_count','MUC5B_H2_count','rs35705950_T_count']], on='pid')
      .merge(pcs[['pid','PC1','PC2','PC3','PC4','PC5']], on='pid')
      .dropna(subset=['age','ever_smoked','peto_exposed','centre']))
df['tq_count'] = df['genotype']
df['tq_bin']   = (df['genotype'] > 0).astype(int)
dm = df[df['centre'] != 16].copy()   # regression dataset (N=798)

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

m_joint_adj = fit(f'case~H2_count+H3_count+rs35705950_T_count+{covs}')
o2,l2,h2,p2 = fmt(m_joint_adj, 'H2_count')
o3,l3,h3,p3 = fmt(m_joint_adj, 'H3_count')
or_,l_,h_,p_ = fmt(m_joint_adj, 'rs35705950_T_count')
print(f"Joint+rs35705950 (H1 ref) | H2: OR {o2:.2f}({l2:.2f}-{h2:.2f}) p={p2:.2e}")
print(f"Joint+rs35705950 (H1 ref) | H3: OR {o3:.2f}({l3:.2f}-{h3:.2f}) p={p3:.2e}")
print(f"Joint+rs35705950 (H1 ref) | rs35705950: OR {or_:.2f}({l_:.2f}-{h_:.2f}) p={p_:.2e}")

# Attenuation on log-OR scale
import math
m_joint = fit(f'case~H2_count+H3_count+{covs}')
o2u,_,_,_ = fmt(m_joint,'H2_count'); o3u,_,_,_ = fmt(m_joint,'H3_count')
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
print('TASK 3 — TaqMan Sensitivity Analysis')
print('='*70)
print(f'All N={len(dm)} in regression dataset have TaqMan genotype')
print(f'tq_count range: {dm.tq_count.min():.0f}-{dm.tq_count.max():.0f}')

m_base    = m1a  # H3 unadjusted, N=798
m_tq_bin  = fit(f'case~H3_count+tq_bin+{covs}')
m_tq_cnt  = fit(f'case~H3_count+tq_count+{covs}')

o_b,l_b,h_b,p_b = fmt(m_base,'H3_count')
o_t1,l_t1,h_t1,p_t1 = fmt(m_tq_bin,'H3_count');  r1_o,r1_l,r1_h,r1_p = fmt(m_tq_bin,'tq_bin')
o_t2,l_t2,h_t2,p_t2 = fmt(m_tq_cnt,'H3_count');  r2_o,r2_l,r2_h,r2_p = fmt(m_tq_cnt,'tq_count')

print(f"H3 unadjusted (N=798):                  H3 OR {o_b:.2f}({l_b:.2f}-{h_b:.2f}) p={p_b:.2e}")
print(f"H3 + TaqMan binary (N=798):             H3 OR {o_t1:.2f}({l_t1:.2f}-{h_t1:.2f}) p={p_t1:.2e} | tq_bin OR {r1_o:.2f}({r1_l:.2f}-{r1_h:.2f}) p={r1_p:.2e}")
print(f"H3 + TaqMan additive (N=798):           H3 OR {o_t2:.2f}({l_t2:.2f}-{h_t2:.2f}) p={p_t2:.2e} | tq_cnt OR {r2_o:.2f}({r2_l:.2f}-{r2_h:.2f}) p={r2_p:.2e}")
# Compare with imputed adjustment (model 1b)
o_imp,l_imp,h_imp,p_imp = fmt(m1b,'H3_count')
print(f"H3 + imputed dosage adj (N=798):        H3 OR {o_imp:.2f}({l_imp:.2f}-{h_imp:.2f}) p={p_imp:.2e}")

print()
print('='*70)
print('TASK 4 — Bootstrap CIs for RERI')
print('='*70)

def reri_point(m, ht, et):
    it = f'{ht}:{et}'
    b = m.params
    lor11 = b[ht]+b[et]+b.get(it,0)
    return np.exp(lor11)-np.exp(b[ht])-np.exp(b[et])+1

results_boot = []
for hc,ec,hl,el in [('H3_count','peto_exposed','H3','Asbestos'),
                     ('H3_count','ever_smoked','H3','Smoking'),
                     ('MUC5B_H1_count','peto_exposed','MUC5B_H1','Asbestos'),
                     ('MUC5B_H1_count','ever_smoked','MUC5B_H1','Smoking')]:
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
print('TASK 5 — Hardy-Weinberg Equilibrium in Controls')
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
print('TASK 6 — Interaction Power Calculations')
print('='*70)

def interaction_power_sim(n_cases, n_controls, p_hap, p_exp,
                          or_hap, or_exp, or_int, alpha=0.05, n_sim=2000, seed=42):
    rng = np.random.default_rng(seed)
    rejected = 0
    for _ in range(n_sim):
        n = n_cases + n_controls
        hap = rng.binomial(1, p_hap, n)
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

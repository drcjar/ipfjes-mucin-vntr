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
# # Verification notebook — every manuscript number, recomputed from source
#
# This notebook re-derives each figure quoted in
# `../MUC5_VNTR_haplogroup_paper_draft.md` directly from the source data and
# compares it against the value printed in the manuscript. Any row flagged
# `FAIL` needs attention. The manuscript values were produced by the `02`/`03`
# scripts' formatted output, so a correct re-run should reproduce them exactly.
#
# Run top-to-bottom in a fresh kernel. Heavy resampling/simulation cells
# (bootstrap CIs, power) are gated behind `RUN_SLOW` (default off).

# %%
import math, warnings
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
warnings.filterwarnings("ignore")
RUN_SLOW = False  # set True to also verify bootstrap RERI CIs and power simulations

ROOT = "."

# Collects every check so we can summarise pass/fail at the end.
_checks = []
def record(section, item, manuscript, computed, ok):
    _checks.append(dict(section=section, item=item,
                        manuscript=manuscript, computed=computed,
                        result="PASS" if ok else "FAIL"))

def chk_num(section, item, manuscript, computed, dp=2):
    """Match when computed rounds to the manuscript value at `dp` decimals."""
    ok = round(float(computed), dp) == round(float(manuscript), dp)
    record(section, item, manuscript, round(float(computed), dp), ok)

def chk_int(section, item, manuscript, computed):
    ok = int(computed) == int(manuscript)
    record(section, item, manuscript, int(computed), ok)

def chk_p(section, item, manuscript, computed):
    """p-values: match on log scale (same order + mantissa)."""
    c, m = float(computed), float(manuscript)
    ok = (c > 0 and m > 0 and abs(math.log10(c) - math.log10(m)) < 0.12)
    record(section, item, f"{m:.2e}", f"{c:.2e}", ok)

def chk_tol(section, item, manuscript, computed, tol):
    """Absolute-tolerance match, for approximate ('~') values and rounding boundaries."""
    ok = abs(float(computed) - float(manuscript)) <= tol
    record(section, item, manuscript, round(float(computed), 2), ok)

# %% [markdown]
# ## 0. Load and merge (identical to `02_primary_analysis.py`)

# %%
pcs = pd.read_csv(f"{ROOT}/code/ipfjes_pca.eigenvec", sep=r"\s+")
pcs["pid"] = pcs["IID"].astype(str)
clin = pd.read_csv(f"{ROOT}/flat_data_genotype_subset.csv")
clin["pid"] = clin["participant_id"].astype(str)
hap = pd.read_csv(f"{ROOT}/muc5_haplogroups.csv")
hap["pid"] = hap["IID"].str.split("_").str[0]

df = (clin
      .merge(hap[["pid", "H1_count", "H2_count", "H3_count",
                  "MUC5B_H1_count", "MUC5B_H2_count", "rs35705950_T_count"]], on="pid")
      .merge(pcs[["pid", "PC1", "PC2", "PC3", "PC4", "PC5"]], on="pid")
      .dropna(subset=["age", "ever_smoked", "peto_exposed", "centre"]))
df["tq_count"] = df["genotype"]
df["tq_bin"] = (df["genotype"] > 0).astype(int)
dm = df[df["centre"] != 16].copy()           # regression dataset
cases, controls = df[df.case == 1], df[df.case == 0]

chk_int("Cohort", "Descriptive N", 829, len(df))
chk_int("Cohort", "Cases", 428, len(cases))
chk_int("Cohort", "Controls", 401, len(controls))
chk_int("Cohort", "Regression N (excl. centre 16)", 798, len(dm))
print(f"Descriptive N={len(df)} (cases={len(cases)}, controls={len(controls)}); regression N={len(dm)}")

# %%
covs = "age + ever_smoked + PC1+PC2+PC3+PC4+PC5 + C(centre)"
def fit(f, d=dm):
    try:    return smf.logit(f, data=d).fit(disp=0)
    except: return smf.logit(f, data=d).fit(method="bfgs", disp=0)
def orci(m, t):
    b, se, p = m.params[t], m.bse[t], m.pvalues[t]
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), p

# %% [markdown]
# ## 1. Table 1 — participant characteristics

# %%
age_p = stats.mannwhitneyu(cases.age, controls.age).pvalue
chk_num("Table 1", "Age p", 0.004, age_p, dp=3)

es_c, es_ct = int(cases.ever_smoked.sum()), int(controls.ever_smoked.sum())
chk_int("Table 1", "Ever-smoked cases n", 324, es_c)
chk_int("Table 1", "Ever-smoked controls n", 276, es_ct)
chk_num("Table 1", "Ever-smoked cases %", 75.7, 100*es_c/len(cases), dp=1)
chk_num("Table 1", "Ever-smoked controls %", 68.8, 100*es_ct/len(controls), dp=1)

ab_c, ab_ct = int(cases.peto_exposed.sum()), int(controls.peto_exposed.sum())
chk_num("Table 1", "Asbestos cases %", 66.8, 100*ab_c/len(cases), dp=1)
chk_num("Table 1", "Asbestos controls %", 63.1, 100*ab_ct/len(controls), dp=1)

for col, lab, mc, mct in [("H1_count","MUC5AC H1",341,327), ("H2_count","MUC5AC H2",114,167),
                          ("H3_count","MUC5AC H3",221,141), ("MUC5B_H1_count","MUC5B H1",88,109)]:
    chk_int("Table 1", f"{lab} carrier cases n", mc, int((cases[col] > 0).sum()))
    chk_int("Table 1", f"{lab} carrier controls n", mct, int((controls[col] > 0).sum()))

rs_c  = int((cases.rs35705950_T_count > 0).sum())
rs_ct = int((controls.rs35705950_T_count > 0).sum())
chk_num("Table 1", "rs35705950 T-carrier cases %", 58.2, 100*rs_c/len(cases), dp=1)
chk_num("Table 1", "rs35705950 T-carrier controls %", 22.9, 100*rs_ct/len(controls), dp=1)
chk_int("Table 1", "TT cases n", 27, int((cases.rs35705950_T_count == 2).sum()))
chk_int("Table 1", "TT controls n", 4, int((controls.rs35705950_T_count == 2).sum()))

# %% [markdown]
# ## 2. Haplogroup frequencies and LD with rs35705950

# %%
# Frequencies are quoted approximately, on the FULL haplogroups file (N=2912),
# not the 829-participant analysis sample — compare with a ~1 percentage-point tolerance.
chk_tol("Freqs", "MUC5AC H1 freq % (full cohort)", 58, 100*hap.H1_count.mean()/2, tol=1.0)
chk_tol("Freqs", "MUC5AC H3 freq % (full cohort)", 23, 100*hap.H3_count.mean()/2, tol=1.0)
chk_tol("Freqs", "MUC5B H1 freq % (full cohort)", 12.5, 100*hap.MUC5B_H1_count.mean()/2, tol=1.0)

# Discordant subset among IPFJES participants (N=829) — the long-read sequencing target.
# NB: the full genotyped file (N=2912) is mostly external controls, so discordance is
# computed on the IPFJES analytical sample (df), not the full file.
hd = (df.H3_count > 0) != (df.rs35705950_T_count > 0)
chk_int("Discordant", "H3/rs35705950 discordant n (IPFJES, N=829)", 265, int(hd.sum()))
chk_num("Discordant", "H3/rs35705950 discordant %", 32, 100*hd.mean(), dp=0)
chk_int("Discordant", "discordant cases (n)",    138, int(hd[df.case==1].sum()))
chk_int("Discordant", "discordant controls (n)", 127, int(hd[df.case==0].sum()))
hd_mb = (df.MUC5B_H1_count > 0) != (df.rs35705950_T_count > 0)
chk_int("Discordant", "MUC5B H1/rs35705950 discordant n", 408, int(hd_mb.sum()))
chk_num("Discordant", "MUC5B H1/rs35705950 discordant %", 49, 100*hd_mb.mean(), dp=0)

# Supplementary Table S4 — per-tSNP counted-allele frequencies (covers code/06_tsnp_table.py)
_dos = pd.read_csv(f"{ROOT}/muc5_tsnps_dosage_pheno.csv")
_mb  = pd.read_csv(f"{ROOT}/muc5b_h1h2_tsnp.raw", sep=r"\s+")
s4_freq = {
    "chr11:1165501:T:A_T": 0.769, "chr11:1195214:G:C_C": 0.573, "chr11:1193830:G:C_G": 0.579,
    "chr11:1194354:A:G_A": 0.579, "chr11:1195265:G:C_C": 0.578, "chr11:1195858:T:C_T": 0.580,
    "chr11:1182048:A:T_T": 0.579, "chr11:1202343:T:C_C": 0.808, "chr11:1202859:A:G_G": 0.809,
    "chr11:1202914:T:C_C": 0.808, "chr11:1178216:C:G_C": 0.774, "chr11:1178023:A:G_A": 0.769,
    "chr11:1177688:A:T_A": 0.769, "chr11:1177356:T:C_T": 0.769, "chr11:1177269:T:C_T": 0.769,
    "chr11:1197482:C:T_C": 0.771, "chr11:1219991:G:T_G": 0.849, "chr11:1244757:C:G_C": 0.875,
}
for _cc, _exp in s4_freq.items():
    _src = _dos if _cc in _dos.columns else _mb
    chk_num("Table S4", f"{_cc} freq", _exp, _src[_cc].mean()/2, dp=3)

r_h3, p_h3 = stats.pearsonr(df.H3_count, df.rs35705950_T_count)
chk_num("LD", "r(H3, rs35705950)", 0.339, r_h3, dp=3)
chk_num("LD", "r²(H3, rs35705950)", 0.115, r_h3**2, dp=3)
chk_p("LD", "p of r(H3, rs35705950)", 9.5e-24, p_h3)
h3c = df.H3_count > 0
chk_num("LD", "% H3 carriers also T-carriers", 60.5,
        100*int((df.loc[h3c, "rs35705950_T_count"] > 0).sum())/int(h3c.sum()), dp=1)
chk_num("LD", "r(MUC5AC H1, rs35705950)", -0.113,
        np.corrcoef(df.H1_count, df.rs35705950_T_count)[0,1], dp=3)
chk_num("LD", "r(MUC5B H1, rs35705950)", -0.104,
        np.corrcoef(df.MUC5B_H1_count, df.rs35705950_T_count)[0,1], dp=3)
chk_num("LD", "r(MUC5AC H2, rs35705950)", -0.225,
        np.corrcoef(df.H2_count, df.rs35705950_T_count)[0,1], dp=3)

# %% [markdown]
# ## 3. Table 2 — haplogroup associations (regression N=798)

# %%
def check_model(section, formula, terms):
    m = fit(formula)
    for term, (mo, ml, mh, mp) in terms.items():
        o, l, h, p = orci(m, term)
        chk_num(section, f"{term} OR", mo, o); chk_num(section, f"{term} CI-lo", ml, l)
        chk_num(section, f"{term} CI-hi", mh, h); chk_p(section, f"{term} p", mp, p)
    return m

check_model("Table 2 / 1a H3 unadj", f"case~H3_count+{covs}",
            {"H3_count": (1.69, 1.31, 2.17, 4.2e-5)})
check_model("Table 2 / 1b H3+rs", f"case~H3_count+rs35705950_T_count+{covs}",
            {"H3_count": (1.17, 0.89, 1.54, 0.27),
             "rs35705950_T_count": (3.75, 2.74, 5.13, 1.9e-16)})
check_model("Table 2 / joint (H1 ref)", f"case~H2_count+H3_count+{covs}",
            {"H2_count": (0.62, 0.47, 0.82, 7.6e-4),
             "H3_count": (1.50, 1.15, 1.94, 2.4e-3)})
check_model("Table 2 / joint+rs", f"case~H2_count+H3_count+rs35705950_T_count+{covs}",
            {"H2_count": (0.71, 0.53, 0.95, 0.022),
             "H3_count": (1.08, 0.81, 1.44, 0.59),
             "rs35705950_T_count": (3.61, 2.63, 4.96, 2.4e-15)})
check_model("Table 2 / 4a MUC5B H1", f"case~MUC5B_H1_count+{covs}",
            {"MUC5B_H1_count": (0.72, 0.52, 0.99, 0.042)})
check_model("Table 2 / 4b MUC5B H1+rs", f"case~MUC5B_H1_count+rs35705950_T_count+{covs}",
            {"MUC5B_H1_count": (0.81, 0.58, 1.14, 0.23),
             "rs35705950_T_count": (3.91, 2.89, 5.29, 8.2e-19)})
check_model("Table 2 / rs alone", f"case~rs35705950_T_count+{covs}",
            {"rs35705950_T_count": (3.96, 2.93, 5.34, 3.0e-19)})

# Supplementary Table S3 — joint-model reference-category sensitivity (H2 as reference)
check_model("Table S3 / H2 ref", f"case~H1_count+H3_count+{covs}",
            {"H1_count": (1.62, 1.22, 2.14, 7.6e-4),
             "H3_count": (2.42, 1.74, 3.37, 1.6e-7)})
check_model("Table S3 / H2 ref + rs35705950", f"case~H1_count+H3_count+rs35705950_T_count+{covs}",
            {"H1_count": (1.40, 1.05, 1.88, 0.022),
             "H3_count": (1.52, 1.06, 2.17, 0.022)})

# %% [markdown]
# ## 4. Table S5 — interaction ORs and RERI (delta method)

# %%
def reri_delta(m, ht, et):
    it = f"{ht}:{et}"; b = m.params; V = m.cov_params()
    reri = np.exp(b[ht]+b[et]+b.get(it,0)) - np.exp(b[ht]) - np.exp(b[et]) + 1
    g = np.zeros(len(b)); idx = {n:i for i,n in enumerate(b.index)}
    g[idx[ht]] = np.exp(b[ht]+b[et]+b.get(it,0)) - np.exp(b[ht])
    g[idx[et]] = np.exp(b[ht]+b[et]+b.get(it,0)) - np.exp(b[et])
    g[idx[it]] = np.exp(b[ht]+b[et]+b.get(it,0))
    se = np.sqrt(g @ V.values @ g)
    return reri, reri-1.96*se, reri+1.96*se

inter = [("H3_count","peto_exposed","H3xAsbestos",0.78,0.46,1.31,-0.29),
         ("H3_count","ever_smoked","H3xSmoking",0.95,0.53,1.71,0.24),
         ("MUC5B_H1_count","peto_exposed","MUC5B_H1xAsbestos",0.69,0.35,1.37,-0.35),
         ("MUC5B_H1_count","ever_smoked","MUC5B_H1xSmoking",0.93,0.47,1.85,-0.22)]
for hc, ec, lab, mo, ml, mh, mreri in inter:
    m = fit(f"case~{hc}+{ec}+{hc}:{ec}+{covs}")
    o, l, h, p = orci(m, f"{hc}:{ec}")
    chk_num("Table S5", f"{lab} intOR", mo, o); chk_num("Table S5", f"{lab} CI-lo", ml, l)
    chk_num("Table S5", f"{lab} CI-hi", mh, h)
    chk_tol("Table S5", f"{lab} RERI", mreri, reri_delta(m, hc, ec)[0], tol=0.01)

# %% [markdown]
# ## 5. Table S1 — rs35705950-stratified

# %%
dm_neg = dm[(dm.rs35705950_T_count == 0) & (~dm.centre.isin([18]))]
dm_pos = dm[(dm.rs35705950_T_count >= 1) & (~dm.centre.isin([17, 21]))]
chk_int("Table S1", "Non-carrier N", 468, len(dm_neg))
chk_int("Table S1", "Carrier N", 320, len(dm_pos))
s1 = {("H3_count","Non-carrier"):(1.23,468), ("H3_count","Carrier"):(1.10,320),
      ("MUC5B_H1_count","Non-carrier"):(0.71,468), ("MUC5B_H1_count","Carrier"):(0.97,320)}
for (hc, strat), (mo, mn) in s1.items():
    d = dm_neg if strat == "Non-carrier" else dm_pos
    o, l, h, p = orci(fit(f"case~{hc}+{covs}", d), hc)
    chk_num("Table S1", f"{hc} {strat} OR", mo, o)

# %% [markdown]
# ## 6. Table S2 — TaqMan sensitivity

# %%
o, *_ = orci(fit(f"case~H3_count+tq_bin+{covs}"), "H3_count")
chk_num("Table S2", "H3+TaqMan binary OR", 1.19, o)
o, *_ = orci(fit(f"case~H3_count+tq_count+{covs}"), "H3_count")
chk_num("Table S2", "H3+TaqMan additive OR", 1.19, o)

# %% [markdown]
# ## 7. Hardy-Weinberg in controls; TaqMan concordance; attenuation

# %%
from scipy.stats import chi2 as chi2_dist
ctrl = dm[dm.case == 0]
def hwe_p(col):
    cnt = ctrl[col]; n = len(ctrl); p_hat = cnt.mean()/2; q = 1-p_hat
    obs = np.array([(cnt==0).sum(), (cnt==1).sum(), (cnt==2).sum()], float)
    exp = np.array([n*q*q, n*2*p_hat*q, n*p_hat*p_hat])
    return 1 - chi2_dist.cdf(((obs-exp)**2/exp).sum(), df=1)
chk_num("HWE", "MUC5AC H3 p", 0.65, hwe_p("H3_count"), dp=2)
chk_num("HWE", "MUC5AC H1 p", 0.76, hwe_p("H1_count"), dp=2)
chk_num("HWE", "MUC5B H1 p", 0.83, hwe_p("MUC5B_H1_count"), dp=2)

conc = int((df.tq_count.round() == df.rs35705950_T_count).sum())
chk_int("TaqMan", "Concordant calls n", 763, conc)
chk_num("TaqMan", "Concordance %", 92.0, 100*conc/len(df), dp=1)
chk_num("TaqMan", "r(TaqMan, imputed)", 0.86,
        np.corrcoef(df.tq_count, df.rs35705950_T_count)[0,1], dp=2)

o1a, *_ = orci(fit(f"case~H3_count+{covs}"), "H3_count")
o1b, *_ = orci(fit(f"case~H3_count+rs35705950_T_count+{covs}"), "H3_count")
chk_num("Attenuation", "H3 attenuation %", 70, 100*(1-math.log(o1b)/math.log(o1a)), dp=0)
o4a, *_ = orci(fit(f"case~MUC5B_H1_count+{covs}"), "MUC5B_H1_count")
o4b, *_ = orci(fit(f"case~MUC5B_H1_count+rs35705950_T_count+{covs}"), "MUC5B_H1_count")
chk_num("Attenuation", "MUC5B H1 attenuation %", 38, 100*(1-math.log(o4b)/math.log(o4a)), dp=0)

# %% [markdown]
# ## 8. (optional, slow) Bootstrap RERI CIs and interaction power

# %%
if RUN_SLOW:
    def reri_point(m, ht, et):
        b = m.params
        return np.exp(b[ht]+b[et]+b.get(f"{ht}:{et}",0)) - np.exp(b[ht]) - np.exp(b[et]) + 1
    boot_exp = {"H3xAsbestos":(-1.64,0.51), "H3xSmoking":(-0.92,1.14),
                "MUC5B_H1xAsbestos":(-1.23,0.31), "MUC5B_H1xSmoking":(-1.19,0.39)}
    for hc, ec, lab, *_ in inter:
        rng = np.random.default_rng(42); boot = []
        for _ in range(500):
            bd = dm.iloc[rng.integers(0, len(dm), len(dm))]
            try: boot.append(reri_point(smf.logit(f"case~{hc}+{ec}+{hc}:{ec}+{covs}", data=bd).fit(disp=0), hc, ec))
            except: pass
        lo, hi = np.percentile(boot, [2.5, 97.5]); ml, mh = boot_exp[lab]
        chk_num("Bootstrap", f"{lab} boot lo", ml, lo, dp=2)
        chk_num("Bootstrap", f"{lab} boot hi", mh, hi, dp=2)

    # Interaction power (Results / Table S5) — the two scenarios cited in the manuscript
    def _power(nc, nctrl, ph, pe, oh, oe, oi, nsim=2000, seed=42):
        rng = np.random.default_rng(seed); rej = 0; n = nc + nctrl
        for _ in range(nsim):
            h = rng.binomial(1, ph, n); e = rng.binomial(1, pe, n)
            pr = 1/(1+np.exp(-(np.log(oh)*h + np.log(oe)*e + np.log(oi)*h*e)))
            y = rng.binomial(1, pr, n)
            try:
                m = smf.logit("y~h+e+h:e", data=pd.DataFrame({"y": y, "h": h, "e": e})).fit(disp=0)
                if m.pvalues.get("h:e", 1) < 0.05: rej += 1
            except Exception: pass
        return rej / nsim
    _ph3 = dm.H3_count.mean()/2; _ph1 = dm.MUC5B_H1_count.mean()/2; _pa = dm.peto_exposed.mean()
    chk_int("Power", "H3xAsbestos OR1.5 (%)", 20, round(100*_power(397, 401, _ph3, _pa, 1.69, 1.3, 1.5)))
    chk_int("Power", "MUC5B_H1xAsbestos OR0.5 (%)", 34, round(100*_power(397, 401, _ph1, _pa, 0.72, 1.3, 0.5)))
else:
    print("RUN_SLOW=False — skipped bootstrap/power (set RUN_SLOW=True to verify those).")

# %% [markdown]
# ## Summary — any FAIL needs attention

# %%
res = pd.DataFrame(_checks)
n_fail = (res.result == "FAIL").sum()
print(f"{len(res)} checks · {len(res)-n_fail} PASS · {n_fail} FAIL\n")
if n_fail:
    print("FAILURES:")
    print(res[res.result == "FAIL"].to_string(index=False))
else:
    print("✅ All recomputed numbers reconcile with the manuscript.")
res

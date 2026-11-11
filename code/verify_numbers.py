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
# This notebook re-derives the manuscript's tabulated figures (Tables 1, 2 and
# Supplementary Tables S1–S4), the primary and conditioned H2 estimates, and the
# LD / HWE / LRT / clade-ordering numbers in `../MUC5_VNTR_haplogroup_paper_draft.md`
# directly from the source data, and compares each against the printed value. Any
# row flagged `FAIL` needs attention. (Supplementary Table S5, imputed length, and
# S6, the cluster scan, are external analyses verified by their own scripts
# `code/12` and `code/13`, not recomputed here.)
#
# **Analytical sample.** The manuscript reports a single sample of N=798
# (397 cases, 401 controls) throughout: the participants with complete
# genotype, imputed chromosome 11 and covariate data, after excluding one
# cases-only recruitment centre (centre 16, n=31) that caused complete
# separation. Every check below therefore runs on `dm`, not on the wider
# 829-participant merge.
#
# **IPFJES cohort only.** The genotype pipeline is restricted to IPFJES at the
# build step (code/00), so muc5_haplogroups.csv contains only IPFJES participants.
# Allele and haplogroup frequencies are computed in the 401 IPFJES controls, which
# approximate a population reference. No external-cohort control data is used here.
#
# Run top-to-bottom in a fresh kernel.

# %%
import math, warnings
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
warnings.filterwarnings("ignore")

ROOT = "/home/drcjar/academic/ipfjes-omics"

# Collects every check so we can summarise pass/fail at the end.
_checks = []
# Provenance for the current block of checks: which variant column, which model
# formula, which N. Set with `prov(...)` before a group of checks. Supplementary
# Table S2's H2 rows were computed from the wrong tagging SNP for weeks because a
# number carried no record of where it came from; this makes that visible.
_prov = {}
from contextlib import contextmanager
@contextmanager
def prov(**kw):
    """Scope provenance to a block. Sticky state would attach one block's tag and
    model to every later check, which is worse than recording nothing."""
    _prov.clear(); _prov.update(kw)
    try:
        yield
    finally:
        _prov.clear()

def record(section, item, manuscript, computed, ok):
    _checks.append(dict(section=section, item=item,
                        manuscript=manuscript, computed=computed,
                        result="PASS" if ok else "FAIL", **_prov))

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
# The eigenvector file is itself restricted to the 798 analytical participants
# (PCs are computed within-sample), so the merge above already yields N=798.
dm = df[(df["centre"] != 16) & df.ethnicity.astype(str).str.contains("White", case=False, na=False)].copy()   # analytical sample: European-ancestry restricted

# The primary MUC5AC H2 tag, named once so a script cannot quietly substitute a
# panel tag for it as 09_h2_direct_and_lrt.py once did.
NEW_TAG_ID = "chr11:1182279"
# Primary H2 measure: direct tag rs879201833 (r2=0.888), replacing the diplotype complement.
_h2raw = pd.read_csv(f"{ROOT}/muc5ac_h2_tag.raw", sep=r"\s+")
_h2raw["pid"] = _h2raw["IID"].astype(str).str.split("_").str[0]
_h2tcol = [c for c in _h2raw.columns if c.startswith(NEW_TAG_ID)][0]
assert _h2tcol.startswith(NEW_TAG_ID), f"H2 tag resolved to {_h2tcol}, not {NEW_TAG_ID}"
dm = dm.merge(_h2raw[["pid", _h2tcol]], on="pid")
_h2f = dm.loc[dm.case == 0, _h2tcol].mean() / 2
dm["H2t"] = dm[_h2tcol] if _h2f < 0.5 else 2 - dm[_h2tcol]   # oriented to the H2-tracking allele
cases, controls = dm[dm.case == 1], dm[dm.case == 0]

chk_int("Cohort", "Analytical N", 771, len(dm))
chk_int("Cohort", "Cases", 385, len(cases))
chk_int("Cohort", "Controls", 386, len(controls))
chk_int("Cohort", "PCA computed within analytical sample only", 771, len(pcs))
print(f"Analytical N={len(dm)} (cases={len(cases)}, controls={len(controls)})")

# %%
covs = "age + ever_smoked + PC1+PC2+PC3+PC4+PC5 + C(centre)"
def fit(f, d=dm):
    try:    return smf.logit(f, data=d).fit(disp=0)
    except: return smf.logit(f, data=d).fit(method="bfgs", disp=0)
def orci(m, t):
    b, se, p = m.params[t], m.bse[t], m.pvalues[t]
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), p

# %% [markdown]
# ## 1. Table 1 — participant characteristics (N=798)

# %%
age_p = stats.mannwhitneyu(cases.age, controls.age).pvalue
chk_num("Table 1", "Age p", 0.01, age_p, dp=3)

# Medians and IQRs were previously unguarded: the Table 1 age IQR survived the
# 829->798 consolidation as a stale value because only p-values and n(%) cells
# were checked. Every continuous Table 1 cell is now verified.
for _lab, _col, _d, _q in [
        ("Age cases", "age", cases, (72.0, 77.0, 82.0)),
        ("Age controls", "age", controls, (71.0, 76.0, 81.0)),
        ("Pack-years cases", "packyrs", cases[cases.ever_smoked == 1], (9.0, 20.0, 39.0)),
        ("Pack-years controls", "packyrs", controls[controls.ever_smoked == 1], (7.0, 19.0, 35.0))]:
    _v = _d[_col].dropna().quantile([.25, .5, .75]).values
    for _n, _obs, _exp in zip(("Q1", "median", "Q3"), _v, _q):
        chk_num("Table 1", f"{_lab} {_n}", _exp, _obs, dp=1)

py_p = stats.mannwhitneyu(cases[cases.ever_smoked == 1].packyrs.dropna(),
                          controls[controls.ever_smoked == 1].packyrs.dropna()).pvalue
chk_num("Table 1", "Pack-years p", 0.11, py_p, dp=3)

es_c, es_ct = int(cases.ever_smoked.sum()), int(controls.ever_smoked.sum())
chk_int("Table 1", "Ever-smoked cases n", 292, es_c)
chk_int("Table 1", "Ever-smoked controls n", 269, es_ct)
chk_num("Table 1", "Ever-smoked cases %", 75.8, 100*es_c/len(cases), dp=1)
chk_num("Table 1", "Ever-smoked controls %", 69.7, 100*es_ct/len(controls), dp=1)

ab_c, ab_ct = int(cases.peto_exposed.sum()), int(controls.peto_exposed.sum())
chk_int("Table 1", "Asbestos cases n", 260, ab_c)
chk_int("Table 1", "Asbestos controls n", 242, ab_ct)
chk_num("Table 1", "Asbestos cases %", 67.5, 100*ab_c/len(cases), dp=1)
chk_num("Table 1", "Asbestos controls %", 62.7, 100*ab_ct/len(controls), dp=1)
chk_int("Table 1", "Asbestos-exposed total n (Methods)", 502, int(dm.peto_exposed.sum()))

for col, lab, mc, mct in [("H1_count","MUC5AC H1",307,318), ("H2t","MUC5AC H2",99,157),
                          ("H3_count","MUC5AC H3",202,140), ("MUC5B_H1_count","MUC5B H1",79,107)]:
    chk_int("Table 1", f"{lab} carrier cases n", mc, int((cases[col] > 0).sum()))
    chk_int("Table 1", f"{lab} carrier controls n", mct, int((controls[col] > 0).sum()))

rs_c  = int((cases.rs35705950_T_count > 0).sum())
rs_ct = int((controls.rs35705950_T_count > 0).sum())
# Table 1's categorical p-values. Chi-squared on the 2x2 of carrier status by
# case status, matching the test 02_primary_analysis.py uses, but recomputed so
# the printed column is checked rather than assumed.
def _cat_p(nc, nct):
    return stats.chi2_contingency([[nc, len(cases) - nc],
                                   [nct, len(controls) - nct]])[1]
chk_p("Table 1", "Ever-smoked p", 0.066, _cat_p(es_c, es_ct))
chk_p("Table 1", "Asbestos p", 0.182, _cat_p(ab_c, ab_ct))
for _col, _lab, _mp in [("H1_count", "MUC5AC H1", 0.398), ("H2t", "MUC5AC H2", 1.5e-5),
                        ("H3_count", "MUC5AC H3", 8.4e-6), ("MUC5B_H1_count", "MUC5B H1", 0.024)]:
    chk_p("Table 1", f"{_lab} carrier p", _mp,
          _cat_p(int((cases[_col] > 0).sum()), int((controls[_col] > 0).sum())))

chk_int("Table 1", "rs35705950 T-carrier cases n", 230, rs_c)
chk_int("Table 1", "rs35705950 T-carrier controls n", 89, rs_ct)
chk_num("Table 1", "rs35705950 T-carrier cases %", 59.7, 100*rs_c/len(cases), dp=1)
chk_num("Table 1", "rs35705950 T-carrier controls %", 23.1, 100*rs_ct/len(controls), dp=1)
chk_int("Table 1", "TT cases n", 26, int((cases.rs35705950_T_count == 2).sum()))
chk_int("Table 1", "TT controls n", 4, int((controls.rs35705950_T_count == 2).sum()))

_, rs_p, _, _ = stats.chi2_contingency([[rs_c, len(cases)-rs_c], [rs_ct, len(controls)-rs_ct]])
chk_p("Table 1", "rs35705950 T-carrier p", 9.85e-25, rs_p)

# %% [markdown]
# ## 2. Haplogroup frequencies (IPFJES controls) and LD with rs35705950

# %%
# Quoted approximately in the Results, computed in the 401 IPFJES controls —
# compare with a ~1 percentage-point tolerance.
chk_tol("Freqs", "MUC5AC H1 freq % (controls)", 57.12, 100*controls.H1_count.mean()/2, tol=1.0)
chk_tol("Freqs", "MUC5AC H2 freq % (controls)", 22.54, 100*controls.H2t.mean()/2, tol=1.0)
chk_tol("Freqs", "MUC5AC H3 freq % (controls)", 20.34, 100*controls.H3_count.mean()/2, tol=1.0)
chk_tol("Freqs", "MUC5B H2 freq % (controls)", 85.23, 100*controls.MUC5B_H2_count.mean()/2, tol=1.0)
chk_tol("Freqs", "MUC5B H1 freq % (controls)", 14.77, 100*controls.MUC5B_H1_count.mean()/2, tol=1.0)

# Discordant subset — H3-carrier status disagrees with rs35705950-carrier status.
hd = (dm.H3_count > 0) != (dm.rs35705950_T_count > 0)
chk_int("Discordant", "H3/rs35705950 discordant n", 241, int(hd.sum()))
chk_num("Discordant", "H3/rs35705950 discordant %", 31.0, 100*hd.mean(), dp=0)
chk_int("Discordant", "discordant cases (n)", 118, int(hd[dm.case==1].sum()))
chk_int("Discordant", "discordant controls (n)", 123, int(hd[dm.case==0].sum()))
bd = (dm.MUC5B_H1_count > 0) != (dm.rs35705950_T_count > 0)
chk_int("Discordant", "MUC5B H1 discordant n", 383, int(bd.sum()))

# Supplementary Table S1 — per-tSNP counted-allele frequencies in the 386 controls.
_dos = pd.read_csv(f"{ROOT}/muc5_tsnps_dosage_pheno.csv")
_dos["pid"] = _dos["IID"].astype(str).str.split("_").str[0]
_ctrl_pid = set(controls.pid)
_dc = _dos[_dos.pid.isin(_ctrl_pid)]
s3_freq = {
    "chr11:1165501:T:A_T": 0.797, "chr11:1195214:G:C_C": 0.571, "chr11:1193830:G:C_G": 0.574,
    "chr11:1194354:A:G_A": 0.574, "chr11:1195265:G:C_C": 0.573, "chr11:1195858:T:C_T": 0.573,
    "chr11:1182048:A:T_T": 0.574, "chr11:1202343:T:C_C": 0.773, "chr11:1202859:A:G_G": 0.773,
    "chr11:1202914:T:C_C": 0.773, "chr11:1178216:C:G_C": 0.803, "chr11:1178023:A:G_A": 0.797,
    "chr11:1177688:A:T_A": 0.797, "chr11:1177356:T:C_T": 0.797, "chr11:1177269:T:C_T": 0.797,
    "chr11:1197482:C:T_C": 0.799, "chr11:1219991:G:T_G": 0.880,
}
for _cc, _exp in s3_freq.items():
    chk_num("Table S1", f"{_cc} freq", _exp, _dc[_cc].mean()/2, dp=3)
# MUC5B discriminator: counted allele tags H2, taken from the haplogroup calls.
chk_num("Table S1", "chr11:1244757 (MUC5B H2) freq", 0.852, controls.MUC5B_H2_count.mean()/2, dp=3)

r_h3, p_h3 = stats.pearsonr(dm.H3_count, dm.rs35705950_T_count)
chk_num("LD", "r(H3, rs35705950)", 0.363, r_h3, dp=3)
chk_num("LD", "r²(H3, rs35705950)", 0.132, r_h3**2, dp=3)
chk_p("LD", "p of r(H3, rs35705950)", 1.74e-25, p_h3)
h3c = dm.H3_count > 0
chk_num("LD", "% H3 carriers also T-carriers", 61.4,
        100*int((dm.loc[h3c, "rs35705950_T_count"] > 0).sum())/int(h3c.sum()), dp=1)
chk_num("LD", "r(MUC5AC H1, rs35705950)", -0.14,
        np.corrcoef(dm.H1_count, dm.rs35705950_T_count)[0,1], dp=3)
chk_num("LD", "r(MUC5B H1, rs35705950)", -0.113,
        np.corrcoef(dm.MUC5B_H1_count, dm.rs35705950_T_count)[0,1], dp=3)
chk_num("LD", "r(MUC5AC H2, rs35705950)", -0.235,
        np.corrcoef(dm.H2t, dm.rs35705950_T_count)[0,1], dp=3)
# Supplementary Table S6 — LD (r2) in the 401 controls (population reference) for H3 and H2.
_ctl = dm[dm.case == 0]
chk_num("Table S6", "H3 controls r2", 0.1,
        np.corrcoef(_ctl.H3_count, _ctl.rs35705950_T_count)[0, 1]**2, dp=2)
chk_num("Table S6", "H2 controls r2", 0.02,
        np.corrcoef(_ctl.H2t, _ctl.rs35705950_T_count)[0, 1]**2, dp=2)
chk_num("Table S6", "H3 full-sample r2", 0.13,
        np.corrcoef(dm.H3_count, dm.rs35705950_T_count)[0, 1]**2, dp=2)

# %% [markdown]
# ## 3. Table 2 — haplogroup associations (N=798)

# %%
def check_model(section, formula, terms):
    m = fit(formula)
    for term, (mo, ml, mh, mp) in terms.items():
        o, l, h, p = orci(m, term)
        chk_num(section, f"{term} OR", mo, o); chk_num(section, f"{term} CI-lo", ml, l)
        chk_num(section, f"{term} CI-hi", mh, h); chk_p(section, f"{term} p", mp, p)
    return m

check_model("Table 2 / 1a H3 unadj", f"case~H3_count+{covs}",
            {"H3_count": (1.63, 1.26, 2.09, 0.000179)})
check_model("Table 2 / 1b H3+rs", f"case~H3_count+rs35705950_T_count+{covs}",
            {"H3_count": (1.10, 0.83, 1.45, 0.521),
             "rs35705950_T_count": (3.90, 2.83, 5.38, 1.1e-16)})
check_model("Table 2 / joint (H1 ref)", f"case~H2t+H3_count+{covs}",
            {"H2t": (0.58, 0.43, 0.78, 0.000267),
             "H3_count": (1.43, 1.10, 1.86, 0.00806)})
check_model("Table 2 / joint+rs", f"case~H2t+H3_count+rs35705950_T_count+{covs}",
            {"H2t": (0.69, 0.51, 0.94, 0.0171),
             "H3_count": (1.02, 0.76, 1.36, 0.91),
             "rs35705950_T_count": (3.72, 2.68, 5.15, 2.6e-15)})
check_model("Table 2 / 4a MUC5B H1", f"case~MUC5B_H1_count+{covs}",
            {"MUC5B_H1_count": (0.71, 0.52, 0.98, 0.0381)})
check_model("Table 2 / 4b MUC5B H1+rs", f"case~MUC5B_H1_count+rs35705950_T_count+{covs}",
            {"MUC5B_H1_count": (0.82, 0.58, 1.15, 0.252),
             "rs35705950_T_count": (3.98, 2.93, 5.41, 1.2e-18)})
check_model("Table 2 / rs alone", f"case~rs35705950_T_count+{covs}",
            {"rs35705950_T_count": (4.03, 2.97, 5.48, 4.0e-19)})

# %% [markdown]
# ## 4. Table S2 — rs35705950-stratified

# %%
dm_neg = dm[(dm.rs35705950_T_count == 0) & (~dm.centre.isin([18]))]
dm_pos = dm[(dm.rs35705950_T_count >= 1) & (~dm.centre.isin([17, 21]))]
chk_int("Table S2", "Non-carrier N", 448, len(dm_neg))
chk_int("Table S2", "Carrier N", 313, len(dm_pos))
s1 = {("H3_count","Non-carrier"):1.15, ("H3_count","Carrier"):1.00, ("MUC5B_H1_count","Non-carrier"):0.72, ("MUC5B_H1_count","Carrier"):1.03}
with prov(tag="H3/MUC5B H1 counts from muc5_haplogroups.csv",
          model="case~<hg>+age+ever_smoked+PC1-5+C(centre)",
          strata="rs35705950 carrier status"):
  for (hc, strat), mo in s1.items():
    d = dm_neg if strat == "Non-carrier" else dm_pos
    o, l, h, p = orci(fit(f"case~{hc}+{covs}", d), hc)
    chk_num("Table S2", f"{hc} {strat} OR", mo, o)
    # intervals recorded so 07_consistency_checks.py can check the Discussion,
    # which quotes the H3 non-carrier interval directly
    _rec_ci = {("H3_count", "Non-carrier"): (0.79, 1.69), ("H3_count", "Carrier"): (0.64, 1.54),
               ("MUC5B_H1_count", "Non-carrier"): (0.46, 1.11), ("MUC5B_H1_count", "Carrier"): (0.54, 1.95)}
    if (hc, strat) in _rec_ci:
        _cl, _ch = _rec_ci[(hc, strat)]
        chk_num("Table S2", f"{hc} {strat} CI-lo", _cl, l)
        chk_num("Table S2", f"{hc} {strat} CI-hi", _ch, h)
        _rec_p = {("H3_count", "Non-carrier"): 0.46, ("H3_count", "Carrier"): 0.98,
                  ("MUC5B_H1_count", "Non-carrier"): 0.14, ("MUC5B_H1_count", "Carrier"): 0.93}
        chk_p("Table S2", f"{hc} {strat} p", _rec_p[(hc, strat)], p)

# %% [markdown]
# ## 5. Hardy-Weinberg in the IPFJES controls, and attenuation

# %%
from scipy.stats import chi2 as chi2_dist
def hwe_p(counts):
    cnt = np.asarray(counts).round().astype(int); n = len(cnt)
    p_hat = cnt.mean()/2; q = 1 - p_hat
    obs = np.array([(cnt==0).sum(), (cnt==1).sum(), (cnt==2).sum()], float)
    exp = np.array([n*q*q, n*2*p_hat*q, n*p_hat*p_hat])
    return 1 - chi2_dist.cdf(((obs-exp)**2/exp).sum(), df=1)

chk_num("HWE", "MUC5AC H3 p", 0.74, hwe_p(controls.H3_count), dp=2)
chk_num("HWE", "MUC5AC H1 p", 0.54, hwe_p(controls.H1_count), dp=2)
chk_num("HWE", "MUC5B H1 p", 0.57, hwe_p(controls.MUC5B_H1_count), dp=2)

# Methods: "All tSNPs passed Hardy-Weinberg testing in the IPFJES controls
# (minimum p=0.41); genotypes were hard-called from imputed data with no missingness."
# rs879201833 lives in muc5ac_h2_tag.raw, not the main dosage export, and was
# previously omitted from this set even though it is a row in Table S1.
_h2raw = pd.read_csv(f"{ROOT}/muc5ac_h2_tag.raw", sep=r"\s+")
_h2raw["pid"] = _h2raw["IID"].astype(str).str.split("_").str[0]
_h2c = _h2raw[_h2raw.pid.isin(set(controls.pid))]
_h2col = [c for c in _h2raw.columns if c.startswith("chr11:1182279")][0]
# The Results HWE sentence now names all four tags, and the paper's only positive
# finding rests on this one, so it gets its own recorded value instead of only
# entering the tSNP minimum below.
chk_num("HWE", "MUC5AC H2 tag p", 0.45, hwe_p(_h2c[_h2col]), dp=2)
_tsnp_hwe = ([hwe_p(_dc[c]) for c in s3_freq]
             + [hwe_p(controls.MUC5B_H2_count), hwe_p(_h2c[_h2col])])
chk_num("HWE", "minimum tSNP HWE p (controls)", 0.39, min(_tsnp_hwe), dp=2)
chk_int("HWE", "tSNPs with any missing genotype", 0,
        int(sum(_dc[c].isna().any() for c in s3_freq)))

o1a, *_ = orci(fit(f"case~H3_count+{covs}"), "H3_count")
o1b, *_ = orci(fit(f"case~H3_count+rs35705950_T_count+{covs}"), "H3_count")
chk_num("Attenuation", "H3 attenuation %", 81.0, 100*(1-math.log(o1b)/math.log(o1a)), dp=0)
o4a, *_ = orci(fit(f"case~MUC5B_H1_count+{covs}"), "MUC5B_H1_count")
o4b, *_ = orci(fit(f"case~MUC5B_H1_count+rs35705950_T_count+{covs}"), "MUC5B_H1_count")
chk_num("Attenuation", "MUC5B H1 attenuation %", 42.0, 100*(1-math.log(o4b)/math.log(o4a)), dp=0)

# %% [markdown]
# ## 6. Direct H2 tagging, nested-model LRTs, reciprocal stratification
#
# Reported in the Discussion; computed by `09_h2_direct_and_lrt.py`.

# %%
_NEW = _h2tcol   # rs879201833 column, already merged into dm above
_d = dm.merge(_dos[["pid", "chr11:1202343:T:C_C", "chr11:1202859:A:G_G",
                    "chr11:1202914:T:C_C"]], on="pid")   # all three panel H2 tags

def _h2_direct(col, d):
    f = d.loc[d.case == 0, col].mean() / 2
    d = d.assign(H2d=d[col] if f < 0.5 else 2 - d[col])
    m = smf.logit(f"case~H2d+H3_count+rs35705950_T_count+{covs}", data=d).fit(disp=0)
    b, se, p = m.params["H2d"], m.bse["H2d"], m.pvalues["H2d"]
    return np.exp(b), np.exp(b - 1.96*se), np.exp(b + 1.96*se), p

for _lab, _c, _mo, _ml, _mh, _mp in [
        ("rs879201833 H2-adj", _NEW, 0.69, 0.51, 0.94, 0.017),
        ("rs28519516 H2-adj", "chr11:1202343:T:C_C", 0.66, 0.49, 0.90, 0.009),
        ("rs28558973 H2-adj", "chr11:1202859:A:G_G", 0.66, 0.49, 0.90, 0.009),
        ("rs28368633 H2-adj", "chr11:1202914:T:C_C", 0.66, 0.49, 0.90, 0.009)]:
    _o, _l, _h, _p = _h2_direct(_c, _d)
    chk_num("H2 direct", f"{_lab} OR", _mo, _o); chk_num("H2 direct", f"{_lab} CI-lo", _ml, _l)
    chk_num("H2 direct", f"{_lab} CI-hi", _mh, _h); chk_num("H2 direct", f"{_lab} p", _mp, _p, dp=3)

# Discussion headline robustness: H2 tag conditioned on the 8 Kousathanas 11p15.5 signals (code/14).
_ks = pd.read_csv(f"{ROOT}/kous_independent_signals.raw", sep=r"\s+")
_ks["pid"] = _ks["IID"].astype(str).str.split("_").str[0]
_krename = {c: "S_" + c.split(":")[1] for c in _ks.columns if c.startswith("chr11:")}
_ks = _ks.rename(columns=_krename)
_dk = dm.merge(_ks[["pid"] + list(_krename.values())], on="pid")
_mk = smf.logit(f"case~H2t+H3_count+rs35705950_T_count+{'+'.join(_krename.values())}+{covs}",
                data=_dk).fit(disp=0)
chk_num("H2 conditioned on 8 signals", "OR", 0.71, np.exp(_mk.params["H2t"]))
chk_num("H2 conditioned on 8 signals", "p", 0.051, _mk.pvalues["H2t"], dp=3)

_l_rs = fit(f"case~rs35705950_T_count+{covs}").llf
_l_h3 = fit(f"case~H3_count+{covs}").llf
_l_bo = fit(f"case~H3_count+rs35705950_T_count+{covs}").llf
chk_num("LRT", "add H3 to rs35705950 chi2", 0.41, 2*(_l_bo - _l_rs), dp=2)
chk_num("LRT", "add H3 to rs35705950 p", 0.52, stats.chi2.sf(2*(_l_bo - _l_rs), 1), dp=2)
chk_num("LRT", "add rs35705950 to H3 chi2", 77.8, 2*(_l_bo - _l_h3), dp=1)
chk_p("LRT", "add rs35705950 to H3 p", 1.15e-18, stats.chi2.sf(2*(_l_bo - _l_h3), 1))

_neg = dm[dm.H3_count == 0]
chk_int("Stratified", "H3 non-carriers N", 429, len(_neg))
_o, _l, _h, _p = orci(fit(f"case~rs35705950_T_count+{covs}", _neg), "rs35705950_T_count")
chk_num("Stratified", "rs35705950 in H3 non-carriers OR", 3.35, _o)
chk_num("Stratified", "rs35705950 in H3 non-carriers CI-lo", 2.11, _l)
chk_num("Stratified", "rs35705950 in H3 non-carriers CI-hi", 5.32, _h)

# %% [markdown]
# ## 7. Additions reported after peer review
#
# H2 stratified by rs35705950 (Table S2); the H3-carrier stratum; the genotypic
# H3 model; and the minimum detectable residual effect quoted in the Discussion.

# %%
_neg2 = dm[(dm.rs35705950_T_count == 0) & (~dm.centre.isin([18]))]
_pos2 = dm[(dm.rs35705950_T_count >= 1) & (~dm.centre.isin([17, 21]))]
with prov(tag=NEW_TAG_ID, model="case~H2t+H3_count+age+ever_smoked+PC1-5+C(centre)",
          strata="rs35705950 carrier status"):
 for _lab, _d, _mo, _mp in [("Non-carrier", _neg2, 0.69, 0.054), ("Carrier", _pos2, 0.54, 0.040)]:
    _o, _l, _h, _p = orci(fit(f"case~H2t+H3_count+{covs}", _d), "H2t")
    chk_num("Table S2 H2", f"H2 {_lab} OR", _mo, _o)
    chk_num("Table S2 H2", f"H2 {_lab} p", _mp, _p, dp=3)

_pos3 = dm[dm.H3_count > 0]
chk_int("Stratified", "H3 carriers N", 342, len(_pos3))
_o, _l, _h, _p = orci(fit(f"case~rs35705950_T_count+{covs}", _pos3), "rs35705950_T_count")
chk_num("Stratified", "rs35705950 in H3 carriers OR", 4.4, _o)
_mi = fit(f"case~H3_count*rs35705950_T_count+{covs}")
chk_num("Stratified", "H3 x rs35705950 interaction p", 0.42,
        _mi.pvalues["H3_count:rs35705950_T_count"], dp=2)

_dg = dm.assign(H3g=dm.H3_count.round().astype(int))
_mg = smf.logit(f"case~C(H3g)+{covs}", data=_dg).fit(disp=0)
for _t, _m in [("C(H3g)[T.1]", 1.87), ("C(H3g)[T.2]", 1.89)]:
    chk_num("Genotypic", f"H3 {_t} OR", _m, np.exp(_mg.params[_t]))
_ma = fit(f"case~H3_count+{covs}")
chk_num("Genotypic", "LRT genotypic vs additive p", 0.14,
        stats.chi2.sf(2 * (_mg.llf - _ma.llf), 1), dp=2)

# Genotypic H2 (joint model H2 tag + H3 + rs35705950, H1 reference) — reported in Results
_dh2 = dm.assign(H2g=dm.H2t.round().astype(int))
_mgh2 = smf.logit(f"case~C(H2g)+H3_count+rs35705950_T_count+{covs}", data=_dh2).fit(disp=0)
for _t, _m in [("C(H2g)[T.1]", 0.65), ("C(H2g)[T.2]", 0.58)]:
    chk_num("Genotypic", f"H2 {_t} OR", _m, np.exp(_mgh2.params[_t]))
_mah2 = fit(f"case~H2t+H3_count+rs35705950_T_count+{covs}")
chk_num("Genotypic", "H2 LRT genotypic vs additive p", 0.58,
        stats.chi2.sf(2 * (_mgh2.llf - _mah2.llf), 1), dp=2)

_m = fit(f"case~H3_count+rs35705950_T_count+{covs}")
chk_num("Precision", "H3|rs35705950 minimum detectable OR (80% power)", 1.5,
        math.exp((1.96 + 0.8416) * _m.bse["H3_count"]))

# The limitations paragraph pools the adjusted H3 estimate with TOPMed. That figure
# was prose only, which is how the supplementary frequencies went stale, so it is
# recomputed here exactly as code/15 does it: fixed-effect inverse-variance weights
# on log ORs, both SEs recovered from two-decimal 95% intervals, so this and code/15
# cannot disagree in the last decimal.
# The IPFJES triple is NOT hardcoded. It is derived from the fitted model `_m` above
# and rounded the way Table 2 prints it, so a change to the analytical sample moves
# these inputs automatically. Hardcoding 1.10/0.83/1.45 would have let the pooled
# figure keep passing on a stale input while Table 2's own check failed, which is
# precisely the drift this is here to catch. Only the TOPMed triple is a literal:
# it is their published estimate, external data we do not hold.
def _logor(o, lo, hi):
    return math.log(o), (math.log(hi) - math.log(lo)) / (2 * 1.959964)
_IPF_1B = (round(math.exp(_m.params["H3_count"]), 2),
           round(math.exp(_m.params["H3_count"] - 1.959964 * _m.bse["H3_count"]), 2),
           round(math.exp(_m.params["H3_count"] + 1.959964 * _m.bse["H3_count"]), 2))
with prov(tag="H3_count", model="IVW fixed-effect pool: IPFJES 1b H3+rs (fitted, rounded as Table 2) with TOPMed (published)",
          strata="IPFJES + TOPMed"):
    (_b1, _s1), (_b2, _s2) = _logor(*_IPF_1B), _logor(0.94, 0.84, 1.05)
    _w1, _w2 = 1 / _s1**2, 1 / _s2**2
    _bm, _sm = (_w1*_b1 + _w2*_b2) / (_w1 + _w2), math.sqrt(1 / (_w1 + _w2))
    chk_num("TOPMed pooled / 1b H3+rs", "H3_count OR", 0.96, math.exp(_bm))
    chk_num("TOPMed pooled / 1b H3+rs", "H3_count CI-lo", 0.87, math.exp(_bm - 1.959964*_sm))
    chk_num("TOPMed pooled / 1b H3+rs", "H3_count CI-hi", 1.07, math.exp(_bm + 1.959964*_sm))

_eth = clin.merge(pcs[["pid"]], on="pid").dropna(
    subset=["age", "ever_smoked", "peto_exposed", "centre"]).query("centre != 16 and ethnicity.str.contains('White', case=False, na=False)")
_white = int((_eth[[c for c in _eth.columns if c.lower() == "ethnicity"][0]] == "White").sum())
chk_int("Ancestry", "white-ethnicity participants", 771, _white)
chk_num("Ancestry", "white-ethnicity %", 100.0, 100 * _white / len(_eth), dp=1)

# %% [markdown]
#
# %% [markdown]
# ## 9. The non-monotonic clade ordering cited in the Discussion
#
# H1 and H3 are indistinguishable after adjustment while both differ from the
# intermediate H2 — the argument that the residual is not a simple repeat-count
# effect. Guarded here because it is a claim about a PATTERN across three
# estimates, not a single number.

# %%
_mj = fit(f"case~H2_count+H3_count+rs35705950_T_count+{covs}")      # H1 reference
_o, _l, _h, _p = orci(_mj, "H3_count")
chk_num("Clade ordering", "H3 vs H1 (adjusted) OR", 1.02, _o)
chk_num("Clade ordering", "H3 vs H1 (adjusted) p", 0.88, _p, dp=2)
_m2 = fit(f"case~H1_count+H3_count+rs35705950_T_count+{covs}")      # H2 reference
chk_num("Clade ordering", "H1 vs H2 (adjusted) OR", 1.37, np.exp(_m2.params["H1_count"]))
chk_num("Clade ordering", "H3 vs H2 (adjusted) OR", 1.4, np.exp(_m2.params["H3_count"]))

# %% [markdown]
# ## 10. Power of the cluster-wide scan for the known 11p15.5 signals
#
# The scan's negative is quoted in the Results with a power figure; a reviewer asked
# for that figure to be backed up rather than asserted. Recomputed here by the same
# method as `17_cluster_scan_power.py` but written out independently, so the two are
# not a single code path. Published AF_EUR and OR are from Kousathanas et al. (2025)
# Supplementary Table 5. Checked with a tolerance of one percentage point, because
# the standard error is obtained from simulated genotypes and shifts slightly with
# the random stream (73.9-74.4% across four seeds).

# %%
_KOUS = [  # pos, rsid, af_eur, OR, nearest gene
    (1013992, "rs148815783", 0.0073, 2.10, "MUC6"),
    (1057948, "rs12574439",  0.12,   1.21, "MUC2"),
    (1078249, "rs7104590",   0.69,   1.15, "MUC2"),
    (1081750, "rs41359951",  0.018,  1.55, "MUC2"),
    (1208601, "rs548688097", 0.0070, 2.84, "MUC5B"),
    (1214934, "rs12802931",  0.20,   1.71, "MUC5B"),
    (1219775, "rs2672794",   0.61,   1.20, "MUC5B"),
    (1224068, "rs56321310",  0.19,   1.39, "MUC5B"),
]
_M_EFF = 54                       # 11_mucin_cluster_scan.py, Li & Ji (2005)
_ZC = stats.norm.isf((0.05 / _M_EFF) / 2)
_rng = np.random.default_rng(4242)        # deliberately not 17_...'s seed

with prov(tag="simulated genotype at published AF_EUR",
          model=f"case ~ g + {covs}", strata=f"N={len(dm)}"):
    _pw = {}
    for _pos, _rs, _af, _or, _gene in _KOUS:
        _ses = []
        for _ in range(200):
            _g = _rng.binomial(2, _af, size=len(dm)).astype(float)
            if _g.std() == 0:
                continue
            try:
                _m = fit(f"case~g+{covs}", dm.assign(g=_g))
            except Exception:
                continue
            _se = _m.bse.get("g", np.nan)
            if np.isfinite(_se) and _se < 10:
                _ses.append(_se)
        _ncp = np.log(_or) / float(np.mean(_ses))
        _pw[_rs] = stats.norm.sf(_ZC - _ncp) + stats.norm.cdf(-_ZC - _ncp)

    # The figure quoted in the Results: their most strongly associated MUC5B signal.
    chk_tol("Cluster-scan power", "rs12802931 (MUC5B, OR 1.71, AF 0.20) power",
            0.74, _pw["rs12802931"], tol=0.01)
    # And the range quoted for the four MUC2/MUC6 signals, which the scan cannot reach.
    _m26 = [_pw[r] for _, r, _, _, g in _KOUS if g in ("MUC2", "MUC6")]
    chk_tol("Cluster-scan power", "MUC2/MUC6 signals minimum power",
            0.012, min(_m26), tol=0.01)
    chk_tol("Cluster-scan power", "MUC2/MUC6 signals maximum power",
            0.017, max(_m26), tol=0.01)
    # The published effect sizes and frequencies quoted alongside them.
    _or26 = [o for _, _, _, o, g in _KOUS if g in ("MUC2", "MUC6")]
    _af26 = [a for _, _, a, _, g in _KOUS if g in ("MUC2", "MUC6")]
    chk_num("Cluster-scan power", "MUC2/MUC6 published OR minimum", 1.15, min(_or26))
    chk_num("Cluster-scan power", "MUC2/MUC6 published OR maximum", 2.10, max(_or26))
    chk_num("Cluster-scan power", "MUC2/MUC6 published AF minimum", 0.007, min(_af26), dp=3)
    chk_num("Cluster-scan power", "MUC2/MUC6 published AF maximum", 0.69, max(_af26))
    chk_int("Cluster-scan power", "signals tested", 8, len(_KOUS))
    # Percent form and the two published figures quoted beside it, so that
    # 07_consistency_checks.py can compare the manuscript sentence word for word.
    chk_tol("Cluster-scan power", "rs12802931 power pct", 74,
            100 * _pw["rs12802931"], tol=1)
    chk_num("Cluster-scan power", "rs12802931 published OR", 1.71,
            [_o for _, _r, _, _o, _ in _KOUS if _r == "rs12802931"][0])
    chk_num("Cluster-scan power", "rs12802931 published AF", 0.20,
            [_a for _, _r, _a, _, _ in _KOUS if _r == "rs12802931"][0])

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
    print("✅ All recomputed values reconcile with the expected values held here.")

# Export every computed value so 07_consistency_checks.py can compare the numbers
# actually printed in the manuscript against them. Without this, the checks above
# only verify this script's own copy of the numbers: updating a literal here and
# forgetting the .md leaves the two silently decoupled, which is exactly how the
# N=798 Hardy-Weinberg values survived the ancestry re-derivation.
import json
with open(f"{ROOT}/code/verified_numbers.json", "w") as _fh:
    json.dump({"values": {f"{c['section']}|{c['item']}": c["computed"] for c in _checks},
               "provenance": {f"{c['section']}|{c['item']}":
                              {k: v for k, v in c.items()
                               if k in ("tag", "model", "strata")}
                              for c in _checks if "tag" in c}}, _fh,
              indent=1, sort_keys=True, default=str)
print(f"wrote code/verified_numbers.json ({len(_checks)} values) for the manuscript cross-check")

res

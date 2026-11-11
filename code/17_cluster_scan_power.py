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
# # 17. Power of the cluster-wide scan for the known 11p15.5 IPF signals
#
# `11_mucin_cluster_scan.py` finds no mucin structural variant associated with IPF
# independently of rs35705950. That is a negative, and a reviewer asked for it to be
# backed by a power calculation rather than asserted. This supplies it.
#
# Kousathanas et al. (2025) resolve nine conditionally-independent IPF signals at
# 11p15.5 by GCTA-COJO; eight lie in the mucin cluster (the ninth, rs34787245, is
# intronic in TOLLIP and outside it). Their Supplementary Table 5 gives, for each,
# the European allele frequency and odds ratio. The question this answers is: if a
# signal of that size and frequency were present in IPFJES, would the cluster-wide
# scan have found it?
#
# Method, as stated in the manuscript Methods:
#   1. simulate genotypes under Hardy-Weinberg equilibrium at the published AF_EUR;
#   2. fit the scan model in the 771 analytical participants to obtain the standard
#      error of the log odds ratio at that frequency;
#   3. evaluate power analytically at the cluster-wide threshold alpha = 0.05/54,
#      where 54 is the effective number of independent tests from `11_...` (Li & Ji).
#
# The SE comes from fitting against the real outcome and the real covariate matrix,
# so it carries the actual case/control split, centre structure and collinearity of
# this sample rather than a balanced-design approximation. It is a null-model SE,
# which is the right scale for a power calculation at these effect sizes.
#
# Reads only shared root inputs (`flat_data_genotype_subset.csv`,
# `code/ipfjes_pca.eigenvec`) plus the published table below. Writes
# `tables/cluster_scan_power.csv`.

# %%
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = "/home/drcjar/academic/ipfjes-omics"
RNG = np.random.default_rng(20261031)   # fixed: the printed numbers must reproduce
N_SIM = 200                             # SE averaged over simulated genotypes
M_EFF = 54                              # from 11_mucin_cluster_scan.py (Li & Ji 2005)
ALPHA = 0.05 / M_EFF

# %% [markdown]
# ## The published signals
#
# Kousathanas et al. 2025, Supplementary Table 5 (`Kousathanas2025mmc2.docx`), rows
# for locus 11p15.5, transcribed verbatim: variant ID, rsID, AF_EUR, OR. Nearest
# gene is theirs. rs34787245 (chr11:1,289,602, TOLLIP) is excluded as outside the
# cluster, matching `12_h2_conditional_on_locus_signals.py`.

# %%
SIGNALS = pd.DataFrame([
    # pos,      rsid,          af_eur,  or_,  nearest
    (1013992,  "rs148815783", 0.0073, 2.10, "MUC6"),
    (1057948,  "rs12574439",  0.12,   1.21, "MUC2"),
    (1078249,  "rs7104590",   0.69,   1.15, "MUC2"),
    (1081750,  "rs41359951",  0.018,  1.55, "MUC2"),
    (1208601,  "rs548688097", 0.0070, 2.84, "MUC5B"),
    (1214934,  "rs12802931",  0.20,   1.71, "MUC5B"),
    (1219775,  "rs2672794",   0.61,   1.20, "MUC5B"),
    (1224068,  "rs56321310",  0.19,   1.39, "MUC5B"),
], columns=["pos", "rsid", "af_eur", "or_pub", "nearest_gene"])

# %% [markdown]
# ## The analytical sample and the scan model
#
# Identical construction to `11_mucin_cluster_scan.py`: N=771, centre 16 excluded,
# self-reported white ethnicity, PCs computed within the sample.

# %%
clin = pd.read_csv(f"{ROOT}/flat_data_genotype_subset.csv")
clin["pid"] = clin.participant_id.astype(str)
pcs = pd.read_csv(f"{ROOT}/code/ipfjes_pca.eigenvec", sep=r"\s+")
pcs["pid"] = pcs.IID.astype(str)

dm = (clin
      .merge(pcs[["pid", "PC1", "PC2", "PC3", "PC4", "PC5"]], on="pid")
      .dropna(subset=["age", "ever_smoked", "centre", "case"]))
dm = dm[(dm.centre != 16)
        & dm.ethnicity.astype(str).str.contains("White", case=False, na=False)].copy()
COVS = "age + ever_smoked + PC1+PC2+PC3+PC4+PC5 + C(centre)"
print(f"N={len(dm)} (cases={int(dm.case.sum())}, controls={int((dm.case == 0).sum())}); "
      f"M_eff={M_EFF}, cluster-wide alpha={ALPHA:.2e}")

# %% [markdown]
# ## Power for each signal

# %%
Z_CRIT = stats.norm.isf(ALPHA / 2)


def se_at_frequency(af):
    """Mean SE of the log OR for a variant at this allele frequency in this sample."""
    ses = []
    for _ in range(N_SIM):
        g = RNG.binomial(2, af, size=len(dm)).astype(float)
        if g.std() == 0:                      # monomorphic draw at very low AF
            continue
        d = dm.assign(g=g)
        try:
            m = smf.logit(f"case ~ g + {COVS}", data=d).fit(disp=0)
        except Exception:
            continue
        se = m.bse.get("g", np.nan)
        if np.isfinite(se) and se < 10:        # discard non-converged separations
            ses.append(se)
    return float(np.mean(ses)), len(ses)


rows = []
for _, s in SIGNALS.iterrows():
    se, n_ok = se_at_frequency(s.af_eur)
    ncp = np.log(s.or_pub) / se
    power = stats.norm.sf(Z_CRIT - ncp) + stats.norm.cdf(-Z_CRIT - ncp)
    rows.append(dict(rsid=s.rsid, pos=s.pos, nearest_gene=s.nearest_gene,
                     af_eur=s.af_eur, or_pub=s.or_pub, se_logor=round(se, 4),
                     n_sim_used=n_ok, power=round(power, 4)))

P = pd.DataFrame(rows)
print(P.to_string(index=False))

# %% [markdown]
# ## The two figures quoted in the manuscript

# %%
muc5b = P[P.nearest_gene == "MUC5B"]
muc2_6 = P[P.nearest_gene.isin(["MUC2", "MUC6"])]

# "Strongest MUC5B signal" is the best-powered one, i.e. the most favourable case
# for the scan: quoting the largest published OR instead would overstate power,
# since that variant (OR 2.84) sits at AF 0.007 where power is negligible.
best_muc5b = muc5b.loc[muc5b.power.idxmax()]
print(f"\nBest-powered MUC5B signal: {best_muc5b.rsid} "
      f"(OR {best_muc5b.or_pub}, AF {best_muc5b.af_eur}) -> power {best_muc5b.power:.0%}")
print(f"MUC2/MUC6 signals: power {muc2_6.power.min():.1%} to {muc2_6.power.max():.1%}; "
      f"OR {muc2_6.or_pub.min():.2f} to {muc2_6.or_pub.max():.2f}; "
      f"AF {muc2_6.af_eur.min()} to {muc2_6.af_eur.max()}")
print(f"All eight: power {P.power.min():.1%} to {P.power.max():.1%}")

P.to_csv(f"{ROOT}/tables/cluster_scan_power.csv", index=False)
print("\nwrote tables/cluster_scan_power.csv")

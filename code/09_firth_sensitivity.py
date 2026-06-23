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

"""09_firth_sensitivity.py — Firth penalised-logistic sensitivity analysis.

Reviewer concern: the primary models drop one cases-only centre (centre 16, n=31) to
avoid complete separation, reducing the regression sample from 829 to 798. Excluding a
whole centre removes cases non-randomly. The standard remedy for separation is Firth's
penalised likelihood (Jeffreys prior), which is finite under separation and lets us
RETAIN all 829 participants including centre 16. This script refits the primary MUC5AC H3
models with Firth on the full N=829 and reports whether the conclusions hold.

It also reports the discordant subset (H3-carrier status != rs35705950-carrier status),
broken down by case/control — the long-read sequencing target.

Run: uv run --python .venv python3 code/09_firth_sensitivity.py
"""
import pandas as pd, numpy as np, warnings
from scipy import stats
warnings.filterwarnings('ignore')


def firth_fit(X, y, max_iter=500, tol=1e-9):
    """Firth penalised logistic regression (Jeffreys-prior) via penalised IRLS.
    Finite under complete separation. Returns (beta, Wald se)."""
    n, p = X.shape
    beta = np.zeros(p)
    for _ in range(max_iter):
        pi = 1.0 / (1.0 + np.exp(-(X @ beta)))
        W = pi * (1 - pi)
        I = (X.T * W) @ X
        Iinv = np.linalg.pinv(I)
        h = W * np.einsum('ij,jk,ik->i', X, Iinv, X)   # hat-matrix diagonal
        U = X.T @ (y - pi + h * (0.5 - pi))             # penalised score
        delta = Iinv @ U
        beta = beta + delta
        if np.max(np.abs(delta)) < tol:
            break
    return beta, np.sqrt(np.diag(Iinv))


def firth_or(design_df, term, y):
    cols = ['intercept'] + list(design_df.columns)
    X = np.column_stack([np.ones(len(design_df)), design_df.values.astype(float)])
    beta, se = firth_fit(X, y.astype(float))
    j = cols.index(term); b, s = beta[j], se[j]
    return np.exp(b), np.exp(b - 1.96 * s), np.exp(b + 1.96 * s), 2 * (1 - stats.norm.cdf(abs(b / s)))

ROOT = '.'
pcs  = pd.read_csv(f'{ROOT}/code/ipfjes_pca.eigenvec', sep=r'\s+'); pcs['pid'] = pcs['IID'].astype(str)
clin = pd.read_csv(f'{ROOT}/flat_data_genotype_subset.csv');        clin['pid'] = clin['participant_id'].astype(str)
hap  = pd.read_csv(f'{ROOT}/muc5_haplogroups.csv');                 hap['pid'] = hap['IID'].str.split('_').str[0]

df = (clin
      .merge(hap[['pid','H3_count','rs35705950_T_count','MUC5B_H1_count']], on='pid')
      .merge(pcs[['pid','PC1','PC2','PC3','PC4','PC5']], on='pid')
      .dropna(subset=['age','ever_smoked','peto_exposed','centre']))
n_case = int(df.case.sum()); n_ctrl = int((df.case == 0).sum())
print(f"Full IPFJES analytical sample N={len(df)} (cases={n_case}, controls={n_ctrl})")
c16 = df[df.centre == 16]
print(f"Centre 16 (cases-only, excluded from standard regression): n={len(c16)}, cases={int(c16.case.sum())}")

# --- Discordant subset: H3-carrier status != rs35705950-carrier status (the long-read target) ---
disc = (df.H3_count > 0) != (df.rs35705950_T_count > 0)
dc, dctrl = int(disc[df.case == 1].sum()), int(disc[df.case == 0].sum())
print("\nDiscordant H3/rs35705950 (long-read target):")
print(f"  total   {int(disc.sum())}/{len(df)} ({100*disc.mean():.0f}%)")
print(f"  cases   {dc}/{n_case} ({100*dc/n_case:.0f}%)")
print(f"  controls {dctrl}/{n_ctrl} ({100*dctrl/n_ctrl:.0f}%)")

# --- Firth models on the FULL N=829 (centre as fixed-effect dummies, centre 16 retained) ---
cen  = pd.get_dummies(df['centre'].astype(int), prefix='c', drop_first=True).astype(float)
base = pd.concat([df[['age','ever_smoked','PC1','PC2','PC3','PC4','PC5']].astype(float), cen], axis=1)
y    = df['case'].astype(int).values

print("\nFirth penalised logistic, full N=829 incl. centre 16 (Wald CIs; cf. primary statsmodels N=798):")
X1a = pd.concat([df[['H3_count']].astype(float), base], axis=1)
o, l, h, p = firth_or(X1a, 'H3_count', y)
print(f"  1a  H3 unadjusted   : OR {o:.2f} ({l:.2f}-{h:.2f}) p={p:.2e}   [primary N=798: 1.69 (1.31-2.17), 4.2e-05]")
X1b = pd.concat([df[['H3_count','rs35705950_T_count']].astype(float), base], axis=1)
o, l, h, p = firth_or(X1b, 'H3_count', y)
print(f"  1b  H3 | rs35705950 : OR {o:.2f} ({l:.2f}-{h:.2f}) p={p:.2f}   [primary N=798: 1.17 (0.89-1.54), 0.27]")
o, l, h, p = firth_or(X1b, 'rs35705950_T_count', y)
print(f"  1b  rs35705950 | H3 : OR {o:.2f} ({l:.2f}-{h:.2f}) p={p:.2e}   [primary N=798: 3.75 (2.74-5.13)]")

# Sanity: Firth on the centre-excluded sample (N=798) should reproduce the standard-logit result
d798 = df[df.centre != 16]
cen8 = pd.get_dummies(d798['centre'].astype(int), prefix='c', drop_first=True).astype(float)
base8 = pd.concat([d798[['age','ever_smoked','PC1','PC2','PC3','PC4','PC5']].astype(float).reset_index(drop=True),
                   cen8.reset_index(drop=True)], axis=1)
y8 = d798['case'].astype(int).values
X8 = pd.concat([d798[['H3_count']].astype(float).reset_index(drop=True), base8], axis=1)
o, l, h, p = firth_or(X8, 'H3_count', y8)
print(f"\nSanity — Firth on N=798 (centre 16 excluded): H3 unadjusted OR {o:.2f} ({l:.2f}-{h:.2f}) p={p:.2e}")
print("  (should sit on top of the standard-logit 1.69; confirms Firth adds nothing where there is no separation)")
print("DONE")

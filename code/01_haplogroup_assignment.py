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

import pandas as pd
import numpy as np

df = pd.read_csv('./muc5_tsnps_dosage_pheno.csv')

# ── MUC5AC concordance SNPs ────────────────────────────────────────────────
# rs28542750: T allele = non-H3 (H3 r²=1.0); H3 count = 2 - T_dosage
# rs769768817: C allele = H1  (H1 r²=0.899); H1 count = C_dosage
df['H3_count'] = 2 - df['chr11:1165501:T:A_T']
df['H1_count'] = df['chr11:1195214:G:C_C']
df['H2_count'] = (2 - df['H3_count'] - df['H1_count']).clip(lower=0)

# ── Cross-check: H1 tag-set average (each counted allele = H1 allele) ─────
h1_tag_cols = [
    'chr11:1193830:G:C_G',  # rs2075842
    'chr11:1194354:A:G_A',  # rs1132433
    'chr11:1195265:G:C_C',  # rs1132434
    'chr11:1195858:T:C_T',  # rs28652890
    'chr11:1182048:A:T_T',  # rs879136008
]
df['H1_tagset_avg'] = df[h1_tag_cols].mean(axis=1)   # scale 0–2

# ── Cross-check: H3 tag-set average (counted allele = non-H3; flip) ────────
h3_tag_cols = [
    'chr11:1178216:C:G_C',  # rs36154966
    'chr11:1178023:A:G_A',  # rs1004828576
    'chr11:1177688:A:T_A',  # rs940158763
    'chr11:1177356:T:C_T',  # rs36151150
    'chr11:1177269:T:C_T',  # rs36132281
    'chr11:1197482:C:T_C',  # rs35779873
]
df['H3_tagset_avg'] = 2 - df[h3_tag_cols].mean(axis=1)   # scale 0–2

# Cross-check: correlation between concordance count and tag-set averages
r_h1 = df['H1_count'].corr(df['H1_tagset_avg'])
r_h3 = df['H3_count'].corr(df['H3_tagset_avg'])
print(f"Concordance ↔ tag-set r (H1): {r_h1:.4f}")
print(f"Concordance ↔ tag-set r (H3): {r_h3:.4f}")

# ── H2 complement-allele average (freq~0.19, may be H2-enriched) ──────────
h2_tag_cols = [
    'chr11:1202343:T:C_C',  # rs28519516
    'chr11:1202859:A:G_G',  # rs28558973
    'chr11:1202914:T:C_C',  # rs28368633
]
df['H2_tagset_complement_avg'] = 2 - df[h2_tag_cols].mean(axis=1)

# ── Dominant (carrier) flags ───────────────────────────────────────────────
df['has_H3'] = (df['H3_count'] > 0).astype(int)
df['has_H1'] = (df['H1_count'] > 0).astype(int)
df['has_H2'] = (df['H2_count'] > 0).astype(int)

# ── Diplotype label ────────────────────────────────────────────────────────
def diplotype_label(row):
    h1 = int(round(row['H1_count']))
    h2 = int(round(row['H2_count']))
    h3 = int(round(row['H3_count']))
    alleles = ['H1'] * h1 + ['H2'] * h2 + ['H3'] * h3
    if len(alleles) != 2:
        return 'uncertain'
    return '/'.join(sorted(alleles))

df['MUC5AC_diplotype'] = df.apply(diplotype_label, axis=1)

# ── MUC5B: rs35705950 kept as covariate (T = risk allele; G = ref non-risk)
df['rs35705950_T_count'] = 2 - df['chr11:1219991:G:T_G']   # count of T (risk) alleles

# ── MUC5B H1/H2 tSNP ─────────────────────────────────────────────────────
# Position GRCh38 1244757, variant chr11:1244757:C:G
# C = ref, G = alt; per Plender H2 r²=1 at this position (H2 ~82% in EUR)
# We'll load this separately and merge
import subprocess
result = subprocess.run(
    ['./plink2', '--bfile', 'chr11', '--snp', 'chr11:1244757:C:G',
     '--export', 'A', '--out', 'muc5b_h1h2_tsnp'],
    capture_output=True, text=True, cwd='.'
)
print(result.stdout[-200:] if result.stdout else "")
print(result.stderr[-200:] if result.stderr else "")

muc5b = pd.read_csv('./muc5b_h1h2_tsnp.raw', sep=r'\s+')
muc5b['IID_bare'] = muc5b['IID'].str.split('_').str[0].astype(int)
muc5b_col = [c for c in muc5b.columns if c.startswith('chr11:1244757')][0]
print(f"\nMUC5B H1/H2 tSNP column: {muc5b_col}")
print(f"Counted allele freq: {muc5b[muc5b_col].mean()/2:.3f}")

# Merge into main df
fam = pd.read_csv('./ipfjes_cases_and_controls.fam', sep=r'\s+',
                  header=None, names=['FID','IID','PAT','MAT','SEX','PHENOTYPE'])
df['IID_bare'] = df['IID'].str.split('_').str[0].astype(int)
df = df.merge(muc5b[['IID_bare', muc5b_col]], on='IID_bare', how='left')
# Per brief: H2 r²=1 at 1244757 — identify which allele based on freq vs Plender (H2~82%)
# counted allele freq > 0.5 suggests it marks the common H2 haplotype
muc5b_counted_freq = df[muc5b_col].mean() / 2
if muc5b_counted_freq > 0.5:
    df['MUC5B_H2_count'] = df[muc5b_col]            # counted allele = H2
    df['MUC5B_H1_count'] = 2 - df[muc5b_col]
    muc5b_h2_allele = muc5b_col.split('_')[-1]
    print(f"  → {muc5b_h2_allele} allele (freq={muc5b_counted_freq:.3f}) assigned as H2 (common)")
else:
    df['MUC5B_H2_count'] = 2 - df[muc5b_col]       # complement = H2
    df['MUC5B_H1_count'] = df[muc5b_col]
    muc5b_h2_allele = 'complement'
    print(f"  → complement (freq={1-muc5b_counted_freq:.3f}) assigned as H2 (common)")

df['MUC5B_diplotype'] = df.apply(
    lambda r: '/'.join(sorted(['H1']*int(round(r['MUC5B_H1_count'])) +
                               ['H2']*int(round(r['MUC5B_H2_count'])))), axis=1)

# ── Summary statistics ────────────────────────────────────────────────────
print("\n── MUC5AC haplotype frequencies (concordance SNPs) ──")
h1_f = df['H1_count'].mean() / 2
h2_f = df['H2_count'].mean() / 2
h3_f = df['H3_count'].mean() / 2
print(f"  H1: {h1_f:.3f}  (Plender EUR: most common, expected ~57%)")
print(f"  H2: {h2_f:.3f}  (Plender EUR: intermediate)")
print(f"  H3: {h3_f:.3f}  (Plender EUR: ~21%)")

print("\n── MUC5AC diplotype counts ──")
print(df['MUC5AC_diplotype'].value_counts().to_string())

print("\n── MUC5B haplotype frequencies ──")
h2b_f = df['MUC5B_H2_count'].mean() / 2
h1b_f = df['MUC5B_H1_count'].mean() / 2
print(f"  H1: {h1b_f:.3f}  (Plender EUR: ~18%)")
print(f"  H2: {h2b_f:.3f}  (Plender EUR: ~82%)")

print("\n── MUC5B diplotype counts ──")
print(df['MUC5B_diplotype'].value_counts().to_string())

# ── Tag-set cross-check ────────────────────────────────────────────────────
print("\n── H2 tag-set complement mean vs H2_count (concordance) ──")
r_h2 = df['H2_count'].corr(df['H2_tagset_complement_avg'])
print(f"  r = {r_h2:.4f}")
print(f"  H2 tag complement avg: mean={df['H2_tagset_complement_avg'].mean():.3f}")

# ── Save ─────────────────────────────────────────────────────────────────
out_cols = ['FID', 'IID', 'pheno',
            'H1_count', 'H2_count', 'H3_count',
            'H1_tagset_avg', 'H3_tagset_avg', 'H2_tagset_complement_avg',
            'has_H1', 'has_H2', 'has_H3', 'MUC5AC_diplotype',
            'rs35705950_T_count', 'MUC5B_H1_count', 'MUC5B_H2_count', 'MUC5B_diplotype']
df[out_cols].to_csv('./muc5_haplogroups.csv', index=False)
print("\nSaved: muc5_haplogroups.csv")

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

"""06_tsnp_table.py — Supplementary Table S4: haplogroup-tagging SNP panel.

Emits a markdown table of the tagging SNPs used for haplogroup assignment, with
the counted (effect) allele, its frequency and the minor-allele frequency (MAF)
in the genotyped panel (N=2,912), the tagged haplogroup, and Plender tagging r².
The MAF column documents the strand-ambiguity argument (Methods): palindromic
C/G and A/T tSNPs are resolved by frequency, which is unambiguous when the MAF
is away from 0.5 (e.g. the MUC5B C/G discriminator, MAF 0.125).
"""
import pandas as pd, re, warnings
warnings.filterwarnings("ignore")
R = "."
ad  = pd.read_csv(f"{R}/muc5_allele_direction.csv")
dos = pd.read_csv(f"{R}/muc5_tsnps_dosage_pheno.csv")
mb  = pd.read_csv(f"{R}/muc5b_h1h2_tsnp.raw", sep=r"\s+")

def row(rsid, col, hap, r2):
    pos, ref, alt, counted = re.match(r"chr11:(\d+):([ACGT]+):([ACGT]+)_([ACGT]+)", col).groups()
    src = dos if col in dos.columns else mb
    freq = src[col].mean() / 2
    maf = min(freq, 1 - freq)
    ambiguous = "Y" if {ref, alt} in ({"C", "G"}, {"A", "T"}) else ""
    return f"| {rsid} | chr11:{pos} | {ref}/{alt} | {counted} | {freq:.3f} | {maf:.3f} | {ambiguous} | {hap} | {r2} |"

print("| tSNP | Position (GRCh38) | Alleles | Counted | Freq | MAF | Palindromic | Tags | Plender r² |")
print("|------|-------------------|:--:|:--:|:--:|:--:|:--:|------|:--:|")
for _, r in ad.iterrows():
    m = re.search(r"r2=([\d.]+)", r["evidence"])
    print(row(r["rsID"], r["column"], r["haplogroup"], m.group(1) if m else "—"))
print(row("(MUC5B disc.)", "chr11:1244757:C:G_C", "MUC5B H2 (H1=complement)", "1.0"))

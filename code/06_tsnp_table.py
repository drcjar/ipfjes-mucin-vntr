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

"""06_tsnp_table.py — Supplementary Table S2: haplogroup-tagging SNP panel.

Emits a markdown table of the tagging SNPs used for haplogroup assignment, with
the counted (effect) allele, its frequency and the minor-allele frequency (MAF)
in the 401 IPFJES controls, the tagged haplogroup, and Plender tagging r².

Frequencies are computed in the 401 IPFJES controls, which are the population
reference used throughout the manuscript (cases are enriched for risk haplotypes).
The MAF column documents the strand-ambiguity argument (Methods): palindromic
C/G and A/T tSNPs are resolved by frequency, which is unambiguous when the MAF
is away from 0.5 (e.g. the MUC5B C/G discriminator, MAF 0.146).

The r² values come from the structured `source_r2` column of muc5_allele_direction.csv
(with `source_table`/`source_hap` recording WHICH Plender table each was read from).
They used to be regex-scraped out of the free-text `evidence` column, which meant the
cited table name was never machine-read — that is how rs35779873's r² was mis-attributed
to Table S4 when it actually comes from Table S6. `10_provenance_audit.py` now validates
every (source_table, source_hap, source_r2) triple against ../mmc2.xlsx directly, and
diffs this script's output against the table pasted into the manuscript.
"""
import pandas as pd, re, warnings
warnings.filterwarnings("ignore")
R = "/home/drcjar/ipfjes-omics"
ad  = pd.read_csv(f"{R}/muc5_allele_direction.csv")
dos = pd.read_csv(f"{R}/muc5_tsnps_dosage_pheno.csv")
mb  = pd.read_csv(f"{R}/muc5b_h1h2_tsnp.raw", sep=r"\s+")
h2  = pd.read_csv(f"{R}/muc5ac_h2_tag.raw", sep=r"\s+")   # direct MUC5AC H2 tag (step 3b)

# Restrict to the 401 IPFJES controls (the manuscript's population reference).
_clin = pd.read_csv(f"{R}/flat_data_genotype_subset.csv")
_clin["pid"] = _clin["participant_id"].astype(str)
_pcs = pd.read_csv(f"{R}/code/ipfjes_pca.eigenvec", sep=r"\s+")
_pcs["pid"] = _pcs["IID"].astype(str)
_ctrl = set(_clin.merge(_pcs[["pid"]], on="pid")
                 .dropna(subset=["age", "ever_smoked", "peto_exposed", "centre"])
                 .query("case == 0 and centre != 16")["pid"])
for _df in (dos, mb, h2):
    _df["pid"] = _df["IID"].astype(str).str.split("_").str[0]
dos = dos[dos.pid.isin(_ctrl)]
mb  = mb[mb.pid.isin(_ctrl)]
h2  = h2[h2.pid.isin(_ctrl)]
print(f"<!-- frequencies computed in {len(dos)} IPFJES controls -->")

def row(rsid, col, hap, r2):
    pos, ref, alt, counted = re.match(r"chr11:(\d+):([ACGT]+):([ACGT]+)_([ACGT]+)", col).groups()
    src = dos if col in dos.columns else (mb if col in mb.columns else h2)
    freq = src[col].mean() / 2
    maf = min(freq, 1 - freq)
    return f"| {rsid} | chr11:{pos} | {ref}/{alt} | {counted} | {freq:.3f} | {maf:.3f} | {hap} | {r2} |"

print("| tSNP | Position | Alleles | Counted | Freq | MAF | Tags | r² |")
print("|--------------|---------|:---:|:---:|:--:|:--:|-----------|:--:|")
for _, r in ad.iterrows():
    hap = "*MUC5B* promoter" if r["haplogroup"] == "MUC5B promoter" else r["haplogroup"].replace(" concordance", " conc.")
    # trailing zeros stripped but at least one decimal kept: 1.000 -> "1.0", 0.899 -> "0.899"
    r2 = "n/a" if pd.isna(r["source_r2"]) else f'{r["source_r2"]:.3f}'.rstrip("0").ljust(3, "0")
    print(row(r["rsID"], r["column"], hap, r2))
print(row("MUC5B discriminator", "chr11:1244757:C:G_C", "MUC5B H2", "1.0"))

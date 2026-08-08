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
# # 15. LD between MUC5AC/MUC5B haplogroups and rs35705950 (Supplementary Table S7)
#
# Requested by a collaborator: how does the LD between each haplogroup and the
# MUC5B promoter variant rs35705950 in IPFJES compare with a reference panel?
#
# Each haplogroup is represented by its concordance tSNP (the tag used to assign
# it): H1 rs769768817, H2 rs879201833, H3 rs28542750, and the MUC5B H1/H2
# discriminator chr11:1,244,757. LD is the Pearson r (and r²) of the tag allele
# count with the rs35705950 T-allele count, computed in the 401 IPFJES controls
# (the population reference, cases being enriched for risk haplotypes) and in the
# full 798.
#
# **Reference column.** Plender et al. do not tabulate haplogroup–rs35705950 LD.
# The reference is the r² between rs35705950 and each tSNP in the 1000 Genomes
# European panel (phase 3) — the panel Plender used for their locus-wide LD
# analysis — obtained live from the Ensembl LD API and recorded below (the API
# does not return pairs whose r² is below ~0.05, so H1/H2/MUC5B come back empty,
# i.e. negligible). Query template:
#   https://rest.ensembl.org/ld/human/pairwise/rs35705950/<tSNP>?population_name=1000GENOMES:phase_3:EUR
# REF_1000G_EUR = {"H3 (rs28542750)": 0.104,        # returned
#                  "H1 (rs769768817)": None,         # empty -> r2 < ~0.05 (negligible)
#                  "H2 (rs879201833)": None,
#                  "MUC5B H1 (chr11:1244757)": None}

# %%
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = "/home/drcjar/ipfjes-omics"

clin = pd.read_csv(f"{ROOT}/flat_data_genotype_subset.csv"); clin["pid"] = clin.participant_id.astype(str)
hap = pd.read_csv(f"{ROOT}/muc5_haplogroups.csv"); hap["pid"] = hap.IID.str.split("_").str[0]
raw = pd.read_csv(f"{ROOT}/muc5ac_h2_tag.raw", sep=r"\s+"); raw["pid"] = raw.IID.astype(str).str.split("_").str[0]
H2COL = [c for c in raw.columns if c.startswith("chr11:1182279")][0]     # rs879201833, the primary H2 tag

d = (clin
     .merge(hap[["pid", "H1_count", "H3_count", "MUC5B_H1_count", "rs35705950_T_count"]], on="pid")
     .merge(raw[["pid", H2COL]], on="pid")
     .dropna(subset=["age", "ever_smoked", "peto_exposed", "centre"]))
d = d[d.centre != 16].copy()
_f = d.loc[d.case == 0, H2COL].mean() / 2
d["H2t"] = d[H2COL] if _f < 0.5 else 2 - d[H2COL]

# (label, dosage column, 1000G-EUR reference r2 from the Ensembl LD API; None = below
#  the API's returned-LD threshold ~0.05, i.e. negligible)
HAPS = [("MUC5AC H1 (rs769768817)", "H1_count", None),
        ("MUC5AC H2 (rs879201833)", "H2t", None),
        ("MUC5AC H3 (rs28542750)", "H3_count", 0.104),
        ("MUC5B H1 (chr11:1,244,757)", "MUC5B_H1_count", None)]

controls = d[d.case == 0]
print(f"{'Haplogroup (tSNP)':30s} {'controls r (r2)':>18s} {'full r (r2)':>16s}   1000G EUR r2")
for lab, col, ref in HAPS:
    rc = np.corrcoef(controls[col], controls.rs35705950_T_count)[0, 1]
    rf = np.corrcoef(d[col], d.rs35705950_T_count)[0, 1]
    ref_s = f"{ref:.2f}" if ref is not None else "<0.05"
    print(f"{lab:30s} {rc:+.2f} ({rc**2:.2f})     {rf:+.2f} ({rf**2:.2f})     {ref_s}")
print("\nOnly MUC5AC H3 has appreciable LD; IPFJES-controls r2 (0.09) matches the 1000G EUR")
print("reference (0.10). Full-sample r2 is larger because cases are enriched for risk haplotypes.")

# Reproducibility — exact commands

Every genotype-handling command used downstream of the imputed data, in order, **verbatim
from the PLINK `.log` files** (so they are guaranteed to match what was run). Together with the
Python scripts (`01`–`09`, see `README.md`) and the three verification layers, this reproduces
every number in the manuscript.

## Software
| Tool | Version | Used for |
|------|---------|----------|
| PLINK | v1.90b6.26 (2 Apr 2022) | single-variant PED/MAP extraction (`--recode`) |
| PLINK2 | v2.00a3.7LM (24 Oct 2022) | LD-pruning, PCA, dosage export (`--export A`), HWE/missingness |
| Python | 3.12 (`.venv`) | analysis: `numpy`, `scipy`, `pandas`, `statsmodels`, `pandas-plink` |

Binaries are invoked as `./plink` / `./plink2` from the repo root (not on `PATH`).

## Upstream (NOT reproduced here)
Genotyping (Affymetrix UK Biobank Axiom, GRCh37), TOPMed imputation (GRCh38), and imputation QC
(variants with r²<0.5 or minor-allele-count ≤3 removed) were performed by the consortium pipeline
(Chin et al., *Eur Respir J* 2026). The inputs below are the **products** of that pipeline.

## Inputs
- `ipfjes_cases_and_controls.{bed,bim,fam}` — typed-array fileset, 862,300 variants × 2,912 samples (for PCA).
- `chr11.{bed,bim,fam}` — TOPMed-imputed chr11, 14,563,388 variants × 2,912 samples (for tSNP/promoter extraction).
- `tsnps_imputed_ids.txt` — 17 variant IDs: the 16 MUC5AC tSNPs + rs35705950 (the *MUC5B* promoter). NB the MUC5B H1/H2 discriminator (`chr11:1244757`) is **not** in this file — it is extracted in its own call (step 3).
- `controls_imputed.txt` — the 2,481 control sample IDs (for control-only HWE).

## Commands (in order)

### 1. Genetic principal components (from typed-array data)
```bash
# 1a. LD-prune  (862,298 autosomal variants -> 335,255 retained)
./plink2 --bfile ipfjes_cases_and_controls \
         --indep-pairwise 500 50 0.2 \
         --out pruned_for_pca

# 1b. PCA on the pruned set  -> ipfjes_pca.eigenvec / .eigenval  (first 5 PCs used in all models)
./plink2 --bfile ipfjes_cases_and_controls \
         --extract pruned_for_pca.prune.in \
         --pca 10 \
         --out ipfjes_pca
```

### 2. MUC5AC tagging-SNP dosages (from imputed chr11)
```bash
# 16 MUC5AC tSNPs + rs35705950 -> additive (0/1/2) dosage matrix muc5_tsnps_dosage.raw
./plink2 --bfile chr11 \
         --extract tsnps_imputed_ids.txt \
         --export A \
         --out muc5_tsnps_dosage
```

### 3. MUC5B H1/H2 discriminator (from imputed chr11)
```bash
# chr11:1,244,757 C/G (GRCh38) -> muc5b_h1h2_tsnp.raw   (also called inside 01_haplogroup_assignment.py)
# This tSNP is in Supplementary Table S4; its control HWE is the MUC5B H1 haplogroup test
# (p=0.83, Results), not the step-4 batch — so "all 17 tSNPs passed HWE" in the paper holds.
./plink2 --bfile chr11 \
         --snp chr11:1244757:C:G \
         --export A \
         --out muc5b_h1h2_tsnp
```

### 4. Tag-SNP QC — Hardy-Weinberg and missingness in controls
```bash
# HWE (chi-squared) + per-variant missingness for the 16 MUC5AC tSNPs + rs35705950, restricted to the 2,481 controls
./plink2 --bfile chr11 \
         --extract tsnps_imputed_ids.txt \
         --keep controls_imputed.txt \
         --hardy \
         --missing variant-only \
         --out muc5_qc
```

### 5. rs35705950 (MUC5B promoter) genotype (imputed)
```bash
# chr11:1,219,991 G/T (GRCh38) = rs35705950  -> muc5b.ped / .map
./plink --bfile chr11 \
        --snp chr11:1219991:G:T \
        --recode \
        --out muc5b
```

## Downstream (Python)
First, `00_build_dosage_pheno.py` assembles `muc5_tsnps_dosage_pheno.csv` from the step-2
`muc5_tsnps_dosage.raw` and the case/control phenotype in `ipfjes_cases_and_controls.fam` (the
`chr11` fileset is unphenotyped, so its PHENOTYPE column is all -9). Then the `.raw` / `.eigenvec`
outputs feed the analysis scripts; run them in order per `README.md`
(`00_build_dosage_pheno.py` → `09_firth_sensitivity.py`), then the verification layers
(`verify_numbers.py`, `07_foundational_audit.py`, `08_consistency_checks.py`). TaqMan rs35705950
genotypes and clinical covariates come from `flat_data_genotype_subset.csv` (provenance in `README.md`).

> The `.raw`, `.ped`, and `.eigenvec` products are individual-level and are **gitignored** (IPFJES
> governance); the commands above regenerate them from the source filesets.

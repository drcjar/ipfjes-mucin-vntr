# Reproducibility — exact commands

Every genotype-handling command used downstream of the imputed data, in order, reconciled
against the PLINK `.log` files. Together with the Python scripts (`00`–`14`, see `README.md`)
and the four verification layers, this reproduces every number in the manuscript.

## Software
| Tool | Version | Used for |
|------|---------|----------|
| PLINK | v1.90b6.26 (2 Apr 2022) | single-variant PED/MAP extraction (`--recode`) |
| PLINK2 | v2.00a3.7LM (24 Oct 2022) | LD-pruning, PCA, dosage export (`--export A`), HWE/missingness |
| Python | 3.12 (`.venv`) | analysis: `numpy`, `scipy`, `pandas`, `statsmodels`, `pandas-plink` |

Binaries are invoked as `./plink` / `./plink2` from the repo root (not on `PATH`).

## Upstream (NOT reproduced here)
Genotyping (Biobank Affymetrix Axiom array, GRCh37), TOPMed imputation (GRCh38), and imputation QC
(variants with r²<0.5 or minor-allele-count ≤3 removed) were performed by the consortium pipeline
(Chin et al., *Eur Respir J* 2026). The inputs below are the **products** of that pipeline.

## Inputs
- `ipfjes_cases_and_controls.{bed,bim,fam}` — typed-array fileset, 862,300 variants (for PCA; restricted to the IPFJES analytical sample via `--keep`).
- `chr11.{bed,bim,fam}` — TOPMed-imputed chr11, 14,563,388 variants (for tSNP/promoter extraction; the build step restricts to IPFJES).
- `tsnps_imputed_ids.txt` — 17 variant IDs: the 16 MUC5AC tSNPs + rs35705950 (the *MUC5B* promoter). NB the MUC5B H1/H2 discriminator (`chr11:1244757`) is **not** in this file — it is extracted in its own call (step 3).
- `flat_data_genotype_subset.csv` — clinical covariates; defines the 798-participant analytical
  sample and the 401 IPFJES controls used for all frequency and Hardy-Weinberg reporting.
- `mmc2.xlsx` — Plender et al. (2024) supplementary workbook (sheets S1–S21). The **primary source**
  for every tagging r² we quote: sheet S4 (92 MUC5AC tSNPs) and sheet S6 (eQTL correlations, the
  source of rs35779873's H3 r²=0.805). `08_provenance_audit.py` reads it directly and fails if any
  citation in `muc5_allele_direction.csv` or the manuscript names the wrong sheet.
- `muc5_tsnps.txt` — the human-readable rsID candidate list (documentation only; PLINK reads the
  position-based `tsnps_imputed_ids.txt`). One candidate, rs1015856541, is commented out: it is in
  neither Plender table and was not present in the imputed fileset.

## Commands (in order)

### 1. Genetic principal components (from typed-array data)
Both steps are restricted to the 798 analytical participants via `--keep`, so that no external
control sample contributes to the ancestry adjustment. `keep798.txt` is a two-column FID/IID list
of the analytical sample. It contains participant identifiers, so it is **not committed**; build it
first with:

```python
# uv run --python .venv python3 - <<'EOF'
import pandas as pd
clin = pd.read_csv("flat_data_genotype_subset.csv"); clin["pid"] = clin.participant_id.astype(str)
hap  = pd.read_csv("muc5_haplogroups.csv");          hap["pid"]  = hap.IID.str.split("_").str[0]
d = (clin.merge(hap[["pid"]], on="pid")
         .dropna(subset=["age", "ever_smoked", "peto_exposed", "centre"]))
d = d[d.centre != 16]                                  # cases-only centre
pd.DataFrame({"FID": d.pid, "IID": d.pid}).to_csv(
    "keep798.txt", sep=" ", header=False, index=False)  # 798 rows
# EOF
```

Note this bootstraps from the haplogroup calls (step 2/3 output) rather than from the eigenvector
file, so it does not depend on the PCA it is used to produce.

```bash
# 1a. LD-prune within the 798  (862,298 autosomal variants -> 314,996 retained)
./plink2 --bfile ipfjes_cases_and_controls \
         --keep keep798.txt \
         --indep-pairwise 500 50 0.2 \
         --out pruned_for_pca

# 1b. PCA on the pruned set  -> ipfjes_pca.eigenvec / .eigenval  (first 5 PCs used in all models)
./plink2 --bfile ipfjes_cases_and_controls \
         --keep keep798.txt \
         --extract pruned_for_pca.prune.in \
         --pca 10 \
         --out ipfjes_pca
```

**Note on the PCA sample.** The eigenvectors are computed **within the 798 IPFJES analytical participants** via `--keep`, so the ancestry adjustment uses only the analytical sample. An earlier version computed the PCs on a wider fileset and subset the eigenvectors afterwards; that version is superseded. Refitting under the within-sample PCs changed no substantive result — the largest
movement in any reported estimate was in the third significant figure — but the eigenvalue
spectrum flattens appreciably (PC1 9.05 → 5.71, PC2 5.05 → 2.47), consistent with the leading
axes of the old PCA reflecting the IPFJES-versus-external-control contrast rather than
structure within the analytical sample.

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
# This tSNP is given in Plender's MAIN TEXT (H1 r2=0.0026 vs H2 r2=1), not their Table S4.
# Its control HWE is the MUC5B H1 haplogroup test (p=0.83, Results); it is one of the 18
# variants covered by the minimum-HWE figure reported in step 4.
./plink2 --bfile chr11 \
         --snp chr11:1244757:C:G \
         --export A \
         --out muc5b_h1h2_tsnp
```

### 3b. Direct MUC5AC H2 tag (from imputed chr11)
```bash
# rs879201833, chr11:1,182,279 A/G (GRCh38). Plender Table S4: H1 r2=0.662, H2 r2=0.888,
# H3 r2=0.074 — the best direct H2 tag available in our imputed data. Used by
# 09_h2_direct_and_lrt.py to test whether the residual H2 association is an artefact of
# assigning H2 by subtraction. It is not: see that script.
# The three H2 tSNPs already in the panel (rs28519516, rs28558973, rs28368633; H2 r2=0.862)
# come from step 2 and need no separate extraction.
./plink2 --bfile chr11 \
         --snp chr11:1182279:A:G \
         --export A \
         --out muc5ac_h2_tag
```

### 3c. Kousathanas independent 11p15.5 signals (for H2 conditioning)
The eight conditionally-independent 11p15.5 IPF signals from Kousathanas et al. (2025),
Supplementary Table S5 — MUC6 `chr11:1013992:C:T`; MUC2 `chr11:1057948:G:C`,
`chr11:1078249:C:T`, `chr11:1081750:C:T`; MUC5B `chr11:1208601:C:T`, `chr11:1214934:A:G`,
`chr11:1219775:T:C`, `chr11:1224068:C:T`. Used by `12_h2_conditional_on_locus_signals.py` to
test whether the MUC5AC H2 residual survives adjustment for every independent signal at the
locus, not just rs35705950. The variant IDs are listed in `kous_independent_signals.txt`.
```bash
./plink2 --bfile chr11 \
         --extract kous_independent_signals.txt \
         --export A \
         --out kous_independent_signals
```

### 4. Tag-SNP QC — Hardy-Weinberg in the IPFJES controls

Not a PLINK step. Tag-SNP Hardy-Weinberg equilibrium is computed in the **401 IPFJES
controls** directly from the exported dosages, in `verify_numbers.py` (section 5), which
reports the minimum p across all 18 variants (0.41). Controls are used as the population
reference for every frequency and HWE figure in the manuscript.

All frequency and Hardy-Weinberg reporting uses the 401 IPFJES controls (the population reference).

Genotypes are hard-called from the imputed fileset, so per-variant missingness is zero by
construction and a call-rate threshold is not informative.

### 5. rs35705950 (MUC5B promoter) genotype (imputed)
```bash
# chr11:1,219,991 G/T (GRCh38) = rs35705950  -> muc5b.ped / .map
./plink --bfile chr11 \
        --snp chr11:1219991:G:T \
        --recode \
        --out muc5b
```

### 6. Imputed VNTR repeat length (ImputeSV server, Bai et al. panel)

Repeat lengths were imputed by the ImputeSV server (`https://yanglab.westlake.edu.cn/ImputeSV/`)
rather than locally, so this step is a data transfer, not a command we run. IPFJES governance
(17/EM/0021) permits it; only pseudonymised genotypes leave the machine. Analysed by
`10_imputed_vntr_length.py`.

```bash
# Build the submission VCF: 798 analytical participants, whole chr11, polymorphic
# biallelic SNPs only, participant IDs replaced by IPFJES_NNNN pseudonyms.
#   rename798.txt   -- 4-column PLINK --update-ids file (old FID/IID -> pseudonym)
#   imputesv_id_map.tsv -- the pseudonym <-> participant_id map. NEVER leaves this machine
#                          and is gitignored, as are all files below.
./plink2 --bfile chr11 --keep keep798_chr11.txt \
         --snps-only just-acgt --max-alleles 2 --mac 1 \
         --make-bed --out tmp_chr11_full          # 1,241,995 of 14,563,388 variants

./plink2 --bfile tmp_chr11_full --update-ids rename798.txt \
         --output-chr chrM --export vcf id-paste=iid bgz \
         --out imputesv_chr11_full                # 109 MB, chr-prefixed GRCh38, GT only
```

Server settings: **phasing = SHAPEIT5** (our export is unphased; PLINK 1 binary cannot store
phase), **pipeline = VNTR** (the SV panel does not contain the mucin coding repeats), **build =
GRCh38**. The server runs BCFtools QC → SHAPEIT5 v5.0.1 → Minimac4 v4.1.6 and returns
password-protected `chr11.VNTR.GRCh38.noSNPs.zip` (`.dose.vcf.gz` with GT/DS/GP/HDS, plus a
`.sites_only.vcf.gz` carrying the per-allele `R2`, and an md5 manifest — both verified OK).

> **Do not submit a sub-chromosomal window.** Minimac chunks the whole chromosome at 20 Mb and
> aborts on the first chunk containing no target variants. An initial submission of
> chr11:1–6 Mb imputed chunk 1 successfully and then died on chunk `20,000,001–40,000,000`
> ("not enough target variants are available to impute this chunk"), producing no output.
> Dropping monomorphic sites (`--mac 1`) removes 12.3 M uninformative variants and keeps the
> whole-chromosome file to ~109 MB.

Two analyses read the returned dose VCFs (`chr11.{VNTR,SV}.GRCh38.noSNPs/…dose.vcf.gz`, both
gitignored):
- `10_imputed_vntr_length.py` — MUC5AC exon-31 repeat length vs IPF (Supplementary Table S6).
- `11_mucin_cluster_scan.py` — agnostic conditional scan of all 109 common, well-imputed VNTR/SV
  alleles across the mucin cluster (MUC6–MUC5B, chr11:1.01–1.263 Mb), marginally and conditional
  on rs35705950, with a Li & Ji (2005) effective-tests threshold (Supplementary Table S7). Gene
  coordinates are GRCh38 (Ensembl); pseudonyms map back through `imputesv_id_map.tsv` (local).

## External figures cited in the Discussion

One set of figures in the Discussion is not computed from IPFJES data and so cannot be
reproduced by `verify_numbers.py`. It is recorded here instead.

### Bai et al. 2026 (imputed VNTR repeat number in UK Biobank)

```bash
# UK Biobank VNTR GWAS, phecode 502 "Postinflammatory pulmonary fibrosis"
#   347 cases, N=349,660, GLM.  Index of all 2,624 traits:
#     https://lab-storage.oss-cn-hangzhou.aliyuncs.com/pheweb_data/VNTR_data/UKB_EUR_VNTR_associations.csv
curl -sO https://lab-storage.oss-cn-hangzhou.aliyuncs.com/pheweb_data/VNTR_data/Binary/502_PheCode.VNTR.v1.0.tsv.gz

# The MUC5AC exon-31 loci with the 24 bp MUC5AC repeat unit:
#   11:1184917-1187696  BETA  0.0053  SE 0.0117  P 0.65
#   11:1189153-1191907  BETA -0.0082  SE 0.0415  P 0.84
# Also present, and null: two further exon-31 loci (96 bp and 48 bp periods), and
#   11:1255547-1256136  period 59  P 0.85 — the MUC5B intron-37 minisatellite of
#   Desseyn et al. 1999, not otherwise discussed in the manuscript.
# Note the VNTR panel lives under pheweb_data/VNTR_data/, NOT pheweb_data/SV_data/;
# the SV panel contains only smaller indels at this locus and not the coding repeats.
zcat 502_PheCode.VNTR.v1.0.tsv.gz \
  | awk -F'\t' 'NR==1 || (($3==11) && $4>1150000 && $4<1270000)'
```

## Downstream (Python)
First, `00_build_dosage_pheno.py` assembles `muc5_tsnps_dosage_pheno.csv` from the step-2
`muc5_tsnps_dosage.raw` and the case/control phenotype in `ipfjes_cases_and_controls.fam` (the
`chr11` fileset is unphenotyped, so its PHENOTYPE column is all -9). Then the `.raw` / `.eigenvec`
outputs feed the analysis scripts; run them in order per `README.md`
(`00_build_dosage_pheno.py` onward), then the verification layers
(`verify_numbers.py`, `06_foundational_audit.py`, `07_consistency_checks.py`, `08_provenance_audit.py`). Clinical covariates come from `flat_data_genotype_subset.csv` (provenance in `README.md`).

> The `.raw`, `.ped`, and `.eigenvec` products are individual-level and are **gitignored** (IPFJES
> governance); the commands above regenerate them from the source filesets.

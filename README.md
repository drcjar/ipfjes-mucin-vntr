# MUC5AC and MUC5B VNTR haplogroups and IPF — IPFJES tagging-SNP analysis

Analysis code and manuscript source for:

> **MUC5AC and MUC5B VNTR haplogroups are not separable from the promoter variant
> rs35705950 in idiopathic pulmonary fibrosis using array-based data: a tagging-SNP
> analysis of the IPF Job Exposures Study (IPFJES).**

A tagging-SNP test of whether MUC5AC (H1/H2/H3, length-defined) and MUC5B (H1/H2,
composition-defined) VNTR haplogroups associate with IPF independently of the *MUC5B*
promoter variant rs35705950, in the IPFJES occupational case–control cohort.

- 📄 Manuscript source: [`MUC5_VNTR_haplogroup_paper_draft.md`](MUC5_VNTR_haplogroup_paper_draft.md)
- 🔬 Exact genotype-handling commands: [`code/REPRODUCE.md`](code/REPRODUCE.md)
- 🧪 Code/data inventory and run instructions: [`code/README.md`](code/README.md)
- 🧪 Preprint: _(medRxiv DOI to be added)_

## Data availability

**No individual-level data are in this repository.** The IPFJES genotype and phenotype data are controlled-access under IPFJES research governance (REC 17/EM/0021) and cannot be redistributed. Third-party
publication PDFs are likewise not included (copyright).

This repository contains **only** analysis code, the manuscript source, and
aggregate / variant-level outputs (e.g. tagging-SNP lists, per-variant QC, PCA
eigenvalues). Scripts that consume the controlled individual-level files (e.g.
`muc5_haplogroups.csv`, `flat_data_genotype_subset.csv`, the genotype filesets) will
not run without access to that data, which must be obtained separately and placed at
the repository root.

## Reproducing the analysis

The pipeline runs in numbered order from the source genotype filesets:

| Step | Script | Purpose |
|------|--------|---------|
| 00 | `code/00_build_dosage_pheno.py` | Build tagging-SNP dosage + phenotype table |
| 01 | `code/01_haplogroup_assignment.py` | Assign MUC5AC/MUC5B haplogroups from tSNPs |
| 02 | `code/02_primary_analysis.py` | Primary logistic-regression associations |
| 03 | `code/03_supplementary_analyses.py` | Stratified, interaction, bootstrap, power |
| 05 | `code/05_h2ref_sensitivity.py` | Reference-category sensitivity |
| 06 | `code/06_tsnp_table.py` | Tagging-SNP table |
| 07 | `code/07_foundational_audit.py` | Label/orientation audit |
| 08 | `code/08_consistency_checks.py` | Manuscript internal-consistency checks |
| 10 | `code/10_provenance_audit.py` | Plender r² provenance audit |
| 11 | `code/11_h2_direct_and_lrt.py` | Direct H2 tag + likelihood-ratio tests |
| 12 | `code/12_imputed_vntr_length.py` | Imputed VNTR repeat length (Table S5) |
| 13 | `code/13_mucin_cluster_scan.py` | Mucin-cluster conditional scan (Table S6) |
| 14 | `code/14_h2_conditional_on_locus_signals.py` | H2 conditioned on locus signals |
| 15 | `code/15_ld_haplogroup_rs35705950.py` | Haplogroup–rs35705950 LD (Table S7) |
| 16 | `code/16_ancestry_and_promoter_gxe.py` | Ancestry sensitivity + gene×env (Table S4) |
| —  | `code/verify_numbers.py` | Recompute every number cited in the manuscript |

- Exact PLINK / PLINK2 commands used to extract genotypes: [`code/REPRODUCE.md`](code/REPRODUCE.md)
- Code/data inventory and run instructions: [`code/README.md`](code/README.md)
- Paths in the scripts are relative to the repository root.

### Software

- PLINK v1.90b6.26 and PLINK2 v2.00a3.7 (genotype extraction / QC / PCA)
- Python 3.12 with `pandas-plink`, `pandas`, `numpy`, `statsmodels`, `matplotlib`
- Build the manuscript PDF/DOCX with `./build_paper.sh` (requires `pandoc`)

## Citation

If you use this code, please cite the manuscript above (preprint DOI to follow).

## Licence

Code is released under the GNU General Public License v3.0 — see [`LICENSE`](LICENSE).

# MUC5AC/MUC5B VNTR haplogroup analysis — code & data provenance

Analysis pipeline for the manuscript `../MUC5_VNTR_haplogroup_paper_draft.md`.
This `code/` directory is the **source of truth for every number** in the paper.
Exact PLINK/PLINK2 commands (PCA, tSNP extraction, QC) are in **`REPRODUCE.md`**.

## Pipeline (run in order, from the repo root, via the project venv)

```bash
uv run --python .venv python3 code/00_build_dosage_pheno.py      # .raw + .fam phenotype -> muc5_tsnps_dosage_pheno.csv
uv run --python .venv python3 code/01_haplogroup_assignment.py   # tSNP dosages -> haplogroup calls
uv run --python .venv python3 code/02_primary_analysis.py        # Table 1, Table 2, LD, interactions
uv run --python .venv python3 code/03_supplementary_analyses.py  # Tables S1/S2, HWE, bootstrap, power
```

Before sending the manuscript anywhere, run all four verification layers; `10_provenance_audit.py`
exits non-zero on failure, so it can gate a build:

```bash
uv run --python .venv python3 code/verify_numbers.py
uv run --python .venv python3 code/07_foundational_audit.py
uv run --python .venv python3 code/08_consistency_checks.py
uv run --python .venv python3 code/10_provenance_audit.py
```

## Verification & sensitivity

- `verify_numbers.ipynb` (+ paired `verify_numbers.py`) — **number audit**: recomputes every
  manuscript number from source data (192 checks (no flags; the bootstrap/power section was removed)).
- `07_foundational_audit.py` — **assignment/orientation audit** (11 checks): re-derives every
  haplogroup-tag, allele-orientation, and case-coding decision and checks it against independent
  references (gnomAD, Plender, the OEM paper). Catches a mislabel the number audit would not — it
  trusts the labels.
- `05_h2ref_sensitivity.py` — MUC5AC joint model under H2/H3 reference (Supplementary Table S3).
- `06_tsnp_table.py` — per-tSNP allele frequencies/MAFs (Supplementary Table S3).
- `11_h2_direct_and_lrt.py` — three post-review analyses reported in the Discussion:
  refits the joint model with four direct MUC5AC H2 tagging SNPs (all give a larger
  residual than the subtraction estimate in Table 2, so the measurement-error
  explanation for the surviving H2 signal fails); nested-model likelihood-ratio tests
  for the conditional asymmetry; and rs35705950 estimated within H3-carrier strata.
  Needs `muc5ac_h2_tag.raw` from REPRODUCE.md step 3b.
- `08_consistency_checks.py` — **internal-consistency audit** (re-runnable): checks the manuscript
  prose against itself — additive cohort sums, fraction↔% agreement, Table 1 cell arithmetic, CI
  ordering/bracketing of every OR, OR-quote drift, and CI-vs-p significance. A third layer that
  catches edit-drift neither the number audit nor the foundational audit would (they don't compare
  statements to each other). The cohort sums are hard-coded to the current data freeze.
- `10_provenance_audit.py` — 67 checks, exits non-zero on
  failure): the only layer that opens the cited document. Reads `../mmc2.xlsx` directly and
  checks that every Plender r² we quote is in the table we say it is in, at the position we say,
  for the haplogroup we say; that our panel counts match Plender's own (92 MUC5AC tSNPs = 35 H1 /
  5 H2 / 52 H3); and that the tagging-SNP table pasted into the manuscript (now Supplementary Table S3) still matches
  `06_tsnp_table.py`'s output cell-for-cell.

**Why layer 4 exists.** Layers 1–3 check numbers against *our* data, labels against *our* logic,
and prose against *other prose*. None of them could see a wrong **citation**, because provenance
lived in a free-text `evidence` string in `../muc5_allele_direction.csv` from which
`06_tsnp_table.py` regex-scraped only the number — the table name was never read by any program,
so nothing could contradict it. It was wrong for rs35779873 (r²=0.805 is from Table S6, not S4)
and the error reached the manuscript. Provenance is now three structured columns
(`source_table`, `source_hap`, `source_r2`) that are machine-checked against the workbook.
**If you add a tSNP or cite a new r², fill those columns and re-run layer 4 — do not put the
source in prose only.**

The tSNP→haplogroup assignment is sealed by **Plender Table S4** (`../mmc2.xlsx`, sheet S4;
the three numeric columns are r² with H1/H2/H3): rs28542750→H3 r²=1.0, rs769768817→H1 r²=0.899,
rs2075842→H1 r²=0.923 — matching the code.

**Known false alarm (flagged by reviewers twice — NOT an error here):** Plender's *main text*
(p.1709) writes "H1 vs. H2/H3 tSNP: rs28542750, H3 vs. H1/H2 tSNP: rs769768817", which *transposes*
the two labels relative to their own Table S4. We follow Table S4 (the quantitative source). Allele
frequency confirms it independently: rs28542750 MAF 0.197 in the 401 IPFJES controls ≈ the ~21% European H3 frequency, and a
tag with r²=1.0 must match its haplogroup's frequency — so rs28542750 tags H3, not H1 (~50%).
`07_foundational_audit.py` enforces this (SWAP GUARD: H1 > H3).

**Second known false alarm (same class):** comparing Plender's S4 and S6 *by position* suggests the
H1 concordance tag is worthless. S6 lists rs28415193 (G/A) at chr11:1,195,214 with H1 r²=0.008; S4
lists rs769768817 at the same position with H1 r²=0.899. These are distinct variants at one
coordinate — the variant extracted here is `chr11:1195214:G:C` (G/C), matching rs769768817, not
rs28415193's G/A. Related: **rs35779873 is not in Table S4 at all**; its H3 r²=0.805 comes from
Table S6, where Plender designate it a MUC5AC H3 tSNP (as they also do in their PheWAS, Table S15).
`../muc5_allele_direction.csv` and the manuscript's Table S3 footnote both record that provenance.
Of Plender's 92 MUC5AC tSNPs (35 H1 / 5 H2 / 52 H3) we use 17; of their 1 MUC5B tSNP we use 1.

## Inputs and their provenance

| File | Provenance | In git? |
|------|-----------|---------|
| `../flat_data_genotype_subset.csv` | Clinical covariates (age, smoking, pack-years, Peto asbestos exposure, centre, ethnicity) + participant IDs for the IPFJES analytical sample. | **No** — individual-level participant data under IPFJES governance (ethics ref 17/EM/0021); gitignored, never commit. |
| `../muc5_tsnps_dosage_pheno.csv` | **Built by `00_build_dosage_pheno.py`**: the `plink2 --export A` dosages (`muc5_tsnps_dosage.raw`) joined to the case/control phenotype from `ipfjes_cases_and_controls.fam` (the imputed `chr11` fileset is unphenotyped). | No (gitignored, large). |
| `../muc5_haplogroups.csv` | Output of `01_haplogroup_assignment.py` (MUC5AC H1/H2/H3 + MUC5B H1/H2 counts, rs35705950 T count). | **No** — per-participant genotype-derived counts; gitignored (governance). |
| `ipfjes_pca.eigenvec` | Genetic principal components from LD-pruned typed array data (`plink2 --pca`). | **No** — per-participant PCs; gitignored (governance). |
| `../chr11.bed/.bim/.fam`, `./plink2` | Per-chromosome **TOPMed-imputed** genotypes (GRCh38, ~14.5M variants) + PLINK2, used by step 01. Also the candidate LD reference for any future conditional/fine-mapping work. | No (gitignored, large). |
| `../mmc2.xlsx` | Plender et al. 2024 supplementary tables (sheets S1–S21); **S4 = tSNP r²** (source of haplogroup assignment), S6 = allele designations. | **No** — copyright; gitignored. |

## Notes

- Step 01 calls `plink2` to extract the MUC5B H1/H2 tSNP (`chr11:1244757:C:G`) from `chr11`.
- All reported analyses, descriptive and regression alike, use the single N=798 analytical
  sample (one cases-only centre, centre 16, excluded). Principal components are computed
  within those 798 IPFJES participants.
- The gene x environment / RERI analyses ARE reported, as Supplementary Table S5; `03` Task 3b emits them exactly as published (carrier-coded, delta-method CIs).

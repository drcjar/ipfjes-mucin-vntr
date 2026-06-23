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

"""00_build_dosage_pheno.py — assemble muc5_tsnps_dosage_pheno.csv.

This is the previously-undocumented step that builds the input read by
01_haplogroup_assignment.py, 06_tsnp_table.py, and verify_numbers.py. It merges:

  - muc5_tsnps_dosage.raw  : PLINK2 `--export A` dosages (REPRODUCE.md step 2):
                             FID/IID + the 16 MUC5AC tSNPs and rs35705950, as
                             additive 0/1/2 counts. The chr11 fileset is
                             unphenotyped, so its PHENOTYPE column is all -9.
  - ipfjes_cases_and_controls.fam : the real case/control coding in column 6
                             (1 = control, 2 = case, -9 = unset).

The chr11 sample IDs are "<id>_<id>"; the .fam uses the bare "<id>", so the
phenotype is matched on the bare ID. Output columns: FID, IID, PAT, MAT, SEX,
the 17 dosage columns (PLINK order), then `pheno`.

Run: uv run --python .venv python3 code/00_build_dosage_pheno.py
"""
import pandas as pd

ROOT = "."

# dosages (PLINK2 --export A); drop the unphenotyped PHENOTYPE column (all -9)
dos = pd.read_csv(f"{ROOT}/muc5_tsnps_dosage.raw", sep=r"\s+").drop(columns=["PHENOTYPE"])

# real phenotype from the merged-analysis fileset (col 1 = IID, col 5 = phenotype)
fam = pd.read_csv(f"{ROOT}/ipfjes_cases_and_controls.fam", sep=r"\s+", header=None)
pmap = dict(zip(fam[1].astype(str), fam[5].astype(int)))

bare = dos["IID"].astype(str).str.split("_").str[0]
missing = sorted(set(bare) - set(pmap))
assert not missing, f"{len(missing)} IDs have no phenotype in the .fam, e.g. {missing[:3]}"
dos["pheno"] = bare.map(pmap).astype(int)

dos.to_csv(f"{ROOT}/muc5_tsnps_dosage_pheno.csv", index=False)
print(f"wrote muc5_tsnps_dosage_pheno.csv: {len(dos)} rows, {dos.shape[1]} cols, "
      f"pheno counts {dict(sorted(dos.pheno.value_counts().items()))}")

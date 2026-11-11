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
# # 12. Imputed MUC5AC VNTR repeat length (ImputeSV / Bai et al. panel)
#
# The manuscript argues, from the tagging structure alone, that an SNP-imputed
# repeat length "is a function of the same haplotype the tSNPs index and inherits
# the same linkage disequilibrium with rs35705950, so it cannot resolve the
# identifiability problem". This script tests that claim directly, using repeat
# lengths imputed into the 798 by the ImputeSV server (Bai et al. 2026 VNTR panel,
# Minimac4 v4.1.6 after SHAPEIT5 phasing). See REPRODUCE.md step 6 for the
# submission and the governance note.
#
# Three questions:
#
# 1. **Does the imputed length reproduce the clade tag?** If yes, the tags and the
#    imputed length are the same information, and imputation buys no identifiability.
# 2. **Is repeat length associated with IPF independently of rs35705950?** This is
#    the conditional test Bai et al.'s UK Biobank analysis does not perform, in an
#    adjudicated IPF case-control series.
# 3. **Does the MUC5AC H2 residual track measured length?** Reported here but *not*
#    interpreted: the composite length is built from alleles that themselves mark
#    H2 membership, so adjusting H2 for it is partly circular. See the note in §4.

# %%
import gzip
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = "/home/drcjar/academic/ipfjes-omics"
VCF = f"{ROOT}/chr11.VNTR.GRCh38.noSNPs/chr11.VNTR.GRCh38.noSNPs.dose.vcf.gz"

# MUC5AC exon 31 (the 24 bp coding repeat); MUC5B intron-37 minisatellite at 1,255,550.
AC_EXON31 = (1_184_000, 1_192_500)
REGION = (1_140_000, 1_280_000)


def read_dosages(lo, hi, min_r2=0.0, min_af=0.0):
    """Return (metadata, DS matrix, sample ids) for VNTR alleles in a window."""
    rows, samples = [], None
    with gzip.open(VCF, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                continue
            pos = int(f[1])
            if not (lo < pos < hi):
                continue
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            r2, af = float(info.get("R2", 0)), float(info.get("AF", 0))
            if r2 < min_r2 or af <= min_af:
                continue
            ix = f[8].split(":").index("DS")
            ds = np.array([float(s.split(":")[ix]) for s in f[9:]])
            rows.append((f[2], pos, len(f[4]) - len(f[3]), af, r2, ds))
    meta = pd.DataFrame([r[:5] for r in rows],
                        columns=["id", "pos", "delta_bp", "AF", "R2"]).set_index("id")
    dos = pd.DataFrame({r[0]: r[5] for r in rows}, index=samples)
    return meta, dos, samples


# %% [markdown]
# ## Assemble the analysis frame
#
# Pseudonyms are mapped back through `imputesv_id_map.tsv` (local, gitignored).

# %%
meta, dos, samples = read_dosages(*REGION, min_r2=0.80, min_af=0.01)
meta_q = meta[meta.AF < 0.99]
print(f"{len(meta_q)} well-imputed common VNTR alleles (R2>=0.80, AF 0.01-0.99) in the mucin region")

idmap = pd.read_csv(f"{ROOT}/imputesv_id_map.tsv", sep="\t", dtype=str).set_index("new")["pid"]
dos.index = [idmap[s] for s in samples]
dos.index.name = "pid"

hap = pd.read_csv(f"{ROOT}/muc5_haplogroups.csv")
hap["pid"] = hap.IID.str.split("_").str[0]
clin = pd.read_csv(f"{ROOT}/flat_data_genotype_subset.csv")
clin["pid"] = clin.participant_id.astype(str)
pcs = pd.read_csv(f"{ROOT}/code/ipfjes_pca.eigenvec", sep=r"\s+")
pcs["pid"] = pcs.IID.astype(str)

dm = (clin
      .merge(hap[["pid", "H1_count", "H2_count", "H3_count",
                  "MUC5B_H1_count", "rs35705950_T_count"]], on="pid")
      .merge(pcs[["pid", "PC1", "PC2", "PC3", "PC4", "PC5"]], on="pid")
      .merge(dos.reset_index(), on="pid")
      .dropna(subset=["age", "ever_smoked", "peto_exposed", "centre"]))
dm = dm[(dm.centre != 16) & dm.ethnicity.astype(str).str.contains("White", case=False, na=False)].copy()

# composite repeat length: total inserted bp across MUC5AC exon 31, vs GRCh38
ac_meta, ac_dos, ac_samples = read_dosages(*AC_EXON31, min_r2=0.30, min_af=0.0)
ac_dos.index = [idmap[s] for s in ac_samples]
ac_len = (ac_dos * ac_meta.delta_bp).sum(axis=1).rename("ac_len_bp")
dm = dm.merge(ac_len.rename_axis("pid").reset_index(), on="pid")
dm["ac_len_kb"] = dm.ac_len_bp / 1000
dm["n24mer"] = dm.ac_len_bp / 24.0          # 24 bp MUC5AC repeat unit
dm["len_hiq_kb"] = (192 * dm["chr11_HQA241TR_6551"] + 750 * dm["chr11_HQA241TR_6559"]) / 1000

COVS = "age + ever_smoked + PC1+PC2+PC3+PC4+PC5 + C(centre)"
print(f"N={len(dm)} (cases={int(dm.case.sum())}, controls={int((dm.case == 0).sum())})")


def fit(f, d=None):
    d = dm if d is None else d
    try:
        return smf.logit(f, data=d).fit(disp=0)
    except Exception:
        return smf.logit(f, data=d).fit(method="bfgs", disp=0)


def orci(m, t):
    b, se, p = m.params[t], m.bse[t], m.pvalues[t]
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), p


# %% [markdown]
# ## 1. Does the imputed length reproduce the clade tag?
#
# `chr11_HQA241TR_6559` is a 751 bp insertion in MUC5AC exon 31, the length variant
# that tracked H3 in the HGSVC3 assemblies.

# %%
L751, L192 = "chr11_HQA241TR_6559", "chr11_HQA241TR_6551"
LAB = {L192: "MUC5AC ex31 +192bp @1,184,994", L751: "MUC5AC ex31 +751bp @1,188,533",
       "chr11_HQA241TR_6596": "MUC5B int37 -118bp @1,255,550",
       "chr11_HQA241TR_6598": "MUC5B int37 -236bp @1,255,550",
       "chr11_HQA241TR_6599": "MUC5B int37 -354bp @1,255,550"}
tgt = ["H1_count", "H2_count", "H3_count", "MUC5B_H1_count", "rs35705950_T_count"]
print(f"{'imputed VNTR allele':32s} {'R2':>6s} " + " ".join(f"{t.replace('_count',''):>12s}" for t in tgt))
for v, lab in LAB.items():
    r = [np.corrcoef(dm[v], dm[t])[0, 1] for t in tgt]
    print(f"{lab:32s} {meta.loc[v,'R2']:6.3f} " + " ".join(f"{x:12.3f}" for x in r))

print("\nSame variable? Imputed length vs the H3 clade tag:")
for lab, term, extra in [("751bp insertion       ", L751, ""),
                         ("751bp + rs35705950    ", L751, "rs35705950_T_count+"),
                         ("H3 clade tag          ", "H3_count", ""),
                         ("H3 tag + rs35705950   ", "H3_count", "rs35705950_T_count+")]:
    o, l, h, p = orci(fit(f"case~{term}+{extra}{COVS}"), term)
    print(f"  {lab} OR {o:.2f} ({l:.2f}-{h:.2f})  p={p:.2g}")

l_h3, l_len = fit(f"case~H3_count+{COVS}").llf, fit(f"case~{L751}+{COVS}").llf
l_both = fit(f"case~H3_count+{L751}+{COVS}").llf
for lab, base in [("add 751bp length to H3 clade", l_h3), ("add H3 clade to 751bp length", l_len)]:
    c = 2 * (l_both - base)
    print(f"  LRT {lab:30s} chi2(1)={c:6.2f}  p={stats.chi2.sf(c, 1):.2f}")

# %% [markdown]
# ## 2. Is repeat length associated with IPF independently of rs35705950?
#
# The conditional test Bai et al. do not perform. Three length definitions, to show
# the answer does not depend on how rare alleles are handled.

# %%
print("Total MUC5AC exon-31 inserted length per person (bp, diploid, vs GRCh38):")
print("  " + dm.ac_len_bp.describe(percentiles=[.05, .25, .5, .75, .95]).round(1)
      .to_string().replace("\n", "\n  "))
print("\nCorrelation of composite length with clade counts and rs35705950:")
for t in ["H1_count", "H2_count", "H3_count", "rs35705950_T_count"]:
    print(f"  {t:20s} r = {np.corrcoef(dm.ac_len_bp, dm[t])[0,1]:+.3f}")

print("\nComposite repeat length vs IPF:")
for lab, v, unit in [("all alleles R2>=0.30  ", "ac_len_kb", "per kb"),
                     ("2 common alleles only ", "len_hiq_kb", "per kb"),
                     ("as 24-bp repeat copies", "n24mer", "per copy")]:
    o1, l1, h1, p1 = orci(fit(f"case~{v}+{COVS}"), v)
    o2, l2, h2, p2 = orci(fit(f"case~{v}+rs35705950_T_count+{COVS}"), v)
    print(f"  {lab}  unadj OR {o1:.3f} ({l1:.3f}-{h1:.3f}) p={p1:.3g}   "
          f"adj OR {o2:.3f} ({l2:.3f}-{h2:.3f}) p={p2:.2g}  [{unit}]")
print(f"\n  r between the two length definitions: {np.corrcoef(dm.ac_len_kb, dm.len_hiq_kb)[0,1]:.4f}")
print(f"  participants carrying any rare (AF<0.01) exon-31 allele: "
      f"{int((np.abs(dm.ac_len_kb - dm.len_hiq_kb) > 0.05).sum())}/{len(dm)}")

# %% [markdown]
# ## 3. MUC5B loci
#
# Plender et al.'s expanded MUC5B haplotypes occur only in individuals of African
# descent; this sample is 96.6% white, so the large MUC5B central-exon insertion is
# essentially absent and poorly imputed. The intron-37 minisatellite of Desseyn et
# al. 1999 is, by contrast, common and well imputed.

# %%
mb, _, _ = read_dosages(1_240_000, 1_260_000, min_r2=0.0, min_af=0.0)
print(mb[mb.AF > 0].to_string())

# %% [markdown]
# ## 4. The H2 residual against measured length — reported, not interpreted
#
# The composite length regresses on the clade counts with R²=0.88
# (`length ≈ 200*H2 + 910*H3`), giving VIF ≈ 9 when H3 and length are entered
# together. More importantly the +192 bp allele marks "not H1", so it carries H2
# membership: adjusting H2 for a length built from it is partly circular. The
# numbers are printed for completeness; the manuscript does not rest on them.

# %%
m = smf.ols("ac_len_bp ~ H2_count + H3_count", data=dm).fit()
print(f"ac_len_bp = {m.params['Intercept']:.1f} + {m.params['H2_count']:.1f}*H2 "
      f"+ {m.params['H3_count']:.1f}*H3   R^2={m.rsquared:.4f}")
for c in ["H2_count", "H3_count", "ac_len_kb"]:
    others = [x for x in ["H2_count", "H3_count", "ac_len_kb"] if x != c]
    r2 = smf.ols(f"{c} ~ " + " + ".join(others), data=dm).fit().rsquared
    print(f"  VIF({c:10s}) = {1/max(1-r2,1e-12):>6.1f}")

for lab, f_ in [("+ rs35705950 (published)         ", f"case~H2_count+H3_count+rs35705950_T_count+{COVS}"),
                ("+ rs35705950 + composite length  ", f"case~H2_count+H3_count+rs35705950_T_count+ac_len_kb+{COVS}"),
                ("+ rs35705950 + length, H3 removed", f"case~H2_count+rs35705950_T_count+ac_len_kb+{COVS}")]:
    o, l, h, p = orci(fit(f_), "H2_count")
    print(f"  H2 {lab} OR {o:.2f} ({l:.2f}-{h:.2f})  p={p:.3f}")
print("DONE")

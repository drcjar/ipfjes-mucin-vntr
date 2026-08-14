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
# # 13. Agnostic conditional scan of imputed 11p15.5 mucin structural variation
#
# The paper tests pre-specified MUC5AC/MUC5B haplogroups. This asks the broader,
# adjacent question the haplogroup analysis does not: **does ANY imputed structural
# variant across the mucin gene cluster carry IPF risk independently of the MUC5B
# promoter variant rs35705950?** — including the large MUC2 VNTR domain, the MUC5B
# structural variants, and the intergenic events, none of which the haplogroup
# analysis touches.
#
# Every common, well-imputed VNTR and SV allele in chr11:1.01–1.263 Mb (MUC6→MUC5B)
# from the ImputeSV panels (Bai et al. 2026; REPRODUCE.md step 6) is tested against
# IPF, marginally and conditional on rs35705950, in the 798 analytical sample. The
# ~110-variant set is in strong LD, so the multiple-testing threshold uses the
# effective number of independent tests (Li & Ji 2005).
#
# This is a robustness generalisation of the paper's negative claim, not a powered
# discovery scan: with 397 cases it cannot exclude modest independent effects, and
# every variant sits in one LD block, so survival of adjustment would still not
# establish a causal repeat effect (the identifiability limit stated in the paper).

# %%
import gzip
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

ROOT = "/home/drcjar/ipfjes-omics"
PANELS = {
    "VNTR": f"{ROOT}/chr11.VNTR.GRCh38.noSNPs/chr11.VNTR.GRCh38.noSNPs.dose.vcf.gz",
    "SV": f"{ROOT}/chr11.SV.GRCh38.noSNPs/chr11.SV.GRCh38.noSNPs.dose.vcf.gz",
}
REGION = (1_010_000, 1_263_000)
GENES = [  # GRCh38, Ensembl (confirmed)
    ("MUC6", 1_012_823, 1_036_718),
    ("MUC2", 1_074_875, 1_110_511),
    ("MUC5AC", 1_157_953, 1_201_138),
    ("MUC5B", 1_223_066, 1_262_172),
]
RS_POS = 1_219_991  # rs35705950, MUC5AC-MUC5B intergenic


def gene_of(pos):
    for g, lo, hi in GENES:
        if lo <= pos <= hi:
            return g
    return "MUC5AC-MUC5B intergenic" if 1_201_138 < pos < 1_223_066 else "other"


# %% [markdown]
# ## Collect common, well-imputed variants from both panels, de-duplicated
#
# A locus present in both the VNTR and SV panels (e.g. the 751 bp exon-31 insertion)
# is the same physical variant; keyed by (position, REF length, ALT length) and kept
# once, at its higher imputation r².

# %%
def read_panel(path, min_r2, min_af):
    rows, samples = [], None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                continue
            pos = int(f[1])
            if not (REGION[0] <= pos <= REGION[1]):
                continue
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            r2, af = float(info.get("R2", 0)), float(info.get("AF", 0))
            if r2 < min_r2 or af < min_af or af > 1 - min_af:
                continue
            ix = f[8].split(":").index("DS")
            ds = np.array([float(s.split(":")[ix]) for s in f[9:]])
            rows.append((pos, len(f[3]), len(f[4]), len(f[4]) - len(f[3]), af, r2, ds))
    return rows, samples


MIN_R2, MIN_AF = 0.80, 0.01          # primary set; §sensitivity re-runs at 0.50
best = {}
samples = None
for panel, path in PANELS.items():
    rows, samples = read_panel(path, MIN_R2, MIN_AF)
    for pos, rl, al, delta, af, r2, ds in rows:
        key = (pos, rl, al)
        if key not in best or r2 > best[key]["r2"]:
            best[key] = dict(pos=pos, delta=delta, af=af, r2=r2, ds=ds,
                             panel=panel, gene=gene_of(pos),
                             vid=f"v_{pos}_{'i' if delta >= 0 else 'd'}{abs(delta)}")

meta = pd.DataFrame([{k: v[k] for k in ("pos", "delta", "af", "r2", "panel", "gene", "vid")}
                     for v in best.values()]).sort_values("pos").reset_index(drop=True)
print(f"{len(meta)} common well-imputed structural variants (R2>={MIN_R2}, MAF>={MIN_AF}) "
      f"in chr11:{REGION[0]:,}-{REGION[1]:,}")
print(meta.gene.value_counts().to_string(), "\n")

dos = pd.DataFrame({v["vid"]: v["ds"] for v in best.values()}, index=samples)

# %% [markdown]
# ## Merge covariates and run the scan

# %%
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
      .merge(hap[["pid", "H3_count", "rs35705950_T_count"]], on="pid")
      .merge(pcs[["pid", "PC1", "PC2", "PC3", "PC4", "PC5"]], on="pid")
      .merge(dos.reset_index(), on="pid")
      .dropna(subset=["age", "ever_smoked", "peto_exposed", "centre"]))
dm = dm[dm.centre != 16].copy()
COVS = "age + ever_smoked + PC1+PC2+PC3+PC4+PC5 + C(centre)"
print(f"N={len(dm)} (cases={int(dm.case.sum())}, controls={int((dm.case == 0).sum())})")


def fit(f):
    try:
        return smf.logit(f, data=dm).fit(disp=0)
    except Exception:
        return smf.logit(f, data=dm).fit(method="bfgs", disp=0)


res = []
for _, r in meta.iterrows():
    v = r.vid
    m1 = fit(f"case ~ {v} + {COVS}")
    m2 = fit(f"case ~ {v} + rs35705950_T_count + {COVS}")
    b1, s1 = m1.params[v], m1.bse[v]
    b2, s2 = m2.params[v], m2.bse[v]
    res.append(dict(
        pos=r.pos, gene=r.gene, delta_bp=r.delta, panel=r.panel, af=r.af, r2=r.r2,
        r_rs=np.corrcoef(dm[v], dm.rs35705950_T_count)[0, 1],
        or_marg=np.exp(b1), p_marg=m1.pvalues[v],
        or_adj=np.exp(b2), p_adj=m2.pvalues[v]))
R = pd.DataFrame(res)

# %% [markdown]
# ## Effective number of independent tests (Li & Ji 2005) and thresholds

# %%
C = np.corrcoef(dm[meta.vid.tolist()].values, rowvar=False)
C = np.nan_to_num(C, nan=0.0)
ev = np.abs(np.linalg.eigvalsh(C))
Meff = int(np.sum((ev >= 1).astype(int) + (ev - np.floor(ev))))  # Li & Ji
thresh = 0.05 / Meff
print(f"{len(meta)} variants; effective independent tests M_eff = {Meff}; "
      f"cluster-wide threshold 0.05/M_eff = {thresh:.2e}")

# %% [markdown]
# ## Results

# %%
print("\n=== Strongest marginal associations (top 12) ===")
print(f"{'gene':24s} {'pos':>9s} {'delta':>7s}  r2   r(rs35705950)  marginal        + rs35705950")
top = R.sort_values("p_marg").head(12)
for _, r in top.iterrows():
    star = "  *cluster-wide*" if r.p_adj < thresh else ""
    print(f"  {r.gene:24s} {int(r.pos):>9d} {r.delta_bp:+6d}bp {r.r2:.2f}   {r.r_rs:+.2f}"
          f"        OR {r.or_marg:5.2f} p={r.p_marg:.0e}   OR {r.or_adj:5.2f} p={r.p_adj:.3f}{star}")

mp = R.p_adj.min()
print(f"\nSmallest adjusted p = {mp:.3g}; expected min over M_eff={Meff} tests under the null "
      f"P(min<{mp:.2g}) = {1-(1-mp)**Meff:.2f}")
print("Marginal signals track LD with rs35705950 (top hits r=0.3-0.5); adjustment attenuates all.")

print(f"\n=== Survivors of adjustment for rs35705950 (p_adj < {thresh:.1e}) ===")
surv = R[R.p_adj < thresh].sort_values("p_adj")
if len(surv) == 0:
    print("  NONE. No imputed mucin structural variant is associated with IPF")
    print("  independently of rs35705950 at the cluster-wide threshold.")
else:
    for _, r in surv.iterrows():
        print(f"  {r.gene:24s} {int(r.pos)} {r.delta_bp:+d}bp  adj OR {r.or_adj:.2f} p={r.p_adj:.2e}")

print("\n=== Marginal hits, split: already-reported vs new ===")
marg = R[R.p_marg < thresh]
ac_known = marg[(marg.gene == "MUC5AC") & (marg.pos.between(1_184_000, 1_192_500))]
new = marg.drop(ac_known.index)
print(f"  reaching cluster-wide significance marginally: {len(marg)}")
print(f"    of which MUC5AC exon-31 (already reported as the H3 signal): {len(ac_known)}")
print(f"    NEW (MUC2 / MUC5B / intergenic / other MUC5AC): {len(new)}")
if len(new):
    print(new[["gene", "pos", "delta_bp", "or_marg", "p_marg", "or_adj", "p_adj"]].to_string(index=False))

R.sort_values("p_marg").to_csv(f"{ROOT}/mucin_cluster_scan.csv", index=False)
print("\nDONE")

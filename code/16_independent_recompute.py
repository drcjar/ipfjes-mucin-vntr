"""16_independent_recompute.py: re-derive Table 2 by a deliberately different route.

`verify_numbers.py` proves the manuscript matches what this pipeline computes. It
cannot prove the pipeline is right: a wrong assumption shared by the analysis and
the audit agrees with itself. Supplementary Table S2's H2 rows were wrong for
weeks for exactly that reason, and were only caught because two scripts happened
to disagree.

This script makes that redundancy deliberate. It differs from the main pipeline
at three points:

  1. Haplogroup counts are derived here from the raw PLINK dosage exports, not
     read from `muc5_haplogroups.csv`. A bug in 01_haplogroup_assignment.py would
     therefore show up rather than propagate.
  2. Models are fitted with statsmodels GLM on a design matrix built by hand,
     not through the formula API, so a formula or patsy quirk cannot be shared.
  3. Centre is dummy-coded explicitly rather than by C(centre).

It then compares its own estimates with `code/verified_numbers.json`. Any
disagreement means one of the two routes is wrong, and both must be inspected.
Exits non-zero on disagreement.
"""
import json, sys, warnings
import numpy as np, pandas as pd, statsmodels.api as sm

warnings.filterwarnings("ignore")
R = "/home/drcjar/ipfjes-omics"

clin = pd.read_csv(f"{R}/flat_data_genotype_subset.csv"); clin["pid"] = clin.participant_id.astype(str)
pcs  = pd.read_csv(f"{R}/code/ipfjes_pca.eigenvec", sep=r"\s+"); pcs["pid"] = pcs.IID.astype(str)

def raw(fn):
    d = pd.read_csv(f"{R}/{fn}", sep=r"\s+"); d["pid"] = d.IID.astype(str).str.split("_").str[0]; return d
tsnp, h2raw, mbraw = raw("muc5_tsnps_dosage.raw"), raw("muc5ac_h2_tag.raw"), raw("muc5b_h1h2_tsnp.raw")

H3C  = "chr11:1165501:T:A_T"          # counted T = NON-H3, so H3 = 2 - dosage
RSC  = "chr11:1219991:G:T_G"          # counted G = non-risk, so T count = 2 - dosage
H2C  = [c for c in h2raw.columns  if c.startswith("chr11:1182279")][0]
MBC  = [c for c in mbraw.columns  if c.startswith("chr11:1244757")][0]

d = (clin.merge(tsnp[["pid", H3C, RSC]], on="pid")
         .merge(h2raw[["pid", H2C]], on="pid")
         .merge(mbraw[["pid", MBC]], on="pid")
         .merge(pcs[["pid", "PC1", "PC2", "PC3", "PC4", "PC5"]], on="pid")
         .dropna(subset=["age", "ever_smoked", "peto_exposed", "centre"]))
d = d[(d.centre != 16) & d.ethnicity.astype(str).str.contains("White", case=False, na=False)].copy()

# Orient each tag toward the haplogroup it names, from control frequency, exactly
# as the Methods describe but computed here from scratch.
ctl = d[d.case == 0]
d["H3"] = 2 - d[H3C] if ctl[H3C].mean() / 2 > 0.5 else d[H3C]
d["T"]  = 2 - d[RSC] if ctl[RSC].mean() / 2 > 0.5 else d[RSC]
d["H2"] = 2 - d[H2C] if ctl[H2C].mean() / 2 > 0.5 else d[H2C]
d["B1"] = 2 - d[MBC] if ctl[MBC].mean() / 2 > 0.5 else d[MBC]   # counted allele tags MUC5B H2

CENTRE = pd.get_dummies(d.centre.astype(int), prefix="c", drop_first=True).astype(float)
BASE   = pd.concat([d[["age", "ever_smoked", "PC1", "PC2", "PC3", "PC4", "PC5"]].astype(float),
                    CENTRE], axis=1)

def fit(terms):
    X = sm.add_constant(pd.concat([d[terms].astype(float), BASE], axis=1), has_constant="add")
    m = sm.GLM(d.case.astype(float), X, family=sm.families.Binomial()).fit()
    return {t: (np.exp(m.params[t]),
                np.exp(m.params[t] - 1.959964 * m.bse[t]),
                np.exp(m.params[t] + 1.959964 * m.bse[t])) for t in terms}

MODELS = {
    "Table 2 / 1a H3 unadj":    (["H3"],             {"H3": "H3_count"}),
    "Table 2 / 1b H3+rs":       (["H3", "T"],        {"H3": "H3_count", "T": "rs35705950_T_count"}),
    "Table 2 / joint (H1 ref)": (["H2", "H3"],       {"H2": "H2t", "H3": "H3_count"}),
    "Table 2 / joint+rs":       (["H2", "H3", "T"],  {"H2": "H2t", "H3": "H3_count", "T": "rs35705950_T_count"}),
    "Table 2 / 4a MUC5B H1":    (["B1"],             {"B1": "MUC5B_H1_count"}),
    "Table 2 / 4b MUC5B H1+rs": (["B1", "T"],        {"B1": "MUC5B_H1_count", "T": "rs35705950_T_count"}),
    "Table 2 / rs alone":       (["T"],              {"T": "rs35705950_T_count"}),
}

V = json.load(open(f"{R}/code/verified_numbers.json"))["values"]   # {values, provenance}
print(f"Independent re-derivation of Table 2  (N={len(d)}, {int(d.case.sum())} cases, "
      f"{int((d.case == 0).sum())} controls)\n")
print(f"{'model':<26}{'term':<22}{'this script':>26}{'pipeline':>26}   ")
bad = 0
for sec, (terms, names) in MODELS.items():
    got = fit(terms)
    for t in terms:
        o, lo, hi = got[t]
        mine = f"{o:.2f} ({lo:.2f}-{hi:.2f})"
        keys = [f"{sec}|{names[t]} {s}" for s in ("OR", "CI-lo", "CI-hi")]
        if not all(k in V for k in keys):
            print(f"{sec:<26}{names[t]:<22}{mine:>26}{'no pipeline value':>26}   ?"); bad += 1; continue
        theirs = f"{float(V[keys[0]]):.2f} ({float(V[keys[1]]):.2f}-{float(V[keys[2]]):.2f})"
        agree = all(abs(a - float(V[k])) <= 0.0051
                    for a, k in zip((o, lo, hi), keys))
        bad += 0 if agree else 1
        print(f"{sec:<26}{names[t]:<22}{mine:>26}{theirs:>26}   {'OK' if agree else 'DIFFERS'}")

n = sum(len(v[0]) for v in MODELS.values())
print(f"\n{n} estimates re-derived independently · {n-bad} agree · {bad} differ")
if bad:
    print("A disagreement means one of the two routes is wrong. Inspect both before trusting either.")
    sys.exit(1)
print("Both routes agree, so the Table 2 estimates do not rest on a single code path.")

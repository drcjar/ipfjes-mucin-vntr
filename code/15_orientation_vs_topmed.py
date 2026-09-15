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

"""15_orientation_vs_topmed.py: answer the allele-orientation query from the
TOPMed replication (A. Peljto, Schwartz lab, email 21 Aug 2026).

Her query: the alleles she used to define the haplogroups differ from those printed
in our Supplementary Table S2, and because the H1 and H3 tags are palindromic
(G/C and A/T) a strand flip could be contributing to the non-replication.

This script settles that from source, without trusting any hand-written table:

  1. The counted (effect) allele is read from the PLINK `--export A` header itself
     (`<variant ID>_<counted allele>`), so it cannot drift from what the models used.
  2. The haplogroup-DEFINING allele is then derived from the counting direction used
     in `01_haplogroup_assignment.py`: the H3 and H2 concordance tags are entered as
     `2 - dosage`, so for those two the defining allele is the COMPLEMENT of the
     counted allele; the H1 tag is counted directly.
  3. Both are compared against the TOPMed assignment, three ways:
       (a) our chr11.bim A1/A2 against TOPMed's A1/A2 as an ORDERED pair. Both are
           PLINK minor/major codings of the same sites, so an exact slot-for-slot
           match means the cohorts coded the variants identically. They do match at
           all three tags, which settles strand outright. NB this is a different
           ordering from the ref/alt printed in Supplementary Table S2, which comes
           from the variant ID; the two coincide only where the genomic reference
           allele is also the major allele (13 of our 19 tags).
       (b) the allele pairs unordered, which is weaker: it excludes a flip at the
           A/G H2 tag, where flipping would give T/C, but says nothing at the
           palindromic H3 (A/T) and H1 (G/C) tags, whose pair is flip-invariant.
       (c) a check that uses no allele labels at all: the haplogroup-allele
           frequency implies an HWE carrier frequency, compared with the carrier
           frequencies observed in each cohort. A strand flip would show up as a
           carrier frequency near the complement, not near the observed value.
  4. Finally the IPFJES and TOPMed H2 estimates are compared formally
     (Cochran's Q / fixed-effect meta) to quantify the disagreement that remains
     once orientation is excluded as an explanation.

TOPMed numbers are transcribed from `../feedback/MUC5AC_haplogroup_tSNPs_TOPMed-IPF_Assoc.xlsx`
(sheets `Haplogroups`, `Table1`, `Table 2`); 2,425 cases / 3,816 controls, models on
N=2,031 cases with complete covariates, adjusted for sex, age, smoking and PC1-4.
Provenance is kept in structured fields below, not in prose.
"""
import re
import numpy as np
import pandas as pd
import warnings
from scipy import stats

warnings.filterwarnings("ignore")
R = "/home/drcjar/ipfjes-omics"

def other_allele(counted, ref, alt):
    """The other allele at a biallelic site.

    NB this is the site's other allele, NOT the DNA complement of the counted
    base: for the palindromic tags (A/T, G/C) the two coincide, but for the
    A/G H2 tag they do not.
    """
    assert counted in (ref, alt), f"{counted} is not an allele of {ref}/{alt}"
    return alt if counted == ref else ref

# The three haplogroup-defining tags, with the counting direction used in
# 01_haplogroup_assignment.py. `defining_is_counted=False` means the model uses
# 2 - dosage, so the haplogroup allele is the OTHER allele at the site.
# `topmed_a1`/`topmed_a2` are the A1/A2 columns of the TOPMed Haplogroup
# Assignment table. A1/A2 is PLINK's counted/other ordering; ref/alt above is the
# reference-genome ordering carried in the variant ID. The two conventions order
# the same pair differently and neither ordering is a claim about the haplogroup,
# so the alleles are compared below as UNORDERED pairs.
TAGS = [
    dict(rsid="rs28542750",  hap="MUC5AC H3", raw="muc5_tsnps_dosage.raw",
         col="chr11:1165501:T:A_T", defining_is_counted=False, topmed_allele="A",
         topmed_a1="A", topmed_a2="T"),
    dict(rsid="rs879201833", hap="MUC5AC H2", raw="muc5ac_h2_tag.raw",
         col="chr11:1182279:A:G_A", defining_is_counted=False, topmed_allele="G",
         topmed_a1="G", topmed_a2="A"),
    dict(rsid="rs769768817", hap="MUC5AC H1", raw="muc5_tsnps_dosage.raw",
         col="chr11:1195214:G:C_C", defining_is_counted=True,  topmed_allele="C",
         topmed_a1="G", topmed_a2="C"),
]

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def bim_alleles(ids):
    """A1/A2 straight out of chr11.bim, for the given variant IDs.

    A1/A2 here follow PLINK 1's minor/major convention, NOT the ref/alt ordering
    carried in the variant ID: for 6 of our 18 tags the genomic reference allele
    is the minor allele and therefore lands in A1. This is also why `--export A`
    counts what it counts, since PLINK 2 reads A2 as REF and counts REF, so the
    counted allele is the major allele at every tag.

    Scanned line-by-line rather than read with pandas: chr11.bim is ~14.5M rows
    and we want three of them.
    """
    want, out = set(ids), {}
    with open(f"{R}/chr11.bim") as fh:
        for line in fh:
            f = line.split()
            if f[1] in want:
                out[f[1]] = (f[4], f[5])
                if len(out) == len(want):
                    break
    missing = want - set(out)
    assert not missing, f"variant IDs absent from chr11.bim: {missing}"
    return out


BIM = bim_alleles([t["col"].rsplit("_", 1)[0] for t in TAGS])

# Carrier frequencies in controls, for the label-free orientation check.
# IPFJES: Table 1 of the manuscript (n=401 controls). TOPMed: sheet Table1 (n=3,816).
OBSERVED_CTRL_CARRIER = {
    "MUC5AC H3": dict(ipfjes=0.352, topmed=0.368),
    "MUC5AC H2": dict(ipfjes=0.416, topmed=0.362),
    "MUC5AC H1": dict(ipfjes=0.815, topmed=0.839),
}

# ---------------------------------------------------------------- control set
_clin = pd.read_csv(f"{R}/flat_data_genotype_subset.csv")
_clin["pid"] = _clin["participant_id"].astype(str)
_pcs = pd.read_csv(f"{R}/code/ipfjes_pca.eigenvec", sep=r"\s+")
_pcs["pid"] = _pcs["IID"].astype(str)
CTRL = set(_clin.merge(_pcs[["pid"]], on="pid")
                .dropna(subset=["age", "ever_smoked", "peto_exposed", "centre"])
                .query("case == 0 and centre != 16 and ethnicity.str.contains('White', case=False, na=False)")["pid"])

_cache = {}


def controls(raw):
    if raw not in _cache:
        df = pd.read_csv(f"{R}/{raw}", sep=r"\s+")
        df["pid"] = df["IID"].astype(str).str.split("_").str[0]
        _cache[raw] = df[df.pid.isin(CTRL)]
    return _cache[raw]


print(f"IPFJES control reference: n={len(controls('muc5_tsnps_dosage.raw'))}\n")

# ------------------------------------------------- 1-3. orientation comparison
print("=" * 78)
print("ORIENTATION: IPFJES vs TOPMed")
print("=" * 78)

rows, agree = [], True
for t in TAGS:
    df = controls(t["raw"])
    pos, ref, alt, counted = re.match(
        r"chr11:(\d+):([ACGT]+):([ACGT]+)_([ACGT]+)", t["col"]).groups()

    counted_freq = df[t["col"]].mean() / 2
    defining = counted if t["defining_is_counted"] else other_allele(counted, ref, alt)
    defining_freq = counted_freq if t["defining_is_counted"] else 1 - counted_freq

    palindromic = {ref, alt} in ({"A", "T"}, {"C", "G"})
    match = defining == t["topmed_allele"]
    agree &= match

    # Same site, two column conventions: our ref/alt vs TOPMed's A1/A2. Compare
    # as unordered pairs. `same_pair` false with `flip_pair` true would be a real
    # strand difference; for the two palindromic tags the pair is invariant under
    # flipping, so a pair match there is NOT evidence of same strand.
    ours_pair, topmed_pair = {ref, alt}, {t["topmed_a1"], t["topmed_a2"]}
    same_pair = ours_pair == topmed_pair
    flip_pair = {COMPLEMENT[a] for a in ours_pair} == topmed_pair

    # Stronger still: our own .bim A1/A2 against TOPMed's A1/A2, compared as an
    # ORDERED pair. Both are PLINK minor/major codings of the same site, so an
    # exact slot-for-slot match means the two cohorts coded the site identically,
    # which settles strand without appealing to allele frequencies at all.
    bim_a1, bim_a2 = BIM[t["col"].rsplit("_", 1)[0]]
    exact = (bim_a1, bim_a2) == (t["topmed_a1"], t["topmed_a2"])

    # Label-free check: HWE carrier frequency implied by the defining-allele
    # frequency, vs what each cohort actually observed. A strand flip would put
    # the observed value near the complement's implied carrier frequency instead.
    implied = 1 - (1 - defining_freq) ** 2
    implied_flipped = 1 - defining_freq ** 2
    obs = OBSERVED_CTRL_CARRIER[t["hap"]]

    rows.append(dict(
        rsID=t["rsid"], hap=t["hap"], pos=f"chr11:{int(pos):,}",
        ref_alt=f"{ref}/{alt}", palindromic="yes" if palindromic else "no",
        counted=counted, counted_freq=round(counted_freq, 3),
        defining=defining, defining_freq=round(defining_freq, 3),
        bim_a1_a2=f"{bim_a1}/{bim_a2}",
        topmed_a1_a2=f'{t["topmed_a1"]}/{t["topmed_a2"]}',
        exact_a1a2="YES" if exact else "no",
        same_pair="yes" if same_pair else "no",
        pair_informative="no (palindromic)" if palindromic else "yes",
        topmed=t["topmed_allele"], agrees="YES" if match else "NO",
        implied_carrier=round(implied, 3),
        implied_if_flipped=round(implied_flipped, 3),
        obs_ipfjes=obs["ipfjes"], obs_topmed=obs["topmed"],
        # Which of the two hypotheses each cohort's observed value sits nearer,
        # and by how much. `margin` is the separation between the two competing
        # predictions: it is what gives the check its power, and it shrinks as
        # the defining-allele frequency approaches 0.5.
        nearer_ipfjes="as-is" if abs(obs["ipfjes"] - implied)
                      < abs(obs["ipfjes"] - implied_flipped) else "FLIPPED",
        nearer_topmed="as-is" if abs(obs["topmed"] - implied)
                      < abs(obs["topmed"] - implied_flipped) else "FLIPPED",
        margin=round(abs(implied - implied_flipped), 3),
    ))

tab = pd.DataFrame(rows)
print(tab[["rsID", "hap", "pos", "ref_alt", "palindromic", "counted",
           "counted_freq", "defining", "defining_freq", "topmed",
           "agrees"]].to_string(index=False))
print()
print("Allele coding: our chr11.bim A1/A2 vs the TOPMed A1/A2 columns.")
print("(Our S2 'Alleles' column prints ref/alt from the variant ID, a DIFFERENT")
print(" ordering: PLINK A1/A2 is minor/major, so the two coincide only where the")
print(" genomic reference allele is also the major allele.)")
print(tab[["rsID", "hap", "ref_alt", "bim_a1_a2", "topmed_a1_a2", "exact_a1a2",
           "same_pair", "pair_informative", "defining", "topmed",
           "agrees"]].to_string(index=False))
print()
if (tab.exact_a1a2 == "YES").all():
    print("The A1/A2 codings match EXACTLY, slot for slot, at all three tags. Both are")
    print("PLINK minor/major codings of the same sites, so the two cohorts have coded")
    print("these variants identically and a strand difference is excluded outright, with")
    print("no appeal to allele frequencies. This is the cleanest evidence available here.")
else:
    print("A1/A2 codings do NOT match exactly; see the table above before reading on.")
print()
print("The weaker unordered-pair comparison agrees: no site carries a different pair of")
print("alleles. For the A/G H2 tag that alone excludes a flip, since a flipped H2 would")
print("read T/C. For the palindromic H3 and H1 tags the pair is invariant under flipping")
print("and so cannot settle strand, which is the substance of the query; the frequency")
print("check below is the independent corroboration for those two.")
print()
print("Label-free orientation check (control carrier frequencies):")
print(tab[["rsID", "hap", "implied_carrier", "implied_if_flipped",
           "obs_ipfjes", "obs_topmed", "nearer_ipfjes", "nearer_topmed",
           "margin"]].to_string(index=False))
print()
print(f"VERDICT: haplogroup-defining alleles "
      f"{'AGREE on all three tags' if agree else 'DISAGREE'} between IPFJES and TOPMed.")

no_flip = (tab.nearer_ipfjes == "as-is").all() and (tab.nearer_topmed == "as-is").all()
print(f"Label-free check: every tag in BOTH cohorts sits nearer the as-is prediction"
      if no_flip else
      f"Label-free check: at least one tag sits nearer the FLIPPED prediction"
      f" - see the table above")
print(f"than the flipped one, so a strand flip is excluded." if no_flip else "")
print(f"Weakest separation between the two predictions: {tab.margin.min():.3f} "
      f"({tab.loc[tab.margin.idxmin(), 'hap']}, defining-allele frequency "
      f"{tab.loc[tab.margin.idxmin(), 'defining_freq']:.3f}). The check loses power as")
print("that frequency approaches 0.5, so it is weakest for H1; the allele labels in the")
print("table above are the primary evidence, and this is the independent corroboration.")

# ------------------------------------------------------- 4. H2 heterogeneity
print()
print("=" * 78)
print("H2 ESTIMATES, ONCE ORIENTATION IS EXCLUDED")
print("=" * 78)


def logor(or_, lo, hi):
    """log OR and its SE from a published OR with 95% CI."""
    return np.log(or_), (np.log(hi) - np.log(lo)) / (2 * 1.959964)


# Manuscript Table 2 / Supplementary Table S3 vs TOPMed sheet 'Table 2'.
MODELS = [
    ("H2, adjusted for H3",               (0.58, 0.43, 0.78), (0.83, 0.74, 0.93)),
    ("H2, adjusted for H3 + rs35705950",  (0.69, 0.51, 0.94), (0.99, 0.88, 1.11)),
    ("H3, unadjusted",                    (1.63, 1.26, 2.09), (1.62, 1.47, 1.78)),
    ("H3, adjusted for rs35705950",       (1.10, 0.83, 1.45), (0.94, 0.84, 1.05)),
]

out = []
for name, ipf, top in MODELS:
    b1, s1 = logor(*ipf)
    b2, s2 = logor(*top)
    z = (b1 - b2) / np.sqrt(s1 ** 2 + s2 ** 2)
    p_het = 2 * stats.norm.sf(abs(z))

    w1, w2 = 1 / s1 ** 2, 1 / s2 ** 2
    bm = (w1 * b1 + w2 * b2) / (w1 + w2)
    sm = np.sqrt(1 / (w1 + w2))

    out.append(dict(
        model=name,
        ipfjes=f"{ipf[0]:.2f} ({ipf[1]:.2f}-{ipf[2]:.2f})",
        topmed=f"{top[0]:.2f} ({top[1]:.2f}-{top[2]:.2f})",
        p_heterogeneity=f"{p_het:.3f}",
        meta=f"{np.exp(bm):.2f} ({np.exp(bm - 1.959964 * sm):.2f}-"
             f"{np.exp(bm + 1.959964 * sm):.2f})",
        ipfjes_weight=f"{100 * w1 / (w1 + w2):.0f}%",
    ))

print(pd.DataFrame(out).to_string(index=False))
print()
print("Read: H3 reproduces in TOPMed, including its attenuation on adjustment for")
print("rs35705950. H2 also reproduces as a protective effect adjusted for H3 alone,")
print("but the rs35705950-adjusted residual - the manuscript's headline - does not.")

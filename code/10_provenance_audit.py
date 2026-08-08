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

"""10_provenance_audit.py — external primary-source seal (fourth audit layer).

WHY THIS EXISTS
---------------
The other three layers never open the primary source:

  verify_numbers.py        numbers vs. OUR data      (recomputes from genotypes)
  07_foundational_audit.py labels/orientation        (internally consistent?)
  08_consistency_checks.py prose vs. prose           (does the paper contradict itself?)

None of them checked a CITATION against the cited document — even though that document,
Plender et al.'s supplement, sits in the repo as ../mmc2.xlsx. 07_foundational_audit.py
said so explicitly in its closing note: "the Plender Table S4 r2 values remain the one
external primary-source seal still to confirm."

That gap let a real error through: rs35779873's tagging r²=0.805 was recorded and published
as coming from Plender Supplementary Table S4. It is not in Table S4 at all — the value is
from Table S6. The mechanism was that provenance lived in a free-text `evidence` string
("Table S4 H3 r2=0.805") from which 06_tsnp_table.py regex-scraped only the NUMBER. The
table name was never read by any program, so nothing could ever contradict it.

WHAT THIS CHECKS
----------------
  A. Every (source_table, source_hap, source_r2) triple in ../muc5_allele_direction.csv
     resolves in ../mmc2.xlsx: the rsID is present in the cited sheet, at the cited GRCh38
     position, with the cited r² for the cited haplogroup (to 3 dp).
  B. Cross-sheet agreement: where an rsID appears in both S4 and S6, their r² values agree.
     The two known position-level quirks are asserted explicitly rather than tolerated.
  C. Panel-level counts we state in the manuscript match Plender's Table S4: 92 MUC5AC
     tSNPs, 35 tagging H1 / 5 H2 / 52 H3 at r²>=0.8, maxima 0.923 / 0.888 / 1.0.
  D. Table-drift: 06_tsnp_table.py's output matches the Supplementary Table S4 block pasted
     into the manuscript, cell for cell. (The published table was a hand-paste; it had
     already drifted — rs28542750's palindromic flag.)
  E. Citation sweep: every "Table S<n>" cited next to an r² in the manuscript and the
     technical methods names a table that actually contains that value.

Exit status is non-zero if any check FLAGs, so this can gate a build.

Run: uv run --python .venv python3 code/10_provenance_audit.py
"""
import re
import subprocess
import sys
import warnings

import pandas as pd

warnings.filterwarnings("ignore")
R = "/home/drcjar/ipfjes-omics"

rows = []


def audit(item, found, expected, ok):
    rows.append((item, str(found), str(expected), "OK" if ok else "FLAG"))
    return ok


# ---------------------------------------------------------------- load sources
def sheet(name, header):
    """Plender supplementary sheet, indexed by rsID, r² columns renamed H1/H2/H3."""
    d = pd.read_excel(f"{R}/mmc2.xlsx", name, header=header)
    d = d[d.iloc[:, 0].astype(str).str.startswith("rs")].copy()
    d = d.rename(columns={c: re.sub(r".*(H[123]).*", r"\1", c) for c in d.columns})
    d = d.rename(columns={d.columns[0]: "rsID"})
    d["pos"] = pd.to_numeric(d["Location (GRCh38)"], errors="coerce").astype("Int64")
    for h in ("H1", "H2", "H3"):
        d[h] = pd.to_numeric(d[h], errors="coerce")
    return d.set_index("rsID")


S4 = sheet("S4", 2)   # "Tagging SNPs (tSNPs) for the three common human haplogroups of MUC5AC"
S6 = sheet("S6", 4)   # "Squared correlations between MUC5AC haplogroups and eQTLs for severe asthma"
SHEETS = {"S4": S4, "S6": S6}

ad = pd.read_csv(f"{R}/muc5_allele_direction.csv")

# ------------------------------------------- A. every cited r² resolves in the cited table
for _, r in ad.iterrows():
    rsid, tbl = r["rsID"], r["source_table"]
    if tbl == "none":                      # rs35705950: the promoter variant, not a tag
        audit(f"{rsid} provenance", "no tagging claim", "not a haplogroup tag",
              pd.isna(r["source_r2"]) and pd.isna(r["source_hap"]))
        continue

    ok_sheet = audit(f"{rsid} cited in {tbl}", tbl in SHEETS and rsid in SHEETS[tbl].index,
                     f"present in Plender {tbl}", tbl in SHEETS and rsid in SHEETS[tbl].index)
    if not ok_sheet:
        continue

    src = SHEETS[tbl].loc[rsid]
    pos = int(re.match(r"chr11:(\d+):", r["column"]).group(1))
    audit(f"{rsid} position", f"chr11:{pos}", f"{tbl}: chr11:{int(src['pos'])}", pos == int(src["pos"]))
    audit(f"{rsid} {r['source_hap']} r2", f"{r['source_r2']:.3f}", f"{tbl}: {src[r['source_hap']]:.3f}",
          abs(float(r["source_r2"]) - float(src[r["source_hap"]])) < 5e-4)

# ------------------------------------------------------- B. cross-sheet agreement + quirks
both = [x for x in ad.rsID if x in S4.index and x in S6.index]
audit("S4/S6 r2 agree where both list an rsID", f"{len(both)} rsIDs compared", "identical to 3dp",
      all(abs(S4.loc[x, h] - S6.loc[x, h]) < 5e-4 for x in both for h in ("H1", "H2", "H3")))

# Quirk 1 — rs35779873 is an H3 tSNP that Plender report in S6, NOT S4. This is the error
# this audit exists to catch; assert the underlying fact so a future "tidy-up" can't silently
# re-attribute it to S4.
audit("rs35779873 absent from S4 (r2 is S6's)", "rs35779873" in S4.index, "absent from S4",
      "rs35779873" not in S4.index and "rs35779873" in S6.index)

# Quirk 2 — S4 and S6 list DIFFERENT variants at chr11:1,195,214. Comparing the sheets by
# position (rather than rsID) makes the H1 concordance tag look worthless (r²=0.008). The
# variant we actually extracted is G/C, matching S4's rs769768817, not S6's G/A rs28415193.
s6_at = S6[S6.pos == 1195214].index.tolist()
audit("chr11:1195214 — S4 vs S6 list different rsIDs", f"S4: rs769768817 / S6: {s6_at}",
      "distinct variants, not a contradiction",
      S4.loc["rs769768817", "pos"] == 1195214 and s6_at == ["rs28415193"]
      and ad.loc[ad.rsID == "rs769768817", "column"].iloc[0].startswith("chr11:1195214:G:C"))

# ------------------------------------------------ C. panel counts vs Plender's own Table S4
# Plender main text: "we discovered 35 tSNPs for H1 (max r2 = 0.92), 5 tSNPs for H2
# (max r2 = 0.89), and 52 tSNPs for H3 (max r2 = 1, Table S4)".
audit("Plender S4 total MUC5AC tSNPs", len(S4), "92 (manuscript Methods)", len(S4) == 92)
for h, n, mx in (("H1", 35, 0.923), ("H2", 5, 0.888), ("H3", 52, 1.0)):
    audit(f"Plender S4 {h} tags (r2>=0.8)", int((S4[h] >= 0.8).sum()), f"{n} (main text)",
          int((S4[h] >= 0.8).sum()) == n)
    audit(f"Plender S4 {h} max r2", f"{S4[h].max():.3f}", f"{mx:.3f} (main text)",
          abs(S4[h].max() - mx) < 5e-4)

tags = ad[ad.source_table != "none"]
audit("MUC5AC tSNPs used", len(tags), "17 (manuscript Methods)", len(tags) == 17)
audit("MUC5AC tSNPs are a subset of Plender's", sum(x in S4.index or x in S6.index for x in tags.rsID),
      f"all {len(tags)}", all(x in S4.index or x in S6.index for x in tags.rsID))

# ------------------------------------------- D. manuscript Table S4 vs 06_tsnp_table.py output
def cells(lines):
    out = {}
    for ln in lines:
        ln = ln.strip()
        if ln.startswith("|") and "chr11:" in ln:
            c = [x.strip() for x in ln.strip("|").split("|")]
            out[c[0]] = c
    return out


def s3_block(path):
    """The Supplementary Table S2 rows only — other tables may also carry chr11: IDs."""
    keep, out = False, []
    for ln in open(path):
        if ln.startswith("**Supplementary Table S2"):
            keep = True
        elif ln.startswith("**Supplementary Table S") and keep:
            break
        if keep:
            out.append(ln)
    return out


gen = cells(subprocess.run([sys.executable, f"{R}/code/06_tsnp_table.py"],
                           capture_output=True, text=True, check=True).stdout.splitlines())
paper = cells(s3_block(f"{R}/MUC5_VNTR_haplogroup_paper_draft.md"))
drift = [f"{k}:{i}" for k, g in gen.items()
         for i, (a, b) in enumerate(zip(g, paper.get(k, [""] * len(g)))) if a != b]
audit("Manuscript Table S2 == 06_tsnp_table.py output", f"{len(gen)} rows, {len(drift)} drifted cells",
      "0 drifted cells", set(gen) == set(paper) and not drift)
if drift:
    print("DRIFTED CELLS:", ", ".join(drift), "\n")

# --------------------------------------------------- E. citation sweep across prose documents
# Any "Table S<n>" appearing within 120 chars of an "r²=<value>" must name a table that
# really contains that value (in either sheet, for any haplogroup).
allvals = set()
for d in SHEETS.values():
    for h in ("H1", "H2", "H3"):
        allvals |= {round(v, 3) for v in d[h].dropna()}

for doc in ("MUC5_VNTR_haplogroup_paper_draft.md", "MUC5_VNTR_technical_methods.md"):
    text = open(f"{R}/{doc}").read()
    bad = []
    for m in re.finditer(r"r²\s*=\s*(\d?\.\d+)", text):
        val = round(float(m.group(1)), 3)
        window = text[max(0, m.start() - 120):m.end() + 120]
        cited = set(re.findall(r"Table (S\d+)", window)) & set(SHEETS)
        if not cited:
            continue                                   # no table cited nearby: nothing to check
        if not any(round(SHEETS[t][h], 3).eq(val).any() for t in cited for h in ("H1", "H2", "H3")):
            # tolerate values that are ours, not Plender's (e.g. MUC5B main-text r², LD r²)
            if val in allvals:
                bad.append(f"r²={val} cited near {sorted(cited)} but only in the other table")
    audit(f"{doc}: r² citations name the right table", f"{len(bad)} suspect", "0 suspect", not bad)
    for b in bad:
        print("  SUSPECT:", b)

# ------------------------------------------------------------------------------- report
w = max(len(r[0]) for r in rows)
print(f"{'PROVENANCE CHECK':{w}}  {'FOUND':34}  {'EXPECTED (primary source)':40}  VERDICT")
for it, f, e, v in rows:
    print(f"{it:{w}}  {f:34}  {e:40}  {v}")
nflag = sum(v == "FLAG" for *_, v in rows)
print(f"\n{len(rows)} provenance checks · {len(rows)-nflag} OK · {nflag} FLAG")
print("Primary source: ../mmc2.xlsx (Plender et al. 2024 supplement), read directly — "
      "no value here is taken on trust from our own prose.")
sys.exit(1 if nflag else 0)

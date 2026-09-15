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

"""07_consistency_checks.py — internal-consistency audit of the manuscript prose.

Third verification layer, complementing:
  - verify_numbers.py        : manuscript numbers vs the code that produced them
  - 06_foundational_audit.py : haplogroup labels / allele orientations vs sources
This script checks the manuscript's INTERNAL self-consistency — statements compared
to each other — which neither of the above catches. It flags:
  1. additive cohort breakdowns that don't sum
  2. a percentage that doesn't match its own fraction (n/N)
  3. Table 1 cells whose n (x%) disagrees with the column N (385 cases / 386 controls)
  4. a confidence interval that is reversed or does not bracket its OR
  5. the same OR quoted with two different CIs (edit-drift)
  6. a CI and p-value that disagree on significance

Run: uv run --python .venv python3 code/07_consistency_checks.py
(The cohort sums in check 1 are hard-coded to the current data freeze — update them
if the cohort numbers change.)
"""
import os, re
from collections import defaultdict

M = "/home/drcjar/ipfjes-omics/MUC5_VNTR_haplogroup_paper_draft.md"
txt = open(M).read()
issues = []
def flag(cls, msg): issues.append((cls, msg))

# 1. Additive cohort relationships (must sum AND the target number must be present)
# The paper reports a single N=771 IPFJES analytical sample (European-ancestry restricted);
# the genotype pipeline
# is restricted to IPFJES (code/00), so only IPFJES cohort relationships are checked.
sums = {
    "analytical 771 = 385 cases + 386 controls":       (385 + 386, 771),
    "798 = 829 merged - 31 (cases-only centre)":       (829 - 31, 798),
    "771 = 798 less 27 non-white-ethnicity":           (798 - 27, 771),
    "829 merged = 902 with clinical data - 73":        (902 - 73, 829),
}
for label, (got, want) in sums.items():
    present = re.search(r'\b' + f"{want:,}".replace(',', '[,]?') + r'\b', txt) is not None
    if got != want:
        flag("ADDITIVE", f"{label}: computes {got}, expected {want}")
    elif not present:
        flag("ADDITIVE", f"{label}: holds, but {want} no longer appears in the text")

# 2. fraction <-> percentage (both orders)
for m in re.finditer(r'(\d[\d,]*)\s*(?:/|of)\s*(\d[\d,]*)[^.\n]{0,45}?(\d+(?:\.\d+)?)\s*%', txt):
    a, b, p = int(m.group(1).replace(',', '')), int(m.group(2).replace(',', '')), float(m.group(3))
    if b and abs(100 * a / b - p) > 0.7:
        flag("FRACTION%", f"'{m.group(0).strip()[:55]}' -> {a}/{b}={100*a/b:.1f}% vs stated {p}%")
for m in re.finditer(r'(\d+(?:\.\d+)?)\s*%[^.\n]{0,45}?(\d[\d,]*)\s*/\s*(\d[\d,]*)', txt):
    p, a, b = float(m.group(1)), int(m.group(2).replace(',', '')), int(m.group(3).replace(',', ''))
    if b and abs(100 * a / b - p) > 0.7:
        flag("FRACTION%", f"'{m.group(0).strip()[:55]}' -> {a}/{b}={100*a/b:.1f}% vs stated {p}%")

# 3. Table 1 cells: "n (x%)" against cases N=385 / controls N=386
for line in txt.splitlines():
    if line.startswith("|") and "%" in line:
        cells = [c.strip() for c in line.split("|")]
        for idx, N, who in [(2, 385, "case"), (3, 386, "control")]:
            if len(cells) > idx:
                mm = re.match(r'(\d[\d,]*)\s*\((\d+(?:\.\d+)?)%\)', cells[idx])
                if mm:
                    a, p = int(mm.group(1).replace(',', '')), float(mm.group(2))
                    if abs(100 * a / N - p) > 0.7:
                        flag("TABLE1", f"{who} cell '{cells[idx]}' -> {a}/{N}={100*a/N:.1f}% vs {p}%")

# 4-6. OR / CI / p consistency
seen = defaultdict(set)
pat = re.compile(r'OR\s*(\d+\.\d+)[\s,(]+95%\s*CI\s*(\d+\.\d+)\s*[–-]\s*(\d+\.\d+)(?:[^.\n]{0,30}?p\s*=\s*([^\s,;)]+))?')
for m in pat.finditer(txt):
    orr, lo, hi = float(m.group(1)), float(m.group(2)), float(m.group(3))
    if lo > hi:
        flag("CI", f"reversed CI (lo>hi): OR {orr} (95% CI {lo}-{hi})")
    elif not (lo <= orr <= hi):
        flag("CI", f"CI does not bracket OR: OR {orr} (95% CI {lo}-{hi})")
    seen[m.group(1)].add((m.group(2), m.group(3)))
    praw = m.group(4)
    if praw:
        psig = ('×10' in praw or 'e-' in praw.lower()) or bool(re.match(r'[\d.]+$', praw) and float(praw) < 0.05)
        if (not (lo < 1 < hi)) != psig:
            flag("SIG", f"OR {orr} (CI {lo}-{hi}) and p={praw} disagree on significance")
for orr, cis in seen.items():
    if len(cis) > 1:
        flag("OR-DRIFT", f"OR {orr} quoted with different CIs {sorted(cis)} (confirm different models)")

# report

# ---------------------------------------------------------------------------
# 6. Manuscript numbers vs recomputed values.
#
# verify_numbers.py compares each recomputed value against a literal held inside
# that script. That verifies the computation, but NOT that the manuscript still
# quotes the same number: update a literal there and forget the .md and the two
# decouple silently. That is how the N=798 Hardy-Weinberg values (minimum p=0.41,
# range 0.65-0.83) survived the ancestry re-derivation into a N=771 paper.
#
# This section closes the loop: it reads the numbers actually printed in the
# manuscript and compares them with code/verified_numbers.json, which
# verify_numbers.py writes on every run.
# ---------------------------------------------------------------------------
import json, os
VN = os.path.join(os.path.dirname(M), "code", "verified_numbers.json")
if not os.path.exists(VN):
    VN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "verified_numbers.json")

if not os.path.exists(VN):
    flag("XCHECK", "code/verified_numbers.json not found - run verify_numbers.py first")
else:
    _vn = json.load(open(VN))
    # verify_numbers.py exports {"values": ..., "provenance": ...}
    V    = _vn["values"]
    PROV = _vn.get("provenance", {})

    # Provenance assertions: a recorded value must have come from the variant the
    # manuscript says it did. Supplementary Table S1's H2 rows were computed from
    # rs28519516, a secondary panel tag, while the table note promised the primary
    # rs879201833; nothing could see that because the number carried no origin.
    for _k, _p in PROV.items():
        if _k.startswith("Table S1 H2|") and not str(_p.get("tag", "")).startswith("chr11:1182279"):
            flag("PROV", f"{_k} computed from {_p.get('tag')}, not the primary H2 tag chr11:1182279")
        if "H2" in _k and "H3_count" not in _k and str(_p.get("tag", "")).startswith("chr11:1202343"):
            flag("PROV", f"{_k} computed from rs28519516, a secondary panel tag")

    def near(a, b, tol=0.0051):
        try:
            return abs(float(a) - float(b)) <= tol
        except (TypeError, ValueError):
            return str(a) == str(b)

    _SUP = str.maketrans("\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b", "0123456789-")

    def pval(cell):
        """Parse a p-value as the manuscript writes it: 0.017, or 1.8x10^-4 set
        in superscripts."""
        cell = cell.strip()
        m = re.match(r"([\d.]+)\s*\u00d7\s*10([\u2070\u00b9\u00b2\u00b3\u2074-\u2079\u207b]+)$", cell)
        if m:
            return float(m.group(1)) * 10 ** int(m.group(2).translate(_SUP))
        try:
            return float(cell)
        except ValueError:
            return None

    def cmp_p(label, cell, key):
        """p-values span orders of magnitude, so compare on the log scale, with an
        absolute fallback for the two-decimal values."""
        if key not in V:
            flag("XCHECK", f"{label}: no computed value keyed '{key}'"); return
        got, want = pval(cell), None
        try: want = float(V[key])
        except (TypeError, ValueError): pass
        if got is None or want is None:
            flag("XCHECK", f"{label}: unparsed p-value {cell!r} or {V[key]!r}"); return
        import math
        if got <= 0 or want <= 0:
            ok_ = abs(got - want) <= 0.0006
        elif min(got, want) < 0.01:
            # below 0.01 the manuscript uses scientific notation, so compare on the
            # log scale only. An absolute tolerance here would let 8.4e-6 and
            # 1.9e-6 match, since they differ by less than 0.0051.
            ok_ = abs(math.log10(got) - math.log10(want)) < 0.012
        else:
            # two-decimal values: allow rounding, e.g. 0.054 printed for 0.0538
            ok_ = abs(got - want) <= 0.0051 or abs(math.log10(got) - math.log10(want)) < 0.012
        if not ok_:
            flag("XCHECK", f"{label}: manuscript p={cell} vs computed {V[key]}")

    def cmp(label, found, key, tol=0.0051):
        """Compare a number read out of the manuscript against the computed value."""
        if key not in V:
            flag("XCHECK", f"{label}: no computed value keyed '{key}'"); return
        if not near(found, V[key], tol):
            flag("XCHECK", f"{label}: manuscript {found} vs computed {V[key]}")

    # --- cohort sizes, wherever they are stated ---
    for pat, key, lab in [
            (r"(\d+) European-ancestry IPF Job Exposures Study participants", "Cohort|Analytical N", "abstract N"),
            (r"European-ancestry IPF Job Exposures Study participants \((\d+) cases", "Cohort|Cases", "abstract cases"),
            (r"IPF Job Exposures Study participants \(\d+ cases, (\d+) controls\)", "Cohort|Controls", "abstract controls"),
            (r"analytical sample of (\d+) participants", "Cohort|Analytical N", "Methods N"),
            (r"Cases \(n=(\d+)\)", "Cohort|Cases", "Table 1 cases"),
            (r"Controls \(n=(\d+)\)", "Cohort|Controls", "Table 1 controls"),
            (r"regression N=(\d+)", "Cohort|Analytical N", "Methods regression N"),
            (r"the (\d+)-participant analytical sample", "Cohort|Analytical N", "Table 1 footnote N")]:
        m = re.search(pat, txt)
        if m: cmp(lab, m.group(1), key, tol=0)
        else: flag("XCHECK", f"{lab}: pattern not found in manuscript")

    # --- Hardy-Weinberg, the pair that drifted ---
    m = re.search(r"Hardy-Weinberg testing in controls \(minimum p=([\d.]+)\)", txt)
    if m: cmp("minimum tSNP HWE p", m.group(1), "HWE|minimum tSNP HWE p (controls)")
    else: flag("XCHECK", "minimum tSNP HWE p: sentence not found")

    m = re.search(r"Hardy-Weinberg equilibrium in controls \(p=([\d.]+)[-\u2013]([\d.]+)\)", txt)
    if m:
        hwe = [V[k] for k in ("HWE|MUC5AC H3 p", "HWE|MUC5AC H1 p", "HWE|MUC5B H1 p",
                              "HWE|MUC5AC H2 tag p") if k in V]
        if hwe:
            if not near(m.group(1), min(hwe)): flag("XCHECK", f"HWE range low: manuscript {m.group(1)} vs computed {min(hwe)}")
            if not near(m.group(2), max(hwe)): flag("XCHECK", f"HWE range high: manuscript {m.group(2)} vs computed {max(hwe)}")
    else:
        flag("XCHECK", "haplogroup HWE range: sentence not found")

    # --- linkage disequilibrium figures quoted in Results ---
    m = re.search(r"correlated with rs35705950 T count \(r=([\d.]+)", txt)
    if m: cmp("r(H3, rs35705950)", m.group(1), "LD|r(H3, rs35705950)", tol=0.0011)
    m = re.search(r"among H3 carriers, ([\d.]+)% also carried the T allele", txt)
    if m: cmp("% H3 carriers also T-carriers", m.group(1), "LD|% H3 carriers also T-carriers", tol=0.11)
    m = re.search(r"discordant in (\d+) of \d+ participants", txt)
    if m: cmp("H3/T discordant n", m.group(1), "Discordant|H3/rs35705950 discordant n", tol=0)

    # --- Discussion prose: figures quoted in running text, previously unchecked ---
    for pat, key, tol in [
            (r"Adjusting rs35705950 for H3 barely changes it \(OR ([\d.]+) to",
             "Table 2 / rs alone|rs35705950_T_count OR", 0.0051),
            (r"barely changes it \(OR [\d.]+ to ([\d.]+)\)",
             "Table 2 / 1b H3+rs|rs35705950_T_count OR", 0.0051),
            (r"adjusting H3 for rs35705950 abolishes it \(([\d.]+) to",
             "Table 2 / 1a H3 unadj|H3_count OR", 0.0051),
            (r"abolishes it \([\d.]+ to ([\d.]+)\)",
             "Table 2 / 1b H3+rs|H3_count OR", 0.0051),
            (r"does not improve fit \(\u03c7\u00b2\u2081=([\d.]+)",
             "LRT|add H3 to rs35705950 chi2", 0.051),
            (r"the reverse improves it decisively \(\u03c7\u00b2\u2081=([\d.]+)",
             "LRT|add rs35705950 to H3 chi2", 0.051),
            (r"adjustment for rs35705950 did not abolish it \(([\d.]+) to",
             r"Table 2 / joint \(H1 ref\)|H2t OR", 0.0051),
            (r"did not abolish it \([\d.]+ to ([\d.]+)\)",
             "Table 2 / joint+rs|H2t OR", 0.0051),
            (r"no longer nominally significant \(([\d.]+), p=",
             "H2 conditioned on 8 signals|OR", 0.0051),
            (r"no longer nominally significant \([\d.]+, p=([\d.]+)\)",
             "H2 conditioned on 8 signals|p", 0.0051),
            (r"The \*MUC5B\* H1 inverse association \(([\d.]+)\)",
             "Table 2 / 4a MUC5B H1|MUC5B_H1_count OR", 0.0051),
            (r"attenuates on adjustment \(([\d.]+), p=",
             "Table 2 / 4b MUC5B H1+rs|MUC5B_H1_count OR", 0.0051),
            (r"the adjusted \*MUC5AC\* H3 interval, ([\d.]+)\u2013",
             "Table 2 / 1b H3+rs|H3_count CI-lo", 0.0051),
            (r"H3 interval, [\d.]+\u2013([\d.]+), remains compatible",
             "Table 2 / 1b H3+rs|H3_count CI-hi", 0.0051),
            # IPFJES pooled with TOPMed's published estimate; recomputed in verify_numbers.py
            (r"pooled with TOPMed it narrows to ([\d.]+) \(",
             "TOPMed pooled / 1b H3+rs|H3_count OR", 0.0051),
            (r"pooled with TOPMed it narrows to [\d.]+ \(([\d.]+)–",
             "TOPMed pooled / 1b H3+rs|H3_count CI-lo", 0.0051),
            (r"pooled with TOPMed it narrows to [\d.]+ \([\d.]+–([\d.]+)\)",
             "TOPMed pooled / 1b H3+rs|H3_count CI-hi", 0.0051)]:
        key = key.replace("\\(", "(").replace("\\)", ")")
        m = re.search(pat, txt)
        if m is None:
            flag("XCHECK", f"Discussion: pattern not found for '{key}'")
        else:
            cmp(f"Discussion value for {key}", m.group(1), key, tol=tol)

    # --- Table 1, every cell ---
    # Two rows are labelled "H1 carrier, n (%)", one MUC5AC and one MUC5B, so the
    # sub-heading rows are tracked to tell them apart.
    T1_COUNT = {"Ever-smoked, n (%)": "Ever-smoked", "Asbestos-exposed (Peto), n (%)": "Asbestos",
                "T carrier, n (%)": "rs35705950 T-carrier", "TT homozygote, n (%)": "TT"}
    T1_IQR = {"Age, median (IQR), years": "Age", "Pack-years (smokers), median (IQR)": "Pack-years"}
    sect, in_t1 = None, False
    for line in txt.split("\n"):
        if line.startswith("| Characteristic |"): in_t1 = True; continue
        if in_t1 and not line.startswith("|"): in_t1 = False
        if not in_t1 or not line.startswith("|"): continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 4: continue
        lab = cells[0].replace("*", "").replace("_", "")
        if "haplogroup" in lab or "promoter" in lab:
            sect = "MUC5AC" if "MUC5AC" in lab else ("MUC5B" if "MUC5B" in lab else "rs")
            continue
        if lab in T1_IQR:                       # "77.0 (72.0-82.0)"
            stem = T1_IQR[lab]
            for col, who in ((1, "cases"), (2, "controls")):
                m = re.match(r"([\d.]+) \(([\d.]+)[-\u2013]([\d.]+)\)", cells[col])
                if not m: flag("XCHECK", f"Table 1 '{lab}' {who}: unparsed cell {cells[col]!r}"); continue
                for val, part in ((m.group(2), "Q1"), (m.group(1), "median"), (m.group(3), "Q3")):
                    cmp(f"Table 1 '{lab}' {who} {part}", val, f"Table 1|{stem} {who} {part}", tol=0.051)
            cmp_p(f"Table 1 '{lab}' p", cells[3], f"Table 1|{stem} p")
            continue
        stem = T1_COUNT.get(lab)
        if stem is None and lab.endswith("carrier, n (%)") and sect in ("MUC5AC", "MUC5B"):
            stem = f"{sect} {lab.split()[0]} carrier"
        if stem is None: continue
        for col, who in ((1, "cases"), (2, "controls")):     # "292 (75.8%)"
            m = re.match(r"(\d+) \(([\d.]+)%\)", cells[col])
            if not m: flag("XCHECK", f"Table 1 '{lab}' {who}: unparsed cell {cells[col]!r}"); continue
            cmp(f"Table 1 '{lab}' {who} n", m.group(1), f"Table 1|{stem} {who} n", tol=0)
            pk = f"Table 1|{stem} {who} %"
            if pk in V: cmp(f"Table 1 '{lab}' {who} %", m.group(2), pk, tol=0.051)
        _t1p = {"Ever-smoked": "Table 1|Ever-smoked p", "Asbestos": "Table 1|Asbestos p",
                "rs35705950 T-carrier": "Table 1|rs35705950 T-carrier p",
                "MUC5AC H1 carrier": "Table 1|MUC5AC H1 carrier p",
                "MUC5AC H2 carrier": "Table 1|MUC5AC H2 carrier p",
                "MUC5AC H3 carrier": "Table 1|MUC5AC H3 carrier p",
                "MUC5B H1 carrier": "Table 1|MUC5B H1 carrier p"}
        if stem in _t1p: cmp_p(f"Table 1 '{lab}' p", cells[3], _t1p[stem])

    # --- Supplementary Table S1: stratified estimates ---
    S1 = {("MUC5AC H3", "Non-carrier"): "Table S1|H3_count Non-carrier OR",
          ("MUC5AC H3", "Carrier"):     "Table S1|H3_count Carrier OR",
          ("MUC5B H1", "Non-carrier"):  "Table S1|MUC5B_H1_count Non-carrier OR",
          ("MUC5B H1", "Carrier"):      "Table S1|MUC5B_H1_count Carrier OR",
          ("MUC5AC H2", "Non-carrier"): "Table S1 H2|H2 Non-carrier OR",
          ("MUC5AC H2", "Carrier"):     "Table S1 H2|H2 Carrier OR"}
    for line in txt.split("\n"):
        if not line.startswith("| MUC5"): continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 6: continue
        hg, stratum = cells[0].replace("*", ""), cells[1].split(" (")[0]
        key = S1.get((hg, stratum))
        if key: cmp(f"Table S1 {hg} {stratum} OR", cells[3], key)
        if stratum in ("Non-carrier", "Carrier"):
            cmp(f"Table S1 {hg} {stratum} N", cells[2], f"Table S1|{stratum} N", tol=0)
            _s1p = {"MUC5AC H3": "Table S1|H3_count", "MUC5B H1": "Table S1|MUC5B_H1_count",
                    "MUC5AC H2": "Table S1 H2|H2"}
            if hg in _s1p:
                cmp_p(f"Table S1 {hg} {stratum} p", cells[5], f"{_s1p[hg]} {stratum} p")

    # --- Supplementary Table S3: H2 across tagging measures ---
    for line in txt.split("\n"):
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 7 or not cells[0].startswith("rs"): continue
        rsid = cells[0].split()[0]
        m = re.match(r"([\d.]+) \(([\d.]+)[-\u2013]([\d.]+)\)", cells[5])
        if not m: continue
        for val, part in ((m.group(1), "OR"), (m.group(2), "CI-lo"), (m.group(3), "CI-hi")):
            k = f"H2 direct|{rsid} H2-adj {part}"
            if k in V: cmp(f"Table S3 {rsid} adjusted {part}", val, k)
        cmp_p(f"Table S3 {rsid} adjusted p", cells[6], f"H2 direct|{rsid} H2-adj p")

    # --- Supplementary Table S6: LD with rs35705950 ---
    S6 = {"MUC5AC H3": ("Table S6|H3 controls r2", "Table S6|H3 full-sample r2"),
          "MUC5AC H2": ("Table S6|H2 controls r2", None)}
    for line in txt.split("\n"):
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 4: continue
        hg = cells[0].replace("*", "").split(" (")[0]
        if hg not in S6: continue
        for col, key in ((1, S6[hg][0]), (2, S6[hg][1])):
            if key is None: continue
            m = re.search(r"\(([\d.]+)\)", cells[col])
            if m: cmp(f"Table S6 {hg} col{col} r2", m.group(1), key, tol=0.0051)

    # --- Table 2, every odds ratio and interval ---
    T2 = {"MUC5AC H3, unadjusted": "Table 2 / 1a H3 unadj|H3_count",
          "MUC5AC H3, adjusted for rs35705950": "Table 2 / 1b H3+rs|H3_count",
          "MUC5AC H2, adjusted for H3": "Table 2 / joint (H1 ref)|H2t",
          "MUC5AC H2, adjusted for H3 + rs35705950": "Table 2 / joint+rs|H2t",
          "MUC5B H1, unadjusted": "Table 2 / 4a MUC5B H1|MUC5B_H1_count",
          "MUC5B H1, adjusted for rs35705950": "Table 2 / 4b MUC5B H1+rs|MUC5B_H1_count",
          "rs35705950 alone (reference)": "Table 2 / rs alone|rs35705950_T_count"}
    for line in txt.split("\n"):
        if not line.startswith("|"): continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 5: continue
        label = cells[0].replace("*", "")
        if label not in T2: continue
        stem = T2[label]
        cmp(f"Table 2 '{label}' OR", cells[2], stem + " OR")
        ci = re.match(r"([\d.]+)[-\u2013]([\d.]+)", cells[3])
        if ci:
            cmp(f"Table 2 '{label}' CI-lo", ci.group(1), stem + " CI-lo")
            cmp(f"Table 2 '{label}' CI-hi", ci.group(2), stem + " CI-hi")
        cmp_p(f"Table 2 '{label}' p", cells[4], stem + " p")


# ---------------------------------------------------------------------------
# 7. Reference integrity: every citation has an entry, every entry is cited.
#    Found by hand twice in one session (Okuda orphaned by a Discussion cut,
#    Abrami by a Methods trim), so it is worth doing automatically.
# ---------------------------------------------------------------------------
_body, _sep, _rest = txt.partition("## References")
if not _sep:
    flag("REFS", "no '## References' heading found")
else:
    _reflist = _rest.split("## Supplementary methods")[0]
    _entries = [l.strip() for l in _reflist.split("\n")
                if l.strip() and not l.startswith(("#", "-", "*", "|"))]
    _surnames = []
    for _e in _entries:
        _m = re.match(r"([A-Z][A-Za-z\-\u2019]+)", _e)
        if _m: _surnames.append((_m.group(1), _e))
    # citations look like (Author, 2024), (Author et al., 2024), (A and B, 2012)
    # "Li and Ji, 2005" and "Li & Ji, 2005" both cite Li; the lookbehinds stop the
    # second author being read as a citation in its own right.
    _CITE = (r"(?<!and )(?<!& )([A-Z][A-Za-z\-\u2019]+)"
             r"(?: (?:and|&) [A-Z][A-Za-z]+)?(?: et al\.)?,\s*(?:19|20)\d\d")
    _cited = set(re.findall(_CITE, _body))
    # the supplementary methods cite things too, and those are legitimate
    _cited |= set(re.findall(_CITE, _rest.split("## Supplementary methods")[-1]))
    for _s, _e in _surnames:
        if _s not in _cited:
            flag("REFS", f"reference listed but never cited: {_e[:70]}")
    for _c in sorted(_cited):
        if _c not in [s for s, _ in _surnames]:
            flag("REFS", f"cited but no reference entry: {_c}")

# ---------------------------------------------------------------------------
# 8. Mechanism-word lint. Not a correctness check: it lists the sentences whose
#    wording has repeatedly named the wrong mechanism (LD vs power vs
#    measurement, clade vs promoter variant, replication vs reproduction), so a
#    reread can be targeted rather than exhaustive. Reported separately from
#    issues, and never fails the run.
# ---------------------------------------------------------------------------
_MECH = re.compile(r"\b(independent(?:ly)?|separabl[ey]|separate|co-inherited|"
                   r"resolve[sd]?|measurement|power(?:ed)?|replicat\w+|reproduc\w+)\b", re.I)
_mech_hits = []
for _para in _body.split("\n\n"):
    if _para.startswith(("|", "#", "*Analysis code", "*Support statement")): continue
    for _s in re.split(r"(?<=[.!?])\s+", _para.strip()):
        if len(_s) > 40 and _MECH.search(_s):
            _mech_hits.append(re.sub(r"\s+", " ", _s))

# ---------------------------------------------------------------------------
# 9. Submission readiness.
# ---------------------------------------------------------------------------
for _need, _lab in [(r"\*Support statement:", "support/funding statement"),
                    (r"[Cc]ode and data availability", "code and data availability statement"),
                    (r"[Cc]onflict of interest|[Cc]ompeting interests", "conflict-of-interest statement")]:
    if not re.search(_need, txt): flag("SUBMIT", f"missing: {_lab}")
_defined = sorted(set(re.findall(r"\*\*Supplementary Table S(\d)", txt)), key=int)
_citedST = sorted({x for pair in re.findall(r"Supplementary Table S(\d)|\(Table S(\d)", _body)
                   for x in pair if x}, key=int)
for _n in _defined:
    if _n not in _citedST: flag("SUBMIT", f"Supplementary Table S{_n} is never cited in the main text")
for _n in _citedST:
    if _n not in _defined: flag("SUBMIT", f"Supplementary Table S{_n} cited but not defined")
if _defined and [int(x) for x in _defined] != list(range(1, len(_defined) + 1)):
    flag("SUBMIT", f"supplementary table numbering not contiguous: {_defined}")
for _stale in ("Word count", "Abstract word count", "TODO", "FIXME", "XXX"):
    if _stale in txt: flag("SUBMIT", f"stale working note in manuscript: '{_stale}'")

print(f"Internal-consistency audit of {M.split('/')[-1]}")
if issues:
    print(f"\n{len(issues)} potential issue(s):")
    for cls, msg in issues:
        print(f"  [{cls}] {msg}")
else:
    print("\nOK - no issues (additive sums, fraction<->%, CI ordering/bracketing, OR-quote drift,\n"
          "     CI-vs-p significance, manuscript numbers vs verify_numbers, references,\n"
          "     submission readiness)")

# The mechanism lint is advisory, never an issue. It exists because every wording
# error in this manuscript's history has been a sentence that read like the
# paper's thesis while naming the wrong mechanism: linkage disequilibrium where
# the barrier was power, measurement where it was co-inheritance, MUC5B clades
# described as differing in length, replication conflated with reproduction.
# Printed only when asked, so it does not drown the issue list.
if os.environ.get("MECH_LINT"):
    print(f"\nMechanism-word sentences to reread ({len(_mech_hits)}):")
    for _s in _mech_hits:
        print(f"  - {_s[:150]}")
else:
    print(f"\n({len(_mech_hits)} sentences use mechanism words "
          f"(independent / separable / co-inherited / resolve / measurement / power /\n"
          f" replicate / reproduce). Run with MECH_LINT=1 to list them for a targeted reread.)")

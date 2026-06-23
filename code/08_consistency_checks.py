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

"""08_consistency_checks.py — internal-consistency audit of the manuscript prose.

Third verification layer, complementing:
  - verify_numbers.py        : manuscript numbers vs the code that produced them
  - 07_foundational_audit.py : haplogroup labels / allele orientations vs sources
This script checks the manuscript's INTERNAL self-consistency — statements compared
to each other — which neither of the above catches. It flags:
  1. additive cohort breakdowns that don't sum
  2. a percentage that doesn't match its own fraction (n/N)
  3. Table 1 cells whose n (x%) disagrees with the column N (428 cases / 401 controls)
  4. a confidence interval that is reversed or does not bracket its OR
  5. the same OR quoted with two different CIs (edit-drift)
  6. a CI and p-value that disagree on significance

Run: uv run --python .venv python3 code/08_consistency_checks.py
(The cohort sums in check 1 are hard-coded to the current data freeze — update them
if the cohort numbers change.)
"""
import re
from collections import defaultdict

M = "./MUC5_VNTR_haplogroup_paper_draft.md"
txt = open(M).read()
issues = []
def flag(cls, msg): issues.append((cls, msg))

# 1. Additive cohort relationships (must sum AND the target number must be present)
sums = {
    "genotyped panel 2,912 = 428 cases + 2,481 controls + 3 unset": (428 + 2481 + 3, 2912),
    "controls 2,481 = 401 IPFJES + 2,080 UK Biobank":               (401 + 2080, 2481),
    "analytical 829 = 428 cases + 401 controls":                    (428 + 401, 829),
    "regression 798 = 829 - 31 (cases-only centre)":                (829 - 31, 798),
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

# 3. Table 1 cells: "n (x%)" against cases N=428 / controls N=401
for line in txt.splitlines():
    if line.startswith("|") and "%" in line:
        cells = [c.strip() for c in line.split("|")]
        for idx, N, who in [(2, 428, "case"), (3, 401, "control")]:
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
print(f"Internal-consistency audit of {M.split('/')[-1]}")
if issues:
    print(f"\n{len(issues)} potential issue(s):")
    for cls, msg in issues:
        print(f"  [{cls}] {msg}")
else:
    print("\nOK - no issues (additive sums, fraction<->%, CI ordering/bracketing, OR-quote drift, CI-vs-p significance)")

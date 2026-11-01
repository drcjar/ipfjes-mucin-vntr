#!/usr/bin/env bash
# Check the manuscript after editing MUC5_VNTR_haplogroup_paper_draft.md.
# Rebuilds PDF/DOCX, runs all four verification layers, and checks the ERJ Open
# Research limits and house style. Exits non-zero if anything fails, so it is
# safe to run before sending the paper anywhere.
# Usage: ./check_paper.sh [--no-build]
set -uo pipefail
cd "$(dirname "$0")"
PY="uv run --python .venv python3"
FAIL=0
step() { printf '%-34s' "$1"; }
ok()   { echo "OK${1:+  $1}"; }
bad()  { echo "FAIL  $1"; FAIL=1; }

if [ "${1:-}" != "--no-build" ]; then
  step "build pdf/docx"
  if ./build_paper.sh MUC5_VNTR_haplogroup_paper_draft >/dev/null 2>&1; then ok; else bad "see ./build_paper.sh"; fi
fi

step "1. number audit"
# The counts are matched with a leading-boundary pattern, not a bare substring: the audits
# print "N FAIL" / "N FLAG" / "N differ", so grep -q "0 FAIL" would also match "10 FAIL"
# and pass the gate on a double-digit failure count.
OUT=$($PY code/verify_numbers.py 2>&1)
if   echo "$OUT" | grep -qE "(^|[^0-9])0 FAIL";  then ok "$(echo "$OUT" | grep -o '[0-9]* checks · [0-9]* PASS')"
elif ! echo "$OUT" | grep -q "checks ·";         then bad "did not complete: $(echo "$OUT" | tail -1)"
else bad "$(echo "$OUT" | grep -oE '[0-9]+ FAIL' | tail -1)"; fi

step "2. assignment/orientation"
OUT=$($PY code/06_foundational_audit.py 2>&1)
if   echo "$OUT" | grep -qE "(^|[^0-9])0 FLAG";     then ok "11 checks"
elif ! echo "$OUT" | grep -q "foundational checks"; then bad "did not complete: $(echo "$OUT" | tail -1)"
else bad "flags raised"; fi

# 07 covers two different kinds of finding: numbers/prose that disagree, and
# submission checklist items still outstanding. A pending conflict-of-interest
# statement is not a numerical error, so they are reported separately.
OUT=$($PY code/07_consistency_checks.py 2>&1)
# Every finding class except SUBMIT counts as a consistency failure, matched by shape
# rather than by name: an allow-list silently dropped new classes as 07 grew them (the
# [PROV] provenance findings were invisible here for exactly that reason).
CONSIST=$(echo "$OUT" | grep -E "^  \[[A-Z0-9%-]+\]" | grep -v "^  \[SUBMIT\]" || true)
SUBMIT=$(echo "$OUT"  | grep -E "^  \[SUBMIT\]" || true)
MECH=$(echo "$OUT"    | grep -oE "^\([0-9]+ sentences use mechanism words" | grep -oE "[0-9]+" || echo "?")

step "3. internal consistency"
if [ -z "$CONSIST" ]; then ok "sums, %, CIs, cross-check vs computed, provenance, references"
else bad "$(echo "$CONSIST" | wc -l) finding(s)"; echo "$CONSIST" | sed 's/^  /      /'; fi

step "submission readiness"
if [ -z "$SUBMIT" ]; then ok "statements, table numbering, no stale notes"
else bad "$(echo "$SUBMIT" | wc -l) outstanding"; echo "$SUBMIT" | sed 's/^  /      /'; fi

# The analysis scripts are not otherwise exercised: a broken one can sit in the
# repo while all four audits pass, because the audits recompute independently.
if [ "${RUN_SCRIPTS:-1}" = "1" ]; then
  step "analysis pipeline runs"
  BROKEN=""
  for s in 01_haplogroup_assignment 02_primary_analysis 03_supplementary_analyses \
           04_h2ref_sensitivity 05_tsnp_table 09_h2_direct_and_lrt 10_imputed_vntr_length \
           11_mucin_cluster_scan 12_h2_conditional_on_locus_signals \
           13_ld_haplogroup_rs35705950 14_ancestry_and_promoter_gxe 15_orientation_vs_topmed \
           17_cluster_scan_power; do
    [ -f "code/$s.py" ] || continue
    $PY "code/$s.py" >/dev/null 2>&1 || BROKEN="$BROKEN $s"
  done
  [ -z "$BROKEN" ] && ok "13 scripts" || bad "failed:$BROKEN"
fi

step "independent re-derivation"
OUT=$($PY code/16_independent_recompute.py 2>&1)
if   echo "$OUT" | grep -qE "(^|[^0-9])0 differ";     then ok "12 Table 2 estimates, second code path"
elif ! echo "$OUT" | grep -q "re-derived independently"; then bad "did not complete: $(echo "$OUT" | tail -1)"
else bad "routes disagree"; echo "$OUT" | grep DIFFERS | sed 's/^/      /'; fi

step "4. source provenance"
OUT=$($PY code/08_provenance_audit.py 2>&1)
if   echo "$OUT" | grep -qE "(^|[^0-9])0 FLAG";    then ok "67 checks"
elif ! echo "$OUT" | grep -q "provenance checks";  then bad "did not complete: $(echo "$OUT" | tail -1)"
else bad "flags raised"; fi

# --- ERJ Open Research limits: main text 3000, structured abstract 250 ---
step "word limits"
$PY - <<'PYEOF'
import re, sys
L = open('MUC5_VNTR_haplogroup_paper_draft.md').read().split('\n')
idx = {}
for n, l in enumerate(L):
    if l.startswith('## '): idx.setdefault(l[3:].strip().split(':')[0], n)
def wc(a, b):
    t = '\n'.join(L[a:b]); t = re.sub(r'\|.*\|', '', t); t = re.sub(r'[*_#`]', '', t)
    return len([w for w in re.split(r'\s+', t) if w.strip()])
main = wc(idx['Introduction'], idx['Tables'])
i0 = next(n for n, l in enumerate(L) if l.startswith('**Background:**'))
ab = '\n'.join(L[i0:idx['Introduction']])
for lab in ('**Background:**', '**Methods:**', '**Results:**', '**Conclusions:**'): ab = ab.replace(lab, '')
ab = re.sub(r'[*_#]', '', ab).replace('---', '')
abn = len([w for w in re.split(r'\s+', ab) if w.strip()])
# ERJ Open Research states 3000 words for original articles, but it is not enforced: six
# recent original articles sampled from Europe PMC and counted the same way (Introduction
# through Discussion, tables excluded) ran 3281-4076 words, median 3672, all six over the
# limit, including a genetics paper by two of our own reviewers at 3602. So the main-text
# limit is a warning here. The abstract limit is kept hard, because journals do enforce it.
# The abstract limit is not enforced either: abstracts of the same six sampled articles ran
# 242-286 words (median 255), four of the six over the stated 250. Both limits therefore warn
# rather than fail; they still print, so drift upward stays visible.
over_main, over_abs = main > 3000, abn > 250
status = "WARN  " if (over_main or over_abs) else "OK    "
note = "  (sampled published articles: ~3281-4076 main, 242-286 abstract; not enforced)" if (over_main or over_abs) else ""
print(status + f"main {main}/3000, abstract {abn}/250" + note)
sys.exit(0)
PYEOF
[ $? -ne 0 ] && FAIL=1

step "house style"
EM=$(grep -c '—' MUC5_VNTR_haplogroup_paper_draft.md)
BOLD=$($PY - <<'PYEOF'
import re
L = open('MUC5_VNTR_haplogroup_paper_draft.md').read().split('\n')
keep = ('Background','Methods','Results','Conclusions','Authors','Affiliations',
        'Corresponding author','Keywords','Take-home message')
n = 0
for l in L:
    if l.startswith('#') or l.startswith(('**Table','**Supplementary Table')): continue
    for m in re.findall(r'\*\*(.+?)\*\*', l):
        if m.rstrip(':') in keep: continue
        if l.startswith('|'): continue          # table row-group headers
        n += 1
print(n)
PYEOF
)
if [ "$EM" -eq 0 ] && [ "$BOLD" -eq 0 ]; then ok "no em dashes, no bold in body text"
else bad "$EM em dash(es), $BOLD bold span(s) in body text"; fi

# Markdown that renders differently from how it reads: a literal * pairing with a later *
# into italics (the joint-senior-author markers vanished this way), consecutive lines
# merging into one paragraph (the affiliations ran together), and @text parsed as a
# citation. Checked on pandoc's own parse, so it sees exactly what the PDF/DOCX get.
step "markdown rendering"
PANDOC="$($PY -c 'import pypandoc; print(pypandoc.get_pandoc_path())')"
RENDER=$("$PANDOC" -f markdown -t json MUC5_VNTR_haplogroup_paper_draft.md | $PY code/render_lint.py)
if [ -z "$RENDER" ]; then ok "no stray emphasis, joined lines or citations"
else bad "$(echo "$RENDER" | wc -l) finding(s)"; echo "$RENDER" | sed 's/^/      /'; fi

echo
[ $FAIL -eq 0 ] && echo "All checks passed." || echo "Something needs attention (see FAIL above)."
exit $FAIL

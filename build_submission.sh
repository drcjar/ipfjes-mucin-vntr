#!/usr/bin/env bash
# Build the ERJ Open Research submission set from the single master markdown.
#
# The master file keeps the supplement inline, because check_paper.sh verifies
# the supplementary tables in place. Submission needs them separated: ERJ counts
# figures and tables against a combined limit of 8, and accepts an online
# supplement as its own file. So we split at build time, not in the source.
#
# Outputs into submission/:
#   main_document.docx    Title through the statements, Tables 1 and 2 only
#   online_supplement.docx  Supplementary methods and Tables S1-S6
set -euo pipefail
cd "$(dirname "$0")"
SRC=MUC5_VNTR_haplogroup_paper_draft.md
OUT=submission
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$OUT"
PANDOC="$(uv run --python .venv python3 -c 'import pypandoc; print(pypandoc.get_pandoc_path())')"

uv run --python .venv python3 - "$SRC" "$TMP" <<'PYEOF'
import sys, pathlib
src, tmp = sys.argv[1], pathlib.Path(sys.argv[2])
L = open(src).read().split('\n')
cut = next(n for n, l in enumerate(L) if l.startswith('## Supplementary methods'))
# Drop the trailing '---' rule that separated the two halves.
main = L[:cut]
while main and main[-1].strip() in ('', '---'):
    main.pop()
title = L[0][2:].strip()
supp = ['# Online supplement', '', '*%s*' % title, ''] + L[cut:]
(tmp / 'main_document.md').write_text('\n'.join(main) + '\n')
(tmp / 'online_supplement.md').write_text('\n'.join(supp) + '\n')
print('main %d lines, supplement %d lines' % (len(main), len(supp)))
PYEOF

for f in main_document online_supplement; do
  "$PANDOC" "$TMP/$f.md" -o "$OUT/$f.docx" --resource-path=.
done
cp cover_letter_erj.docx STROBE_checklist.docx "$OUT/"
echo "built: $OUT/main_document.docx  $OUT/online_supplement.docx  (+ cover letter, STROBE)"

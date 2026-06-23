#!/usr/bin/env bash
# Regenerate manuscript PDF (for reading/checking) and DOCX (for Google Docs + Paperpile editing)
# from the markdown source. PDF/DOCX are build artifacts — regenerate on demand, don't hand-edit.
# Usage: ./build_paper.sh [basename-without-extension]   (default: MUC5_VNTR_haplogroup_paper_draft)
set -euo pipefail
cd "$(dirname "$0")"
BASE="${1:-MUC5_VNTR_haplogroup_paper_draft}"
PANDOC="$(uv run --python .venv python3 -c 'import pypandoc; print(pypandoc.get_pandoc_path())')"

# DOCX — clean Word heading/table styles for Google Docs import + Paperpile citation insertion.
"$PANDOC" "$BASE.md" -o "$BASE.docx" --resource-path=.

# PDF — DejaVu Serif gives full Unicode coverage so superscript p-values (e.g. ×10⁻²⁴) render.
"$PANDOC" "$BASE.md" -o "$BASE.pdf" --pdf-engine=xelatex \
  -V mainfont="DejaVu Serif" -V geometry:margin=2.5cm -V fontsize=11pt --resource-path=.

echo "built: $BASE.pdf  $BASE.docx"

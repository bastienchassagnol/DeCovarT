#!/usr/bin/env bash
# Render the review Quarto document, then shrink the PDF with Ghostscript.
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
cd "$root"
stem=retour_figure_1_and_4
rendered="${stem}_render.pdf"
compressed="${stem}.pdf"
preset="${PDFSETTINGS:-/ebook}"

quarto render "${stem}.qmd"

if ! command -v gs >/dev/null 2>&1; then
  echo "Ghostscript (gs) is required; install with: sudo apt install ghostscript" >&2
  exit 1
fi

mv -f "$compressed" "$rendered"
gs -sDEVICE=pdfwrite \
  -dCompatibilityLevel=1.4 \
  -dPDFSETTINGS="${preset}" \
  -dNOPAUSE -dQUIET -dBATCH \
  -sOutputFile="$compressed" \
  "$rendered"

ls -lh "$compressed" "$rendered"

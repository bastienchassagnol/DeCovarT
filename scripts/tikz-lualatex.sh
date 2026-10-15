#!/usr/bin/env bash
# LuaLaTeX for the pandoc-ext/diagram TikZ engine. The filter calls this
# script with pdflatex's arguments (-interaction=nonstopmode,
# -output-directory DIR, file.tex); LuaLaTeX accepts the same interface.
# After compiling, Ghostscript turns every glyph into an outline, so the SVG
# that Inkscape makes from the PDF does not depend on TeX fonts being
# installed in the reader's browser (quarto-dev/quarto-cli#13113).
# Needs TeX Live (lualatex), Ghostscript and Inkscape.
set -euo pipefail

lualatex "$@"

outdir="."
tex=""
while [ $# -gt 0 ]; do
  case "$1" in
    -output-directory) outdir="$2"; shift 2 ;;
    -*) shift ;;
    *) tex="$1"; shift ;;
  esac
done
pdf="$outdir/$(basename "${tex%.tex}").pdf"
gs -q -dNOPAUSE -dBATCH -dNoOutputFonts -sDEVICE=pdfwrite \
  -o "$pdf.outlined" "$pdf"
mv "$pdf.outlined" "$pdf"

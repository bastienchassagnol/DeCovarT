#!/usr/bin/env bash
# Pandoc's gfm writer inserts a space after fence openers (``` mermaid).
# GitHub requires ```mermaid with no space.
set -euo pipefail
readme="${1:-README.md}"
sed -i \
  -e 's/^``` mermaid$/```mermaid/' \
  -e 's/^``` sh$/```sh/' \
  -e 's/^``` html$/```html/' \
  -e 's/^``` r$/```r/' \
  "$readme"

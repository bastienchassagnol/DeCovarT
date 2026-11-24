---
description: How to find, draw, or generate a figure
alwaysApply: true
---
# Figure sources

When a figure is needed, try these in order. Stop at the first one that works.

## 1. Reuse a licensed image

Search the web for a figure of the exact concept. Use only public domain or
Creative Commons (CC BY or CC BY-SA). Prefer Wikimedia Commons, or an
open-access paper that states the licence on the article. Skip paywalled
journal figures and any image with no stated licence.

For each candidate, give the page URL, the direct image URL, the licence, the
author and a one-line attribution. Do not download anything until the user
picks one.

## 2. Draw it (code)

Follow the repository's figure scripts and seeds.

- **Data plots:** R (`ggplot2`) or OJS/D3 in Quarto, per `data-visualisation`
  (SVG first; a PDF for LaTeX is made from the SVG).
- **Diagrams, DAGs, plates, schematics, outlines:** TikZ, `pgfplots` or
  `curve2e` under LuaLaTeX, per `latex-scientific` (and `graphical-models`
  for probabilistic models). No PSTricks.

## 3. Generate it

If neither a licensed image nor a code-drawn figure will do, generate the
image with the session's image generator. Say that the picture is generated,
and do not present it as a published figure.

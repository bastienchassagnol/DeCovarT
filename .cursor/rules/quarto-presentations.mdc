---
description: Quarto authoring — labels, callouts, code display, tables, Mermaid, slides routing
globs: "**/*.{qmd,Rmd},**/_quarto.yml,**/_publish.yml"
alwaysApply: false
paths:
  - "**/*.{qmd,Rmd}"
  - "**/_quarto.yml"
  - "**/_publish.yml"
---
# Quarto: documents, websites and reveal.js

Reproducible execution; cite claims at point of use; native Quarto/Pandoc
features before custom JS.

Route by output type before major edits:

- **Slides** (`format: revealjs`): **quarto-revealjs** (theme, layout, motion,
  export, OJS/D3/React in decks).
- **Website** (`project.type: website`): **quarto-website**.
- **Book** (`project.type: book`): **quarto-book**.

## Code display

- Documents, books and vignettes fold code by default: `code-fold: true`,
  `code-summary: "Show the code"` and `code-tools: true` under
  `format: html:` in `_quarto.yml` (once per project). A chunk that must
  stay open sets `#| code-fold: false`. Slides are not folded.

## Markdown and cross-references

- Labels: `#fig-…`, `#tbl-…`, `#eq-…`, `{#sec-…}`; cite mid-sentence (`@fig-…`).
- Lists: standard Markdown only — no Pandoc example lists, task checkboxes, or
  unicode step markers (Visual Editor safe).
- Package names in backticks; citations `@key` unchanged.
- Callouts for real warnings only: `::: {#nte-… .callout-warning}` with
  `@nte-…` references. `echo: false` when readers need output only.
- Emoji: `:shortcode:` with `from: markdown+emoji` when needed.

## Tables

- Static output (PDF, manuscript, static HTML, slides): `tinytable::tt()`.
- `reactable` only when HTML interactivity (sort, filter, search) helps, with
  `pagination = FALSE` on slides.
- Hand-written LaTeX: `booktabs` (`latex-scientific`). No `knitr::kable`.

## Mermaid

Executable chunks:

````markdown
```{mermaid}
%%| label: fig-example
%%| fig-cap: "Caption."
flowchart TD
  A --> B
```
````

`eval: true` in site config or chunk. Inspect `pre.mermaid` and the
quarto-diagram JS before client-side fallbacks. Mermaid theme `neutral`
(`%%{init: {'theme': 'neutral'}}%%`) unless the project sets one. Mermaid is
for workflows; probabilistic models use `graphical-models`.

## Pages that run R

Set seeds in the page; end it with `sessionInfo()` in a collapsed callout
(quarto-website skill, `references/web-layout.md`).

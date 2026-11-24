---
description: Plate notation for probabilistic models in Quarto TikZ and LaTeX
globs: "**/*.{qmd,Rmd,tex},**/_quarto.yml"
alwaysApply: false
paths:
  - "**/*.{qmd,Rmd,tex}"
  - "**/_quarto.yml"
---
# Probabilistic graphical models

Use **classical plate notation**, implemented with the TikZ **`bayesnet`**
library ([tikz-bayesnet](https://github.com/jluttine/tikz-bayesnet)), for every
generative DAG in LaTeX and Quarto. One shared preamble; same node styles in PDF
and SVG. Workflows (pipelines, decision paths) stay **Mermaid**. Primary **causal
DAG** figures use causal-inference notation when causality is the subject, not
this generative template.

**Terminology:** PGM / directed graphical model; plate = indexed repetition, not
automatic independence — nested plates and the sampling design decide dependence.
Arrows are parent-to-child dependencies in the model spec, not causal claims
unless the text says so.

## Symbols (native `bayesnet`)

| Style | Shape | Use |
| --- | --- | --- |
| `latent` | unshaded circle | latent stochastic variable or unknown parameter |
| `obs` | grey circle | observed stochastic variable |
| `det` | unshaded diamond | deterministic quantity ($\mu_i$, $\eta$, link) |
| `const` | text or framed square | fixed input; `\tikzset{fixed/.style={const, draw, …}}` when framed |
| `factor` | filled square | factor graph only |
| `\plate{name}{(fit)(list)}{label}` | solid rounded rectangle | indexed subgraph; index at bottom right |
| `note` (extension) | text beside node | equation or distribution; not a random variable |

Monochrome by default; colour may highlight a meaningful layer only if the
article already defines it — do not replace shape with colour. Equations and
likelihoods in `note` nodes, never inside stochastic/deterministic nodes.
Include every enclosed node and in-plate `note` in the plate fit list. One
diagram per float or slide.

**Layout:** top-down generative hierarchy when possible — globals on top, group
latents, observations inside inner plates; minimal crossings; consistent
`node distance`.

## Shared TikZ preamble

Quarto: `diagram.engine.tikz.header-includes` in `_quarto.yml`. LaTeX: macros
preamble once. The first five lines are the `latex-scientific` font block
(STIX Two Math, TeX Gyre Heros; `\symbf`, not `bm`), repeated so the preamble
can be pasted whole ([quarto-cli#13113](https://github.com/quarto-dev/quarto-cli/issues/13113)).

```latex
\usepackage{amsmath}
\usepackage{fontspec}
\usepackage{unicode-math}
\setmathfont{STIX Two Math}
\setmainfont{TeX Gyre Heros}[Ligatures=CommonOff] % figures; \setsansfont in a manuscript
\usepackage{tikz}
\usetikzlibrary{bayesnet,positioning,fit,calc,arrows.meta}
\tikzset{
  fixed/.style={const, draw, minimum size=7mm, inner sep=2pt},
  note/.style={font=\footnotesize, align=left},
  dependency/.style={-{Stealth[length=2mm]}, semithick}
}
```

Use `\draw[dependency]` or `\draw[->]`, not `\edge` (nullfont `;` warnings);
do not redefine `latent`, `obs`, `det`, or `\plate`. Only the outermost plate
is the independent sampling unit; inner plates index repeated measures within
it.

## Quarto (`diagram` filter)

[`pandoc-ext/diagram`](https://github.com/pandoc-ext/diagram): `{.tikz}` block
with `%%| filename:` and `%%| alt:`. Tall decks: `.fit-fig` when CSS fixes
image width to `auto`.

````markdown
```{.tikz}
%%| filename: hier-gaussian
%%| alt: Global beta and sigma generate mu_i from x_i, then y_i, over i.
\begin{tikzpicture}[node distance=12mm]
  \node[latent] (beta) {$\beta$};
  \node[latent, right=of beta] (sigma) {$\sigma$};
  \node[det, below=of beta] (mu) {$\mu_i$};
  \node[const, left=of mu] (x) {$x_i$};
  \node[obs, below=of mu] (y) {$y_i$};
  \node[note, right=of mu] (mean) {$\mu_i=x_i^\top\beta$};
  \node[note, right=of y] (lik) {$y_i\sim\mathcal{N}(\mu_i,\sigma^2)$};
  \draw[dependency] (beta) -- (mu);
  \draw[dependency] (x) -- (mu);
  \draw[dependency] (mu) -- (y);
  \draw[dependency] (sigma) -- (y);
  \plate{obs}{(x)(mu)(y)(mean)(lik)}{$i=1,\ldots,n$}
\end{tikzpicture}
```
````

Filter builds `standalone`; default binary name is `pdflatex` but passes
`-interaction=nonstopmode` and `-output-directory` — **point `TIKZ_BIN` at a
wrapper** that runs **LuaLaTeX** (XeLaTeX acceptable), then outlines fonts.
Web output is the portable SVG only (no editable variant):

```bash
#!/usr/bin/env bash
set -euo pipefail
lualatex "$@"
outdir="."; tex=""
while [ $# -gt 0 ]; do   # filter passes: -interaction=..., -output-directory DIR, file.tex
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
```

The filter converts the PDF to SVG; Ghostscript `-dNoOutputFonts` keeps the
labels portable. Set `TIKZ_BIN` in `_environment.local`. Cache key is source only — clear cache
after preamble edits.

## LaTeX articles

Keep each `tikzpicture` in its own `standalone` source and follow the TikZ
SVG-first workflow in `latex-scientific` (editable SVG, portable SVG, PDF from
the portable SVG).

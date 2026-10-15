# Probabilistic graphical models

Draw a DAG or a generative model in TikZ, with plate notation. A workflow that
is not a model (a pipeline, a decision path) stays Mermaid in Quarto.

Symbols follow the
[RevBayes graphical-model convention](https://revbayes.github.io/tutorials/intro/getting_started#probabilistic-graphical-models).
Define the styles and the `\plate` macro once, not on each picture.

| Style | Shape | Use |
| --- | --- | --- |
| `const` | square, white | fixed quantity: a covariate, a design value, a parameter treated as fixed |
| `stoch` | blue circle | latent random variable |
| `obs` | grey circle | observed (clamped) random variable |
| `det` | dotted circle | quantity computed from its parents ($\eta$, $\mu$) |
| `var` | grey ellipse | measured variable in a causal graph |
| `plate` | dashed rounded rectangle | a repeated subgraph, labelled with its index |
| `note` | small text beside a node | an equation or a word; not a random variable |

- Put equations in `note` nodes. Do not write an equation inside a node.
- A nested plate repeats inside the outer one. The outer plate is the
  independent unit. The inner plate is its correlated replicates. Only the
  outer plate counts as replication.
- Include a `note` in the plate's fit list when the note sits inside the plate.
- One diagram per slide or per float. Do not invent a second palette.

The macro is `\plate{name}{(fit)(list)}{label}`. The fit list is the node
names the rectangle must enclose, each in parentheses. The label is the index,
for example `$i = 1,\ldots,n$`.

```latex
\plate{units}{(uj)(yj)}{independent units $j$}
\plate{reps}{(yi)}{correlated replicates $i$}
```

Shared preamble fragment. In Quarto this belongs under
`diagram.engine.tikz.header-includes` in `_quarto.yml`. In a LaTeX article it
belongs in the macros block of the preamble. Load `amsmath`, `amssymb`,
`amsfonts`, and `bm` as well: node labels use that mathematics
([quarto-cli#13113](https://github.com/quarto-dev/quarto-cli/issues/13113)).

```latex
\usetikzlibrary{arrows.meta,positioning,fit,backgrounds,calc,shapes.geometric}
\definecolor{ink}{HTML}{1F2937}
\definecolor{latfill}{HTML}{DAE9FB}
\definecolor{latdraw}{HTML}{4477B5}
\definecolor{obsfill}{HTML}{C9CDD3}
\definecolor{obsdraw}{HTML}{747B86}
\definecolor{platedraw}{HTML}{98A2B3}
\tikzset{
  every node/.style={text=ink},
  const/.style={rectangle, draw=ink, fill=white, minimum size=9mm, line width=0.8pt},
  stoch/.style={circle, draw=latdraw, fill=latfill, minimum size=11mm, line width=0.9pt},
  obs/.style={circle, draw=obsdraw, fill=obsfill, minimum size=11mm, line width=0.9pt},
  det/.style={circle, draw=ink, densely dotted, fill=white, minimum size=11mm, line width=0.9pt},
  var/.style={ellipse, draw=obsdraw, fill=obsfill, minimum height=9mm, inner xsep=2mm, line width=0.9pt},
  plate/.style={draw=platedraw, dashed, rounded corners=5pt, inner sep=4mm, line width=0.8pt},
  platelabel/.style={anchor=north east, font=\small\itshape, text=obsdraw, inner sep=1mm},
  note/.style={font=\small, text=ink!80},
  edge/.style={-{Stealth[length=2.6mm]}, line width=0.8pt, draw=ink}
}
\newcommand{\plate}[3]{%
  \begin{scope}[on background layer]
    \node[plate, fit=#2] (#1) {};
  \end{scope}
  \node[platelabel] (#1-lab) at (#1.south east) {#3};}
```

## Quarto

Use the [`diagram` filter](https://github.com/pandoc-ext/diagram). A diagram is
a `tikz` block with a filename and alt text. `%%| width:` has no effect when
the deck CSS sets image width to `auto`; add `.fit-fig` to a slide that must
shrink a tall diagram.

````markdown
```{.tikz}
%%| filename: dag-example
%%| alt: A covariate x determines the mean mu, which generates the observed y, repeated over i.
\node[const] (x) {$x$};
\node[det, right=of x] (mu) {$\mu$};
\node[obs, right=of mu] (y) {$y$};
\node[note, below=2mm of mu] (n) {$\mu = g^{-1}(x\beta)$};
\draw[edge] (x) -- (mu);
\draw[edge] (mu) -- (y);
\plate{p}{(x)(mu)(y)(n)}{$i = 1,\ldots,n$}
```
````

The filter writes a `standalone` document and compiles it. Its fallback binary
is named `pdflatex`, and it always passes `-interaction=nonstopmode` and
`-output-directory`. Those flags are the ordinary TeX Live interface: `xelatex`
and `lualatex` accept them. pdfLaTeX is not required by TikZ or by the plate
styles above. Point `TIKZ_BIN` at a wrapper whose first command is the engine
you want. Prefer LuaLaTeX, which is the engine for articles; XeLaTeX is
acceptable. Keep the wrapper's remaining steps.

Inkscape turns the PDF into an SVG and can leave the labels as live text in a
TeX font the browser does not have. After the engine finishes, run Ghostscript
with `-dNoOutputFonts` so every glyph becomes an outline, then replace the PDF
with that file. That outline step is what makes the SVG portable. Changing the
engine does not remove it. The wrapper needs TeX, Ghostscript, and Inkscape.

```bash
#!/usr/bin/env bash
# TIKZ_BIN. The filter calls this with pdflatex's arguments.
set -euo pipefail
lualatex "$@"
# parse -output-directory and the .tex path the same way, then:
gs -q -dNOPAUSE -dBATCH -dNoOutputFonts -sDEVICE=pdfwrite \
  -o "$pdf.outlined" "$pdf"
mv "$pdf.outlined" "$pdf"
```

Set `TIKZ_BIN` to the wrapper's absolute path in `_environment.local`. The
diagram cache key is the picture source only. After editing the shared styles,
delete the cache directory and render again.

## LaTeX articles

Use the same node styles and the same `\plate` macro. The picture is a
`tikzpicture` in the article, compiled by LuaLaTeX (XeLaTeX acceptable), so it
stays a vector in the PDF. The Ghostscript outline step is only for the Quarto
SVG path.

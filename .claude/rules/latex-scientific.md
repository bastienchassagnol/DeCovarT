---
description: LuaLaTeX articles — preamble, maths, floats, TikZ figures (SVG first), curve2e, minted, audit
globs: "**/*.{tex,cls,sty,bib},**/article/**/*"
alwaysApply: false
paths:
  - "**/*.{tex,cls,sty,bib}"
  - "**/article/**"
---
# Modern LaTeX (scientific articles)

Semantic markup, modular sources, maintainable packages. References:
[latex-advice](https://github.com/dspinellis/latex-advice),
[TeXtidote](https://github.com/sylvainhalle/textidote). Prose tone:
`scientific-writing` (deslop skill).

## Engine and preamble layout

- Compile with **LuaLaTeX** (XeLaTeX acceptable). Do **not** use `inputenc` or
  `fontenc` with LuaLaTeX — use **`fontspec`** and UTF-8 source directly.
- Guard legacy pdfLaTeX accidentally:

```latex
\usepackage{iftex}
\ifPDFTeX
  \PackageError{project}{Use LuaLaTeX or XeLaTeX}{}
\fi
\usepackage{fontspec}
```

**Two-block preamble habit** for new PDFs:

1. **Packages block** — load with options once (`amsmath`, `geometry`, `xcolor`,
   `tikz`, `pgfplots`, `tcolorbox`, `siunitx`, `microtype`, `cleveref`,
   `enumitem`, `booktabs`, `graphicx`, **`fancyhdr`**, **`titlesec`**, etc.).
2. **Macros block** — `\definecolor{...}`, `\DeclareMathOperator`, custom maths
   shortcuts, `\newtcolorbox` callout styles.

Large projects:

```text
project/
├── main.tex
├── preamble/
├── chapters/
├── figures/
├── bibliography/
└── appendices/
```

`\include` chapters; `\input` snippets; `\includeonly` while debugging.

## Document structure and spacing

`\frontmatter` / `\mainmatter` / `\backmatter` for theses. Prefer `\medskip` or
`\setlength{\parskip}{1ex}` over decorative `\\`.

## Mathematics

- Default display: **`equation`** (numbered, `\label` + `\cref`), not bare `$$`.
- No `eqnarray`. Use `align`, `split`, `cases`, `mathtools`.
- Operators: `\DeclareMathOperator{\Var}{var}` etc., not `\text{var}`.
- Subscripts upright where semantic: `x_{\mathrm{max}}`.

## Text

`\textbf`, `\textit`, `\emph` — not `\bf`/`\it`. `csquotes` + `\enquote{}`.

## Figures and tables

- Float specifier **`[tbp]`**; avoid habitual `[H]`.
- Inside floats: **`\centering`**, not `center` environment.
- **`caption`** with `singlelinecheck=off` when captions contain lists.
- Tables: `booktabs`, `longtable`, `tabularx`.

## Headers and callouts

**`fancyhdr`** for informative running headers/footers (chapter/section titles).

**`tcolorbox`** with `[most]` (loads `enhanced` skins). Callout pattern:

```latex
\usepackage[most]{tcolorbox}
\newtcolorbox{notebox}[1][]{
  colback=mainblue!5!white, colframe=mainblue!80!black,
  fonttitle=\bfseries, coltitle=white, colbacktitle=mainblue!90!black,
  enhanced,
  attach boxed title to top left={yshift=-2mm, xshift=5mm},
  boxed title style={sharp corners},
  title={\textcolor{white}{Note}}, #1}
```

## Code listings

Prefer **`minted`** + Pygments over **`listings`** for R/Python — UTF-8 and
GitHub-like highlighting:

```latex
\usepackage{minted}
\setminted{breaklines,fontsize=\footnotesize}
% build: -shell-escape; style=github-dark or monokai as appropriate
\begin{minted}{r}
# code
\end{minted}
```

Use `listings` only when `minted` is unavailable; then shared `literate` table for
UTF-8 in comments.

## Graphics

Use only packages that compile in one LuaLaTeX run on Ubuntu TeX Live:
TikZ, pgfplots, and `curve2e`. Do not use PSTricks. Its pictures need
`latex` and `dvips`, and TikZ or `curve2e` already draw those pictures.

- **TikZ** (`arrows.meta`, `positioning`, `shapes.geometric`, `calc`, `fit`)
  for DAGs, pipelines, and box diagrams.
- **`pgfplots`** (`\pgfplotsset{compat=1.18}`) for function and data plots
  that must live in TeX.
- Data plots from R: SVG from `svglite`, converted to PDF, then
  `\includegraphics` (`data-visualisation`).
- **`curve2e`** (Beccari; extends `pict2e`; `texdoc curve2e-manual`) for
  stroked or filled outlines and Bézier curves in the same LuaLaTeX run.
  Load **`xcolor` before `curve2e`**.
  - Author from nodes and outgoing directions:
    `\Curve(x,y)<dx,dy>...(x,y)<dx,dy>`. `\Curve*` fills. A second
    direction at a node is a cusp. Quadratic authoring is `\Qurve`
    (same syntax). `\CurveBetween` / `\CbezierBetween` draw one cubic
    and can set the distances to the control points.
  - Trace a path whose control points you already have (SVG `C`, or an
    SVG `Q` lifted to a cubic) with pict2e `\moveto`, `\lineto`,
    `\curveto`, `\closepath`, then `\fillpath` or `\strokepath`.
    `\Curve` and `\Qurve` recompute controls from tangents; they are the
    wrong tool for a trace. PDF has no quadratic operator, so lift `Q`
    to a cubic: each control sits two-thirds of the way from an endpoint
    to the quadratic control.
  - The `picture` origin is the bottom left. Flip SVG *y* before you
    emit coordinates. `\fillpath` is non-zero winding (PDF `f`). An SVG
    `evenodd` path needs `f*` (redefine `\pIIe@fill@op` locally), or
    holes fill in.
  - Efficiency: pict2e appends each segment of one path with `\xdef`,
    so one path gets slower as the square of its own length. Splitting
    the drawing into the original paths avoids that. A trace of this
    kind (about 6,000 segments over about 60 paths) compiled with
    LuaLaTeX in under two seconds. Do not build one giant path of tens
    of thousands of segments; include the SVG or PDF instead. Hobby
    splines are TikZ-only; `curve2e` has no Hobby interface.
- Transparent page: class `standalone`, and no `\pagecolor`. Delete a
  full-canvas white rectangle. After any PDF-to-SVG conversion, delete a
  converter-added white page rectangle.

### Fonts (default)

STIX Two Math for equations, TeX Gyre Heros for diagram labels (including
`bayesnet` figures, see `graphical-models`). Both ship with TeX Live and are
installed as system fonts, so Inkscape can edit Heros text.

```latex
\usepackage{amsmath}               % before unicode-math
\usepackage{fontspec}
\usepackage{unicode-math}          % no amssymb or bm; use \symbf
\setmathfont{STIX Two Math}
\setmainfont{TeX Gyre Heros}[Ligatures=CommonOff] % figure sources only
```

In a figure source, set Heros as the main font: a TikZ `font=\footnotesize`
key replaces `\sffamily`, so a sans-serif style alone leaves Latin Modern in
some labels. `Ligatures=CommonOff` keeps "fi" and "fl" as two letters: ligature
glyphs carry no Unicode, so Inkscape drops them from the SVG text. In a
manuscript, keep the body font and use
`\setsansfont{TeX Gyre Heros}` with `font=\sffamily` in TikZ styles.

### TikZ figures: SVG first, PDF last

Scope: scientific figures drawn in TikZ (diagrams, DAGs, plates,
schematics). R plots follow `data-visualisation`; Quarto `{.tikz}` blocks
get the portable SVG only (`graphical-models`).

Inkscape is the vector target. Each figure gets an **editable** SVG (labels
as live Heros text, maths as exact outlines or TexText objects) and a
**portable** SVG (all glyphs as outlines). The PDF is rendered from the
portable SVG, never from the editable one.

1. Write the figure as its own `standalone` source, `figures/name.tex`, with
   the font block. The `.tex` file stays authoritative for layout.
2. Compile once: `lualatex -interaction=nonstopmode -halt-on-error name.tex`.
3. Rename the maths font in the PDF to a same-length name that is not
   installed, so Inkscape outlines maths and keeps labels as text:

   ```bash
   qpdf --qdf --object-streams=disable name.pdf name-qdf.pdf
   sed 's/STIXTwoMath-Regular/STIXTwoMathOutlined/g' name-qdf.pdf > name-src.pdf
   ```

4. Editable SVG:
   `inkscape name-src.pdf --pages=1 --pdf-font-strategy=draw-missing --export-plain-svg -o name-edit.svg`.
   Heros labels are `<text>`; maths glyphs are paths with the TeX shapes.
5. Equations you expect to rewrite: in Inkscape, replace them with
   **TexText** objects (Extensions → Text → TexText; LuaLaTeX; preamble file
   `figures/textext-preamble.tex` holding the font block, with
   `\setmainfont` swapped for `\setsansfont`). Put them on their own layer.
   TikZ output never becomes TexText by itself.
6. Portable SVG:
   `inkscape name-edit.svg --export-text-to-path --export-plain-svg -o name.svg`.
7. Final PDF from the portable SVG:
   `inkscape name.svg --export-type=pdf -o name-final.pdf`. Check with
   `pdffonts name-final.pdf` (no fonts listed), then
   `\includegraphics{figures/name-final.pdf}` in the manuscript.
8. Commit `name.tex`, `name-edit.svg` and `name.svg`. The PDFs are build
   products. Regenerating step 4 overwrites manual edits, so make lasting
   changes in the `.tex` file or on the TexText layer.

Why not `dvisvgm`: tested with dvisvgm 3.6 and Inkscape 1.4, Inkscape ignores
the embedded WOFF2 fonts, and unicode-math glyphs land on wrong code points.
Why not keep maths as text (`--pdf-font-strategy=keep`): glyphs are mapped
through ToUnicode and some are wrong (σ became ς). `pdftocairo -svg` and
`pdf2svg` outline everything, so nothing stays editable.

TexText, per user, no sudo (Inkscape runs extensions with `/usr/bin/python3`):

```bash
# from the unpacked TexText-Linux-<version>.tar.gz release
uv run --no-project --python /usr/bin/python3 \
  python setup.py --inkscape-executable /usr/bin/inkscape
# installs into ~/.config/inkscape/extensions/textext/
```

## Acronyms and glossaries

Minimise acronyms in prose. When a manuscript needs recurring abbreviations, use
the [`glossaries`](https://ctan.org/pkg/glossaries) package rather than a
hand-maintained table.

- Preamble: `\usepackage[acronym,toc]{glossaries}` then `\makeglossaries`.
- Terms: `\newglossaryentry{alr}{name=additive log-ratio, description={…}}`
  (lowercase singular `name`).
- Acronyms: `\newacronym{alr}{ALR}{additive log-ratio}`.
- Cite in prose with `\gls{}` / `\Gls{}` / `\glspl{}`; `\acrshort{}`,
  `\acrlong{}`, `\acrfull{}` for acronym forms.
- Print with `\printglossary[type=\acronymtype]` and `\printglossaries`;
  retitle with `\printglossary[title=…, toctitle=…]` when the journal requires.
- **LuaLaTeX builds:** run `makeglossaries` (or configure `latexmk` custom
  dependency) so acronym lists stay current—not a single-pass compile.
- Keep entries consistent with `inst/WORDLIST` (R packages) and notation in
  TikZ figures.

## Bibliography

`.bib` exported from Zotero with native citation keys
(`documentation-bibliography`), tidied with `bibtex-tidy`. `cleveref` for
`\cref{eq:...}` / `\Cref{...}`.

## Anti-patterns

`$$` in `.tex` sources, `eqnarray`, deprecated font switches, hard-coded visual
spacing, decorative `\\`, unlabelled maths, floats without `\centering`.

## Pre-submission audit

Before treating a manuscript as complete, review `.tex` and `preamble/*.tex`
against this rule and produce a short report:

1. Engine: LuaLaTeX or XeLaTeX; `fontspec`; no `inputenc`/`fontenc` unless
   pdfLaTeX is intentional (error).
2. Maths: no `eqnarray`; numbered `equation`/`align` with labels; flag bare
   `$$` in `.tex` sources (not Quarto exports).
3. Floats: `[tbp]`; `\centering` in floats; flag habitual `[H]`.
4. References: `cleveref` preferred over manual `Eq.~\ref`.
5. Code: `minted` + `-shell-escape`; if `listings`, UTF-8 `literate` or accent risk.
6. Graphics: R plots as SVG-derived PDF; TikZ figures SVG-first (editable +
   portable SVG, PDF from the portable one); STIX Two Math for maths, TeX
   Gyre Heros for labels; no PSTricks.
7. Packages: `siunitx`, `booktabs`, `tcolorbox[most]` callouts when used.
8. Bibliography: citekeys resolvable; no invented `\cite{}`.
9. Glossaries: `\gls` used consistently when `glossaries` is loaded.

Output template:

```markdown
## LaTeX modern audit — <main.tex>

### Pass
- …

### Should fix (methodology / portability)
- [file:line] issue → recommended change

### Optional polish
- …

### Build notes
- Engine: LuaLaTeX; shell-escape for minted; makeglossaries / latexmk …
```

Prioritise compilation, portability and journal markup. Optional: TeXtidote
or `chktex` if installed.

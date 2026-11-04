# Figure sources

When a figure is needed, try these in order. Stop at the first one that works.

## 1. Reuse a licensed image

Search the web for a figure of the exact concept. Use only public domain or
Creative Commons (CC BY or CC BY-SA). Prefer Wikimedia Commons, or an
open-access paper that states the licence on the article. Skip paywalled
journal figures and any image with no stated licence.

For each candidate, give:

- the page URL
- the direct image URL
- the licence
- the author
- a one-line attribution

Do not download anything until the user picks one.

## 2. Draw it in R or Python

If no licensed image exists, draw the figure in R or Python when the picture
is a plot, diagram, or teaching sketch the code can produce. Follow the
repository's figure scripts and seeds.

## 3. Generate it

If neither a licensed image nor an R or Python drawing will do, generate the
image with the session's image generator. Say that the picture is generated,
and do not present it as a published figure.

## Archive replaced figures

When a figure file (`.png`, `.svg`, `.jpg`, `.jpeg`) is replaced by a TikZ
drawing or by a figure script, keep the original in `archives/figures/` at
the repository root:

- Move it with `git mv`, or recover a file already deleted with
  `git show <commit>^:<path> > archives/figures/<name>`.
- Keep the original file name.
- Add one line to `archives/figures/README.md`: original path, the figure and
  script that replace it, the commit, and the reason.

Do not archive outputs that are regenerated from a source kept in the
repository: Mermaid or TikZ diagrams, `_diagram-cache/` SVGs, and the outputs
of figure scripts. No vignette or chapter links to `archives/`, and the folder
is in `.Rbuildignore`.

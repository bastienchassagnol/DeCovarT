# Links to the paper

- [Overleaf repository, shared with Anais Baudot](https://www.overleaf.com/project/697398ed149e81845cb5208e)
- [Hal Archive](https://arxiv.org/pdf/2309.09557)

# Compiling the LaTeX article

Source file: `article/main.tex` (bibliography: `article/decovart_library.bib`).

From the repository root, generate `article/main.pdf` with:

```sh
cd article
latexmk -pdf -interaction=nonstopmode -file-line-error -synctex=1 main.tex
```

In Cursor / VS Code with **LaTeX Workshop**, use the default recipe **latexmk** (configured in `.vscode/settings.json`), or run **“LaTeX Workshop: Build with recipe”**.

In **LaTeX Workshop**, auxiliaries are removed automatically after a **successful** build (`latex-workshop.latex.autoClean.run`: `onSucceeded`). The PDF is kept. Manual clean: Command Palette → **“LaTeX Workshop: Clean up auxiliary files”**.
# Github actions

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/bastienchassagnol/DeCovarT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bastienchassagnol/DeCovarT/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->
  

<!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/bastienchassagnol/DeCovarT/graph/badge.svg)](https://app.codecov.io/gh/bastienchassagnol/DeCovarT)
  <!-- badges: end -->

# Project structure



```html
<iframe
  src="output/package_network/decovart_function_network.html"
  title="DeCovarT function call graph"
  width="100%"
  height="720"
  style="border: none;"
></iframe>
```

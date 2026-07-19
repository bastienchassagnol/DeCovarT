

<!-- README.md is generated from README.qmd. Please edit that file. -->

# DeCovarT <a href="https://bastienchassagnol.github.io/DeCovarT/"><img src="man/figures/logo.svg" alt="DeCovarT logo" align="right" height="139"/></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/bastienchassagnol/DeCovarT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bastienchassagnol/DeCovarT/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/bastienchassagnol/DeCovarT/graph/badge.svg)](https://app.codecov.io/gh/bastienchassagnol/DeCovarT)
<!-- badges: end -->

DeCovarT is an R package for bulk transcriptomic deconvolution that accounts for
gene–gene covariance structure in purified reference populations.

## Links to the paper

- [Overleaf repository (shared with Anaïs Baudot)](https://www.overleaf.com/project/697398ed149e81845cb5208e)
- [HAL / arXiv archive](https://arxiv.org/pdf/2309.09557)

## The generative model of DeCovarT

Constrained DeCovarT maps unconstrained coordinates
$\boldsymbol{\rho}_{i}\in\mathbb{R}^{J-1}$ to cellular proportions
$\boldsymbol{p}_{i}$ on the open simplex via a regularised soft-max
$\boldsymbol{\psi}$, then forms the bulk observation as a Gaussian convolution
of cell-type profiles
$\boldsymbol{x}_{j}\sim\mathcal{N}_{G}(\boldsymbol{\mu}_{j},\boldsymbol{\Sigma}_{j})$:

$$
\boldsymbol{y}_{i}\,|\,\cdot
\sim\mathcal{N}_{G}\!\Bigl(
  \boldsymbol{\mu}\boldsymbol{p}_{i},\;
  \sum_{j=1}^{J}p_{ji}^{2}\boldsymbol{\Sigma}_{j}
\Bigr).
$$

```mermaid
%%{init: {
  "theme": "base",
  "themeVariables": {
    "primaryColor": "#f5f7fa",
    "primaryBorderColor": "#334155",
    "primaryTextColor": "#0f172a",
    "lineColor": "#475569",
    "fontSize": "14px"
  },
  "flowchart": { "curve": "basis" }
}}%%
flowchart TB
  subgraph plateJ["j = 1,...,J"]
    direction TB
    mu["μ_j"]
    Sj["Σ_j"]
    xj(("x_j"))
    mu --> xj
    Sj --> xj
  end

  subgraph plateI["i = 1,...,N"]
    direction TB
    rho(("ρ_i"))
    p["p_i = ψ(ρ_i)"]
    y(["y_i"])
    rho -->|"ψ soft-max"| p
    p -.-> y
  end

  xj -.->|"y_i = μ p_i"| y

  classDef param fill:#ffffff,stroke:#0f172a,stroke-width:1.5px;
  classDef latent fill:#ffffff,stroke:#0f172a,stroke-width:1.5px;
  classDef det fill:#ffffff,stroke:#0f172a,stroke-width:1.5px,stroke-dasharray:5 3;
  classDef obs fill:#cbd5e1,stroke:#0f172a,stroke-width:1.5px,stroke-dasharray:5 3;

  class mu,Sj param;
  class rho,xj latent;
  class p det;
  class y obs;
```

Forward map
$\boldsymbol{\psi}:\boldsymbol{\rho}\mapsto\boldsymbol{p}$
($C^{2}$-diffeomorphism):

$$
p_{j}=\frac{e^{\rho_{j}}}{\sum_{k<J}e^{\rho_{k}}+1}\ (j<J),\qquad
p_{J}=\frac{1}{\sum_{k<J}e^{\rho_{k}}+1},\qquad
\boldsymbol{\psi}^{-1}(\boldsymbol{p})=\bigl(\ln(p_{j}/p_{J})\bigr)_{j<J}.
$$

## Compiling the LaTeX article

Source file: `article/main.tex` (bibliography: `article/decovart_library.bib`).
DAG panels in Fig.~DAG-model are drawn with [`tikz-bayesnet`](https://ctan.org/pkg/tikz-bayesnet).

From the repository root, generate `article/main.pdf` with:

```sh
cd article
latexmk -pdf -interaction=nonstopmode -file-line-error -synctex=1 main.tex
```

In Cursor / VS Code with **LaTeX Workshop**, use the default recipe **latexmk**
(configured in `.vscode/settings.json`), or run **“LaTeX Workshop: Build with recipe”**.

Auxiliaries are removed automatically after a successful build
(`latex-workshop.latex.autoClean.run`: `onSucceeded`). Manual clean: Command Palette
→ **“LaTeX Workshop: Clean up auxiliary files”**.

## Regenerating this README

Edit `README.qmd` (and optionally `man/figures/decovart_constrained_dag.mmd`), then run:

```sh
quarto render README.qmd
./scripts/fix_readme_fences.sh README.md
```

This overwrites `README.md` (GitHub Flavoured Markdown).
The fence-fix script removes the space Pandoc inserts after code-fence
language tags (needed for GitHub Mermaid rendering).

## Project structure

Package function call graph (DeCovarT → DeCovarT edges only):

```html
<iframe
  src="output/package_network/decovart_function_network.html"
  title="DeCovarT function call graph"
  width="100%"
  height="720"
  style="border: none;"
></iframe>
```

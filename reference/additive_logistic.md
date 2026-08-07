# Additive logistic transform (unconstrained coordinates to the simplex)

Implements the reparametrisation
\\\boldsymbol{\psi}:\boldsymbol{\rho}\mapsto\boldsymbol{p}\\ used in the
article, sending unconstrained coordinates
\\\boldsymbol{\rho}\in\mathbb{R}^{J-1}\\ to cellular proportions
\\\boldsymbol{p}\in\Delta^{J-1}\\. This is the *additive logistic
transform* of Aitchison, i.e. the inverse additive log-ratio map
(\\\mathrm{alr}^{-1}\\), equivalently a softmax with the last category
\\J\\ pinned as reference (\\\rho_J\equiv 0\\).

Recovers the unconstrained additive log-ratio coordinates
\\\rho_j=\ln(p_j/p_J)\\ for \\j=1,\ldots,J-1\\, with the last part
\\p_J\\ as reference. This is Aitchison's additive log-ratio
(\\\mathrm{alr}\\) transform, equivalently the multinomial-logit link
with reference category \\J\\ (see
[`compositions::alr()`](https://rdrr.io/pkg/compositions/man/alr.html)
and
[`vignette("softmax-alr-derivatives", package = "DeCovarT")`](https://bastienchassagnol.github.io/DeCovarT/articles/softmax-alr-derivatives.md)).

## Usage

``` r
additive_logistic(rho)

additive_log_ratio(p)
```

## Arguments

- rho:

  Numeric vector \\\boldsymbol{\rho}\in\mathbb{R}^{J-1}\\ of
  unconstrained additive log-ratio coordinates (reference cell type
  \\J\\).

- p:

  Numeric vector \\\boldsymbol{p}\\ on the open simplex.

## Value

Numeric vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\ on the unit simplex
(\\\mathbf{1}^{\mathsf{T}}\boldsymbol{p}=1\\,
\\\boldsymbol{p}\ge\mathbf{0}\\).

Numeric vector \\\boldsymbol{\rho}\in\mathbb{R}^{J-1}\\.

## Details

With \\A=\sum\_{k=1}^{J-1}\mathrm{e}^{\rho_k}+1\\, \$\$
p_j=\frac{\mathrm{e}^{\rho_j}}{A}\quad(j\<J),\qquad p_J=\frac{1}{A}.
\$\$ Equivalently,
\\\boldsymbol{p}=\boldsymbol{\psi}(\boldsymbol{\rho})\\ with
\\\boldsymbol{\psi}(\boldsymbol{\rho})\propto
(\mathrm{e}^{\rho_1},\ldots,\mathrm{e}^{\rho\_{J-1}},1)^{\mathsf{T}}\\.
Jacobians and Hessians of both maps are derived in the package vignette
[`vignette("softmax-alr-derivatives", package = "DeCovarT")`](https://bastienchassagnol.github.io/DeCovarT/articles/softmax-alr-derivatives.md).
See also
[`compositions::alrInv()`](https://rdrr.io/pkg/compositions/man/alr.html).

## See also

The inverse map (additive log-ratio) is documented as
`additive_log_ratio()` on this help page.

## Examples

``` r
rho <- c(0.2, -0.5)
p <- additive_logistic(rho)
sum(p)
#> [1] 1
additive_log_ratio(p)
#> [1]  0.2 -0.5
```

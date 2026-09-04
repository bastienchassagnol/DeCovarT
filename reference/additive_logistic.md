# Additive logistic transform (unconstrained coordinates to the simplex)

Inverse additive log-ratio map
\\\boldsymbol{\psi}:\boldsymbol{\rho}\mapsto\boldsymbol{p}\\ of
Aitchison (\\\mathrm{alr}^{-1}\\): a softmax with the last category
\\J\\ pinned as reference (\\\rho_J\equiv 0\\). Solvers and
[`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
use the isometric log-ratio chart
[`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md)
instead. This ALR helper is retained for the vignette appendix and for
reference-invariance checks against
[`vcov_alr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_alr_delta.md).

Recovers the unconstrained additive log-ratio coordinates
\\\rho_j=\ln(p_j/p_J)\\ for \\j=1,\ldots,J-1\\, with the last part
\\p_J\\ as reference. This is Aitchison's additive log-ratio
(\\\mathrm{alr}\\) transform, equivalently the multinomial-logit link
with reference category \\J\\ (see
[`compositions::alr()`](https://rdrr.io/pkg/compositions/man/alr.html)
and
[`vignette("theory-decovart-generative-model", package = "DeCovarT")`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.md)).

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
[`vignette("theory-decovart-generative-model", package = "DeCovarT")`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.md).
See also
[`compositions::alrInv()`](https://rdrr.io/pkg/compositions/man/alr.html).

## Numerical stability

Evaluated through the log-sum-exp shift
\\p_i\propto\exp(\tilde{\rho}\_i-\max_k\tilde{\rho}\_k)\\ with
\\\tilde{\boldsymbol{\rho}}=(\rho_1,\ldots,\rho\_{J-1},0)\\, which is
algebraically identical to the ratio above but never forms
\\\mathrm{e}^{\rho_i}\\ directly. The naive quotient overflows to `NaN`
for \\\rho_i\gtrsim 710\\, which unconstrained ascent can reach when the
MLE approaches a simplex face; the shifted form returns an exact zero
instead and keeps \\\boldsymbol{\Sigma}(\boldsymbol{p})\\ positive
definite.

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
# Stable far out in ALR space, where exp(rho) would overflow.
additive_logistic(c(800, 1200))
#> [1] 1.91517e-174  1.00000e+00  0.00000e+00
```

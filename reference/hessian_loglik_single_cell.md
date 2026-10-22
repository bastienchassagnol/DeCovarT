# Observed Hessian of the single-cell log-likelihood in \\\boldsymbol{p}\\

Second derivatives of
[`loglik_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_single_cell.md)
with respect to \\\boldsymbol{p}\\. Newton-Raphson uses this observed
Hessian. The expected information in the article, \\I\_{jk}=c^2\phi^{-1}
\tilde{\boldsymbol{\mu}}\_{\cdot j}^{\top}\boldsymbol{W}
\tilde{\boldsymbol{\mu}}\_{\cdot k} +\tfrac{1}{2}\mathrm{tr}(
\boldsymbol{W}\boldsymbol{\Sigma}\_j\boldsymbol{W}\boldsymbol{\Sigma}\_k)\\,
is \\\mathbb{E}\[-\mathbf{H}\]\\ and is not what
[`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html)
computes, so it is not the matrix passed to the solver.

## Usage

``` r
hessian_loglik_single_cell(
  p,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale
)
```

## Arguments

- p:

  Composition of length \\J\\. It does not have to sum to 1 (numerical
  derivatives step off the simplex).

- y:

  Bulk vector of length \\G\\.

- mean_signature_matrix:

  Reference mean, \\G\times J\\.

- Sigma:

  Cell-level covariances, \\G\times G\times J\\.

- mean_scale:

  Positive plug-in \\c=an\\.

- variance_scale:

  Positive plug-in \\\phi=a^2 n\\.

## Value

Symmetric numeric \\J\times J\\ matrix.

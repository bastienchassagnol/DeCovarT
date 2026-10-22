# Cached Cholesky factorisation of \\\boldsymbol{\Sigma}(\boldsymbol{p})\\

Assembles
[`.compute_global_variance()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-compute_global_variance.md)
and factorises it once via
[`base::chol()`](https://rdrr.io/r/base/chol.html), returning the
matrix, its upper Cholesky factor, the log-determinant (from the
factor's diagonal) and the precision
\\\boldsymbol{\Omega}(\boldsymbol{p})\\ (from the same factor, via
[`base::chol2inv()`](https://rdrr.io/r/base/chol2inv.html)). The
covariance is treated as dense: no sparsity or low-rank structure is
assumed.

## Usage

``` r
.sigma_p_factorisation(p, Sigma)
```

## Arguments

- p:

  Numeric vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\.

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

## Value

A list with `matrix` (\\\boldsymbol{\Sigma}(\boldsymbol{p})\\), `chol`
(upper-triangular Cholesky factor), `log_det`
(\\\log\det\boldsymbol{\Sigma}(\boldsymbol{p})\\) and `inverse`
(\\\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\\).

## Details

Newton-type solvers evaluate the log-likelihood, score and Hessian at
the same trial \\\boldsymbol{p}\\ within one iteration. The shared
single-slot cache (keyed on \\p_j^2\\ and `Sigma`) factorises
\\\boldsymbol{\Sigma}(\boldsymbol{p})\\ once (\\O(G^{3})\\) and returns
the cached result to the later calls.

## Examples

``` r
p <- c(0.6, 0.4)
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
names(.sigma_p_factorisation(p, Sigma))
#> [1] "matrix"  "chol"    "log_det" "inverse"
```

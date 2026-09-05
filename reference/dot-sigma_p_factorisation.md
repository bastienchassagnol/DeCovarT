# Cached Cholesky factorisation of \\\boldsymbol{\Sigma}(\boldsymbol{p})\\

Assembles
[`.compute_global_variance()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-compute_global_variance.md)
and factorises it once via
[`base::chol()`](https://rdrr.io/r/base/chol.html), returning the matrix
itself together with its Cholesky factor, log-determinant (via the
factor's diagonal, not
[`base::det()`](https://rdrr.io/r/base/det.html)), and inverse (via
[`qr.solve()`](https://rdrr.io/r/base/qr.html) after a uniform spectral
shift if the mixture is not numerically SPD).

## Usage

``` r
.sigma_p_factorisation(p, Sigma, backend = NULL)
```

## Arguments

- p:

  Numeric vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\.

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

- backend:

  Optional
  [`new_decovart_covariance()`](https://bastienchassagnol.github.io/DeCovarT/reference/new_decovart_covariance.md)
  object declaring an exploitable covariance structure (block, band,
  sparse or diagonal-plus-low-rank). When `NULL` (default) the universal
  dense Cholesky is used and cached, exactly as before. When supplied,
  the structured operators of that backend supply `log_det` and the
  `solve` / `quadform` closures; `inverse` is then materialised only for
  drop-in compatibility with callers that still expect an explicit
  precision, and the operator path (`solve`) should be preferred to keep
  the structural speed-up (a dense `chol2inv` would undo it).

## Value

A list with elements: `matrix` (\\\boldsymbol{\Sigma}(\boldsymbol{p})\\
itself, or `NULL` for factored backends), `chol` (upper-triangular
Cholesky factor, or `NULL`), `log_det`
(\\\log\det\boldsymbol{\Sigma}(\boldsymbol{p})\\) and `inverse`
(\\\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\\). When `backend` is
supplied the list additionally carries `solve` and `quadform` closures.

## Details

[`stats::optim()`](https://rdrr.io/r/stats/optim.html),
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) and
[`marqLevAlg::marqLevAlg()`](https://rdrr.io/pkg/marqLevAlg/man/marqLevAlg.html)
each treat the log-likelihood, gradient and Hessian as three independent
callback functions, but all Newton-type solvers evaluate them at the
SAME trial point \\\boldsymbol{p}\\ within one iteration (the Hessian
chain rule in
[`hessian_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_constrained.md)
even re-evaluates the unconstrained gradient a second time internally).
Without sharing work, one (log-lik, gradient, Hessian) triple pays for
the \\O(G^{3})\\ assembly-and-factorisation of
\\\boldsymbol{\Sigma}(\boldsymbol{p})\\ up to four times over; this
single-slot cache, keyed on exact equality of `p` and `Sigma`, means
only the first of those calls actually factorises, and the rest simply
return the cached result in \\O(G^{2})\\ (the cost of the equality
check). Profiling on a 38-gene / 3-cell-type scenario showed the
redundancy

## See also

[`new_decovart_covariance()`](https://bastienchassagnol.github.io/DeCovarT/reference/new_decovart_covariance.md)

## Examples

``` r
p <- c(0.6, 0.4)
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
names(.sigma_p_factorisation(p, Sigma))
#> [1] "matrix"  "chol"    "log_det" "inverse"
```

# Cached Cholesky factorisation of \\\boldsymbol{\Sigma}(\boldsymbol{p})\\

Assembles
[`.compute_global_variance()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-compute_global_variance.md)
and factorises it once via
[`base::chol()`](https://rdrr.io/r/base/chol.html), returning the matrix
itself together with its Cholesky factor, log-determinant (via the
factor's diagonal, not
[`base::det()`](https://rdrr.io/r/base/det.html)), and inverse (via
[`base::chol2inv()`](https://rdrr.io/r/base/chol2inv.html), which reuses
the factor rather than an independent
[`base::solve()`](https://rdrr.io/r/base/solve.html)).

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

A list with elements: `matrix` (\\\boldsymbol{\Sigma}(\boldsymbol{p})\\
itself), `chol` (upper-triangular Cholesky factor), `log_det`
(\\\log\det\boldsymbol{\Sigma}(\boldsymbol{p})\\) and `inverse`
(\\\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\\).

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

# Restricted maximum likelihood with fixed cellular ratios

Maximises the DeCovarT log-likelihood over the simplex while holding a
subset of coordinates at prescribed values. Writing \\A\\ for the
constrained index set and \\s=\sum\_{j\in A}c_j\\, the free block is
reparametrised as \$\$ \boldsymbol{p}\_{A^{c}} =
(1-s)\\\boldsymbol{\psi}(\boldsymbol{\rho}), \qquad
\boldsymbol{\rho}\in\mathbb{R}^{\|A^{c}\|-1}, \$\$ with
\\\boldsymbol{\psi}\\ the additive logistic map
([`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)).
The constrained coordinates are *substituted* rather than pushed through
a logarithm, so a null such as \\p_j=0\\ is representable exactly: this
is what makes boundary likelihood-ratio tests computable (see
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md)).

## Usage

``` r
restricted_mle_decovart(
  y,
  mean_signature_matrix,
  Sigma,
  fixed,
  epsilon = 10^-8,
  itmax = 500L
)
```

## Arguments

- y:

  Numeric vector (or one-column matrix)
  \\\boldsymbol{y}\in\mathbb{R}^{G}\\.

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (plug-in means).

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

- fixed:

  Named or integer-indexed numeric vector of constrained ratios, e.g.
  `c(celltype_1 = 0)`. Names are matched against
  `colnames(mean_signature_matrix)`.

- epsilon, itmax:

  Relative convergence tolerance and iteration budget passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

A list with `coefficients` (the restricted \\J\\ vector), `loglik`,
`fixed`, and `convergence`.

## Details

Zero mass on \\A\\ keeps \\\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
p_j^{2} \boldsymbol{\Sigma}\_j\\ positive definite as long as at least
one coordinate remains strictly positive. Optimisation uses
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) (`"BFGS"`) with
the analytic score obtained by chaining
[`gradient_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_unconstrained.md)
through \\(1-s)\mathbf{J}\_{\boldsymbol{\psi}}\\.

## See also

[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`profile_loglik_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_decovart.md),
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md),
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md),
[`profile_loglik_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_decovart.md),
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md)

## Examples

``` r
mu <- matrix(c(20, 22, 18, 22, 20, 24), nrow = 2)
colnames(mu) <- paste0("ct", 1:3)
Sigma <- array(c(diag(2), diag(2), diag(2)), dim = c(2, 2, 3))
y <- drop(mu %*% c(0.5, 0.5, 0))
restricted_mle_decovart(y, mu, Sigma, fixed = c(ct3 = 0))$coefficients
#> ct1 ct2 ct3 
#> 0.5 0.5 0.0 
```

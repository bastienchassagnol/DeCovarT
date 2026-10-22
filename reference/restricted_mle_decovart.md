# Maximum likelihood on the closed simplex, optionally on a face

Maximises the DeCovarT log-likelihood, summed over the replicate bulk
columns of `y`, over the closed simplex, while holding a (possibly
empty) subset of coordinates at prescribed values. It is the fitting
engine of
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md)
and
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md).
Writing \\A\\ for the constrained index set and \\s=\sum\_{j\in A}c_j\\,
the free block is reparametrised as \$\$ \boldsymbol{p}\_{A^{c}} =
(1-s)\\\boldsymbol{\psi}(\boldsymbol{z}), \qquad
\boldsymbol{z}\in\mathbb{R}^{\|A^{c}\|-1}, \$\$ with
\\\boldsymbol{\psi}\\ the isometric logistic map
([`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md))
on a Helmert basis of the free face. The constrained coordinates are
*substituted* rather than pushed through a logarithm, so a null such as
\\p_j=0\\ is representable exactly: this is what makes boundary
likelihood-ratio tests computable (see
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md)).

## Usage

``` r
restricted_mle_decovart(
  y,
  mean_signature_matrix,
  Sigma,
  fixed = NULL,
  epsilon = 10^-8,
  itmax = 500L
)
```

## Arguments

- y:

  Numeric vector of length \\G\\, or \\G\times N\\ matrix of replicate
  bulk profiles that share one composition; their log-likelihoods are
  summed (the pooled, population-level fit).

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (plug-in means).

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

- fixed:

  Numeric vector of constrained ratios named by cell types, e.g.
  `c(celltype_1 = 0)`; the names must be columns of
  `mean_signature_matrix`. `NULL` (default) fits the unrestricted MLE.

- epsilon, itmax:

  Relative convergence tolerance and iteration budget passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

## Value

A list with `coefficients` (the fitted \\J\\ vector), `loglik` (summed
over the columns of `y`), `fixed`, and `convergence`.

## Details

Zero mass on \\A\\ keeps \\\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
p_j^{2} \boldsymbol{\Sigma}\_j\\ positive definite as long as at least
one coordinate remains strictly positive. The free face is fitted by
Newton-Raphson with the PORT trust region of
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), from its
barycentre, with the pooled analytic score and Hessian chained through
\\(1-s)\mathbf{J}\_{\boldsymbol{\psi}}\\. The bounded first step
matters: with several columns the pooled score is large, and an
unbounded quasi-Newton step can land on a vertex, where the likelihood
is flat, and stop there.

Without `fixed`, the interior ILR optimum is compared with the
restricted fit on the face where it stops (coordinates below `1e-3` set
to zero), and the larger likelihood is kept: the ILR chart cannot
represent an exact zero, so this is how a boundary MLE is reached.

## See also

[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_single_cell.md),
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_decovart.md),
[`confint_ilr_backtransform()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_ilr_backtransform.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`confint_profile_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_single_cell.md),
[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md),
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md),
[`fit_firth_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_firth_decovart.md),
[`fit_firth_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_firth_single_cell.md),
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`lrt_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_single_cell.md),
[`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md),
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md),
[`reference_bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_single_cell.md),
[`restricted_mle_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_single_cell.md)

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

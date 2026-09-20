# Log-likelihood profile in the ILR coordinate \\\rho\in\mathbb{R}^{J-1}\\

For the bivariate toy (\\J=2\\) the free coordinate is scalar. The
profile evaluates
[`loglik_multivariate_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate_constrained.md)
on a grid of \\\rho\\ (Helmert ILR). A white diamond marks
[`isometric_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md)\\(p^{\star})\\;
a red circle with a black stroke marks the numerical MLE (grid argmax of
\\\ell\\).

## Usage

``` r
plot_bulk_loglik_ilr_profile(
  true_theta,
  grid = 400L,
  y = NULL,
  rho_lim = NULL,
  expected = FALSE
)
```

## Arguments

- true_theta:

  List with `p`, `mu`, `sigma` for \\G=2\\.

- grid:

  Number of \\\rho\\ evaluation points.

- y:

  Optional bulk observation. Default is \\\mu p^{\star}\\. Ignored when
  `expected = TRUE`.

- rho_lim:

  Length-2 range for \\\rho\\. Default spans -4 to 10, expanded if
  needed to include the ILR image of \\p^{\star}\\ (unbalanced
  compositions sit at large positive \\\rho\\).

- expected:

  If `TRUE`, evaluate the expected log-likelihood under
  \\Y\sim\mathcal{N}(\mu p^{\star},\Sigma(p^{\star}))\\.

## Value

A `ggplot` (likelihood versus \\\rho\\ on a log10 y-axis).

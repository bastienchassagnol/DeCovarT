# Log-likelihood profile in the ILR coordinate \\\rho\in\mathbb{R}^{J-1}\\

For the bivariate toy (\\J=2\\) the free coordinate is scalar. The
profile evaluates
[`loglik_multivariate_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate_constrained.md)
on a grid of \\\rho\\ (Helmert ILR). The MLE for \\y=\mu p^{\star}\\ is
marked at
[`isometric_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md)\\(p^{\star})\\.

## Usage

``` r
plot_bulk_loglik_ilr_profile(true_theta, grid = 400L, y = NULL, rho_lim = NULL)
```

## Arguments

- true_theta:

  List with `p`, `mu`, `sigma` for \\G=2\\.

- grid:

  Number of \\\rho\\ evaluation points.

- y:

  Optional bulk observation. Default is \\\mu p^{\star}\\.

- rho_lim:

  Length-2 range for \\\rho\\. Default spans -4 to 10, expanded if
  needed to include the ILR image of \\p^{\star}\\ (unbalanced
  compositions sit at large positive \\\rho\\).

## Value

A `ggplot` (likelihood versus \\\rho\\ on a log10 y-axis).

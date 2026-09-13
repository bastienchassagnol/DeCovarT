# Contour of the bulk log-likelihood on hypothesised ratios

For \\J=2\\ the axes are \\(p_1,p_2)\\ on the unit square (the dashed
line is the simplex). For \\J=3\\ the axes are additive log-ratio
coordinates \\\rho_1=\ln(p_1/p_3)\\, \\\rho_2=\ln(p_2/p_3)\\
([`additive_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)).
The true simulation proportions (MLE for \\y=\mu p^{\star}\\) are
marked.

## Usage

``` r
plot_bulk_loglik_surface_p(true_theta, grid = 50L, y = NULL)
```

## Arguments

- true_theta:

  List with `p`, `mu`, `sigma` for \\G=2\\.

- grid:

  Length of the lattice per axis.

- y:

  Optional bulk observation. Default is \\\mu p^{\star}\\.

## Value

A `ggplot`.

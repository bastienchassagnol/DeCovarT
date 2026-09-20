# Contour of the bulk log-likelihood on hypothesised ratios

For \\J=2\\ the axes are \\(p_1,p_2)\\ on the unit square (the dashed
line is the simplex). For \\J=3\\ the axes are additive log-ratio
coordinates \\\rho_1=\ln(p_1/p_3)\\, \\\rho_2=\ln(p_2/p_3)\\
([`additive_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)).
White diamonds mark the true simulation proportions; red circles with a
black stroke mark the numerical MLE (grid argmax of \\\ell\\). For a
single observation \\y=\mu p^{\star}\\ these two typically diverge
unless \\p^{\star}\\ is equi-balanced. Set `expected = TRUE` to plot
\\\mathbb{E}\_{Y\mid p^{\star}}\[\ell(p;Y)\]\\ instead.

## Usage

``` r
plot_bulk_loglik_surface_p(true_theta, grid = 50L, y = NULL, expected = FALSE)
```

## Arguments

- true_theta:

  List with `p`, `mu`, `sigma` for \\G=2\\.

- grid:

  Length of the lattice per axis.

- y:

  Optional bulk observation. Default is \\\mu p^{\star}\\. Ignored when
  `expected = TRUE`.

- expected:

  If `TRUE`, evaluate the expected log-likelihood under
  \\Y\sim\mathcal{N}(\mu p^{\star},\Sigma(p^{\star}))\\.

## Value

A `ggplot`.

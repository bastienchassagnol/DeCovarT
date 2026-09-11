# Contour of the bulk log-likelihood on a \\(p_1,p_2)\\ lattice

Axes are hypothesised cell-type ratios, not gene expression. The true
simulation proportions (MLE for \\y=\mu p^{\star}\\) are marked; the
dashed line is the simplex \\p_1+p_2=1\\.

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

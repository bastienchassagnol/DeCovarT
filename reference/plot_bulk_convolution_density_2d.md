# 2-D density of the bulk convolution \\y\sim N(\mu p,\Sigma(p))\\

Centroids of the purified components are marked; 95% ellipses are
omitted so the convolution density fills the panel.

## Usage

``` r
plot_bulk_convolution_density_2d(true_theta, n = 1200L)
```

## Arguments

- true_theta:

  List with `p`, `mu`, `sigma` for \\G=2\\.

- n:

  Draws **per cell type**.

## Value

A `ggplot`.

# 2-D density of purified Gaussians

A kernel-density raster fills the shared gene-axis window so panels have
no inner white frame. Cell-type 1 is a red circle; cell-type 2 is a
green triangle. Outlines are exact 95% Gaussian ellipses for the known
means and covariances
([`gaussian_confidence_ellipse()`](https://bastienchassagnol.github.io/DeCovarT/reference/gaussian_confidence_ellipse.md)):
\\(x-\mu)^{\mathsf{T}}\Sigma^{-1}(x-\mu)\le\chi^{2}\_{2,0.95}\\.

## Usage

``` r
plot_purified_density_2d(true_theta, n = 800L)
```

## Arguments

- true_theta:

  List with `p`, `mu`, `sigma` for \\G=2\\.

- n:

  Draws **per cell type**.

## Value

A `ggplot`.

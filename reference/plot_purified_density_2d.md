# 2-D density of purified Gaussians

An exact mixture raster \\\sum_j
p_j\\\varphi(x;\boldsymbol{\mu}\_j,\boldsymbol{\Sigma}\_j)\\ fills the
shared gene-axis window, so rare types contribute mass in proportion to
the scenario's composition (not an equal-sized kernel density per cell
type). This is **not** the bulk convolution
\\\mathcal{N}(\boldsymbol{\mu}\boldsymbol{p},\boldsymbol{\Sigma}(\boldsymbol{p}))\\.
Cell-type 1 is a red circle; cell-type 2 is a green triangle. Outlines
are exact 95% Gaussian ellipses for the known means and covariances
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

  Unused for the raster (kept for API compatibility with
  [`plot_bulk_convolution_density_2d()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_bulk_convolution_density_2d.md));
  mixing weights come from `true_theta$p`.

## Value

A `ggplot`.

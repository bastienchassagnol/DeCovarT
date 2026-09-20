# rgl surface of the bulk log-likelihood

For \\J=2\\ the lattice is \\(p_1,p_2)\\. For \\J\ge 3\\ it is the ALR
plane \\(\rho_1,\rho_2)\\.

## Usage

``` r
plot_bulk_loglik_rgl(
  true_theta,
  grid = 40L,
  y = NULL,
  title = "Bulk log-likelihood"
)
```

## Arguments

- true_theta:

  List with `p`, `mu`, `sigma` for \\G=2\\.

- grid:

  Length of the lattice per axis.

- y:

  Optional bulk observation. Default is \\\mu p^{\star}\\. Ignored when
  `expected = TRUE`.

- title:

  Plot title (bold).

## Value

A `ggplot` raster of the `rgl` snapshot.

## Details

Opens an `rgl` window, draws
[`rgl::persp3d()`](https://dmurdoch.github.io/rgl/dev/reference/persp3d.html)
of
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md),
and marks the true composition (pale sphere) and the numerical MLE (red
sphere with a dark halo). The snapshot is a ggplot raster; a **single**
ggplot2 legend is attached once per PDF page (including 2-by-2 books),
not per subplot.

## Examples

``` r
if (FALSE) { # interactive() && requireNamespace("rgl", quietly = TRUE) && requireNamespace("png", quietly = TRUE)
mu <- matrix(c(20, 22, 22, 20), 2)
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
th <- list(p = c(0.5, 0.5), mu = mu, sigma = Sigma)
plot_bulk_loglik_rgl(th, grid = 20L)
}
```

# rgl surface of the bulk log-likelihood on a \\(p_1,p_2)\\ lattice

Opens an `rgl` window, draws
[`rgl::persp3d()`](https://dmurdoch.github.io/rgl/dev/reference/persp3d.html)
of
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)
versus hypothesised ratios, and marks the MLE (true simulation
proportions for \\y=\mu p^{\star}\\) with a sphere. Returns a ggplot
snapshot suitable for a PDF page.

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

  Optional bulk observation. Default is \\\mu p^{\star}\\.

- title:

  Plot title (bold).

## Value

A `ggplot` raster of the `rgl` snapshot.

## Examples

``` r
if (FALSE) { # interactive() && requireNamespace("rgl", quietly = TRUE) && requireNamespace("png", quietly = TRUE)
mu <- matrix(c(20, 22, 22, 20), 2)
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
th <- list(p = c(0.5, 0.5), mu = mu, sigma = Sigma)
plot_bulk_loglik_rgl(th, grid = 20L)
}
```

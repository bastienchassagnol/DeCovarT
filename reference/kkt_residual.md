# Simplex KKT / projected-score residual

\\\\\Pi\_\Delta(p + \alpha g) - p\\\_2 / \alpha\\ for maximisation, with
\\\Pi\_\Delta\\ the Euclidean projection onto the probability simplex
(sort-and-threshold algorithm). It is zero exactly at a
Karush-Kuhn-Tucker point of the log-likelihood on the closed simplex,
interior or on a face, which makes it a solver-independent check of an
estimate.

## Usage

``` r
kkt_residual(p, grad, alpha = 1)
```

## Arguments

- p:

  Composition on the simplex.

- grad:

  Score (gradient of the log-likelihood) at `p`, same length.

- alpha:

  Positive step; a non-finite or non-positive value is replaced by 1.

## Value

A non-negative number, or `NA` when the lengths differ.

## Examples

``` r
# Interior stationary point: the score is parallel to the ones vector.
kkt_residual(c(0.3, 0.7), c(2, 2))
#> [1] 2.775558e-16
# Vertex with an outward score: still a KKT point.
kkt_residual(c(1, 0), c(1, -5))
#> [1] 0
kkt_residual(c(0.5, 0.5), c(1, -1))
#> [1] 0.7071068
```

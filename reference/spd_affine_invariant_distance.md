# Affine-invariant Riemannian (AIRM) distance between two SPD matrices

AIRM is the **affine-invariant Riemannian metric** on the cone of
symmetric positive-definite matrices:
\\d_R(A,B)=\lVert\log(A^{-1/2}BA^{-1/2})\rVert_F\\. Do **not** use the
Frobenius \\\lVert A-B\rVert_F\\: that Euclidean chord leaves the cone,
is not inversion-invariant, and suffers the swelling effect.

## Usage

``` r
spd_affine_invariant_distance(a, b)
```

## Arguments

- a, b:

  Symmetric positive-definite matrices of equal size.

## Value

Non-negative scalar.

## See also

[`compute_average_riemannian()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_riemannian.md)

## Examples

``` r
# Two bivariate Gaussians: identity vs a correlated SPD covariance.
a <- diag(2)
b <- matrix(c(2, 0.5, 0.5, 1), nrow = 2)
spd_affine_invariant_distance(a, a)
#> [1] 0
d <- spd_affine_invariant_distance(a, b)
d > 0
#> [1] TRUE
```

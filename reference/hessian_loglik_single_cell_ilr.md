# Observed Hessian of the single-cell log-likelihood in ILR coordinates

Second-order chain rule through
[`jacobian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_isometric_logistic.md)
and
[`hessian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_isometric_logistic.md),
as for
[`hessian_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_constrained.md).

## Usage

``` r
hessian_loglik_single_cell_ilr(
  z,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale,
  V = NULL
)
```

## Arguments

- z:

  ILR coordinate of length \\J-1\\.

- y:

  Bulk vector of length \\G\\.

- mean_signature_matrix:

  Reference mean, \\G\times J\\.

- Sigma:

  Cell-level covariances, \\G\times G\times J\\.

- mean_scale:

  Positive plug-in \\c=an\\.

- variance_scale:

  Positive plug-in \\\phi=a^2 n\\.

- V:

  Optional ILR basis. See
  [`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md).

## Value

Symmetric numeric \\(J-1)\times(J-1)\\ matrix.

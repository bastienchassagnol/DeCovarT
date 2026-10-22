# Single-cell log-likelihood in ILR coordinates

Single-cell log-likelihood in ILR coordinates

## Usage

``` r
loglik_single_cell_ilr(
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

Scalar log-likelihood.

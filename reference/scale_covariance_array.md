# Scale component covariances, keeping precision zeros

Replaces each \\\Sigma_j\\ by \\s\Sigma_j\\. Then
\\\Omega_j(s)=\Omega_j/s\\, so the **support** of the precision
(structural zeros from the graph) is unchanged.

## Usage

``` r
scale_covariance_array(sigma, scale)
```

## Arguments

- sigma:

  \\G\times G\times J\\ covariance array.

- scale:

  Positive scalar \\s\\.

## Value

Scaled array.

## See also

[`scale_covariances_to_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/scale_covariances_to_overlap.md)

## Examples

``` r
sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
scaled <- scale_covariance_array(sigma, 4)
scaled[, , 1]
#>      [,1] [,2]
#> [1,]    4    0
#> [2,]    0    4
```

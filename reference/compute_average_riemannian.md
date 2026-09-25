# Mean pairwise AIRM distance of component covariances

Mean pairwise AIRM distance of component covariances

## Usage

``` r
compute_average_riemannian(true_theta, J = NULL)
```

## Arguments

- true_theta:

  List validated by
  [`check_true_theta()`](https://bastienchassagnol.github.io/DeCovarT/reference/check_true_theta.md):
  `mu` (\\G\times J\\), `sigma` (\\G\times G\times J\\), and optionally
  `p` (length \\J\\ or \\J\times N\\). If `p` is missing it is set to
  \\(1/J,\ldots,1/J)\\.

- J:

  Number of cell types. Defaults to the third dimension of `sigma`.

## Value

Scalar average of \\d_R(\Sigma_j,\Sigma\_\ell)\\ over \\j\<\ell\\.

## See also

[`spd_affine_invariant_distance()`](https://bastienchassagnol.github.io/DeCovarT/reference/spd_affine_invariant_distance.md)

## Examples

``` r
a <- diag(2)
b <- matrix(c(2, 0.5, 0.5, 1), nrow = 2)
theta <- list(
  mu = cbind(c(0, 0), c(1, 0)),
  sigma = array(c(a, b), dim = c(2, 2, 2))
)
compute_average_riemannian(theta)
#> [1] 0.8249946
# J = 2 has a single pair, so the average equals the pairwise AIRM.
all.equal(
  compute_average_riemannian(theta),
  spd_affine_invariant_distance(a, b)
)
#> [1] TRUE
```

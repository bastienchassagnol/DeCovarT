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

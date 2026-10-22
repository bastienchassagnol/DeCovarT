# Mean pairwise AIRM distance of component covariances

Average over \\j\<\ell\\ of the **affine-invariant Riemannian metric**
on the cone of symmetric positive-definite matrices,
\\d_R(\Sigma_j,\Sigma\_\ell)=\lVert\log(\Sigma_j^{-1/2}\Sigma\_\ell
\Sigma_j^{-1/2})\rVert_F\\, computed from the eigenvalues of
\\R^{-\mathsf{T}}\Sigma\_\ell R^{-1}\\ with
\\\Sigma_j=R^{\mathsf{T}}R\\. The Frobenius chord
\\\lVert\Sigma_j-\Sigma\_\ell\rVert_F\\ is not used: it leaves the cone,
is not inversion-invariant, and suffers the swelling effect. With
\\J=2\\ the average is the pairwise distance.

## Usage

``` r
compute_average_riemannian(true_theta, J = NULL)
```

## Arguments

- true_theta:

  List with the generative parameters (see `.parse_true_theta()`): `mu`
  (\\G\times J\\), `sigma` (\\G\times G\times J\\), and optionally `p`
  (length \\J\\ or \\J\times N\\). If `p` is missing it is set to
  \\(1/J,\ldots,1/J)\\.

- J:

  Number of cell types. Defaults to the third dimension of `sigma`.

## Value

Scalar average of \\d_R(\Sigma_j,\Sigma\_\ell)\\ over \\j\<\ell\\.

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
```

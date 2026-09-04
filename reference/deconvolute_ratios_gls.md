# Generalised least squares with a fixed residual covariance

Fits \\\boldsymbol{y}\approx\boldsymbol{\mu}\boldsymbol{p}\\ by
[`MASS::lm.gls()`](https://rdrr.io/pkg/MASS/man/lm.gls.html) with a
**known** \\G\times G\\ residual covariance `W` that does not depend on
\\\boldsymbol{p}\\, then
[`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md).
Use
[`fixed_gls_covariance()`](https://bastienchassagnol.github.io/DeCovarT/reference/fixed_gls_covariance.md)
for the global-diagonal competitor \\\operatorname{diag}\\\sum_j \bar
p_j^2\boldsymbol{\Sigma}\_j\\\\. This is the mean-only GLS gold standard
for covariance-structure benchmarks, not the DeCovarT convolution
\\\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
p_j^2\boldsymbol{\Sigma}\_j\\.

## Usage

``` r
deconvolute_ratios_gls(y, mean_signature_matrix, W, inverse = TRUE)
```

## Arguments

- y:

  Numeric bulk vector of length \\G\\.

- mean_signature_matrix:

  Numeric matrix \\G\times J\\.

- W:

  Residual covariance (when `inverse = TRUE`) or precision (when
  `inverse = FALSE`), of size \\G\times G\\.

- inverse:

  Passed to
  [`MASS::lm.gls()`](https://rdrr.io/pkg/MASS/man/lm.gls.html); `TRUE`
  means `W` is a covariance matrix.

## Value

Named simplex vector of length \\J\\.

## See also

[`fixed_gls_covariance()`](https://bastienchassagnol.github.io/DeCovarT/reference/fixed_gls_covariance.md),
[`deconvolute_ratios_deconrnaseq()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)

## Examples

``` r
mu <- cbind(ct1 = c(20, 40), ct2 = c(40, 20))
y <- drop(mu %*% c(0.4, 0.6))
Sigma <- array(c(diag(2), 2 * diag(2)), dim = c(2, 2, 2))
deconvolute_ratios_gls(y, mu, fixed_gls_covariance(Sigma))
#> ct1 ct2 
#> 0.4 0.6 
```

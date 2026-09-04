# Fixed residual covariance for a GLS competitor

Evaluates the convolution covariance
\\\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
p_j^2\boldsymbol{\Sigma}\_j\\ at a **known** composition (default the
barycentre \\p_j=1/J\\) and optionally keeps only the diagonal. That
matrix does **not** depend on the unknown \\\boldsymbol{p}\\ at fit
time. Do not copy it into every slice of a DeCovarT tensor: that would
yield \\\\p\\\_2^2\boldsymbol{\Sigma}\_{\mathrm{GLS}}\\, which still
depends on \\p\\.

## Usage

``` r
fixed_gls_covariance(Sigma, p = NULL, diagonal = TRUE)
```

## Arguments

- Sigma:

  Array \\G\times G\times J\\ of cell-type covariances.

- p:

  Optional composition of length \\J\\; default \\(1/J,\ldots,1/J)\\.

- diagonal:

  If `TRUE` (default), return
  \\\operatorname{diag}\\\boldsymbol{\Sigma}(\boldsymbol{p})\\\\.

## Value

A \\G\times G\\ covariance matrix.

## See also

[`deconvolute_ratios_gls()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_gls.md)

## Examples

``` r
Sigma <- array(c(diag(2), 2 * diag(2)), dim = c(2, 2, 2))
fixed_gls_covariance(Sigma)
#>      [,1] [,2]
#> [1,] 0.75 0.00
#> [2,] 0.00 0.75
```

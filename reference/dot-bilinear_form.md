# Evaluate a matrix-induced bilinear form

Computes \\\boldsymbol{x}^{\mathsf{T}}\boldsymbol{A}\boldsymbol{y}\\.
When \\\boldsymbol{A}\\ is symmetric positive definite (as for a
non-degenerate Gaussian covariance or precision), this is the
\\\boldsymbol{A}\\-inner product. Prefer this name over "dot product",
which is reserved for \\\boldsymbol{x}^{\mathsf{T}}\boldsymbol{y}\\.

## Usage

``` r
.bilinear_form(x, A, y = x)
```

## Arguments

- x:

  Numeric vector.

- A:

  Numeric square matrix, compatible with `x` and `y`.

- y:

  Numeric vector of the same length as `x` (default `x`, which yields
  the quadratic form
  \\\boldsymbol{x}^{\mathsf{T}}\boldsymbol{A}\boldsymbol{x}\\).

## Value

Numeric scalar.

## Details

Implementation uses
[`base::crossprod()`](https://rdrr.io/r/base/crossprod.html) as
`drop(crossprod(x, A %*% y))`, which is the standard efficient route to
a bilinear / quadratic form in R (avoids an explicit transpose of
\\\boldsymbol{x}\\ and a temporary outer product).

## Examples

``` r
.bilinear_form(c(1, 2), diag(2), c(3, 4))
#> [1] 11
```

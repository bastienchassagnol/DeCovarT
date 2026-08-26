# Evaluate a matrix-induced inner product

Computes \\\langle\boldsymbol{x},\boldsymbol{y}\rangle\_{\boldsymbol{A}}
=\boldsymbol{x}^{\mathsf{T}}\boldsymbol{A}\boldsymbol{y}\\. When
\\\boldsymbol{A}\\ is symmetric positive definite (as for a
non-degenerate Gaussian covariance or precision), this is the
\\\boldsymbol{A}\\-inner product on \\\mathbb{R}^{p}\\
(<https://en.wikipedia.org/wiki/Inner_product_space#Basic_properties>).
Prefer this name over the generic "bilinear form" when
\\\boldsymbol{A}\\ is SPD; the Euclidean inner product is the special
case \\\boldsymbol{A}=\mathbf{I}\\
(\\\boldsymbol{x}^{\mathsf{T}}\boldsymbol{y}\\).

## Usage

``` r
.inner_product(x, A, y = x)
```

## Arguments

- x:

  Numeric vector.

- A:

  Numeric square matrix, compatible with `x` and `y` (SPD when
  interpreted as an inner-product metric).

- y:

  Numeric vector of the same length as `x` (default `x`, which yields
  the squared \\\boldsymbol{A}\\-norm
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
.inner_product(c(1, 2), diag(2), c(3, 4))
#> [1] 11
```

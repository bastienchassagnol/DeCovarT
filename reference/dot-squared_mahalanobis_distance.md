# Squared Mahalanobis distance

Computes \\(\boldsymbol{x}-\boldsymbol{m})^{\mathsf{T}}
\boldsymbol{\Sigma}^{-1}(\boldsymbol{x}-\boldsymbol{m})\\ for a
symmetric positive-definite covariance (or scatter) matrix
\\\boldsymbol{\Sigma}\\. This is the squared distance; take the square
root for the Mahalanobis distance itself. Solves
\\\boldsymbol{\Sigma}\boldsymbol{z}=\boldsymbol{\delta}\\ instead of
forming \\\boldsymbol{\Sigma}^{-1}\\ explicitly.

## Usage

``` r
.squared_mahalanobis_distance(x, center = numeric(length(x)), covariance)
```

## Arguments

- x:

  Numeric vector.

- center:

  Numeric centre \\\boldsymbol{m}\\ (default the origin).

- covariance:

  Symmetric positive-definite \\G\times G\\ matrix.

## Value

Non-negative numeric scalar.

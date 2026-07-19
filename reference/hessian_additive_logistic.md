# Second derivatives (Hessian tensor) of the additive logistic map

Tensor of mixed partials of
[`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md),
\\\partial^{2}p_i/(\partial\rho_k\partial\rho_j)\\, stored as an array
of size \\(J-1)\times(J-1)\times J\\.

## Usage

``` r
hessian_additive_logistic(rho)
```

## Arguments

- rho:

  Numeric vector \\\boldsymbol{\rho}\in\mathbb{R}^{J-1}\\.

## Value

Numeric array used in the constrained Hessian chain rule.

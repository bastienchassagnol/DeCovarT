# Jacobian \\\mathbf{J}\_{\boldsymbol{\psi}}\\ of the additive logistic map

Returns the Jacobian of
[`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md),
\\\mathbf{J}\_{\boldsymbol{\psi}}\in\mathcal{M}\_{J\times(J-1)}\\ with
entries \\(\mathbf{J}\_{\boldsymbol{\psi}})\_{i,j}=\partial
p_i/\partial\rho_j\\.

## Usage

``` r
jacobian_additive_logistic(rho)
```

## Arguments

- rho:

  Numeric vector \\\boldsymbol{\rho}\in\mathbb{R}^{J-1}\\.

## Value

Numeric matrix of size \\J\times(J-1)\\.

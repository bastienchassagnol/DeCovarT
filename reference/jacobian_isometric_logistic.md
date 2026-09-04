# Jacobian of the isometric logistic map

\$\$ \mathbf{J}\_{\boldsymbol{\psi}}(\boldsymbol{z}) =
\frac{\partial\boldsymbol{p}}{\partial\boldsymbol{z}^{\mathsf{T}}} =
\mathbf{S}(\boldsymbol{p})\mathbf{V}, \$\$ with
\\\mathbf{S}(\boldsymbol{p})=\operatorname{diag}(\boldsymbol{p})
-\boldsymbol{p}\boldsymbol{p}^{\mathsf{T}}\\.

## Usage

``` r
jacobian_isometric_logistic(z, V = NULL)
```

## Arguments

- z:

  Numeric vector \\\boldsymbol{z}\in\mathbb{R}^{J-1}\\.

- V:

  Optional ILR basis; see
  [`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md).

## Value

Numeric matrix \\J\times(J-1)\\.

## Examples

``` r
jacobian_isometric_logistic(c(0.2, -0.5))
#>            [,1]       [,2]
#> [1,]  0.1952773  0.1742416
#> [2,] -0.1704937  0.1313152
#> [3,] -0.0247836 -0.3055568
```

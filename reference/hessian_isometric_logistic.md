# Hessian tensor of the isometric logistic map

Slice \\i\\ is \$\$ \mathbf{H}\_{\psi_i} =
p_i\bigl(\boldsymbol{q}\_i\boldsymbol{q}\_i^{\mathsf{T}}
-\mathbf{C}\_V\bigr), \$\$ where
\\\bar{\boldsymbol{v}}=\mathbf{V}^{\mathsf{T}}\boldsymbol{p}\\,
\\\boldsymbol{q}\_i=\mathbf{V}\_{i\cdot}-\bar{\boldsymbol{v}}\\ and
\\\mathbf{C}\_V=\mathbf{V}^{\mathsf{T}}\mathbf{S}(\boldsymbol{p})
\mathbf{V}\\. Array shape \\(J-1)\times(J-1)\times J\\.

## Usage

``` r
hessian_isometric_logistic(z, V = NULL)
```

## Arguments

- z:

  Numeric vector \\\boldsymbol{z}\in\mathbb{R}^{J-1}\\.

- V:

  Optional ILR basis; see
  [`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md).

## Value

Numeric array used in the constrained Hessian chain rule.

## Examples

``` r
hessian_isometric_logistic(c(0.2, -0.5))
#> , , 1
#> 
#>            [,1]         [,2]
#> [1,] 0.05085575  0.105113732
#> [2,] 0.10511373 -0.009675559
#> 
#> , , 2
#> 
#>            [,1]         [,2]
#> [1,]  0.0713137 -0.106489954
#> [2,] -0.1064900 -0.007291872
#> 
#> , , 3
#> 
#>              [,1]        [,2]
#> [1,] -0.122169445 0.001376222
#> [2,]  0.001376222 0.016967431
#> 
```

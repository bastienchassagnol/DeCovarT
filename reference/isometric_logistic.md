# Isometric logistic transform (ILR coordinates to the simplex)

Inverse isometric log-ratio map
\\\boldsymbol{\psi}:\boldsymbol{z}\mapsto\boldsymbol{p}\\ with a Helmert
basis \\\mathbf{V}\\, \$\$ \boldsymbol{p} =
\operatorname{softmax}(\mathbf{V}\boldsymbol{z}) =
\mathcal{C}\bigl(\exp(\mathbf{V}\boldsymbol{z})\bigr). \$\$ No cell type
is pinned as a reference (unlike
[`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)).
Evaluated with a log-sum-exp shift so large \\\\\boldsymbol{z}\\\\
cannot overflow.

\$\$ \boldsymbol{z} = \mathbf{V}^{\mathsf{T}}\log\boldsymbol{p} =
\mathbf{V}^{\mathsf{T}}\operatorname{clr}(\boldsymbol{p}), \$\$ which is
well-defined because \\\mathbf{V}^{\mathsf{T}}\mathbf{1}=0\\. Requires a
strictly positive composition (open simplex).

## Usage

``` r
isometric_logistic(z, V = NULL)

isometric_log_ratio(p, V = NULL)
```

## Arguments

- z:

  Numeric vector \\\boldsymbol{z}\in\mathbb{R}^{J-1}\\.

- V:

  Optional ILR basis; see `isometric_logistic()`.

- p:

  Numeric vector on the open simplex.

## Value

Numeric vector \\\boldsymbol{p}\\ on the unit simplex.

Numeric vector \\\boldsymbol{z}\in\mathbb{R}^{J-1}\\.

## See also

`isometric_log_ratio()`,
[`helmert_basis()`](https://bastienchassagnol.github.io/DeCovarT/reference/helmert_basis.md),
[`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)

## Examples

``` r
z <- c(0.2, -0.5)
p <- isometric_logistic(z)
sum(p)
#> [1] 1
isometric_log_ratio(p)
#> [1]  0.2 -0.5
isometric_logistic(c(800, 1200))
#> [1] 1 0 0
```

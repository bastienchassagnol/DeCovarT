# Complete a weighted adjacency to an SPD precision

Applies a uniform spectral shift that preserves off-diagonal signs and
support: \$\$ \boldsymbol{\Omega} = \boldsymbol{W} +
\bigl(\lvert\lambda\_{\min}(\boldsymbol{W})\rvert + u\bigr) \mathbf{I}.
\$\$

## Usage

``` r
build_normalised_precision(weighted_adjacency, precision_shift)
```

## Arguments

- weighted_adjacency:

  Symmetric numeric matrix (zero diagonal).

- precision_shift:

  Positive diagonal cushion \\u\\.

## Value

Symmetric positive-definite matrix.

## Examples

``` r
A <- generate_random_network_skeleton(6L, "hub")
W <- assign_iid_signed_weights(A)
build_normalised_precision(W, precision_shift = 0.2)
#>            [,1]       [,2]       [,3]      [,4]      [,5]      [,6]
#> [1,]  0.6242641 -0.3000000 -0.3000000 0.0000000 0.0000000 0.0000000
#> [2,] -0.3000000  0.6242641  0.0000000 0.0000000 0.0000000 0.0000000
#> [3,] -0.3000000  0.0000000  0.6242641 0.0000000 0.0000000 0.0000000
#> [4,]  0.0000000  0.0000000  0.0000000 0.6242641 0.3000000 0.3000000
#> [5,]  0.0000000  0.0000000  0.0000000 0.3000000 0.6242641 0.0000000
#> [6,]  0.0000000  0.0000000  0.0000000 0.3000000 0.0000000 0.6242641
```

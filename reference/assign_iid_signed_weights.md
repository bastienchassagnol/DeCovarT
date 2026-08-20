# Assign i.i.d. signed weights on a binary undirected adjacency

For each undirected edge, draws an inhibitory versus activatory sign and
a common magnitude. By convention a **positive** precision off-diagonal
is an **inhibitory** partial correlation (\\\rho\_{jk\mid
-\\j,k\\}=-\Omega\_{jk}/ \sqrt{\Omega\_{jj}\Omega\_{kk}}\\). Thus
`prop_inhibitory` is the target fraction of edges with
\\\Omega\_{jk}\>0\\ (negative partial correlation). The complementary
fraction is activatory (\\\Omega\_{jk}\<0\\).

## Usage

``` r
assign_iid_signed_weights(
  adjacency_matrix,
  prop_inhibitory = 0.5,
  weight_magnitude = 0.3
)
```

## Arguments

- adjacency_matrix:

  Symmetric binary matrix.

- prop_inhibitory:

  Numeric in \\\[0,1\]\\; expected fraction of inhibitory precision
  edges.

- weight_magnitude:

  Positive magnitude \\v\\ for every signed off-diagonal.

## Value

Symmetric numeric matrix with the same support as `adjacency_matrix` and
zero diagonal.

## Examples

``` r
A <- generate_random_network_skeleton(6L, "hub")
assign_iid_signed_weights(A, prop_inhibitory = 0.5)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  0.0  0.3 -0.3  0.0  0.0  0.0
#> [2,]  0.3  0.0  0.0  0.0  0.0  0.0
#> [3,] -0.3  0.0  0.0  0.0  0.0  0.0
#> [4,]  0.0  0.0  0.0  0.0  0.3 -0.3
#> [5,]  0.0  0.0  0.0  0.3  0.0  0.0
#> [6,]  0.0  0.0  0.0 -0.3  0.0  0.0
```

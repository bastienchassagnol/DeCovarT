# Helmert contrast matrix for isometric log-ratio coordinates

Returns the \\J\times(J-1)\\ Helmert sub-matrix \\\mathbf{V}\\ with \$\$
\mathbf{V}^{\mathsf{T}}\mathbf{V}=\mathbf{I}\_{J-1},\qquad
\mathbf{V}^{\mathsf{T}}\mathbf{1}=\mathbf{0}. \$\$ Column \\k\\
contrasts the first \\k\\ parts against part \\k+1\\. This is the
standard ILR basis used with
[`isometric_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md)
/
[`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md);
any other valid ILR basis is \\\mathbf{V}Q\\ for an orthogonal \\Q\\,
which rotates coordinates but leaves simplex geometry, eigenvalues of
quadratic forms, and \\\mathrm{Var}(\hat{\boldsymbol{p}})\\ unchanged.

## Usage

``` r
helmert_basis(n_parts)
```

## Arguments

- n_parts:

  Integer \\J\ge 2\\ (number of simplex parts).

## Value

Numeric matrix \\J\times(J-1)\\.

## Examples

``` r
V <- helmert_basis(3L)
crossprod(V)
#>              [,1]         [,2]
#> [1,] 1.000000e+00 2.451427e-17
#> [2,] 2.451427e-17 1.000000e+00
drop(crossprod(V, rep(1, 3)))
#> [1] 0 0
```

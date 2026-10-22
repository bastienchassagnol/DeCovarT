# Stage-i scales for one single-cell bulk sample

Generalised method of moments for the two technical scales, at a
composition that is either supplied or taken from the free-scale
DeconRNASeq start \\(c^{(0)},\boldsymbol{p}^{(0)})
=\arg\min\_{c\>0,\\p\in\Delta}\lVert y-c\tilde\mu p\rVert_2^2\\. The
start is floored into the open simplex. At that composition, \\\hat c
=(\boldsymbol{m}^{\top}\boldsymbol{W}\boldsymbol{m})^{-1}
\boldsymbol{m}^{\top}\boldsymbol{W}\boldsymbol{y}\\ and
\\\hat\phi=\boldsymbol{r}^{\top}\boldsymbol{W}\boldsymbol{r}/G\\, with
\\\boldsymbol{W}=\boldsymbol{S}\_W(\boldsymbol{p})^{-1}\\. Then \\\hat
a=\hat\phi/\hat c\\ and \\\hat n=\hat c^2/\hat\phi\\. \\\hat n\\ is a
cell count only when the two platforms share their scale; it is not
rounded.

\\\hat c\\ is not the library-size ratio \\\mathbf{1}^{\top}y/\sum_j p_j
S_j\\.

## Usage

``` r
single_cell_nuisance_scales(y, mean_signature_matrix, Sigma, p = NULL)
```

## Arguments

- y:

  Bulk vector of length \\G\\.

- mean_signature_matrix:

  Reference mean, \\G\times J\\.

- Sigma:

  Cell-level covariances, \\G\times G\times J\\.

- p:

  Optional composition. `NULL` runs the free-scale start.

## Value

List with `p`, `ols_scale` (the least-squares \\c^{(0)}\\, or `NA` when
`p` was supplied), `mean_scale`, `variance_scale`, `platform_scale` and
`n_cells`.

## See also

[`deconvolute_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_single_cell.md),
[`deconvolute_ratios_deconrnaseq()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)

## Examples

``` r
mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
set.seed(1)
y <- 30 * drop(mu %*% c(0.6, 0.4)) + rnorm(4, sd = 0.5)
single_cell_nuisance_scales(y, mu, sigma)
#> $p
#> [1] 0.596913 0.403087
#> 
#> $ols_scale
#> [1] 30.13053
#> 
#> $mean_scale
#> [1] 30.13053
#> 
#> $variance_scale
#> [1] 2.82171
#> 
#> $platform_scale
#> [1] 0.09364953
#> 
#> $n_cells
#> [1] 321.7371
#> 
```

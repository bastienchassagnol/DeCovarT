# Expected Fisher information of unconstrained \\\boldsymbol{p}\\

For the multivariate-normal mean–covariance map of the DeCovarT
convolution,
\\\boldsymbol{y}\sim\mathcal{N}\_{G}(\boldsymbol{\mu}\boldsymbol{p},
\boldsymbol{\Sigma}(\boldsymbol{p}))\\ with
\\\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
p_j^{2}\boldsymbol{\Sigma}\_j\\ and precision
\\\boldsymbol{\Omega}(\boldsymbol{p})=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\\,
the expected Fisher information has entries \$\$ I(\boldsymbol{p})\_{jk}
= \boldsymbol{\mu}\_{\cdot j}^{\top} \boldsymbol{\Omega}(\boldsymbol{p})
\boldsymbol{\mu}\_{\cdot k} + 2 p_j p_k\\ \mathrm{tr}\bigl(
\boldsymbol{\Omega}(\boldsymbol{p})\boldsymbol{\Sigma}\_j
\boldsymbol{\Omega}(\boldsymbol{p})\boldsymbol{\Sigma}\_k \bigr). \$\$
The first summand is the mean contribution (an
\\\boldsymbol{\Omega}\\-inner product of signature columns); the second
is the covariance contribution of the quadratic map
\\\boldsymbol{p}\mapsto\boldsymbol{\Sigma}(\boldsymbol{p})\\. See the
multivariate-normal formula on
<https://en.wikipedia.org/wiki/Fisher_information#Multivariate_normal_distribution>.

## Usage

``` r
expected_fisher_unconstrained(p, mean_signature_matrix, Sigma)
```

## Arguments

- p:

  Numeric proportions on the open simplex.

- mean_signature_matrix:

  Mean signature \\\boldsymbol{\mu}\\ (\\G\times J\\).

- Sigma:

  Cell-type covariances \\G\times G\times J\\.

## Value

Symmetric \\J\times J\\ expected Fisher information matrix
\\I(\boldsymbol{p})\\.

## See also

[`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md),
[`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),
[`.inner_product()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-inner_product.md)

## Examples

``` r
# Two cell types, bivariate Gaussians: I(p) is 2 x 2.
p <- c(0.6, 0.4)
mu <- cbind(c(0, 0), c(3, 0))
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
I <- expected_fisher_unconstrained(p, mu, Sigma)
I
#>          [,1]      [,2]
#> [1,] 5.325444  3.550296
#> [2,] 3.550296 19.674556
```

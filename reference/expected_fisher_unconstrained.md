# Expected Fisher information of unconstrained \\\boldsymbol{p}\\

For the multivariate-normal mean–covariance map of the DeCovarT
convolution,
\\\boldsymbol{y}\sim\mathcal{N}\_{G}(\boldsymbol{\mu}\boldsymbol{p},
\boldsymbol{\Sigma}(\boldsymbol{p}))\\ with
\\\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
p_j^{2}\boldsymbol{\Sigma}\_j\\ and precision
\\\boldsymbol{\Theta}(\boldsymbol{p})=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\\,
the expected Fisher information has entries \$\$ I(\boldsymbol{p})\_{jk}
= \boldsymbol{\mu}\_{\cdot j}^{\top} \boldsymbol{\Theta}(\boldsymbol{p})
\boldsymbol{\mu}\_{\cdot k} + 2 p_j p_k\\ \mathrm{tr}\bigl(
\boldsymbol{\Theta}(\boldsymbol{p})\boldsymbol{\Sigma}\_j
\boldsymbol{\Theta}(\boldsymbol{p})\boldsymbol{\Sigma}\_k \bigr). \$\$
The first summand is the mean contribution (an
\\\boldsymbol{\Theta}\\-inner product of signature columns); the second
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

[`vcov_alr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_alr_delta.md),
[`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),
[`.inner_product()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-inner_product.md)

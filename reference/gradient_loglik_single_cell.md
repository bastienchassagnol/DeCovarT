# Score of the single-cell log-likelihood in \\\boldsymbol{p}\\

With \\\boldsymbol{W}=\boldsymbol{S}\_W(\boldsymbol{p})^{-1}\\ and
\\\boldsymbol{r}=\boldsymbol{y}-c\tilde{\boldsymbol{\mu}}\boldsymbol{p}\\,
\$\$ U_j = \phi^{-1}\bigl( c\\\tilde{\boldsymbol{\mu}}\_{\cdot
j}^{\top}\boldsymbol{W}\boldsymbol{r} +\tfrac{1}{2}
\boldsymbol{r}^{\top}\boldsymbol{W}\boldsymbol{\Sigma}\_j
\boldsymbol{W}\boldsymbol{r} \bigr)
-\tfrac{1}{2}\mathrm{tr}(\boldsymbol{W}\boldsymbol{\Sigma}\_j). \$\$
This is the score in the article. It is not
[`gradient_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_unconstrained.md),
whose trace term carries an extra \\p_j\\ because
\\\partial\boldsymbol{\Sigma}/\partial p_j
=2p_j\boldsymbol{\Sigma}\_j\\.

## Usage

``` r
gradient_loglik_single_cell(
  p,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale
)
```

## Arguments

- p:

  Composition of length \\J\\. It does not have to sum to 1 (numerical
  derivatives step off the simplex).

- y:

  Bulk vector of length \\G\\.

- mean_signature_matrix:

  Reference mean, \\G\times J\\.

- Sigma:

  Cell-level covariances, \\G\times G\times J\\.

- mean_scale:

  Positive plug-in \\c=an\\.

- variance_scale:

  Positive plug-in \\\phi=a^2 n\\.

## Value

Numeric vector of length \\J\\.

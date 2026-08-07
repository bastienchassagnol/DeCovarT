# Hessian \\\mathbf{H}\\ of the unconstrained log-likelihood

Analytic Hessian \\\mathbf{H}\in\mathcal{M}\_{J\times J}\\ with entries
\\\mathbf{H}\_{i,j}=\partial^{2}\ell/(\partial p_i\partial p_j)\\,
matching the matricial formulae of the article (quadratic forms in
\\\boldsymbol{\Theta}\\, \\\boldsymbol{\Sigma}\_i\\,
\\\boldsymbol{\mu}\_{\cdot i}\\ and residual
\\\boldsymbol{r}=\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}\\).

## Usage

``` r
hessian_loglik_unconstrained(p, y, mean_signature_matrix, Sigma)
```

## Arguments

- p:

  Numeric vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\.

- y:

  Numeric vector (or one-column matrix)
  \\\boldsymbol{y}\in\mathbb{R}^{G}\\.

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (plug-in means).

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

## Value

Symmetric numeric matrix \\\mathbf{H}\\.

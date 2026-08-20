# Unconstrained DeCovarT log-likelihood

Evaluates the conditional log-likelihood
\\\ell\_{\boldsymbol{y}\\\|\\\boldsymbol{\zeta}}(\boldsymbol{p})\\ of a
bulk profile under the Gaussian convolution model of the article, \$\$
\boldsymbol{y}\\\|\\(\boldsymbol{\zeta},\boldsymbol{p})
\sim\mathcal{N}\_{G}\\\bigl(\boldsymbol{\mu}\boldsymbol{p},\\
\boldsymbol{\Sigma}(\boldsymbol{p})\bigr), \$\$ with plug-in parameters
\\\boldsymbol{\zeta}=(\boldsymbol{\mu},\\\boldsymbol{\Sigma}\_j\\\_{j=1}^{J})\\
and mixture covariance
\\\boldsymbol{\Sigma}(\boldsymbol{p})=\sum\_{j}p_j^{2}\boldsymbol{\Sigma}\_j\\.

## Usage

``` r
loglik_multivariate(p, y, mean_signature_matrix, Sigma)
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

Scalar log-likelihood value.

## Details

Up to an additive constant independent of \\\boldsymbol{p}\\, \$\$
\ell\_{\boldsymbol{y}\\\|\\\boldsymbol{\zeta}}(\boldsymbol{p}) =
-\log\det\boldsymbol{\Sigma}(\boldsymbol{p}) -\tfrac{1}{2}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\mathsf{T}}
\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}). \$\$ Argument
`mean_signature_matrix` stores the plug-in mean signature
\\\boldsymbol{\mu}\\. Latent sample-specific profiles
\\\boldsymbol{x}\_{\cdot j}\\ are **not** observed; the frequentist
likelihood treats \\\boldsymbol{\mu}\\ as a fixed proxy. Estimating
those latents jointly with \\\boldsymbol{p}\\ requires a Bayesian / MAP
step (see
[`.map_gaussian_convolution()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-map_gaussian_convolution.md)).

## See also

[`gradient_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_unconstrained.md),
[`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)

## Examples

``` r
genes <- paste0("g", 1:2)
cts <- paste0("ct", 1:2)
mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
Sigma <- array(
  c(1, 0, 0, 1, 1, 0, 0, 1),
  dim = c(2, 2, 2),
  dimnames = list(genes, genes, cts)
)
p <- c(0.6, 0.4)
y <- drop(mu %*% p)
loglik_multivariate(p, y, mu, Sigma)
#> [1] 1.307853
```

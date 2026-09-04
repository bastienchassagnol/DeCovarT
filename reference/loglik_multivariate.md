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
-\tfrac{1}{2}\log\det\boldsymbol{\Sigma}(\boldsymbol{p}) -\tfrac{1}{2}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\mathsf{T}}
\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}). \$\$ Both terms carry
the factor \\1/2\\ of the Gaussian log-density. An earlier release used
\\-\log\det\boldsymbol{\Sigma}(\boldsymbol{p})\\, which doubled the
determinant contribution and left the objective inconsistent with
[`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md).
With the factor restored,
\\\mathbb{E}\[-\mathbf{H}\]=I(\boldsymbol{p})\\ exactly.

Computationally this is the same Cholesky-and-backsolve evaluation as
`mvtnorm::dmvnorm(..., log = TRUE)` (Genz and Bretz), omitting only the
additive \\-\tfrac{G}{2}\log(2\pi)\\ that does not depend on
\\\boldsymbol{p}\\. The cached factor from
[`.sigma_p_factorisation()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-sigma_p_factorisation.md)
supplies the upper-triangular \\\boldsymbol{R}\\ with
\\\boldsymbol{\Sigma}(\boldsymbol{p})=\boldsymbol{R}^{\mathsf{T}}\boldsymbol{R}\\;
the Mahalanobis term is then \\\lVert\boldsymbol{R}^{-\mathsf{T}}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})\rVert^{2}\\, obtained by
[`base::backsolve()`](https://rdrr.io/r/base/backsolve.html) without
forming the explicit inverse. The inverse is still cached because the
analytic score and Hessian need \\\boldsymbol{\Theta}(\boldsymbol{p})
=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\\. A QR factorisation of
\\\boldsymbol{\Sigma}(\boldsymbol{p})\\ would be a more expensive route
to the same SPD quantities; Cholesky is the natural factorisation.

Argument `mean_signature_matrix` stores the plug-in mean signature
\\\boldsymbol{\mu}\\. Latent sample-specific profiles
\\\boldsymbol{x}\_{\cdot j}\\ are **not** observed; the frequentist
likelihood treats \\\boldsymbol{\mu}\\ as a fixed proxy. Estimating
those latents jointly with \\\boldsymbol{p}\\ requires a Bayesian / MAP
step (see
[`.map_gaussian_convolution()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-map_gaussian_convolution.md)).

## See also

[`gradient_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_unconstrained.md),
[`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md)

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
#> [1] 0.6539265
```

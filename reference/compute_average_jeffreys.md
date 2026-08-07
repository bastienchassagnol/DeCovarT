# Average pairwise Jeffreys divergence of a Gaussian mixture

Returns the \\p_j p\_{\ell}\\-weighted mean of pairwise Jeffreys
(symmetrised KL) divergences between purified Gaussians \$\$
\overline{J} = \frac{ \sum\_{1\le j\<\ell\le J} p_j p\_{\ell}\\
J(f\_{j},f\_{\ell}) }{ \sum\_{1\le j\<\ell\le J} p_j p\_{\ell} }
\in\[0,\infty), \$\$ with \\f\_{j}=\mathcal{N}\_{G}(
\boldsymbol{\mu}\_{\cdot j},\boldsymbol{\Sigma}\_{j} )\\. If `p` is
omitted it defaults to the equi-balanced vector \\1/J\\, which recovers
the uniform pairwise average
\\2/(J(J-1))\sum\_{j\<\ell}J(f\_{j},f\_{\ell})\\.

## Usage

``` r
compute_average_jeffreys(true_theta, J = NULL)
```

## Arguments

- true_theta:

  List validated by
  [`check_true_theta()`](https://bastienchassagnol.github.io/DeCovarT/reference/check_true_theta.md):
  `mu` (\\G\times J\\), `sigma` (\\G\times G\times J\\), and optionally
  `p` (length \\J\\ or \\J\times N\\). If `p` is missing it is set to
  \\(1/J,\ldots,1/J)\\.

- J:

  Number of cell types. Defaults to the third dimension of `sigma`.

## Value

Scalar average pairwise Jeffreys divergence.

## See also

[`check_true_theta()`](https://bastienchassagnol.github.io/DeCovarT/reference/check_true_theta.md),
[`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md),
[`compute_glmnet_gene_scores()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_glmnet_gene_scores.md)

## Examples

``` r
set.seed(1)
theta <- list(
  mu = cbind(c(0, 0), c(3, 0)),
  sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
)
compute_average_jeffreys(theta)
#> [1] 9
```

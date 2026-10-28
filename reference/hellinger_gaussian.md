# Hellinger distance between two multivariate Gaussians

Closed form through the Bhattacharyya coefficient \$\$ \log\mathrm{BC} =
\tfrac14\log\lvert\boldsymbol{\Sigma}\_j\rvert
+\tfrac14\log\lvert\boldsymbol{\Sigma}\_\ell\rvert
-\tfrac12\log\lvert\bar{\boldsymbol{\Sigma}}\rvert
-\tfrac18\boldsymbol{\delta}^{\top}\bar{\boldsymbol{\Sigma}}^{-1}
\boldsymbol{\delta}, \$\$ with
\\\bar{\boldsymbol{\Sigma}}=(\boldsymbol{\Sigma}\_j+
\boldsymbol{\Sigma}\_\ell)/2\\ and \\\boldsymbol{\delta}=
\boldsymbol{\mu}\_j-\boldsymbol{\mu}\_\ell\\; then
\\H=\sqrt{1-\mathrm{BC}}\in\[0,1\]\\. It is symmetric, unlike the MixSim
misclassification overlap
([`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md)),
and bounded, unlike the Jeffreys divergence
([`compute_average_jeffreys()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_jeffreys.md)).

## Usage

``` r
hellinger_gaussian(mu_j, mu_l, sigma_j, sigma_l)
```

## Arguments

- mu_j, mu_l:

  Mean vectors of length \\G\\.

- sigma_j, sigma_l:

  \\G\times G\\ covariance matrices.

## Value

The Hellinger distance, a number in \\\[0,1\]\\.

## Examples

``` r
hellinger_gaussian(c(0, 0), c(1, 1), diag(2), diag(2))
#> [1] 0.4703182
```

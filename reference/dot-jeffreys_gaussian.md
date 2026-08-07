# Pairwise Jeffreys divergence between two multivariate Gaussians

Symmetrised KL divergence
\\J(f\_{j},f\_{\ell})=D\_{\mathrm{KL}}(f\_{j}\parallel f\_{\ell})+
D\_{\mathrm{KL}}(f\_{\ell}\parallel f\_{j})\\ for
\\f\_{j}=\mathcal{N}\_{G}( \boldsymbol{\mu}\_{\cdot
j},\boldsymbol{\Sigma}\_{j} )\\ and likewise for \\\ell\\, using the
closed form in the feature-selection vignette.

## Usage

``` r
.jeffreys_gaussian(mu_j, mu_l, sigma_j, sigma_l)
```

## Arguments

- mu_j, mu_l:

  Numeric mean vectors of length \\G\\.

- sigma_j, sigma_l:

  Numeric \\G\times G\\ covariances.

## Value

Non-negative scalar Jeffreys divergence.

## References

Kullback S, Leibler RA (1951). "On Information and Sufficiency." *The
Annals of Mathematical Statistics* 22(1), 79–86.
[doi:10.1214/aoms/1177729694](https://doi.org/10.1214/aoms/1177729694) .

Multivariate normal KL closed form:
<https://statproofbook.github.io/P/mvn-kl.html>.

Symmetrised (Jeffreys) divergence:
<https://en.wikipedia.org/wiki/Kullback-Leibler_divergence>.

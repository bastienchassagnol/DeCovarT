# DeCovarT: covariance-aware bulk transcriptomic deconvolution

**DeCovarT** estimates cellular proportions from bulk RNA-seq by
modelling mixtures as Gaussian convolutions of purified cell-type means
and covariances. Notation follows the article: genes \\g=1,\ldots,G\\,
cell types \\j=1,\ldots,J\\, samples \\i=1,\ldots,N\\; bulk
\\\boldsymbol{y}\\, mean signature \\\boldsymbol{\mu}\\, proportions
\\\boldsymbol{p}\\, covariances / precisions
\\\boldsymbol{\Sigma}\_j\\/\\\boldsymbol{\Theta}\_j\\. Proportions live
on the open simplex and are optimised in unconstrained ALR coordinates
\\\boldsymbol{\rho}\in\mathbb{R}^{J-1}\\
([`vignette("softmax-alr-derivatives", package = "DeCovarT")`](https://bastienchassagnol.github.io/DeCovarT/articles/softmax-alr-derivatives.md)).

The frequentist API plugs in \\\boldsymbol{\mu}\\ for unobserved latent
profiles \\\boldsymbol{x}\_{\cdot j}\\; MAP recovery of those latents is
the Bayesian extension.

To learn more, start with the package website vignettes, or in an R
session: `browseVignettes(package = "DeCovarT")` and
[`?DeCovarT::deconvolute_ratios`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md).

Reference manuals (PDF / HTML) can be regenerated with
`source("scripts/auxiliary/generate_package_manual.R")`.

## See also

Useful links:

- <https://github.com/bastienchassagnol/DeCovarT>

- <https://bastienchassagnol.github.io/DeCovarT/>

- Report bugs at <https://github.com/bastienchassagnol/DeCovarT/issues>

## Author

**Maintainer**: Bastien Chassagnol <bastien_chassagnol@laposte.net>
([ORCID](https://orcid.org/0000-0002-8955-2391))

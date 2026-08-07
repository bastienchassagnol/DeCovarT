# Gene scores from multinomial elastic-net cell-type classification

Fits a multinomial elastic net
([`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html))
that predicts cell type from expression features, using only the mean
signature \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\ (no
covariances). Each column \\\boldsymbol{\mu}\_{\cdot j}\\ is expanded
into `n_rep` isotropic Gaussian perturbations so that cross-validation
is well-defined with one prototype per type. Gene scores are the sum
over classes of absolute coefficients at `lambda.min` (intercept
excluded).

This is the `glmnet` screen in the four-score shortlist of the
feature-selection vignette (paired with Jeffreys, MixSim overlap and
DEGs).

## Usage

``` r
compute_glmnet_gene_scores(
  mu,
  alpha = 0.5,
  n_rep = 20L,
  noise_sd = NULL,
  nfolds = NULL,
  ...
)
```

## Arguments

- mu:

  Numeric mean signature \\\boldsymbol{\mu}\\ with genes in rows and
  cell types in columns (\\G\times J\\).

- alpha:

  Elastic-net mixing parameter in \\\[0,1\]\\ (default `0.5`).

- n_rep:

  Number of noisy replicates per cell-type mean (default `20L`).

- noise_sd:

  Isotropic Gaussian noise sd for replicates. Default is `1e-3` times
  the mean absolute column scale of `mu`.

- nfolds:

  Number of CV folds (default `min(10L, n_rep)`).

- ...:

  Additional arguments forwarded to
  [`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html).

## Value

Named numeric vector of length \\G\\ (gene scores; larger means stronger
multinomial signal).

## See also

[`compute_average_jeffreys()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_jeffreys.md),
[`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md)

## Examples

``` r
set.seed(1)
mu <- cbind(c(0, 0, 5), c(5, 0, 0), c(0, 5, 0))
rownames(mu) <- paste0("g", seq_len(nrow(mu)))
scores <- compute_glmnet_gene_scores(mu)
names(scores)[which.max(scores)]
#> [1] "g2"
```

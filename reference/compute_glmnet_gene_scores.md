# Gene scores from multinomial elastic-net cell-type classification

Fits a multinomial (or binomial) elastic net
([`glmnet::glmnet()`](https://rdrr.io/pkg/glmnet/man/glmnet.html)) that
predicts cell type from expression features. Inputs are purified
expression profiles \\\boldsymbol{X}\in\mathcal{M}\_{G\times J\times
N}\\ (genes \\\times\\ cell types \\\times\\ samples) and length-\\J\\
cell-type labels. Variability across samples replaces synthetic
isotropic noise. Gene scores are the sum over classes of absolute
coefficients at a chosen \\\lambda\\ (intercept excluded). For nested /
CV selection of \\\lambda\\, see the experimental
`compute_glmnet_gene_scores_cv()` helper (not shipped in the package
build).

## Usage

``` r
compute_glmnet_gene_scores(
  expression_profiles,
  celltype_labels,
  alpha = 0.5,
  lambda = NULL,
  ...
)
```

## Arguments

- expression_profiles:

  Numeric array \\G\times J\times N\\ of purified profiles.

- celltype_labels:

  Character or factor labels of length \\J\\ (one per cell-type slice).

- alpha:

  Elastic-net mixing parameter in \\\[0,1\]\\ (default `0.5`).

- lambda:

  Optional penalty value at which coefficients are extracted. When
  `NULL`, uses the smallest \\\lambda\\ on the fitted path.

- ...:

  Additional arguments forwarded to
  [`glmnet::glmnet()`](https://rdrr.io/pkg/glmnet/man/glmnet.html).

## Value

Named numeric vector of length \\G\\ (gene scores; larger means stronger
multinomial signal).

## See also

[`compute_average_jeffreys()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_jeffreys.md),
[`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md)

## Examples

``` r
set.seed(1)
G <- 4L
J <- 3L
N <- 12L
profiles <- array(0, dim = c(G, J, N))
for (j in seq_len(J)) {
  profiles[j, j, ] <- 5 + stats::rnorm(N, sd = 0.2)
}
labels <- paste0("ct", seq_len(J))
scores <- compute_glmnet_gene_scores(profiles, labels)
names(scores)[which.max(scores)]
#> [1] "gene_3"
```

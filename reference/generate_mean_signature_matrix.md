# Generate mean profiles from a target Gram matrix

Builds \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\ whose *unit*
columns realise a prescribed Gram matrix \\R\\ (pairwise cosines) via
the symmetric square root, \$\$ \boldsymbol{\mu} = s\\ \boldsymbol{Q}
R^{1/2}, \qquad
\boldsymbol{Q}^{\mathsf{T}}\boldsymbol{Q}=\boldsymbol{I}\_{J}. \$\$ Then
\\\boldsymbol{\mu}^{\mathsf{T}}\boldsymbol{\mu}=s^{2}R\\ whenever \\R\\
has unit diagonal. The default \\R\\ is the equicorrelation matrix
\\(1-\rho)I+\rho\mathbf{1}\mathbf{1}^{\mathsf{T}}\\
([`equicorrelation_gram()`](https://bastienchassagnol.github.io/DeCovarT/reference/equicorrelation_gram.md)).
Supply `target_gram` to set unequal pairwise cosines (for example two
close pairs and a distant background), provided \\R\\ is symmetric
nonnegative-definite.

## Usage

``` r
generate_mean_signature_matrix(
  n_genes,
  n_celltypes,
  mean_scale = 10,
  target_cosine = 0,
  target_gram = NULL,
  seed = NULL,
  nonnegative = FALSE,
  gene_names = NULL,
  celltype_names = NULL
)
```

## Arguments

- n_genes:

  Integer \\G\\; must be at least `n_celltypes`.

- n_celltypes:

  Integer \\J\ge 2\\.

- mean_scale:

  Positive scalar \\s\\ (centroid norms). Default `10`. Hold fixed when
  studying cosine / collinearity alone.

- target_cosine:

  Numeric pairwise cosine \\\rho\\ for the equicorrelation Gram. Ignored
  when `target_gram` is supplied. Must lie in \\\[-1/(J-1),1\]\\.

- target_gram:

  Optional symmetric \\J\times J\\ cosine matrix (unit diagonal).
  Overrides `target_cosine`.

- seed:

  Optional integer. When supplied, the orthonormal frame
  \\\boldsymbol{Q}\\ is drawn from a Gaussian QR; otherwise it is
  deterministic. Ignored when `nonnegative = TRUE`.

- nonnegative:

  If `TRUE`, use a disjoint-support nonnegative frame so that
  \\\boldsymbol{\mu}\ge 0\\ whenever \\R^{1/2}\ge 0\\ (required by
  [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)).
  The Gram is unchanged.

- gene_names:

  Optional character vector of length \\G\\.

- celltype_names:

  Optional character vector of length \\J\\.

## Value

Numeric matrix \\\boldsymbol{\mu}\\ with dimensions \\G\times J\\.

## Details

Realisations are unique up to an orthogonal transformation of the gene
space: replacing \\\boldsymbol{Q}\\ by \\\boldsymbol{Q}U\\ with
\\U^{\mathsf{T}}U=\boldsymbol{I}\_{J}\\ leaves the Gram unchanged. By
default \\\boldsymbol{Q}\\ is a deterministic thin QR frame; pass `seed`
for a Haar-like Gaussian frame. Cholesky \\R=LL^{\mathsf{T}}\\ is an
alternative square root used in Monte Carlo simulation; see the
how-to-build-synthetic-scenarios vignette appendix.

## Examples

``` r
generate_mean_signature_matrix(
  n_genes = 6L,
  n_celltypes = 2L,
  target_cosine = 0.5
)
#>        celltype_1 celltype_2
#> gene_1  -6.607061   1.260622
#> gene_2  -2.187420  -8.163563
#> gene_3  -3.590332  -2.818114
#> gene_4  -3.590332  -2.818114
#> gene_5  -3.590332  -2.818114
#> gene_6  -3.590332  -2.818114
k <- matrix(c(1, 0.98, 0.2, 0.98, 1, 0.2, 0.2, 0.2, 1), 3, 3)
generate_mean_signature_matrix(
  n_genes = 8L,
  n_celltypes = 3L,
  target_gram = k
)
#>        celltype_1 celltype_2 celltype_3
#> gene_1 -3.7888055 -2.2171550   2.058065
#> gene_2 -5.1519186 -6.3537691   1.730545
#> gene_3 -0.7386983 -0.7386983  -8.784618
#> gene_4 -3.4221926 -3.2914485  -1.766356
#> gene_5 -3.4221926 -3.2914485  -1.766356
#> gene_6 -3.4221926 -3.2914485  -1.766356
#> gene_7 -3.4221926 -3.2914485  -1.766356
#> gene_8 -3.4221926 -3.2914485  -1.766356
```

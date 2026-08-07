# Generate mean profiles with a target pairwise cosine

Builds \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\ by blending a
shared unit direction \\\boldsymbol{u}\\ with cell-type-private
orthogonal marker directions \\\boldsymbol{v}\_{j}\\: \$\$
\tilde{\boldsymbol{\mu}}\_{\cdot j} = \sqrt{\rho}\\\boldsymbol{u}
+\sqrt{1-\rho}\\\boldsymbol{v}\_{j}, \qquad \boldsymbol{\mu}\_{\cdot j}
= s\\ \frac{\tilde{\boldsymbol{\mu}}\_{\cdot j}}{
\\\tilde{\boldsymbol{\mu}}\_{\cdot j}\\\_2 }. \$\$ The private vectors
\\\boldsymbol{v}\_{j}\\ are indicator directions on a partition of the
\\G\\ genes (type \\j\\ only) and then \\\ell_2\\-normalised, so
\\\boldsymbol{v}\_{j}^{\mathsf{T}}\boldsymbol{v}\_{k}=0\\ for \\j\neq
k\\. With a shared unit \\\boldsymbol{u}\\, \$\$
\tilde{\boldsymbol{\mu}}\_{\cdot j}^{\mathsf{T}}
\tilde{\boldsymbol{\mu}}\_{\cdot k} = \rho + \sqrt{\rho(1-\rho)}\\
\bigl( \boldsymbol{u}^{\mathsf{T}}\boldsymbol{v}\_{j} +
\boldsymbol{u}^{\mathsf{T}}\boldsymbol{v}\_{k} \bigr) \qquad (j\neq k).
\$\$ After column normalisation the pairwise cosines of
\\\boldsymbol{\mu}\\ therefore track \\\rho\\ closely when the cross
terms \\\boldsymbol{u}^{\mathsf{T}}\boldsymbol{v}\_{j}\\ are small
relative to the leading \\\rho\\ (many genes per block). The global
scale \\s\\ sets column norms (and hence Euclidean separation) without
changing angles: for fixed \\\rho\\, \\\\\boldsymbol{\mu}\_{\cdot
j}-\boldsymbol{\mu}\_{\cdot k}\\\_2 \propto s\\. Prefer dialling
\\\rho\\ when second-order precision weights already control interaction
strength; keep \\s\\ fixed across scenarios that compare mean
collinearity alone (Aliee and Theis 2021) .

## Usage

``` r
generate_mean_signature_matrix(
  n_genes,
  n_celltypes,
  mean_scale = 10,
  target_cosine = 0,
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

  Positive scalar \\s\\ (centroid norms). Default `10`, as in the nine
  factorial scenarios. Hold fixed when studying cosine / collinearity
  alone.

- target_cosine:

  Numeric in \\\[0,1\]\\, the collinearity dial \\\rho\\.

- gene_names:

  Optional character vector of length \\G\\.

- celltype_names:

  Optional character vector of length \\J\\.

## Value

Numeric matrix \\\boldsymbol{\mu}\\ with dimensions \\G\times J\\.

## Details

**Private marker blocks.** Genes are partitioned into \\J\\ nearly equal
contiguous blocks. Type \\j\\'s private direction
\\\boldsymbol{v}\_{j}\\ is the indicator of its block, then
\\\ell_2\\-normalised. Distinct blocks are orthogonal, so type-specific
signal does not leak across columns before the shared component is
added.

**Shared–private blend.** With unit shared direction
\\\boldsymbol{u}=G^{-1/2}\mathbf{1}\\, each column is
\\\sqrt{\rho}\\\boldsymbol{u}+\sqrt{1-\rho}\\\boldsymbol{v}\_{j}\\,
re-normalised, then scaled by \\s\\. Thus \\\rho\\ dials collinearity
while \\s\\ dials Euclidean separation without changing angles.

## Examples

``` r
generate_mean_signature_matrix(
  n_genes = 6L,
  n_celltypes = 2L,
  target_cosine = 0.5
)
#>        celltype_1 celltype_2
#> gene_1   5.334021   2.209424
#> gene_2   5.334021   2.209424
#> gene_3   5.334021   2.209424
#> gene_4   2.209424   5.334021
#> gene_5   2.209424   5.334021
#> gene_6   2.209424   5.334021
```

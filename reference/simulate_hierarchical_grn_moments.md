# Simulate GRN first- and second-order moments

Builds a mean matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\ and
a shared covariance array \\(\boldsymbol{\Sigma}\_j)\_{j}\\ under a
graph-constrained precision model. Means follow the AutoGeneS-inspired
cosine construction of
[`generate_mean_signature_matrix()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_mean_signature_matrix.md)
with target pairwise cosine \\\rho\\. The adjacency \\\boldsymbol{A}\\
is drawn from a random-graph model (igraph: scale-free, stochastic
block, or Watts–Strogatz small-world); i.i.d. signed weights with
inhibitory fraction `prop_inhibitory` form \\\boldsymbol{W}\\; the
precision is completed by a spectral shift; each cell type shares
\\\boldsymbol{\Sigma}\_j=\boldsymbol{\Omega}^{-1}\\.

## Usage

``` r
simulate_hierarchical_grn_moments(
  n_genes,
  n_celltypes = 2L,
  mean_scale = 10,
  target_cosine = 0,
  precision_shift,
  precision_scale,
  prop_inhibitory = 0.5,
  graph_model = c("scale_free", "stochastic_block_model", "small_world"),
  graph_params = list()
)
```

## Arguments

- n_genes:

  Integer; number of genes \\G\\.

- n_celltypes:

  Integer; number of cell types \\J\\ (default 2).

- mean_scale:

  Positive scalar \\s\\ for centroid norms (default `10`, as in the
  nine-scenario grid).

- target_cosine:

  Numeric in \\\[0,1\]\\; target pairwise cosine similarity between
  columns of \\\boldsymbol{\mu}\\.

- precision_shift:

  Diagonal cushion \\u\\ for the spectral shift.

- precision_scale:

  Positive magnitude \\v\\ of signed off-diagonal precision weights.

- prop_inhibitory:

  Numeric in \\\[0,1\]\\; fraction of edges with positive precision
  weight (inhibitory partial correlation). Default `0.5` balances
  inhibitory and activatory edges.

- graph_model:

  One of `"scale_free"`, `"stochastic_block_model"`, `"small_world"`.

- graph_params:

  Named list of generator parameters (see
  [`generate_random_network_skeleton()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_random_network_skeleton.md)).

## Value

Named list with:

- `mean_profiles`: matrix \\\boldsymbol{\mu}\\;

- `covariance_matrices`: array
  \\(\boldsymbol{\Sigma}\_j)\_{j}\in\mathcal{M}\_{G\times G\times J}\\;

- `graph_structure`: `adjacency_matrix`, `weighted_adjacency`, and
  `normalised_precision` \\\boldsymbol{\Omega}\\;

- `objectives`: `mean_abs_cosine` and `sum_euclidean_distance`.

## Examples

``` r
set.seed(42)
moments <- simulate_hierarchical_grn_moments(
  n_genes = 40L,
  n_celltypes = 3L,
  mean_scale = 10,
  target_cosine = 0.1,
  precision_shift = 0.1,
  precision_scale = 0.3,
  prop_inhibitory = 0.5,
  graph_model = "scale_free"
)
str(moments, max.level = 2)
#> List of 4
#>  $ mean_profiles      : num [1:40, 1:3] 2.61 2.61 2.61 2.61 2.61 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ covariance_matrices: num [1:40, 1:40, 1:3] 1.809 1.376 -0.62 0.372 -0.535 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ graph_structure    :List of 3
#>   ..$ adjacency_matrix    : int [1:40, 1:40] 0 1 0 0 0 0 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ weighted_adjacency  : num [1:40, 1:40] 0 -0.3 0 0 0 0 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ normalised_precision: num [1:40, 1:40] 0.967 -0.3 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>  $ objectives         :List of 2
#>   ..$ mean_abs_cosine       : num 0.332
#>   ..$ sum_euclidean_distance: num 34.7
```

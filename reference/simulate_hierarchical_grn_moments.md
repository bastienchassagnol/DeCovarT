# Simulate GRN first- and second-order moments

Builds a mean matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\ and
**cell-type-specific** second-order moments \\(\boldsymbol{\Omega}\_j,
\boldsymbol{\Sigma}\_j=\boldsymbol{\Omega}\_j^{-1})\_{j=1}^{J}\\ under a
graph-constrained precision model. Means follow the AutoGeneS-inspired
cosine construction of
[`generate_mean_signature_matrix()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_mean_signature_matrix.md)
with target pairwise cosine \\\rho\\. For each cell type, an adjacency
is drawn from a random-graph model (or supplied), i.i.d. signed weights
with inhibitory fraction `prop_inhibitory` form \\\boldsymbol{W}\_j\\,
and the precision is completed by a spectral shift. Distinct cell types
receive **independent** precision draws by default (biology rarely
shares one network across types); pass a length-\\J\\ `graph_model` /
`graph_params` or a pre-built `adjacency` list / array for hybrid
designs.

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
  graph_model = c("erdos_renyi", "hub", "star", "scale_free", "stochastic_block_model",
    "small_world"),
  graph_params = list(),
  adjacency = NULL
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

- precision_shift:

  Diagonal cushion \\u\\ for the spectral shift (scalar, or length
  \\J\\).

- precision_scale:

  Positive magnitude \\v\\ of signed off-diagonal precision weights
  (scalar, or length \\J\\).

- prop_inhibitory:

  Numeric in \\\[0,1\]\\; fraction of edges with positive precision
  weight (inhibitory partial correlation). Default `0.5` balances
  inhibitory and activatory edges (scalar, or length \\J\\).

- graph_model:

  One of `"erdos_renyi"`, `"hub"`, `"scale_free"`,
  `"stochastic_block_model"`, `"small_world"`, or a character vector of
  length \\J\\.

- graph_params:

  Named list of generator parameters (shared), or a list of length \\J\\
  of such named lists (see
  [`generate_random_network_skeleton()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_random_network_skeleton.md)).

- adjacency:

  Optional pre-built undirected adjacencies: a list of \\J\\ \\G\times
  G\\ matrices, or a \\G\times G\times J\\ array. When supplied,
  `graph_model` / `graph_params` are ignored for skeleton generation.

## Value

Named list with:

- `mean_profiles`: matrix \\\boldsymbol{\mu}\\;

- `covariance_matrices`: array
  \\(\boldsymbol{\Sigma}\_j)\_{j}\in\mathcal{M}\_{G\times G\times J}\\;

- `precision_matrices`: array
  \\(\boldsymbol{\Omega}\_j)\_{j}\in\mathcal{M}\_{G\times G\times J}\\;

- `graph_structure`: `adjacency_matrices`, `weighted_adjacencies`, and
  `normalised_precision` (all \\G\times G\times J\\);

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
#> List of 5
#>  $ mean_profiles      : num [1:40, 1:3] 2.61 2.61 2.61 2.61 2.61 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ covariance_matrices: num [1:40, 1:40, 1:3] 1.809 -1.376 -0.62 -0.372 -0.535 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ precision_matrices : num [1:40, 1:40, 1:3] 0.967 0.3 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ graph_structure    :List of 3
#>   ..$ adjacency_matrices  : int [1:40, 1:40, 1:3] 0 1 0 0 0 0 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 3
#>   ..$ weighted_adjacencies: num [1:40, 1:40, 1:3] 0 0.3 0 0 0 0 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 3
#>   ..$ normalised_precision: num [1:40, 1:40, 1:3] 0.967 0.3 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 3
#>  $ objectives         :List of 2
#>   ..$ mean_abs_cosine       : num 0.332
#>   ..$ sum_euclidean_distance: num 34.7
```

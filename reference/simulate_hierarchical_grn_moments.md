# Simulate GRN first- and second-order moments

Builds a mean matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\ and
a shared covariance array \\(\boldsymbol{\Sigma}\_j)\_{j}\\ under a
graph-constrained precision model. Means follow the AutoGeneS-inspired
construction of `generate_mean_signature_matrix()` with target pairwise
cosine similarity \\\rho\\ and column scale \\s\\, so that one can dial
collinearity (minimise cosine) against centroid separation (maximise
Euclidean distance). The adjacency \\\boldsymbol{A}\\ is drawn from a
random-graph model; the precision \\\boldsymbol{\Omega}\\ is obtained
via the affine spectral shift of `build_normalised_precision()`, and
each cell type shares
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
  graph_model = c("power_law", "stochastic_block_model"),
  graph_params = list()
)
```

## Arguments

- n_genes:

  Integer; number of genes \\G\\.

- n_celltypes:

  Integer; number of cell types \\J\\ (default 2).

- mean_scale:

  Positive scalar \\s\\ for centroid norms.

- target_cosine:

  Numeric in \\\[0,1\]\\; target pairwise cosine similarity between
  columns of \\\boldsymbol{\mu}\\.

- precision_shift, precision_scale:

  Diagonal shift \\u\\ and off-diagonal scale \\v\\ used to build
  \\\boldsymbol{\Omega}\\.

- graph_model:

  `"power_law"` or `"stochastic_block_model"`.

- graph_params:

  Named list of generator parameters: `power` / `edges_per_node`
  (power-law) or `block_prob` / `p_within` / `p_between` (SBM).

## Value

Named list with:

- `mean_profiles`: matrix \\\boldsymbol{\mu}\\;

- `covariance_matrices`: array
  \\(\boldsymbol{\Sigma}\_j)\_{j}\in\mathcal{M}\_{G\times G\times J}\\;

- `graph_structure`: `adjacency_matrix` and `normalised_precision`
  \\\boldsymbol{\Omega}\\;

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
  graph_model = "power_law",
  graph_params = list(power = 1, edges_per_node = 2)
)
str(moments, max.level = 2)
#> List of 4
#>  $ mean_profiles      : num [1:40, 1:3] 2.57 2.57 2.57 2.57 2.57 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ covariance_matrices: num [1:40, 1:40, 1:3] 1.773 -0.8118 -0.4864 1.1831 0.0232 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ graph_structure    :List of 2
#>   ..$ adjacency_matrix    : num [1:40, 1:40] 0 1 1 0 0 0 1 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ normalised_precision: num [1:40, 1:40] 1.32 0.3 0.3 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>  $ objectives         :List of 2
#>   ..$ mean_abs_cosine       : num 0.414
#>   ..$ sum_euclidean_distance: num 32.5

## Verify positive-definiteness of the shared covariance
eigen_vals <- eigen(
  moments$covariance_matrices[, , 1],
  only.values = TRUE
)$values
stopifnot(all(eigen_vals > 0))
```

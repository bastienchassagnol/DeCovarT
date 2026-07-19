# Simulate hierarchical GRN first- and second-order moments

Builds mean matrices \\\boldsymbol{\mu}\\ and covariance arrays
\\(\boldsymbol{\Sigma}\_j)\\ for parent/child cell populations under a
graph-constrained precision model. Parent means use complementary block
structure; child means are
\\\boldsymbol{\mu}^{(k)}=\boldsymbol{\mu}^{(\mathrm{parent})}+\boldsymbol{\delta}^{(k)}\\
with
\\\boldsymbol{\delta}^{(k)}\sim\mathcal{N}(\mathbf{0},\sigma\_\delta^{2}\mathbf{I})\\.
Marginal variances follow \\\sigma_g^{2}=\mu_g+\mu_g^{\alpha}/L\\;
covariances are obtained from a shared normalised precision
\\\boldsymbol{\Omega}\\ via
\\\boldsymbol{\Sigma}\_k=\mathbf{D}\_k^{1/2}\mathbf{R}\mathbf{D}\_k^{1/2}\\
with \\\mathbf{R}=\mathrm{cov2cor}(\boldsymbol{\Omega}^{-1})\\.

## Usage

``` r
simulate_hierarchical_grn_moments(
  n_expressed_genes,
  mean_lower_expressed,
  mean_upper_expressed,
  mean_lower_background,
  mean_upper_background,
  library_size,
  alpha,
  precision_shift,
  precision_scale,
  child_perturbation_sd,
  graph_model = c("power_law", "stochastic_block_model"),
  graph_params = list()
)
```

## Arguments

- n_expressed_genes:

  Integer; expressed genes per parent block (total dimension
  \\G=2\times\\ `n_expressed_genes`).

- mean_lower_expressed, mean_upper_expressed:

  Uniform bounds for expressed block means.

- mean_lower_background, mean_upper_background:

  Uniform bounds for background block means.

- library_size:

  Positive scalar \\L\\ in the mean–variance law.

- alpha:

  Positive power in \\\mu_g^{\alpha}/L\\ (`2` recovers a classical
  NB-like law).

- precision_shift, precision_scale:

  Diagonal shift \\u\\ and off-diagonal scale \\v\\ used to build
  \\\boldsymbol{\Omega}\\.

- child_perturbation_sd:

  Standard deviation \\\sigma\_\delta\\.

- graph_model:

  `"power_law"` or `"stochastic_block_model"`.

- graph_params:

  Named list of generator parameters (see Details in source for `power`
  / `edges_per_node` or `block_prob` / `p_within` / `p_between`).

## Value

Named list with `parent_parameters`, `child_parameters` (each holding
`mean_profiles` \\\boldsymbol{\mu}\\ and `covariance_matrices`
\\(\boldsymbol{\Sigma}\_j)\\), and `graph_structure`
(`adjacency_matrix`, `normalised_precision` \\\boldsymbol{\Omega}\\).

## Examples

``` r
set.seed(42)
moments <- simulate_hierarchical_grn_moments(
    n_expressed_genes     = 50,
    mean_lower_expressed  = 2,
    mean_upper_expressed  = 6,
    mean_lower_background = 0.1,
    mean_upper_background = 0.5,
    library_size          = 10000,
    alpha                 = 2,
    precision_shift       = 0.1,
    precision_scale       = 0.3,
    child_perturbation_sd = 0.1,
    graph_model           = "power_law",
    graph_params          = list(power = 1, edges_per_node = 2)
)
str(moments, max.level = 2)
#> List of 3
#>  $ parent_parameters:List of 2
#>   ..$ mean_profiles      : num [1:2, 1:100] 5.659 0.35 5.748 0.187 3.145 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ covariance_matrices: num [1:100, 1:100, 1:2] 5.662 -3.114 -0.286 -2.775 2.177 ...
#>   .. ..- attr(*, "dimnames")=List of 3
#>  $ child_parameters :List of 2
#>   ..$ mean_profiles      : num [1:4, 1:100] 5.779 5.459 0.35 0.484 5.853 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ covariance_matrices: num [1:100, 1:100, 1:4] 5.783 -3.175 -0.284 -2.853 2.184 ...
#>   .. ..- attr(*, "dimnames")=List of 3
#>  $ graph_structure  :List of 2
#>   ..$ adjacency_matrix    : num [1:100, 1:100] 0 1 1 1 0 1 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ normalised_precision: num [1:100, 1:100] 1.68 0.3 0.3 0.3 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2

## Verify positive-definiteness of a child covariance
eigen_vals <- eigen(
    moments$child_parameters$covariance_matrices[, , 1],
    only.values = TRUE
)$values
stopifnot(all(eigen_vals > 0))
```

# Recommend a covariance backend from a network topology

Maps a `graph_model` name (see
[`generate_random_network_skeleton()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_random_network_skeleton.md))
to the covariance-factorisation `structure` most likely to help. This
encodes the decision-tree workflow: the dense Cholesky is always a safe
default; the structured backends are refinements when the covariance
(not merely the precision) inherits exploitable structure.

## Usage

``` r
covariance_structure_from_graph_model(graph_model)
```

## Arguments

- graph_model:

  One of `"erdos_renyi"`, `"hub"`, `"star"`, `"scale_free"`,
  `"stochastic_block_model"`, `"small_world"`, `"band"`, `"ar"`.
  Matching is case-insensitive.

## Value

A single string: one of `"dense"`, `"block"`, `"band"`, `"sparse"`,
`"diag_lowrank"`.

## Details

A sparse gene-network **precision** matrix
\\\boldsymbol{\Omega}\_j=\boldsymbol{\Sigma}\_j^{-1}\\ does **not**
imply a sparse **covariance** \\\boldsymbol{\Sigma}\_j\\; the inverse of
a sparse precision is generally dense. The mapping below therefore
expresses modelling intent, and the returned structure only accelerates
the likelihood when \\\boldsymbol{\Sigma}(\boldsymbol{p})\\ itself (or a
cheap representation of it) has that structure. Disconnected gene
modules give an exactly block-diagonal covariance; shared regulatory
programs give a diagonal-plus-low-rank covariance; a small set of bridge
/ housekeeping hubs gives an arrow pattern handled by the Schur
complement.

## See also

[`new_decovart_covariance()`](https://bastienchassagnol.github.io/DeCovarT/reference/new_decovart_covariance.md)

## Examples

``` r
covariance_structure_from_graph_model("stochastic_block_model")
#> [1] "block"
covariance_structure_from_graph_model("erdos_renyi")
#> [1] "sparse"
```

# Sample a symmetric adjacency from a random-graph family

Draws undirected skeletons with igraph: Erdős–Rényi (`erdos_renyi`), hub
/ star modules (`hub`), Barabási–Albert preferential attachment
(`scale_free`), a stochastic block model (`stochastic_block_model`), or
Watts–Strogatz small-world (`small_world`) (Barabási and Albert 1999) .

## Usage

``` r
generate_random_network_skeleton(n_genes, graph_model, graph_params = list())
```

## Arguments

- n_genes:

  Integer \\G\\, number of nodes (genes).

- graph_model:

  One of `"erdos_renyi"`, `"hub"`, `"scale_free"`,
  `"stochastic_block_model"`, `"small_world"` (also `"star"`, mapped to
  `"hub"`). Matching is case-insensitive.

- graph_params:

  Named list of generator parameters:

  `erdos_renyi`

  :   `p` edge probability for
      [`igraph::sample_gnp()`](https://r.igraph.org/reference/sample_gnp.html)
      (default targets average degree about 2)

  `hub`

  :   `n_hubs`: partition nodes into that many groups and connect each
      group's hub to all other members (star modules)

  `scale_free`

  :   `power`, `edges_per_node` (`m` in
      [`igraph::sample_pa()`](https://r.igraph.org/reference/sample_pa.html))

  `stochastic_block_model`

  :   `block_prob`, `p_within`, `p_between`

  `small_world`

  :   `nei`, `p` (rewiring probability) for
      [`igraph::sample_smallworld()`](https://r.igraph.org/reference/sample_smallworld.html)

## Value

Symmetric integer matrix \\G\times G\\ with zero diagonal.

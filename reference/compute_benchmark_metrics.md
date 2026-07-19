# Compute summary metrics for estimated proportions

When a ground truth \\\boldsymbol{p}^{\star}\\ is supplied, scores
compare \\\hat{\boldsymbol{p}}\\ to \\\boldsymbol{p}^{\star}\\.
Otherwise scores compare the reconstituted bulk
\\\hat{\boldsymbol{y}}=\boldsymbol{\mu}\hat{\boldsymbol{p}}\\ to the
observed \\\boldsymbol{y}\\.

## Usage

``` r
compute_benchmark_metrics(
  y,
  mean_signature_matrix,
  estimated_p,
  true_ratios = NULL
)
```

## Arguments

- y:

  Bulk expression vector \\\boldsymbol{y}\in\mathbb{R}^{G}\\ (one
  heterogeneous sample).

- mean_signature_matrix:

  Mean signature \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (columns = cell types).

- estimated_p:

  Estimated proportions \\\hat{\boldsymbol{p}}\in\mathbb{R}^{J}\\.

- true_ratios:

  Optional ground-truth \\\boldsymbol{p}^{\star}\in\mathbb{R}^{J}\\ used
  only for benchmark scores.

## Value

A `tibble` with mse/rmse/mae, optionally \\R^{2}\\ / adjusted \\R^{2}\\,
and Pearson correlation.

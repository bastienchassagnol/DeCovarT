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
  (columns = cell types; plug-in for latent profiles).

- estimated_p:

  Estimated proportions \\\hat{\boldsymbol{p}}\in\mathbb{R}^{J}\\.

- true_ratios:

  Optional ground-truth proportions
  \\\boldsymbol{p}^{\star}\in\mathbb{R}^{J}\\. When supplied, metrics
  compare \\\hat{\boldsymbol{p}}\\ to \\\boldsymbol{p}^{\star}\\;
  otherwise they compare
  \\\hat{\boldsymbol{y}}=\boldsymbol{\mu}\hat{\boldsymbol{p}}\\ to
  \\\boldsymbol{y}\\.

## Value

A `tibble` with mse/rmse/mae, optionally \\R^{2}\\ / adjusted \\R^{2}\\,
and Pearson correlation.

## Examples

``` r
mu <- matrix(c(20, 22, 22, 20), nrow = 2,
             dimnames = list(paste0("g", 1:2), paste0("ct", 1:2)))
y <- drop(mu %*% c(0.4, 0.6))
compute_benchmark_metrics(y, mu, estimated_p = c(0.45, 0.55),
                          true_ratios = c(0.4, 0.6))
#> # A tibble: 1 × 6
#>   model_mse model_rmse model_mae model_coef_determination model_coef_determina…¹
#>       <dbl>      <dbl>     <dbl>                    <dbl>                  <dbl>
#> 1   0.00250     0.0500    0.0500                     0.75                   0.75
#> # ℹ abbreviated name: ¹​model_coef_determination_adjusted
#> # ℹ 1 more variable: model_cor <dbl>
```

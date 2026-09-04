# Compute deconvolution benchmark metrics

Returns three blocks: composition / regression scores, Monte Carlo
summaries (ADEMP-style), and optimisation / runtime diagnostics. When
\\\boldsymbol{p}^{\star}\\ is a matrix \\J\times N\\ (or
\\\hat{\boldsymbol{p}}\\ is), cell-type Pearson correlation and Monte
Carlo summaries are computed **across samples**, not within a single
composition.

## Usage

``` r
compute_benchmark_metrics(
  y,
  mean_signature_matrix,
  estimated_p,
  true_ratios = NULL,
  se = NULL,
  lower = NULL,
  upper = NULL,
  elapsed = NULL,
  memory_bytes = NULL,
  kkt_residual = NULL,
  numerical_converged = NULL,
  loglik_hat = NULL,
  loglik_true = NULL,
  presence_threshold = 1e-04,
  level = 0.95,
  coverage_interval = "wilson",
  algorithm = NA_character_
)
```

## Arguments

- y:

  Numeric vector (or one-column matrix)
  \\\boldsymbol{y}\in\mathbb{R}^{G}\\.

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (plug-in means).

- estimated_p:

  Estimated proportions: a length-\\J\\ vector or a numeric matrix
  \\J\times N\\.

- true_ratios:

  Optional ground truth with the same layout. When omitted,
  reconstitution scores compare
  \\\hat{\boldsymbol{y}}=\boldsymbol{\mu}\hat{\boldsymbol{p}}\\ with
  \\\boldsymbol{y}\\.

- se:

  Optional standard errors matching `estimated_p`.

- lower, upper:

  Optional interval bounds matching `estimated_p` (used for coverage and
  mean width). When omitted but `se` is supplied, Wald intervals at
  `level` are used.

- elapsed, memory_bytes, kkt_residual:

  Optional per-sample runtime diagnostics (recycled to \\N\\).

- numerical_converged:

  Optional logical per sample: the solver returned a finite simplex
  vector.

- loglik_hat, loglik_true:

  Optional per-sample log-likelihoods used for theoretical convergence
  (regret \\\ell(\boldsymbol{p}^{\star})-\ell(\hat{\boldsymbol{p}})\\).

- presence_threshold:

  Threshold \\\varepsilon\\ for presence / F1 / false-positive mass
  (default `1e-4`).

- level:

  Wald coverage level when `lower` / `upper` are omitted (default
  `0.95`).

- coverage_interval:

  Interval for the Monte Carlo coverage *rate* \\\hat\pi\\: `"wilson"`
  (default), `"wald"`, or `"agresti_coull"`. This is not the interval
  for \\p_j\\ itself.

- algorithm:

  Character label stored on each metrics row (typically the solver
  name). Recycled to length one.

## Value

A named list with:

- `regression`: `global` (one row per sample) and `cell_type` (one row
  per cell type);

- `monte_carlo`: one row per cell type;

- `optimisation`: one row per sample.

## Examples

``` r
mu <- matrix(c(20, 22, 22, 20), nrow = 2,
             dimnames = list(paste0("g", 1:2), paste0("ct", 1:2)))
y <- drop(mu %*% c(0.4, 0.6))
compute_benchmark_metrics(y, mu, estimated_p = c(0.45, 0.55),
                          true_ratios = c(0.4, 0.6))
#> $regression
#> $regression$global
#> # A tibble: 1 × 9
#>   sample_id algorithm     tv   rmse angular  sdid maxae reconstitution_mae
#>   <chr>     <chr>      <dbl>  <dbl>   <dbl> <dbl> <dbl>              <dbl>
#> 1 sample_1  NA        0.0500 0.0500  0.0622 0.312  0.05                 NA
#> # ℹ 1 more variable: reconstitution_cor <dbl>
#> 
#> $regression$cell_type
#> # A tibble: 2 × 5
#>   algorithm cell_type pearson presence_f1 false_positive_mass
#>   <chr>     <chr>       <dbl>       <dbl>               <dbl>
#> 1 NA        ct1            NA           1                   0
#> 2 NA        ct2            NA           1                   0
#> 
#> 
#> $monte_carlo
#> # A tibble: 2 × 14
#>   algorithm cell_type    bias empirical_sd mean_model_sd mean_model_se
#>   <chr>     <chr>       <dbl>        <dbl>         <dbl>         <dbl>
#> 1 NA        ct1        0.05             NA            NA            NA
#> 2 NA        ct2       -0.0500           NA            NA            NA
#> # ℹ 8 more variables: se_sd_ratio <dbl>, rmse <dbl>, coverage <dbl>,
#> #   coverage_lower <dbl>, coverage_upper <dbl>, coverage_interval <chr>,
#> #   mean_interval_width <dbl>, mcse_coverage <dbl>
#> 
#> $optimisation
#> # A tibble: 1 × 10
#>   sample_id algorithm elapsed_sec memory_bytes kkt_residual numerical_converged
#>   <chr>     <chr>           <dbl>        <dbl>        <dbl> <lgl>              
#> 1 sample_1  NA                 NA           NA           NA TRUE               
#> # ℹ 4 more variables: theoretical_converged <lgl>, loglik_regret <dbl>,
#> #   ct1 <dbl>, ct2 <dbl>
#> 
```

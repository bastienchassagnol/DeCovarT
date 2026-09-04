# Parallel deconvolution of a bulk expression matrix

For each column \\\boldsymbol{y}\_{\cdot i}\\ of the bulk matrix
\\\boldsymbol{Y}\in\mathcal{M}\_{G\times N}\\, estimates
\\\hat{\boldsymbol{p}}\_{\cdot i}\\ with every supplied algorithm.
Samples are iterated with `furrr` when `cores > 1`; algorithms are
sequential, so workers are never nested. When covariance information is
provided, DeCovarT methods maximise
\\\ell\_{\boldsymbol{y}\\\|\\\boldsymbol{\zeta}}(\boldsymbol{p})\\ under
\\\boldsymbol{y}\\\|\\(\boldsymbol{\zeta},\boldsymbol{p})\sim
\mathcal{N}\_{G}( \boldsymbol{\mu}\boldsymbol{p},
\boldsymbol{\Sigma}(\boldsymbol{p}) )\\.

## Usage

``` r
deconvolute_ratios(
  signature_matrix,
  bulk_expression,
  true_ratios = NULL,
  Sigma = NULL,
  deconvolution_functions = NULL,
  standardise = FALSE,
  scaled = FALSE,
  cores = getOption("mc.cores", 1L),
  coverage_interval = "wilson"
)
```

## Arguments

- signature_matrix:

  Mean signature \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (rownames = genes, colnames = cell types).

- bulk_expression:

  Bulk matrix \\\boldsymbol{Y}\in\mathcal{M}\_{G\times N}\\ (rownames =
  genes, colnames = samples).

- true_ratios:

  Optional ground-truth proportions: a length-\\J\\ numeric vector
  (recycled over samples) or a numeric matrix \\J\times N\\ (or
  \\N\times J\\).

- Sigma:

  Optional array \\(\boldsymbol{\Sigma}\_j)\_{j}\in\mathcal{M}\_{G\times
  G\times J}\\ of cell-type covariances (numeric; off-diagonals may be
  negative).

- deconvolution_functions:

  Named list; each element has `FUN` and optional
  `additional_parameters`.

- standardise:

  If `TRUE`, apply a **gene-wise** affine z-score computed once from
  \\\boldsymbol{\mu}\\ to bulk, means and covariances (see Details).
  Cell-type-wise or sample-wise transforms are not supported.

- scaled:

  Deprecated. `TRUE` (log2 mixing) always errors.

- cores:

  Number of `furrr` workers for the **sample** loop. Defaults to
  `getOption("mc.cores", 1L)`. Use `cores = 1` to stay sequential (CRAN
  examples and nested callers). With `cores > 1`, `.map_samples()` uses
  `furrr_options(seed = TRUE)`, which assigns an independent
  L'Ecuyer-CMRG stream per worker (R 4.1 / `future` streams). Do not
  offset a shared seed by the sample index.

- coverage_interval:

  Passed to
  [`compute_benchmark_metrics()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_benchmark_metrics.md)
  (`"wilson"`, `"wald"`, or `"agresti_coull"`).

## Value

A named list from
[`compute_benchmark_metrics()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_benchmark_metrics.md)
with `regression` (global and cell-type subtables), `monte_carlo`, and
`optimisation` (per-sample elapsed time, memory, KKT residual, and
\\\hat{\boldsymbol{p}}\\). First-generation solvers still call
[`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md);
the three native DeCovarT maps (ILR or \\p/\sum p\\) already lie on the
simplex.

## References

Erkkilä T, Lehmusvaara S, Ruusuvuori P, Visakorpi T, Shmulevich I,
Lähdesmäki H (2010). “Probabilistic Analysis of Gene Expression
Measurements from Heterogeneous Tissues.” *Bioinformatics*, **26**.
[doi:10.1093/bioinformatics/btq406](https://doi.org/10.1093/bioinformatics/btq406)
. Ahn J, Yuan Y, Parmigiani G, Suraokar MB, Diao L, Wistuba II, Wang W
(2013). “DeMix: Deconvolution for Mixed Cancer Transcriptomes Using Raw
Measured Data.” *Bioinformatics (Oxford, England)*, **29**.
[doi:10.1093/bioinformatics/btt301](https://doi.org/10.1093/bioinformatics/btt301)
. Wang Z, Cao S, Morris JS, Ahn J, Liu R, Tyekucheva S, Gao F, Li B, Lu
W, Tang X, Wistuba II, Bowden M, Mucci L, Loda M, Parmigiani G, Holmes
CC, Wang W (2018). “Transcriptome Deconvolution of Heterogeneous Tumor
Samples with Immune Infiltration.” *iScience*, **9**.
[doi:10.1016/j.isci.2018.10.028](https://doi.org/10.1016/j.isci.2018.10.028)
. Anghel CV, Quon G, Haider S, Nguyen F, Deshwar AG, Morris QD, Boutros
PC (2015). “ISOpureR: An R Implementation of a Computational
Purification Algorithm of Mixed Tumour Profiles.” *BMC Bioinformatics*,
**16**.
[doi:10.1186/s12859-015-0597-x](https://doi.org/10.1186/s12859-015-0597-x)
. Ogundijo OE, Wang X (2017). “A Sequential Monte Carlo Approach to Gene
Expression Deconvolution.” *PLOS ONE*, **12**.
[doi:10.1371/journal.pone.0186167](https://doi.org/10.1371/journal.pone.0186167)
. Hafemeister C, Satija R (2019). “Normalization and variance
stabilization of single-cell RNA-seq data using regularized negative
binomial regression.” *Genome Biology*, **20**(1), 296.
[doi:10.1186/s13059-019-1874-1](https://doi.org/10.1186/s13059-019-1874-1)
. Chion M, Leroy A (2023). “A Bayesian Framework for Multivariate
Differential Analysis Accounting for Missing Data.” arXiv.
[doi:10.48550/arxiv.2307.08975](https://doi.org/10.48550/arxiv.2307.08975)
. Newman A, Liu C, Green M, Gentles A, Feng W, Xu Y, Hoang C, Diehn M,
Alizadeh A (2015). “Robust Enumeration of Cell Subsets from Tissue
Expression Profiles.” *Nature methods*, **12**.
[doi:10.1038/nmeth.3337](https://doi.org/10.1038/nmeth.3337) .

## See also

[`deconvolute_ratios_Marquardt_Levenberg()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)

## Examples

``` r
set.seed(1)
genes <- paste0("g", 1:2)
cts <- paste0("ct", 1:2)
mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
Sigma <- array(
  c(1, 0, 0, 1, 1, 0, 0, 1),
  dim = c(2, 2, 2),
  dimnames = list(genes, genes, cts)
)
Y <- simulate_bulk_mixture(mu, Sigma, p = c(0.5, 0.5), n = 2)$Y
deconvolute_ratios(
  signature_matrix = mu,
  bulk_expression = Y,
  true_ratios = c(0.5, 0.5),
  Sigma = Sigma,
  deconvolution_functions = list(
    "nnls" = list(FUN = deconvolute_ratios_nnls)
  ),
  cores = 1
)
#> $regression
#> $regression$global
#> # A tibble: 2 × 9
#>   sample_id algorithm     tv   rmse angular  sdid  maxae reconstitution_mae
#>   <chr>     <chr>      <dbl>  <dbl>   <dbl> <dbl>  <dbl>              <dbl>
#> 1 sample_1  nnls      0.0806 0.0806   0.102 0.399 0.0806                 NA
#> 2 sample_2  nnls      0.220  0.220    0.264 0.634 0.220                  NA
#> # ℹ 1 more variable: reconstitution_cor <dbl>
#> 
#> $regression$cell_type
#> # A tibble: 2 × 5
#>   algorithm cell_type pearson presence_f1 false_positive_mass
#>   <chr>     <chr>       <dbl>       <dbl>               <dbl>
#> 1 nnls      ct1            NA           1                   0
#> 2 nnls      ct2            NA           1                   0
#> 
#> 
#> $monte_carlo
#> # A tibble: 2 × 14
#>   algorithm cell_type    bias empirical_sd mean_model_sd mean_model_se
#>   <chr>     <chr>       <dbl>        <dbl>         <dbl>         <dbl>
#> 1 nnls      ct1        0.0696        0.212            NA            NA
#> 2 nnls      ct2       -0.0696        0.212            NA            NA
#> # ℹ 8 more variables: se_sd_ratio <dbl>, rmse <dbl>, coverage <dbl>,
#> #   coverage_lower <dbl>, coverage_upper <dbl>, coverage_interval <chr>,
#> #   mean_interval_width <dbl>, mcse_coverage <dbl>
#> 
#> $optimisation
#> # A tibble: 2 × 10
#>   sample_id algorithm elapsed_sec memory_bytes kkt_residual numerical_converged
#>   <chr>     <chr>           <dbl>        <dbl>        <dbl> <lgl>              
#> 1 sample_1  nnls                0    537381888        0.444 TRUE               
#> 2 sample_2  nnls                0    541645824        0.152 TRUE               
#> # ℹ 4 more variables: theoretical_converged <lgl>, loglik_regret <dbl>,
#> #   ct1 <dbl>, ct2 <dbl>
#> 
```

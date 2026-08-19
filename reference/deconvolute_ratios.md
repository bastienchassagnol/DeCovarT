# Parallel deconvolution of a bulk expression matrix

For each column \\\boldsymbol{y}\_{\cdot i}\\ of the bulk matrix
\\\boldsymbol{Y}\in\mathcal{M}\_{G\times N}\\, estimates
\\\hat{\boldsymbol{p}}\_{\cdot i}\\ with every supplied algorithm. When
covariance information is provided, DeCovarT methods maximise
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
  scaled = FALSE,
  cores = ifelse(.Platform$OS.type == "unix", getOption("mc.cores",
    parallel::detectCores()), 1)
)
```

## Arguments

- signature_matrix:

  Mean signature \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (rownames = genes, colnames = cell types). Numeric matrix of
  non-negative expression; \\G \ge J\\ is required (the mixture is
  undetermined if \\J \> G\\). Used as a frequentist **plug-in** for
  unobserved latent cell-type profiles \\\boldsymbol{x}\_{\cdot j}\\.

- bulk_expression:

  Bulk matrix \\\boldsymbol{Y}\in\mathcal{M}\_{G\times N}\\: numeric,
  non-negative, gene rownames matching `signature_matrix`.

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

- scaled:

  If `TRUE`, apply a log2 transform before estimation.

- cores:

  Number of parallel workers.

## Value

A `tibble` of estimated \\\hat{\boldsymbol{p}}\\ and metrics, after
[`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md).

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
.

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
#> # A tibble: 2 × 10
#>   model_mse model_rmse model_mae model_coef_determination model_coef_determina…¹
#>       <dbl>      <dbl>     <dbl>                    <dbl>                  <dbl>
#> 1   0.00650     0.0806    0.0806                        0                      0
#> 2   0.0483      0.220     0.220                         0                      0
#> # ℹ abbreviated name: ¹​model_coef_determination_adjusted
#> # ℹ 5 more variables: model_cor <dbl>, ct1 <dbl>, ct2 <dbl>, OMIC_ID <chr>,
#> #   algorithm <chr>
```

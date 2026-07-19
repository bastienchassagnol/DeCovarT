# Parallel deconvolution of a bulk expression matrix

For each column \\\boldsymbol{y}\_{\cdot i}\\ of the bulk matrix
\\\boldsymbol{Y}\in\mathcal{M}\_{G\times N}\\, estimates
\\\hat{\boldsymbol{p}}\_{\cdot i}\\ with every supplied algorithm. When
covariance information is provided, DeCovarT methods maximise
\\\ell\_{\boldsymbol{y}\\\|\\\boldsymbol{\zeta}}(\boldsymbol{p})\\ under
\\\boldsymbol{y}\\\|\\(\boldsymbol{\zeta},\boldsymbol{p})\sim
\mathcal{N}\_{G}(\boldsymbol{\mu}\boldsymbol{p},\boldsymbol{\Sigma}(\boldsymbol{p}))\\.

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
  (rownames = genes, colnames = cell types).

- bulk_expression:

  Bulk matrix \\\boldsymbol{Y}\in\mathcal{M}\_{G\times N}\\.

- true_ratios:

  Optional ground-truth proportions for scoring.

- Sigma:

  Optional array \\(\boldsymbol{\Sigma}\_j)\_{j}\in\mathcal{M}\_{G\times
  G\times J}\\.

- deconvolution_functions:

  Named list; each element has `FUN` and optional
  `additional_parameters`.

- scaled:

  If `TRUE`, apply a log2 transform before estimation.

- cores:

  Number of parallel workers.

## Value

A `tibble` of estimated \\\hat{\boldsymbol{p}}\\ and metrics, after
[`enforce_identifiability()`](https://bastienchassagnol.github.io/DeCovarT/reference/enforce_identifiability.md).

## See also

[`deconvolute_ratios_Marquardt_Levenberg()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)

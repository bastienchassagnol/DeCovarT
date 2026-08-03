# Score mean profiles with AutoGeneS-style objectives

Computes the two objectives used by AutoGeneS for signature quality:
mean absolute pairwise cosine similarity between cell-type columns of
\\\boldsymbol{\mu}\\ (to minimise) and the sum of pairwise Euclidean
distances between those columns (to maximise).

## Usage

``` r
compute_mean_profile_objectives(mean_signature_matrix)
```

## Arguments

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (columns = cell types).

## Value

A named list with `mean_abs_cosine` and `sum_euclidean_distance`.

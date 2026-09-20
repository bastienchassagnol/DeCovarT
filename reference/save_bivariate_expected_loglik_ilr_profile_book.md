# Multi-page PDF of expected ILR log-likelihood profiles (\\J=2\\)

Multi-page PDF of expected ILR log-likelihood profiles (\\J=2\\)

## Usage

``` r
save_bivariate_expected_loglik_ilr_profile_book(
  config,
  theta_tbl,
  file,
  data_rds = NULL
)
```

## Arguments

- config:

  Slim config tibble with `ID`.

- theta_tbl:

  Tibble with `ID` and `true_theta`.

- file:

  Output PDF path.

- data_rds:

  Optional directory; when set, writes `purified_density.rds` (ggplot
  raster and overlay tables).

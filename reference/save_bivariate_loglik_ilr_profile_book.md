# Multi-page PDF of ILR log-likelihood profiles (\\\rho\in\mathbb{R}\\)

Multi-page PDF of ILR log-likelihood profiles (\\\rho\in\mathbb{R}\\)

## Usage

``` r
save_bivariate_loglik_ilr_profile_book(
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

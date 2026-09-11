# Multi-page PDF of bulk convolution densities

Multi-page PDF of bulk convolution densities

## Usage

``` r
save_bivariate_bulk_density_book(
  config,
  theta_tbl,
  file,
  n = 800L,
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

- n:

  Draws per cell type.

- data_rds:

  Optional directory; when set, writes `purified_density.rds` (ggplot
  raster and overlay tables).

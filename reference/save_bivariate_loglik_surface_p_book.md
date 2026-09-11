# Multi-page PDF of bulk log-likelihood surfaces on \\(p_1,p_2)\\

Multi-page PDF of bulk log-likelihood surfaces on \\(p_1,p_2)\\

## Usage

``` r
save_bivariate_loglik_surface_p_book(config, theta_tbl, file, data_rds = NULL)
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

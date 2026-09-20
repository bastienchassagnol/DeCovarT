# Multi-page PDF of expected log-likelihood surfaces

Same layout as
[`save_bivariate_loglik_surface_p_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_bivariate_loglik_surface_p_book.md),
but each panel is \\\mathbb{E}\_{Y\mid p^{\star}}\[\ell(p;Y)\]\\ rather
than one realised bulk column \\y=\mu p^{\star}\\.

## Usage

``` r
save_bivariate_expected_loglik_surface_p_book(
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

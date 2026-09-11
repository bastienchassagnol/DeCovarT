# Multi-page PDF (and optional HTML) of rgl bulk log-likelihood surfaces

Requires `rgl` and `png` for PDF snapshots. When `html_file` is set,
also writes an interactive
[`rgl::rglwidget()`](https://dmurdoch.github.io/rgl/dev/reference/rglwidget.html)
HTML (needs `htmlwidgets`).

## Usage

``` r
save_bivariate_loglik_rgl_book(config, theta_tbl, file, html_file = NULL)
```

## Arguments

- config:

  Slim config tibble with `ID`.

- theta_tbl:

  Tibble with `ID` and `true_theta`.

- file:

  Output PDF path.

- html_file:

  Optional path for a self-contained interactive HTML.

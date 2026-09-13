# Multi-page PDF of clustered algorithm-similarity heatmaps

Twelve pages (CLD \\\times\\ variance \\\times\\ composition). Each page
is a 2-by-2 of the correlation corners, with Ward D2 clustering of
\\1-r\\ and a dendrogram to the right of the tiles, with leaves flush
against the heatmap. One shared colourbar per page.

## Usage

``` r
save_bivariate_similarity_book(artefacts, file, data_rds = NULL)
```

## Arguments

- artefacts:

  List from
  [`read_simulation_artefacts()`](https://bastienchassagnol.github.io/DeCovarT/reference/read_simulation_artefacts.md)
  with `assemble = TRUE`, or the same named pieces.

- file:

  Output PDF path.

- data_rds:

  Optional directory for ggplot `data` RDS files.

## Value

`file`, invisibly.

# Drop duplicated design columns from a scenario-tagged table

Aliases (`scenario_idx`, `rho_ct1` / `rho_ct2`, `proportion_name`,
`centroid`) are removed after copying them onto the canonical names
`proportions` and `centroids` when those are missing.

## Usage

``` r
slim_scenario_table(tbl)
```

## Arguments

- tbl:

  A tibble that may contain redundant aliases.

## Value

The same table without alias columns.

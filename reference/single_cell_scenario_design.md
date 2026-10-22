# Four single-cell simulation scenarios

The scenario grid used to check the single-cell fit. Two compositions,
equal and imbalanced, crossed with equal and unequal cell sizes. Every
scenario is drawn at \\G=10\\ genes, \\n_i=10^{4}\\ cells and platform
scales \\a_i\in\\1, 2\\\\ (same capture efficiency, or twice the bulk
capture of the single-cell platform).

## Usage

``` r
single_cell_scenario_design()
```

## Value

A list with `n_genes`, `n_cells`, `platform_scale` and `scenarios` (four
named lists with `proportions` and `cell_size`).

## See also

[`rescale_single_cell_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/rescale_single_cell_moments.md),
[`simulate_bulk_from_cells()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_from_cells.md)

## Examples

``` r
names(single_cell_scenario_design()$scenarios)
#> [1] "balanced_equal_size"     "balanced_varying_size"  
#> [3] "imbalanced_equal_size"   "imbalanced_varying_size"
```

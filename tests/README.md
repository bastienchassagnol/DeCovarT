# Tests and fixtures

Unit tests live in `tests/testthat/` and mirror `R/` stems
(`test-03_03_DeCovarT.R` for `R/03_03_DeCovarT_estimate_ratios_frequentist.R`).

## Bundled example data (G5.1)

| Object | Where | Used by |
|---|---|---|
| `toy_deconvolution` | `inst/extdata/toy_deconvolution.rds` (copy under `tests/testthat/fixtures/`) | Input-check and `deconvolute_ratios()` smoke tests |
| `bivariate_configuration.rds`, `bivariate_estimation.rds` | `tests/testthat/fixtures/` | Golden-file comparison in `test-03_03_DeCovarT.R` |

Regenerate the bivariate fixtures from the package root after `devtools::load_all()`:

```r
source("tests/testthat/fixtures/make-useful-things.R")
```

## Timing (current suite)

All regular files finish in well under one second on a desktop (slowest:
`test-02_02_generate_synthetic_networks.R`, hierarchical GRN with 30 genes).
The bivariate golden-file test re-runs every solver, including simulated
annealing, and is the heaviest single test.

## Extended tests (G5.10)

Set `DECOVART_EXTENDED_TESTS=true` to run tests wrapped in
`skip_if_not_extended()` (see `tests/testthat/helper-srr.R`). Paper-scale
simulations and any future large downloads belong in the companion
reproducibility repository (G1.5 / G5.11), not in this package.

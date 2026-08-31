# Tests and fixtures

Unit tests live in `tests/testthat/` and mirror `R/` stems
(`test-03_03_DeCovarT.R` for `R/03_03_DeCovarT_estimate_ratios_frequentist.R`).

## Bundled example data (G5.1)

| Object | Where | Used by |
|---|---|---|
| `toy_deconvolution` | `inst/extdata/toy_deconvolution.rds` (copy under `tests/testthat/fixtures/`) | Input-check and `deconvolute_ratios()` smoke tests |

The bivariate toy scenario used by `test-03_03_DeCovarT.R` is built in
memory by `new_bivariate_toy_scenario()` in `tests/testthat/helper.R`
(r-pkgs “create useful_things with a helper”). Construction is cheap, so
nothing is written to the package tree, `getwd()`, or the home filespace.
On the CRAN tarball `scripts/` is excluded; that test then skips.

If a test must touch disk, wrap the write in `withr::local_tempfile()` /
`withr::local_tempdir()` (cleaned when the test exits). Do not leave
files under `tempdir()`.

## Timing (current suite)

All regular files finish in well under one second on a desktop (slowest:
`test-02_02_generate_synthetic_networks.R`, hierarchical GRN with 30 genes).
The bivariate helper re-runs every solver, including simulated
annealing, and is the heaviest single test.

## Extended tests (G5.10)

Set `DECOVART_EXTENDED_TESTS=true` to run tests wrapped in
`skip_if_not_extended()` (see `tests/testthat/helper-srr.R`). Paper-scale
simulations and any future large downloads belong in the companion
reproducibility repository (G1.5 / G5.11), not in this package.

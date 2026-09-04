## Resubmission

This is a resubmission of 2.3.0. The remaining human-review item was
writing under the package directory from
`tests/testthat/fixtures/make-useful-things.R`.

* That script (and the golden RDS files it wrote) has been removed.
  Tests now build the bivariate toy scenario in memory via
  `new_bivariate_toy_scenario()` in `tests/testthat/helper.R`. No
  default output path; nothing is written to the home filespace, the
  package tree, or `getwd()`. Where tests still need a file they use
  `withr::local_tempfile()` / `with_tempfile()`, which delete the
  artefact when the test ends.
* CRAN vignettes (`fig02-bivariate-toy`,
  `theory-decovart-generative-model`) no longer execute chunks that call
  `library(DeCovarT)` / `DeCovarT::` during `R CMD build`. Quarto's
  Windows CLI starts a new R process that cannot see the temporary
  library used to build vignettes (quarto-dev/quarto-r#217). Live
  checks remain in testthat. R-hub Windows pre-installs the package
  (`local::.`), matching GitHub Actions.

Previous 2.3.0 / 2.2.3 tarball notes (still apply): optional heatmap
Suggests (`ComplexHeatmap`, `circlize`, `viridis`); no
`Additional_repositories`; Quarto Markdown tables in CRAN vignettes;
absolute pkgdown URLs for excluded articles; Rd figure widths in pixels;
DESCRIPTION Hunspell wording unchanged (British English / technical
terms).

## Test environments

* Local: Ubuntu 24.04, R 4.5.2 (`devtools::check()` / `rcmdcheck --as-cran`)
* GitHub Actions: Ubuntu / Windows (R-release)
* R-hub: linux, windows (no macOS)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission to CRAN (version 2.3.1). The GitHub series
  already used 2.0.x / 2.2.x / 2.3.x; there is no previous CRAN release.
* Maintainer email: bastien_chassagnol@laposte.net
* Incoming feasibility may report "Possibly misspelled words" in
  DESCRIPTION. These are intentional:
  - British English: optimisers, reparametrisation, methodology
  - Technical terms: Deconvolution, Hessians
  - Mid-word splits reported by Hunspell (tran/ript, wi/de, odology,
    ribed) are false positives on unbroken words.
* Method reference: Chassagnol, Nuel and Becht (2023)
  <doi:10.48550/arXiv.2309.09557>.
* The CRAN tarball ships two Quarto vignettes
  (`fig02-bivariate-toy`, `theory-decovart-generative-model`). Longer articles
  remain on the pkgdown site and are linked with absolute URLs.
* Heatmap packages `ComplexHeatmap`, `circlize`, and `viridis` are in
  Suggests only (optional `plot_correlation_Heatmap()`).

## Reverse dependencies

There are currently no reverse dependencies on CRAN / Bioconductor.

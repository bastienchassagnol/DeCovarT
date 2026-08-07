## Test environments

* Local: Ubuntu 24.04, R 4.5.2
* GitHub Actions: Ubuntu / Windows / macOS (release + devel as configured)
* R-hub: linux (R-devel), windows (R-devel)

## R CMD check results

0 errors | 0 warnings | 0 notes (local `--as-cran`)

* This is a new submission to CRAN (version 2.0.0, matching the GitHub
  release tag `v2.0.0`).
* Method reference: Chassagnol, Nuel and Becht (2023)
  <doi:10.48550/arXiv.2309.09557> (also on HAL).
* The CRAN tarball ships two Quarto vignettes
  (`DeCoVart-use-cases`, `softmax-alr-derivatives`); additional articles
  remain on the pkgdown site only.
* Local check may report an installed-size INFO (~9 Mb after slimming
  vignettes/logos); previously ~61 Mb.

## Reverse dependencies

There are currently no reverse dependencies on CRAN / Bioconductor
(`revdep/` scaffolding in place; `revdep/cran.md` records 0 packages).

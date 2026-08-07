## Test environments

* Local: Ubuntu 24.04, R 4.5.2
* GitHub Actions: Ubuntu / Windows / macOS (release + devel as configured)
* R-hub: linux (R-devel), windows (R-devel)

## R CMD check results

0 errors | 0 warnings | 0 notes (local `--as-cran`)

* This is a new submission to CRAN (not yet on CRAN).
* Method reference: Chassagnol, Nuel and Becht (2023)
  <doi:10.48550/arXiv.2309.09557> (also on HAL).
* Local check reports an INFO (not a NOTE) that the installed size is
  about 61 Mb (`doc` ~52 Mb from vignette HTML/figures; `logo` ~3.6 Mb;
  `help` ~4.8 Mb). Happy to compress logos / slim vignette assets if
  requested.

## Reverse dependencies

There are currently no reverse dependencies on CRAN / Bioconductor
(`revdep/` scaffolding in place; `revdep/cran.md` records 0 packages).

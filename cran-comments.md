## Test environments

* Local: Ubuntu 24.04, R 4.5.2 (`rcmdcheck --as-cran`)
* GitHub Actions: Ubuntu / Windows (R-release)
* R-hub: linux, windows (no macOS)
* win-builder: devel and release (submitted)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission to CRAN (version 2.2.0). The GitHub series
  already used 2.0.x / 2.2.0; there is no previous CRAN release.
* Maintainer email: bastien_chassagnol@laposte.net
* Incoming feasibility also reports:
  Availability using Additional_repositories specification:
    ? ? https://bioconductor.org/packages/release/bioc
  Imports `ComplexHeatmap` from Bioconductor (heatmaps only;
  `plot_correlation_Heatmap()`). Estimators do not call it. `circlize`
  is on CRAN.
* Method reference: Chassagnol, Nuel and Becht (2023)
  <doi:10.48550/arXiv.2309.09557>.
* The CRAN tarball ships two Quarto vignettes
  (`DeCoVart-use-cases`, `softmax-alr-derivatives`); additional articles
  remain on the pkgdown site.

## Reverse dependencies

There are currently no reverse dependencies on CRAN / Bioconductor.

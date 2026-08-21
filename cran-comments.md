## Resubmission

This is a resubmission of 2.2.2. Changes addressing the remaining
incoming NOTE and dependency layout:

* Removed `Additional_repositories` pointing at Bioconductor. That field
  is for packages outside CRAN / Bioconductor / Omegahat / R-Forge;
  listing the official Bioconductor software repo produced the
  availability `? ?` NOTE without helping installability.
* Moved `ComplexHeatmap`, `circlize`, and `viridis` from Imports to
  Suggests. Only `plot_correlation_Heatmap()` needs them; estimators
  do not. The helper checks `requireNamespace()` and tells users how
  to install (CRAN packages via `install.packages()`, `ComplexHeatmap`
  via `BiocManager::install("ComplexHeatmap")`). Examples and tests are
  conditional on those Suggests.

Previous 2.2.2 fixes (still in this tarball): Quarto Markdown tables in
CRAN vignettes (no Windows `tinytable` S4 failure); absolute pkgdown
URLs for excluded articles; Rd figure widths in pixels; DESCRIPTION
Hunspell wording unchanged (British English / technical terms).

## Test environments

* Local: Ubuntu 24.04, R 4.5.2 (`rcmdcheck --as-cran`)
* GitHub Actions: Ubuntu / Windows (R-release)
* R-hub: linux, windows (no macOS)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission to CRAN (version 2.2.3). The GitHub series
  already used 2.0.x / 2.2.x; there is no previous CRAN release.
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
  (`DeCoVart-use-cases`, `softmax-alr-derivatives`). Longer articles
  remain on the pkgdown site and are linked with absolute URLs.
* Heatmap packages `ComplexHeatmap`, `circlize`, and `viridis` are in
  Suggests only (optional `plot_correlation_Heatmap()`).

## Reverse dependencies

There are currently no reverse dependencies on CRAN / Bioconductor.

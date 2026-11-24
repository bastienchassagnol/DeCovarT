---
description: R package layout, CRAN incoming, roxygen, testthat, GHA deps, pkgdown
globs: "**/DESCRIPTION,**/NAMESPACE,**/man/**,**/tests/testthat/**,**/vignettes/**,**/.github/workflows/*pkgdown*,**/.github/workflows/*R-CMD*,**/.github/workflows/*coverage*,**/_pkgdown.yml,**/pkgdown/**"
alwaysApply: false
paths:
  - "**/DESCRIPTION"
  - "**/NAMESPACE"
  - "**/man/**"
  - "**/tests/testthat/**"
  - "**/vignettes/**"
  - "**/.github/workflows/*pkgdown*"
  - "**/.github/workflows/*R-CMD*"
  - "**/.github/workflows/*coverage*"
  - "**/_pkgdown.yml"
  - "**/pkgdown/**"
---
# R package development

Applies when the repo is an R package (`DESCRIPTION`); `R/` code also follows
`r-development`.

Follow the repository's own checks first. Encode
[CRAN Cookbook code issues](https://contributor.r-project.org/cran-cookbook/code_issues.html)
in `R/`, examples, vignettes, and tests. Prefer **`withr`** (in Suggests when
needed) over bare `on.exit()` / `unlink()` for temporary changes.

## Package workflow

Standard layout: `DESCRIPTION`, roxygen-driven `NAMESPACE`, `R/`, `man/`,
`tests/testthat/`, `vignettes/`, `inst/`.

After API or doc changes:

```r
devtools::document()
devtools::load_all()
devtools::test()
devtools::check()
```

Use `pkgcheck::pkgcheck()` and `rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"))`
with `options(repos = c(CRAN = "https://cloud.r-project.org"))` before CRAN
submission—not only local `devtools::check()` or GHA green alone.

Also trigger **R-hub** / win-builder (Windows and Linux as your policy requires)
and read **incoming pretest** logs, not only GitHub Actions.

### CRAN incoming (beyond local check)

- **Relative HTML links** between vignettes (`foo.html#sec-…`) break when a target
  vignette is `.Rbuildignore`d. From any CRAN-shipped vignette, link excluded
  articles with **absolute pkgdown URLs**
  (`https://<user>.github.io/<pkg>/articles/…`).
- DESCRIPTION Hunspell flags: explain British English or domain terms in
  `cran-comments.md`; keep jargon in `inst/WORDLIST` for `spelling` / aspell.
- **Do not set `Additional_repositories` for CRAN or Bioconductor.** CRAN already
  treats Bioconductor software as a standard source; redundant
  `bioconductor.org/packages/release/bioc` triggers incoming availability `? ?`.
  Use `Additional_repositories` only for packages outside CRAN / Bioconductor /
  Omegahat / R-Forge.

### CRAN Cookbook (code)

- **`TRUE` / `FALSE`:** never `T` / `F` as logicals or object names.
- **Seeds:** no hard-coded `set.seed(n)` in package functions; take `seed = NULL`
  and call `set.seed()` only when the caller supplies a value. Tests, examples,
  and vignettes *should* use `withr::with_seed()` for reproducibility.
- **Console output:** no unsuppressable `print()` / `cat()` / `writeLines()` in
  package functions; use `message()` / `warning()` / `stop()`, or gate with
  `verbose`. (`print` / `summary` methods and interactive functions are exempt.)
- **`options()`, `par()`, `setwd()`:** restore state with `withr` or `on.exit`;
  never `options(warn = -1)`—use `suppressWarnings()` on one expression.
- **Writing files:** no default writes to home, package dir, or `getwd()`.
  Examples, vignettes, tests: `withr::local_tempfile()` / `local_tempdir()` /
  `with_tempfile()` / `with_tempdir()` only—not bare `tempfile()` + forgotten
  `unlink()`. Test helpers in `tests/testthat/helper.R` as in-memory
  constructors; no `saveRDS()` into `tests/testthat/fixtures/` during check.
  Committed read-only fixtures under `inst/extdata/` are fine.
- **`.GlobalEnv`:** no `<<-` assignment except documented Shiny exceptions; no
  `rm(list = ls())` in examples or vignettes.
- **`installed.packages()`:** never; use `requireNamespace()` / `require()`.
- **Installing software:** no `install.packages()` / `BiocManager::install()` in
  examples, tests, vignettes, or ordinary functions exercised during check.
- **Cores:** at most **2** in examples, vignettes, tests, and installation code;
  pass `cores = 1L` or `2L` explicitly.
- **Quarto vignettes on Windows:** `quarto::html` starts a fresh R without the
  check-time temp library. CRAN-shipped vignettes must not `library(pkg)` live
  unless `eval: false`. Put live checks in testthat; CI may pre-install the
  package (`local::.` / source install) for optional evaluated chunks.

### Formatting and pre-commit

- **Air** (line width 80, 2 spaces) when configured.
- **pre-commit** when `.pre-commit-config.yaml` exists; no `--no-verify` unless
  the user asks.
- Honour project `.lintr`; do not add global **lintr** by default.

## Roxygen

- Substantive roxygen (title, `@param`, `@return`) → **`@export`** (may add
  `@keywords internal` for advanced helpers).
- **`@noRd` only** when there is no substantive documentation (e.g. `.onLoad`).
- Never `@noRd` to hide docs you already wrote.
- **`@inheritParams`** for shared parameter blocks—do not copy-paste.

## Testing and validation

[R Packages — Designing your test suite](https://r-pkgs.org/testing-design.html):

- **testthat** edition 3 (`Config/testthat/edition: 3`).
- **`usethis::use_test("file_stem")`:** `R/foo.R` → `tests/testthat/test-foo.R`.
- Self-contained files: mostly `test_that()`; minimal top-level scope.
- Behaviour, edges, missing values, regressions—not every trivial line.
- Smoke or correctness test per exported function unless explicitly deferred.
- Style (goodpractice / lintr):
  - `expect_identical()` when exact match matters; `expect_equal()` for tolerance.
  - `expect_named`, `expect_length`, `expect_type` over verbose `expect_true(is.*)`.
  - Split `expect_true(a && b)` into separate expectations.
- Stochastic: `withr::with_seed()`; in-memory `helper.R` constructors over golden
  RDS when cheap; `local_*()` + `withr::defer()` for side effects.
- Interactive: `devtools::test()`, `test_coverage()`, `test_coverage_active_file()`.
- Checks: `devtools::check()` **and** `pkgcheck::pkgcheck()`; optional
  `pkgcheck::checks_to_markdown()`.
- CI: keep `R-CMD-check` and `test-coverage` meaningful (`needs: check|coverage`).

Record `sessioninfo` for scientific outputs when the repo does so.

## GitHub Actions: system libraries

Do **not** add a custom `ci_sysreqs.R` or bash `apt install` list as the default.

Use [`r-lib/actions/setup-r-dependencies@v2`](https://github.com/r-lib/actions/tree/v2/setup-r-dependencies)
**after** `r-lib/actions/setup-r@v2`. It installs `pak`, caches the library, and on
Ubuntu and other supported Linux distros installs system libraries via `pak` /
`rstudio/r-system-requirements`.

Typical pattern for check, coverage, and pkgdown workflows:

```yaml
- uses: r-lib/actions/setup-r-dependencies@v2
  with:
    extra-packages: any::rcmdcheck, local::.
    needs: check|coverage|website
```

Treat `pak` sysreqs as the expected set—not a second parallel apt installer.

## pkgdown

- Initialise GitHub Pages once with `quarto publish gh-pages` when appropriate;
  see [Quarto publish action](https://quarto.org/docs/publishing/github-pages.html#publish-action).
- Prefer **external includes** under `pkgdown/` over raw HTML in `_pkgdown.yml`
  (`template.includes.after_body`, pkgdown `{{#site}}{{root}}{{/site}}` URLs).
- Static assets in `pkgdown/assets/`.
- Prefer **Sass** (`pkgdown/extra.scss`) over ad-hoc CSS
  ([pkgdown customise](https://pkgdown.r-lib.org/articles/customise.html)).
# R and statistical computing

When R is primary, follow the repository first. Defaults below apply when the repo
is silent.

## Package workflow

Standard layout: `DESCRIPTION`, roxygen-driven `NAMESPACE`, `R/`, `man/`,
`tests/testthat/`, `vignettes/`, `inst/`.

Prefer `usethis`, `devtools`, `roxygen2`, `testthat`. Document and export the
public API deliberately; `@inheritParams` for shared parameters; `@noRd` only for
undocumented internals.

After API/doc changes:

```r
devtools::document()
devtools::load_all()
devtools::test()
devtools::check()
```

Run `pkgcheck::pkgcheck()` locally and enforce it in CI with
[pkgcheck-action](https://github.com/ropensci-review-tools/pkgcheck-action)
(`.github/workflows/pkgcheck.yaml`, from `pkgcheck::use_github_action_pkgcheck()`).
Add `--as-cran` checks when the package targets CRAN. pkgcheck does not see
directories in `.Rbuildignore` (such as `scripts/`); guard those with a
pre-commit `pygrep` hook on `pkg:::` and `asNamespace("pkg")`.

### CRAN-quality hygiene (summary)

- `Suggests` + `requireNamespace()` for optional deps; no `install.packages()` in
  examples/tests/vignettes.
- No writes to home, package dir, or `getwd()` by default; `withr` temp files.
- Restore `options()`, `par()`, wd; no `<<-`, no `.GlobalEnv` mutation, no
  package-level mutable environments (caches, flags): compute statelessly or
  pass values explicitly (`withr::with_seed()`, `withr::with_output_sink()`,
  `cli`/`rlang` `.frequency = "once"`).
- Never reach a package's internals with `:::`, `get()` or `asNamespace()`,
  from scripts, vignettes or another package. A function needed outside the
  package is documented and exported (`@export`, with `@examples`). Scripts
  run against the installed package (`library()`), where internals are not
  visible, unlike `devtools::load_all()`
  ([load_all](https://r-pkgs.org/workflow101.html#benefits-of-load_all),
  [R CMD check: R code](https://r-pkgs.org/R-CMD-check.html#r-code)).
- `TRUE`/`FALSE`; seeds only when caller supplies; ≤2 cores in check-time code.
- ASCII-only `R/` strings when required for portability.

### Formatting and pre-commit

- **Air** for R formatting (line width 80, 2 spaces) when configured.
- **pre-commit** when the repo has `.pre-commit-config.yaml`:
  `pre-commit clean` → stage all → `pre-commit run` → re-stage autofixes → commit.
  No `--no-verify` unless the user asks.

Do not add `lintr` to global defaults; honour project `.lintr` if present.

## R coding style

`snake_case`; `<-`; native `|>`; `pkg::fun()` in package code; explicit missing
values; `file.path()`/`fs`; `fixed = TRUE` for literal `grep`; cohesive functions.

**dplyr:** `{{ }}`, `.data[[var]]`, `all_of()`/`any_of()`, `:=` for dynamic names.

**Interfaces:** validate at boundary; forward `...` only with clear semantics;
actionable `stop()`/`warning()` without `paste()` wrappers.

## Performance, testing, reproducibility

Profile first; vectorise; `data.table` or DuckDB when scale demands it.

`testthat` edition 3: behaviour, edges, regressions, controlled RNG (`withr`).

Record seeds and `sessioninfo` for scientific outputs; `renv` when the repo uses
it; `targets` for heavy pipelines when justified.

## Visualisation

Follow the data-visualisation rule (`~/.claude/rules/data-visualisation.md`;
Cursor: `~/.cursor/rules/data-visualisation.mdc`): house theme `theme_slide()`,
colour roles, export, and the caveats checklist. When the repository defines a
theme, use it.

R packages: `ggplot2` for every static figure; `ragg` for PNG export; `scales`
(axis labels, rounding); `forcats` (ordering); `ggrepel` (direct labels);
`ggdist` (distributions, rainclouds); `ggiraph` (interactive, linked hover);
`ComplexHeatmap` only for richly annotated heatmaps; `grid` for multi-panel
layout.

## Parallelism and Shiny

`furrr` for reproducible map-style parallel work (independent RNG streams);
`mirai` for low-latency async when appropriate; no nested parallelism; configurable
workers.

Shiny: modules; no globals; testable domain logic.

## Preferred tools (when justified by repo)

`dplyr`/`tidyr`, `purrr`, `fs`/`here`, `ggplot2`, `tinytable` for static HTML
tables, `checkmate`, `cli`, `glue`/`stringr`, `jsonlite`, `devtools` stack.

Resampling/tuning: **`mlr3`** ecosystem (see statistical-modelling rule).

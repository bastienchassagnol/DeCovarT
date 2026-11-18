---
description: R packages, CRAN hygiene, style, testing, ggplot2 (see data-visualisation rule), parallelism
globs: "**/*.{R,r},**/DESCRIPTION,**/NAMESPACE,**/tests/**/*,**/R/**/*"
alwaysApply: false
paths:
  - "**/*.{R,r}"
  - "**/DESCRIPTION"
  - "**/NAMESPACE"
  - "**/tests/**/*"
  - "**/R/**/*"
---
# R and statistical computing

## Local R layout

This laptop runs Ubuntu 26.04 (Resolute). The system interpreter is the official CRAN apt build, not the Ubuntu archive package and not a copy under `/usr/local`.

- Repository: `resolute-cran40` from `https://cloud.r-project.org/bin/linux/ubuntu`, signed with the Michael Rutter key (`E298A3A825C0D65DFD57CBB651716619E084DAB9`).
- Binary: `/usr/bin/R` and `/usr/bin/Rscript` (R 4.6). Cursor `r.rpath.linux` stays on that binary so the language server uses GNU R.

`rig` is installed in user mode (`~/.local/bin/rig`) for extra R versions. It does not replace `/usr/bin/R`. Use `rig add` and `R-x.y.z` when an `renv.lock` records another R version. Do not run `rig default` over the apt R.

`renv` freezes the package list. Switch to the intended R version before `renv::init()` or `renv::snapshot()`. After a minor-version upgrade, do not reuse the previous user library (for example `~/R/x86_64-pc-linux-gnu-library/4.5`).

The local Cursor terminal is radian at `~/.local/bin/radian`, installed with `uv tool` for this user. Remote SSH sessions keep their own R terminal. On the system R, the editor packages are `languageserver`, `httpgd`, and `jsonlite`.

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

Use `pkgcheck::pkgcheck()` and `--as-cran` checks when the package targets CRAN.

### CRAN-quality hygiene (summary)

- `Suggests` + `requireNamespace()` for optional deps; no `install.packages()` in
  examples/tests/vignettes.
- No writes to home, package dir, or `getwd()` by default; `withr` temp files.
- Restore `options()`, `par()`, wd; no `<<-` / `.GlobalEnv` mutation.
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

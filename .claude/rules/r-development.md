---
description: R versions, project root, coding style, performance, parallelism (not package CRAN/CI)
globs: "**/*.{R,r}"
alwaysApply: false
paths:
  - "**/*.{R,r}"
---
# R and statistical computing

Packages add `r-package-development`. Follow the repository first; the
defaults below apply when it is silent.

## R versions

- Keep the host's system R (CRAN apt build on the laptop, see `software.qmd`;
  the cluster module or build there). Extra versions come from `rig` in user
  mode (`rig add`, `R-x.y.z`) when an `renv.lock` asks for one; never
  `rig default` over the system R.
- `renv`: switch to the intended R version before `renv::init()` or
  `renv::snapshot()`; after a minor upgrade, do not reuse the old user
  library.
- Editor: `r.rpath.linux` points to the host's R binary (User or Remote
  settings, `r-python-project`); the editor needs `languageserver`, `httpgd`
  and `jsonlite`.

## Project root

Non-package repositories (Quarto books, analysis projects, scripts that are
not part of an R package) anchor paths with `here::here()`. `here` wraps
`rprojroot`, caches the root (`.here`, an RStudio project, `DESCRIPTION`, or
the VCS root), and is the function to call from scripts and reports.

`fromhere` (`from_rproj()`, `from_git()`) is for a nested layout where several
markers exist and the root has to be named, such as a package inside a
monorepo. `fs` does not search for a root: use it to list, copy, and create
files (`fs::path()`, `fs::dir_create()`, `fs::dir_ls()`) on a path from `here`.

`rprojroot` criteria (`is_r_package`, `has_file()`) belong in package code
that implements its own root detection. Analysis scripts, including scripts
in a Quarto book, stay on `here`.

## R coding style

`snake_case`; `<-`; native `|>`; `pkg::fun()` in package code; explicit missing
values; `file.path()`/`fs`; `fixed = TRUE` for literal `grep`; cohesive functions.

**dplyr:** `{{ }}`, `.data[[var]]`, `all_of()`/`any_of()`, `:=` for dynamic names.

**Interfaces:** validate at boundary; forward `...` only with clear semantics;
actionable `stop()`/`warning()` without `paste()` wrappers.

## Performance and reproducibility

Profile first; vectorise; `data.table` or DuckDB when scale demands it.

Record seeds and `sessioninfo` for scientific outputs; `renv` when the repo uses
it; `targets` for heavy pipelines when justified.

## Visualisation

`data-visualisation` (house theme, devices, caveats); a repository theme wins.

## Parallelism and Shiny

`furrr` for reproducible map-style parallel work (independent RNG streams);
`mirai` for low-latency async when appropriate; no nested parallelism; configurable
workers.

Shiny: modules; no globals; testable domain logic.

## Preferred tools (when justified by repo)

`dplyr`/`tidyr`, `purrr`, `fs`/`here`, `checkmate`, `cli`, `glue`/`stringr`, `jsonlite`, `devtools` stack.

Resampling/tuning: **`mlr3`** ecosystem (see statistical-modelling rule).
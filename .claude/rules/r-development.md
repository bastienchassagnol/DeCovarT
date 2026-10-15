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

## Visualisation (`ggplot2`)

Before custom hacks: [R Graphics Cookbook](https://r-graphics.org/),
[ggplot2 book](https://ggplot2-book.org/), then current package docs.

Question-first plots; explicit `ggsave()` size/dpi; PDF/SVG for vector output.
When the repository already defines a ggplot theme, use that theme.

Otherwise use `theme_slide()` on every ggplot. It is `theme_minimal()` with the
grid, type, and background fixed:

```r
theme_slide <- function() {
  ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e6e1d8"),
      plot.title = ggplot2::element_text(
        face = "bold", colour = "#0e4d3a", size = 18
      ),
      plot.subtitle = ggplot2::element_text(colour = "#3d4a42", size = 13),
      axis.title = ggplot2::element_text(colour = "#1c2b24"),
      strip.text = ggplot2::element_text(face = "bold", colour = "#0e4d3a"),
      legend.position = "bottom",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}
```

Colour roles, unless the repository already assigns them:

- Series and density curves: `#1d4e89`
- Second series or a marked point: `#c45c26`
- Histogram bars: `#d9d3c7` with white bin outlines
- Guides, reference lines, and annotation text: `#1c2b24`

A multi-panel figure keeps `theme_slide()` and shrinks type so each panel stays
readable (titles about 15, subtitles about 11, axis titles about 12, axis text
about 10). Place the panels with `grid::grid.layout()`, not a new theme. Write
the file with `ragg::agg_png()` (or `ggsave(..., device = ragg::agg_png)`), on a
white background, and state width, height, and resolution in the call.

**Distributions:** `ggdist` rainclouds when shape matters; report binning rule
for histograms.

**Benchmarks:** forest/dot-whisker for summaries; facet metrics — do not map two
primary metrics to colour and size.

**Heatmaps:** `geom_tile()` default; `ComplexHeatmap` only for rich annotations.

**Interactive:** prefer **`ggiraph`** (scales better in HTML than Plotly for many
points). Keep a static archival figure.

## Parallelism and Shiny

`furrr` for reproducible map-style parallel work (independent RNG streams);
`mirai` for low-latency async when appropriate; no nested parallelism; configurable
workers.

Shiny: modules; no globals; testable domain logic.

## Preferred tools (when justified by repo)

`dplyr`/`tidyr`, `purrr`, `fs`/`here`, `ggplot2`, `tinytable` for static HTML
tables, `checkmate`, `cli`, `glue`/`stringr`, `jsonlite`, `devtools` stack.

Resampling/tuning: **`mlr3`** ecosystem (see statistical-modelling rule).

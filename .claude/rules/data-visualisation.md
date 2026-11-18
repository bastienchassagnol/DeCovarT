---
description: Figure design, house ggplot2 theme and the dataviz caveats checklist
globs: "**/*.{R,r,Rmd,qmd,py,tex,ipynb}"
alwaysApply: false
paths:
  - "**/*.{R,r,Rmd,qmd,py,tex,ipynb}"
---
# Data visualisation

Use this rule for any chart, plot or figure: R, Python, Quarto, LaTeX, slides.
Each figure answers one question. When the repository already defines a theme
or palette, use it; the defaults below apply when the repository is silent.
Code is R/`ggplot2`; the principles hold in any language.

Long material lives in the `data-visualisation` skill
(`~/.claude/skills/data-visualisation/references/` after
`agent-sync install-user`; same under `~/.cursor/skills/`): `chart-chooser.md`
(decision tree), `caveats.md` (avoid/prefer snippets), `line-chart-makeover.md`,
`libraries.md` (static → interactive ladder), `colour-palettes.md`,
`fonts.md`, and one file per chart module (`annotated-scatter.md`,
`latent-projection.md`, `treemap.md`, `ternary.md`,
`waterfall-clinical.md`). Read the one you need.

## House style

Use `theme_slide()` on every ggplot. It is `theme_minimal()` with the grid,
type, and background fixed:

```r
theme_slide <- function(base_family = NULL) {
  if (is.null(base_family)) base_family <- slide_family()
  ggplot2::theme_minimal(base_size = 18, base_family = base_family) +
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
about 10). Place the panels with `grid::grid.layout()`, not a new theme.

`slide_family()` is defined in `fonts.md`. It returns the single family in the
project `fonts/` directory, or a generic family that `systemfonts::match_fonts()`
resolves, and it names the family actually used. Pass another `base_family`
when `fonts/` holds more than one family.

Fonts are three steps (`fonts.md`): files in `fonts/`, discovery with
`systemfonts`, rendering with the device. Write a raster with
`ragg::agg_png()` (or `ggsave(..., device = ragg::agg_png)`), on a white
background, and state width, height, and resolution in the call. Vector output
is `svglite::svglite`, with an explicit `ggsave()` size. Do not draw the figure
with `pdf()` or `cairo_pdf`; convert the SVG if a venue wants PDF. Ordinary
text is `theme_slide()`. Markdown, mixed faces, or a text box use `marquee`
(`marquefy_theme()` on the finished theme). Do not use `showtext`, `extrafont`,
or `ggtext`.

## Choose the chart

1. Name the data type: numeric, categoric, categoric + numeric, relational,
   map, or time series.
2. Name the purpose: distribution, correlation, ranking, part of a whole,
   evolution, map, or flow.
3. Follow that branch in `chart-chooser.md`; apply its house overrides.

- **Distributions:** `ggdist` rainclouds when shape matters; report the binning
  rule for histograms.
- **Benchmarks:** forest or dot-whisker plots; facet the metrics; do not map
  two primary metrics to colour and size.
- **Heatmaps:** `geom_tile()` by default; `ComplexHeatmap` only for rich
  annotations; cluster rows and columns.
- **Line charts:** direct labels at line ends, no legend, a title that states
  the finding, context in grey and one or two highlighted series
  (`line-chart-makeover.md`).
- **Annotated scatter:** dashed guides at a pre-specified split, a few labels
  on the extremes (`annotated-scatter.md`).
- **Latent projections:** scaled PCA with variance on the axes
  (`latent-projection.md`). A nonlinear embedding is labelled as such.
- **Compositions:** a treemap when the reader needs the absolute count and the
  share of a parent (`treemap.md`). Exactly three closed parts: a ternary plot
  (`ternary.md`), not a radar chart.
- **Clinical response:** a waterfall of best percent change, coloured by best
  overall response (`waterfall-clinical.md`). It does not show survival.
- **Interactive:** `ggiraph` for ggplots in HTML; OJS + Plot or D3 for
  playgrounds in Quarto pages and decks; React only past that
  (`libraries.md`). Keep a static archival figure.

## Checklist before publishing a figure

Card numbers (#nn) refer to Yan Holtz's *61 dataviz caveats*; details and code
in `caveats.md`.

1. Scales: full axis or a justified log scale (#01); no dual axis (#05), 3D
   (#06), pie (#13) or radial bars (#14); the data or n shown behind every
   summary (#02, #04, #10).
2. Bars: ranked (#50); one colour unless highlighting (#18); lollipop or
   dumbbell when bars are long or paired (#12, #16).
3. Scatter and bubbles: deliberate aspect ratio (#58); overplotting handled
   (#61); size by area (#52); large bubbles at the back (#55).
4. Lines: no spaghetti (#53); dates readable (#59).
5. Colour: classify the mapped variable, build the palette with `colorspace`
   (qualitative, sequential, or diverging), preview the hex list on iWantHue
   (qualitative) or the r-graph-gallery palette finder, then run the ten-rule
   check in `colour-palettes.md`. Holtz cards still apply: sequential for
   numbers (#22), distinct hues for groups (#24), diverging with a light
   neutral midpoint for centred data (#26); few, soft, colour-blind-safe
   colours (#21, #23, #27); no rainbow (#25); dark annotations (#20).
6. Labels: horizontal (#40), direct (#41), rounded (#43), annotated (#44), not
   on the data.
7. Layout: aligned (#28), two fonts (#29), no box (#30), readable sizes (#31,
   nothing below about 14 px on screen), a title that states the finding (#33).
8. Several panels: same colours (#46), shared axes (#47), margins (#48),
   linked hover when interactive (#49); UpSet beyond three sets (#56).
9. Heatmaps: clustered order (#60), diverging palette for centred values.

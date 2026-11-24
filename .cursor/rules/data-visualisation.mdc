---
description: Figure design, house ggplot2 theme, output devices and the dataviz caveats checklist
globs: "**/*.{R,r,Rmd,qmd,tex}"
alwaysApply: false
paths:
  - "**/*.{R,r,Rmd,qmd,tex}"
---
# Data visualisation

Any chart, plot or figure in R, Quarto, LaTeX or slides. Plots are R
(`ggplot2`, `ggdist`, `ggiraph`) or OJS/D3 (React only past that); no Python
plotting. Each figure answers one question. A repository theme or palette
wins over the defaults below.

How-to lives in the `data-visualisation` skill (`references/`: inspiration,
chart chooser, ggplot2 production and the house theme, caveats, colours,
fonts, libraries, chart modules). Read the file you need.

**Workflow:** inspiration → chart choice → ggplot2 with `theme_slide()` →
checklist → `roast-my-chart` skill before publishing. Claude Code enforces
the inspiration step with the `dataviz-nudge` hook (`DATAVIZ_NUDGE=off`
disables it).

## House style and output

- `theme_slide()` on every ggplot (code, colour roles and multi-panel sizes in
  `references/ggplot2-production.md`; `slide_family()` in `fonts.md`).
- Fonts: files in the project `fonts/`, discovery with `systemfonts`, then the
  device. Mixed faces or text boxes use `marquee`; no `showtext`, `extrafont`
  or `ggtext`.
- Devices, SVG first: vector with `svglite::svglite`, raster with `ragg`
  (`ggsave(..., device = ragg::agg_png)`) on a white background, size and
  resolution stated in the call. Never `pdf()` or `cairo_pdf`. A venue or
  LaTeX PDF is made from the SVG:
  `inkscape f.svg --export-text-to-path --export-type=pdf -o f.pdf`.
- Captions define panels, axes, colours and uncertainty; readable type at
  final size; every display item cited in the text.

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

Card numbers (#nn) refer to Yan Holtz's dataviz caveats; details and code in
`caveats.md`.

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

Then run the `roast-my-chart` skill on the finished figure.

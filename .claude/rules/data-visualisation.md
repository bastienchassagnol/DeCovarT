# Data visualisation

Use this rule for any chart, plot or figure: R, Python, Quarto, LaTeX, slides.
Each figure answers one question. When the repository already defines a theme
or palette, use it; the defaults below apply when the repository is silent.
Code is R/`ggplot2`; the principles hold in any language.

References before custom hacks: [R Graphics Cookbook](https://r-graphics.org/),
[ggplot2 book](https://ggplot2-book.org/), then current package docs.

## House style

Use `theme_slide()` on every ggplot. It is `theme_minimal()` with the grid,
type, and background fixed:

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
white background, and state width, height, and resolution in the call. Use
PDF/SVG for vector output, with an explicit `ggsave()` size.

## Chart types

- **Distributions:** `ggdist` rainclouds when shape matters; report the binning
  rule for histograms.
- **Benchmarks:** forest or dot-whisker plots for summaries; facet the metrics;
  do not map two primary metrics to colour and size.
- **Heatmaps:** `geom_tile()` by default; `ComplexHeatmap` only for rich
  annotations.
- **Interactive:** prefer `ggiraph` (scales better in HTML than Plotly for many
  points). Keep a static archival figure.

## Caveats, by theme

Paraphrased from Yan Holtz, *61 dataviz caveats*
([data-to-viz.com](https://www.data-to-viz.com)); card numbers in brackets
(#nn). Maps (#15, #35–39) are left out. Snippets assume `library(ggplot2)` and
show `# avoid` then `# prefer`.

### 1. Honest scales and chart types

- **Full scale (#01).** Do not cut the y-axis of bars or areas: a cut axis
  inflates differences. When effects are multiplicative, change the scale
  (log) rather than the origin.

  ```r
  d <- data.frame(g = c("A", "B"), v = c(99, 80))
  # avoid
  ggplot(d, aes(g, v)) + geom_col() + coord_cartesian(ylim = c(50, 100))
  # prefer
  ggplot(d, aes(g, v)) + geom_col() + scale_y_continuous(limits = c(0, NA))
  # multiplicative effects: a log scale, not a cut axis
  ggplot(diamonds, aes(carat, price)) + geom_point(alpha = 0.1) + scale_y_log10()
  ```

- **No dual y-axis (#05):** the choice of the two scales manufactures or hides
  a correlation. Use two aligned panels, or index both series to a common base.
- **No 3D (#06):** it distorts perception and adds nothing.
- **No pie charts (#13), no radial bar charts (#14):** angles and arcs are hard
  to compare, and outer bars look longer. Use a sorted bar or lollipop chart.
  Circular barplots (#17) only with many bars, a clear pattern and a large
  inner circle.
- **Show the data behind summaries.** Boxplots hide bimodality and sample size
  (#02): add the points or use a raincloud. A mean alone hides n (#10): show
  the points or at least n. Error bars on bars hide the distribution (#04):
  show the points, and always say whether the bar is an SD, SE or CI.
- **Other traps:** follow conventions, no inverted axes (#03); show the wider
  context, not a convenient window (#07); draw arcs, bars and areas from the
  data, never by eye (#08); place time on a true scale so uneven gaps show
  (#09); check for Simpson's paradox by splitting groups (#11); try several
  histogram bin widths (#51).

### 2. Bar charts, lollipops, dumbbells

- **Meaningful order (#50):** rank the bars, unless the categories have a
  natural order.
- **One variable, one colour (#18, #23):** length already encodes the value;
  use a second colour only to highlight a subgroup.
- **Lollipop for one variable (#12):** many long bars create a moiré effect; a
  lollipop carries the same information with less ink.
- **Dumbbell for two groups or nested bars (#16):** side-by-side bars leave the
  reader to estimate the gap; a dumbbell draws it, the paired counterpart of a
  forest plot.
- **Stacked charts are hard to read (#57):** only the bottom series has a fixed
  baseline; use lines or small multiples.

```r
d <- data.frame(item = LETTERS[1:6], v = c(31, 52, 18, 44, 27, 39))
# avoid: unordered, multicoloured bars
ggplot(d, aes(item, v, fill = item)) + geom_col()
# prefer: ranked lollipop, one colour, one highlighted item
d$hi <- d$item == "B"
ggplot(d, aes(v, forcats::fct_reorder(item, v))) +
  geom_segment(aes(x = 0, xend = v, yend = item), colour = "grey70") +
  geom_point(aes(colour = hi), size = 3) +
  scale_colour_manual(values = c(`FALSE` = "#1d4e89", `TRUE` = "#c45c26"),
                      guide = "none") +
  labs(x = NULL, y = NULL)

# dumbbell: two groups per category
w <- data.frame(item = LETTERS[1:5], y1 = c(40, 55, 30, 62, 48),
                y2 = c(52, 50, 41, 70, 60))
ggplot(w, aes(y = forcats::fct_reorder(item, y2))) +
  geom_segment(aes(x = y1, xend = y2, yend = item), colour = "grey70") +
  geom_point(aes(x = y1), colour = "grey55", size = 3) +
  geom_point(aes(x = y2), colour = "#1d4e89", size = 3) +
  labs(x = NULL, y = NULL)
```

### 3. Scatter and bubble plots

- **Shape (#58):** stretching the panel changes how steep and strong a relation
  looks. Default to a square panel; use a 1:1 ratio when both axes share units.
- **Overplotting (#61):** with thousands of points, use small transparent
  points, a 2D density or hexagonal bins.
- **Bubble area, not radius (#52):** mapping a value to the radius inflates
  large bubbles with the square of the value. This applies to dot-cloud plots.
- **Z-index (#55):** draw large bubbles first, so small ones stay visible on
  top.

```r
# avoid: radius scaling, large bubbles drawn last, opaque blob
ggplot(mpg, aes(displ, hwy, size = cyl)) + geom_point() + scale_radius()
# prefer: area scaling, large bubbles at the back, square panel
d <- mpg[order(-mpg$cyl), ]
ggplot(d, aes(displ, hwy, size = cyl)) +
  geom_point(colour = "#1d4e89", alpha = 0.5) +
  scale_size_area(max_size = 8) +
  theme(aspect.ratio = 1)

# overplotting: transparency, hexbins or a 2D density
ggplot(diamonds, aes(carat, price)) + geom_point(size = 0.3, alpha = 0.05)
ggplot(diamonds, aes(carat, price)) + geom_hex(bins = 60)
ggplot(diamonds, aes(carat, price)) + geom_density_2d_filled()
```

### 4. Line charts

- **Spaghetti (#53):** many overlapping lines are unreadable. Grey the context
  and highlight the one or two series that carry the message, or use small
  multiples.
- **Connect the dots (#54):** when x is ordered (time), join the points.
- **Date axes (#59):** ticks at regular, meaningful intervals (months,
  quarters), with the year shown.

```r
e <- economics_long
# avoid
ggplot(e, aes(date, value01, colour = variable)) + geom_line()
# prefer: context in grey, one highlighted series, readable dates
ggplot(e, aes(date, value01, group = variable)) +
  geom_line(colour = "grey80") +
  geom_line(data = subset(e, variable == "unemploy"), colour = "#c45c26",
            linewidth = 1) +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y")
```

### 5. Colour

- **Intentional (#18):** colour only when it adds meaning; one consistent
  scheme.
- **Background (#19):** light, low-saturation backgrounds; softened dark tones,
  not pure black, in dark mode.
- **Contrast (#20):** enough contrast between marks, background and labels;
  write annotations in a dark colour.
- **No pure, flashy colours (#21):** pure red, green and blue are harsh and
  fail colour-blind readers; use slightly desaturated hues.
- **Numbers take a sequential palette (#22), groups take distinct hues (#24):**
  a gradient suggests an order between groups that does not exist; categorical
  hues on numbers force readers back to the legend.
- **Few colours (#23):** beyond a handful, readers cannot match them to the
  legend; highlight what matters and leave the rest grey.
- **No rainbow (#25):** it draws bands that are not in the data. Use a
  perceptually uniform palette (viridis, or a single-hue sequential).
- **Diverging data, diverging palette (#26):** when values fall on both sides of
  a meaningful midpoint (zero, a fold change of 1), use one hue per direction,
  meeting at a light neutral colour at the midpoint. This includes heatmaps.
- **Colour-blind readers (#27):** avoid red against green; blue and orange
  pairs are safe; label lines directly.

```r
z <- expand.grid(x = 1:20, y = 1:20)
z$v <- with(z, sin(x / 3) * cos(y / 4))
# avoid: rainbow on centred data
ggplot(z, aes(x, y, fill = v)) + geom_tile() + scale_fill_gradientn(colours = rainbow(7))
# prefer: diverging, light neutral midpoint at zero
ggplot(z, aes(x, y, fill = v)) + geom_tile() +
  scale_fill_gradient2(low = "#1d4e89", mid = "#f7f7f7", high = "#c45c26",
                       midpoint = 0)
# numbers without a midpoint: sequential
ggplot(z, aes(x, y, fill = abs(v))) + geom_tile() + scale_fill_viridis_c()
```

### 6. Labels, legends and annotation

- **Horizontal labels (#40):** choose a layout whose labels read horizontally,
  such as horizontal bars, rather than rotated text.
- **Direct labels, few legends (#41, #42):** a legend makes the eye travel back
  and forth. Label series directly, in empty space next to what they describe,
  typically the right margin at the end of each line, never on top of the data.
- **Round numbers (#43):** round to what matters, and drop an axis that only
  repeats the value labels.
- **Annotate (#44):** short notes, shaded areas and value labels that point at
  the finding.
- **Do not label every tick (#45):** label round steps, keep minor ticks.

```r
e <- subset(economics_long, variable %in% c("psavert", "uempmed"))
last <- e[e$date == max(e$date), ]
# avoid: legend and rotated labels
ggplot(e, aes(date, value, colour = variable)) + geom_line() +
  theme(axis.text.x = element_text(angle = 45))
# prefer: labels in the right margin, rounded values, an annotation
ggplot(e, aes(date, value, colour = variable)) +
  geom_line() +
  geom_text(data = last, aes(label = variable), hjust = 0, nudge_x = 200) +
  annotate("text", x = as.Date("1980-01-01"), y = 24, label = "Peak saving",
           colour = "#1c2b24", hjust = 0) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  scale_colour_manual(values = c(psavert = "#1d4e89", uempmed = "#c45c26"),
                      guide = "none") +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 70, 5.5, 5.5))
# crowded labels: ggrepel::geom_text_repel(direction = "y", nudge_x = 200)
```

### 7. Typography and layout

- **Align everything with something (#28):** title, subtitle, axis labels and
  notes share one left edge.
- **Two fonts at most (#29):** build the hierarchy with size and weight.
- **No box (#30):** drop the frame and the dense grid; keep a few light
  gridlines.
- **Font size (#31):** text must be readable at final size; nothing meant to be
  read below about 14 px on screen, and the title larger. In print, use the
  multi-panel sizes of the house style.
- **Group what belongs together (#32):** tight spacing within a group, more
  space between groups; labels next to their data.
- **The title says the finding (#33):** "Sales doubled since 2020", not
  "Sales by year"; highlight the finding on the chart.
- **Declutter (#34):** remove every element that carries no information.

```r
# avoid: framed panel, title offset from the y labels
ggplot(mpg, aes(displ, hwy)) + geom_point() + theme_bw() +
  labs(title = "Highway mileage by displacement")
# prefer: no frame, title and subtitle aligned with the y-axis labels
ggplot(mpg, aes(displ, hwy)) + geom_point(colour = "#1d4e89") +
  labs(title = "Larger engines use more fuel",
       subtitle = "Highway miles per gallon by displacement (L)",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(plot.title.position = "plot", panel.grid.minor = element_blank())
```

### 8. Several charts together

- **Consistent colours (#46):** the same group has the same colour in every
  panel; fix it with a named vector.
- **Shared y-axis (#47):** panels that compare one measure share its scale.
- **Margins between charts (#48):** clear gaps between panels and padding
  inside each.
- **Link charts that share an axis (#49):** on hover, highlight the same date
  or group in every panel.
- **Set intersections (#56):** beyond three sets, use an UpSet plot rather than
  a Venn diagram.

```r
pal <- c(compact = "#1d4e89", suv = "#c45c26", midsize = "#7a8b99")
d <- subset(mpg, class %in% names(pal))
# shared scale, fixed colours, space between panels
ggplot(d, aes(displ, hwy, colour = class)) + geom_point() +
  facet_wrap(~year, scales = "fixed") +
  scale_colour_manual(values = pal) +
  theme(panel.spacing = unit(1.5, "lines"))
# linked hover across panels: one data_id in every panel
p1 <- ggplot(d, aes(displ, hwy, data_id = model)) +
  ggiraph::geom_point_interactive()
p2 <- ggplot(d, aes(displ, cty, data_id = model)) +
  ggiraph::geom_point_interactive()
ggiraph::girafe(ggobj = patchwork::wrap_plots(p1, p2))
```

### 9. Heatmaps

- **Reorder rows and columns (#60):** alphabetical or random order looks like
  noise; order by hierarchical clustering so blocks appear.
- **Centred values take a diverging palette (#26):** see theme 5.

```r
m <- scale(as.matrix(mtcars))
ro <- hclust(dist(m))$order
co <- hclust(dist(t(m)))$order
h <- as.data.frame(as.table(m[ro, co]))
ggplot(h, aes(Var2, Var1, fill = Freq)) + geom_tile() +
  scale_fill_gradient2(low = "#1d4e89", mid = "#f7f7f7", high = "#c45c26") +
  labs(x = NULL, y = NULL)
```

## Checklist before publishing a figure

1. Scales: full axis or a justified log scale; no dual axis, 3D, pie or radial
   bars; the data or n shown behind every summary.
2. Bars: ranked; one colour unless highlighting; lollipop or dumbbell when
   bars are long or paired.
3. Scatter and bubbles: deliberate aspect ratio; overplotting handled; size
   by area; large bubbles at the back.
4. Lines: no spaghetti; dates readable.
5. Colour: sequential for numbers, distinct hues for groups, diverging with a
   light neutral midpoint for centred data; few, soft, colour-blind-safe
   colours; dark annotations.
6. Labels: horizontal, direct, rounded, annotated, not on the data.
7. Layout: aligned, two fonts, no box, readable sizes, a title that states the
   finding.
8. Several panels: same colours, shared axes, margins, linked hover when
   interactive.
9. Heatmaps: clustered order, diverging palette for centred values.

# Scenario book: the descriptors of every scenario, drawn twice.
# * scenario_metrics.pdf: a funky heatmap of shapes, one geom and one
#   continuous palette per column group, each column min-max scaled on its
#   own (log10 first when it spans more than two orders of magnitude); a
#   column constant over the design is printed as text.
# * scenario_metrics_table.{html,pdf}: the same values as numbers, a
#   reactable table for the vignettes and a tinytable PDF for the article.
# The mean signature has the same conditioning in every scenario, so
# kappa(mu) and lambda_min(mu'mu) are reported in a note, not as columns.

.scenario_columns <- list(
  design = c(
    n_cells = "Cells per bulk n_i",
    cell_size_ratio = "Cell-size ratio"
  ),
  composition = c(
    h_star = "Shannon H*",
    concentration = "||p||^2"
  ),
  network = c(
    network_density = "Edge density",
    network_mean_degree = "Mean degree",
    hoyer_abs_correlation = "Hoyer sparsity"
  ),
  information = c(
    mixsim_baromega = "Mean overlap",
    mixsim_maxomega = "Max overlap",
    hellinger = "Hellinger",
    f_cov_max = "f_cov max"
  ),
  complexity = c(
    kappa_sigma_p = "kappa Sigma(p)",
    lambda_min_sigma_p = "lambda_min Sigma(p)"
  )
)

# Column headers of the PDF table, in LaTeX.
.scenario_latex_names <- c(
  n_cells = "Cells per bulk $n_i$",
  cell_size_ratio = "Cell-size ratio",
  h_star = "Shannon $H^{\\star}$",
  concentration = "$\\lVert\\boldsymbol{p}\\rVert_2^2$",
  network_density = "Edge density",
  network_mean_degree = "Mean degree",
  hoyer_abs_correlation = "Hoyer sparsity",
  mixsim_baromega = "Mean overlap $\\bar\\omega$",
  mixsim_maxomega = "Max overlap $\\omega_{\\max}$",
  hellinger = "Hellinger",
  f_cov_max = "$\\max f_{\\mathrm{cov}}$",
  kappa_sigma_p = "$\\kappa\\{\\boldsymbol{\\Sigma}(\\boldsymbol{p})\\}$",
  lambda_min_sigma_p = "$\\lambda_{\\min}\\{\\boldsymbol{\\Sigma}(\\boldsymbol{p})\\}$"
)

.scenario_group_titles <- c(
  design = "Design",
  composition = "Composition",
  network = "Networks",
  information = "Statistical information",
  complexity = "Numerical conditioning"
)

.scenario_group_geoms <- c(
  design = "rect",
  composition = "circle",
  network = "funkyrect",
  information = "bar",
  complexity = "circle"
)

.scenario_group_palettes <- lapply(
  c(
    design = "Purples",
    composition = "Greens",
    network = "Greys",
    information = "Blues",
    complexity = "Oranges"
  ),
  function(name) RColorBrewer::brewer.pal(9L, name)[3:8]
)

# Which end of a column is the easiest scenario, for the green and red marks
# of the PDF table: "low" when the smallest value is the best one.
.scenario_best_end <- c(
  n_cells = "high",
  cell_size_ratio = "low",
  h_star = "high",
  concentration = "low",
  network_density = "low",
  network_mean_degree = "low",
  hoyer_abs_correlation = "high",
  mixsim_baromega = "low",
  mixsim_maxomega = "low",
  hellinger = "high",
  f_cov_max = "high",
  kappa_sigma_p = "low",
  lambda_min_sigma_p = "high"
)

# The mean signature is shared by every scenario of a figure; its condition
# number and smallest Gram eigenvalue go in a note. `latex` selects the
# markup of the PDF table.
.signature_note <- function(data, latex = FALSE) {
  value <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      return("not available")
    }
    if (dplyr::n_distinct(signif(x, 3L)) == 1L) {
      return(trimws(.signif_label(x[[1L]])))
    }
    sprintf(
      "%s to %s",
      trimws(.signif_label(min(x))),
      trimws(.signif_label(max(x)))
    )
  }
  same <- all(
    dplyr::n_distinct(signif(data$kappa_mu, 3L)) == 1L,
    dplyr::n_distinct(signif(data$lambda_min_mu, 3L)) == 1L
  )
  note <- sprintf(
    if (latex) {
      "Mean signature: $\\kappa(\\boldsymbol{\\mu})$ = %s, $\\lambda_{\\min}(\\boldsymbol{\\mu}^{\\top}\\boldsymbol{\\mu})$ = %s."
    } else {
      "Mean signature: kappa(mu) = %s, lambda_min(mu'mu) = %s."
    },
    value(data$kappa_mu),
    value(data$lambda_min_mu)
  )
  if (same) {
    sub(
      "Mean signature:",
      "Mean signature, the same in every scenario:",
      note,
      fixed = TRUE
    )
  } else {
    paste(note, "(range over the scenarios)")
  }
}

.data_scenario_metrics <- function(benchmark, layout, figure) {
  descriptors <- dplyr::left_join(
    benchmark$descriptors,
    benchmark$supplementary[c(
      "ID",
      setdiff(names(benchmark$supplementary), names(benchmark$descriptors))
    )],
    by = "ID"
  )
  descriptors <- descriptors[setdiff(
    names(descriptors),
    setdiff(names(layout), "ID")
  )]
  metrics <- dplyr::inner_join(layout, descriptors, by = "ID")
  # Smallest eigenvalue of the Gram matrix of the mean signature: near zero
  # when the means of two types are nearly collinear.
  metrics$lambda_min_mu <- vapply(
    metrics$ID,
    function(id) {
      mu <- benchmark$theta_true[[id]]$mu
      min(eigen(crossprod(mu), symmetric = TRUE, only.values = TRUE)$values)
    },
    numeric(1)
  )
  metrics$row_label <- as.character(metrics$panel)
  metrics$row_group <- as.character(metrics$page)
  metrics[order(metrics$page, metrics$panel), ]
}

# The descriptor columns present and finite in `data`, with their group.
.scenario_metric_columns <- function(data) {
  columns <- purrr::imap_dfr(.scenario_columns, function(cols, group) {
    tibble::tibble(id = names(cols), name = unname(cols), group = group)
  })
  columns[
    columns$id %in%
      names(data) &
      vapply(columns$id, function(id) any(is.finite(data[[id]])), logical(1)),
  ]
}

# Min-max scaling of one column to [0.25, 1], so the minimum stays visible;
# a positive column spanning more than two orders of magnitude is scaled on
# log10.
.scale_to_shape <- function(x) {
  finite <- x[is.finite(x)]
  if (all(finite > 0) && max(finite) / min(finite) > 100) {
    x <- log10(x)
    finite <- log10(finite)
  }
  0.25 + 0.75 * (x - min(finite)) / (max(finite) - min(finite))
}

plot_scenario_metrics <- function(data) {
  columns <- .scenario_metric_columns(data)
  table <- tibble::tibble(id = data$ID, scenario = data$row_label)
  info <- list(tibble::tibble(
    id = "scenario",
    name = "",
    geom = "text",
    group = "label",
    palette = NA_character_,
    width = 7,
    overlay = FALSE,
    legend = FALSE,
    hjust = 0
  ))
  for (k in seq_len(nrow(columns))) {
    id <- columns$id[[k]]
    group <- columns$group[[k]]
    value <- data[[id]]
    constant <- dplyr::n_distinct(signif(value[is.finite(value)], 6)) <= 1L
    if (constant) {
      table[[id]] <- .signif_label(value)
    } else {
      table[[id]] <- .scale_to_shape(value)
    }
    info[[length(info) + 1L]] <- tibble::tibble(
      id = id,
      name = columns$name[[k]],
      geom = if (constant) "text" else .scenario_group_geoms[[group]],
      group = group,
      palette = if (constant) NA_character_ else group,
      width = if (identical(.scenario_group_geoms[[group]], "bar")) {
        2.5
      } else {
        1.4
      },
      overlay = FALSE,
      legend = !constant,
      hjust = 0.5
    )
  }
  column_info <- dplyr::bind_rows(info)
  column_info$options <- lapply(column_info$geom, function(geom) {
    if (identical(geom, "text")) list(size = 2.6) else list()
  })
  groups <- unique(columns$group)
  column_groups <- tibble::tibble(
    group = c("label", groups),
    palette = c("label", groups),
    level1 = c("Scenario", unname(.scenario_group_titles[groups]))
  )
  palettes <- c(
    list(label = c("grey20", "grey20")),
    .scenario_group_palettes[groups]
  )
  # funkyheatmap lays its legends out on one fixed-width row, which crops
  # long titles; its legends are switched off and drawn by
  # .scenario_shape_legend() instead.
  legends <- lapply(c("label", groups), function(palette) {
    list(palette = palette, enabled = FALSE)
  })
  row_info <- tibble::tibble(id = data$ID, group = data$row_group)
  row_groups <- tibble::tibble(
    group = unique(data$row_group),
    level1 = unique(data$row_group)
  )
  heat <- funkyheatmap::funky_heatmap(
    data = table,
    column_info = column_info,
    row_info = row_info,
    row_groups = row_groups,
    column_groups = column_groups,
    palettes = palettes,
    legends = legends,
    scale_column = FALSE,
    add_abc = FALSE,
    position_args = funkyheatmap::position_arguments(
      col_annot_offset = 4,
      col_annot_angle = 40,
      expand_xmin = 1,
      expand_xmax = 6,
      expand_ymax = 2
    )
  )
  page <- patchwork::wrap_plots(
    heat,
    .scenario_shape_legend(unique(
      column_info$group[column_info$geom != "text"]
    )),
    ncol = 1L,
    heights = c(heat$height, 1.4)
  ) +
    patchwork::plot_annotation(
      title = "Descriptors of the simulated scenarios",
      caption = paste(
        "Shape size and colour scale each column on its own, from its minimum to its maximum (log10 first when",
        "the column spans more than two orders of magnitude);\nthe smallest shape is the column minimum.",
        "Columns constant over the design show their value. Numbers: scenario_metrics_table.html / .pdf.",
        paste0("\n", .signature_note(data))
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 15),
        plot.caption = ggplot2::element_text(
          hjust = 0,
          size = 9,
          colour = "grey25"
        )
      )
    )
  page$width <- max(heat$width * 1.1 + 1, 15)
  page$height <- heat$height * 1.05 + 3.5
  list(page)
}

# One key per column group, under the heatmap: its shape and palette from
# the column minimum to the column maximum.
.scenario_shape_legend <- function(groups) {
  keys <- purrr::map_dfr(groups, function(group) {
    tibble::tibble(
      title = .scenario_group_titles[[group]],
      geom = .scenario_group_geoms[[group]],
      k = 1:5,
      # funkyheatmap rectangles fill their cell; only colour varies
      size = if (identical(.scenario_group_geoms[[group]], "rect")) {
        rep(0.8, 5L)
      } else {
        c(0.25, 0.44, 0.62, 0.81, 1)
      },
      colour = grDevices::colorRampPalette(.scenario_group_palettes[[group]])(
        5L
      )
    )
  })
  keys$title <- factor(keys$title, unique(keys$title))
  areas <- keys[keys$geom != "bar", ]
  bars <- keys[keys$geom == "bar", ]
  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = areas,
      ggplot2::aes(
        .data$k,
        0,
        size = .data$size,
        fill = .data$colour,
        shape = .data$geom
      ),
      colour = "black",
      stroke = 0.3
    ) +
    ggplot2::geom_rect(
      data = bars,
      ggplot2::aes(
        xmin = .data$k - 0.4,
        xmax = .data$k + 0.4,
        ymin = -0.35,
        ymax = -0.35 + 0.45 * .data$size,
        fill = .data$colour
      ),
      colour = "black",
      linewidth = 0.2
    ) +
    ggplot2::annotate(
      "text",
      x = 1,
      y = -0.75,
      label = "min",
      size = 3
    ) +
    ggplot2::annotate(
      "text",
      x = 5,
      y = -0.75,
      label = "max",
      size = 3
    ) +
    ggplot2::facet_wrap(~title, nrow = 1L) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_size_area(max_size = 7, guide = "none") +
    ggplot2::scale_shape_manual(
      values = c(circle = 21, rect = 22, funkyrect = 22),
      guide = "none"
    ) +
    ggplot2::coord_cartesian(xlim = c(0.4, 5.6), ylim = c(-0.95, 0.5)) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# Numeric table --------------------------------------------------------------------

.scenario_metric_values <- function(data) {
  columns <- .scenario_metric_columns(data)
  values <- data[c("row_group", "row_label", columns$id)]
  names(values) <- c("Group", "Scenario", columns$name)
  list(values = values, columns = columns)
}

# Conditional marks of the PDF table, on the original body rows: the easiest
# value of each numeric column in green and the hardest in red; a shaded
# background and a font size growing with the value for three columns; thin
# grey rules between the column groups.
.style_scenario_cells <- function(pdf, values, columns, starts) {
  # group_tt() inserts one row before each group start and does not move the
  # styles set earlier, so the marks are set after it on the shifted rows.
  shift <- function(i) i + findInterval(i, starts)
  green <- "#1B7F3B"
  red <- "#B3261E"
  shaded <- c("h_star", "mixsim_baromega", "kappa_sigma_p")
  ramp <- grDevices::colorRamp(
    rev(grDevices::hcl.colors(9L, "Blues 3"))[1:4]
  )
  n <- nrow(values)
  for (k in seq_len(nrow(columns))) {
    id <- columns$id[[k]]
    x <- values[[columns$name[[k]]]]
    j <- k + 1L
    finite <- is.finite(x)
    if (dplyr::n_distinct(signif(x[finite], 3L)) <= 1L) {
      next
    }
    if (id %in% shaded) {
      unit <- (x - min(x[finite])) / diff(range(x[finite]))
      if (all(x[finite] > 0) && max(x[finite]) / min(x[finite]) > 100) {
        y <- log10(x)
        unit <- (y - min(y[finite])) / diff(range(y[finite]))
      }
      unit[!finite] <- 0
      pdf <- tinytable::style_tt(
        pdf,
        i = shift(seq_len(n)),
        j = j,
        background = grDevices::rgb(ramp(unit), maxColorValue = 255),
        fontsize = 0.8 + 0.4 * unit
      )
    }
    if (id %in% names(.scenario_best_end)) {
      high <- identical(unname(.scenario_best_end[[id]]), "high")
      best <- which(
        finite &
          signif(x, 6L) ==
            signif(if (high) max(x[finite]) else min(x[finite]), 6L)
      )
      worst <- which(
        finite &
          signif(x, 6L) ==
            signif(if (high) min(x[finite]) else max(x[finite]), 6L)
      )
      pdf <- tinytable::style_tt(
        pdf,
        i = shift(best),
        j = j,
        color = green,
        bold = TRUE
      )
      pdf <- tinytable::style_tt(
        pdf,
        i = shift(worst),
        j = j,
        color = red,
        bold = TRUE
      )
    }
  }
  first <- 1L + match(unique(columns$group), columns$group)
  tinytable::style_tt(
    pdf,
    i = shift(seq_len(n)),
    j = first,
    line = "l",
    line_color = "grey70",
    line_width = 0.05
  )
}

# reactable HTML (grouped by page, cells tinted by their column rank) and
# tinytable PDF of the descriptors, three significant digits.
write_scenario_metrics_table <- function(data, file_stem) {
  table <- .scenario_metric_values(data)
  values <- table$values
  numeric_cols <- table$columns$name
  tint <- function(value, column) {
    x <- values[[column]]
    if (!is.finite(value) || dplyr::n_distinct(x[is.finite(x)]) <= 1L) {
      return(NULL)
    }
    group <- table$columns$group[table$columns$name == column]
    ramp <- grDevices::colorRamp(.scenario_group_palettes[[group]])
    rank <- (value - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))
    list(
      background = grDevices::rgb(ramp(0.15 + 0.6 * rank), maxColorValue = 255)
    )
  }
  html <- reactable::reactable(
    values,
    groupBy = "Group",
    defaultExpanded = TRUE,
    pagination = FALSE,
    compact = TRUE,
    striped = FALSE,
    columns = c(
      list(Group = reactable::colDef(minWidth = 260)),
      stats::setNames(
        lapply(numeric_cols, function(column) {
          reactable::colDef(
            cell = function(value) .signif_label(value),
            style = function(value) tint(value, column),
            align = "right"
          )
        }),
        numeric_cols
      )
    ),
    columnGroups = lapply(unique(table$columns$group), function(group) {
      reactable::colGroup(
        name = .scenario_group_titles[[group]],
        columns = table$columns$name[table$columns$group == group]
      )
    })
  )
  html <- htmlwidgets::appendContent(
    html,
    htmltools::tags$p(.signature_note(data), style = "font-size: 0.85em;")
  )
  dir.create(dirname(file_stem), recursive = TRUE, showWarnings = FALSE)
  htmlwidgets::saveWidget(
    html,
    normalizePath(paste0(file_stem, ".html"), mustWork = FALSE),
    selfcontained = TRUE,
    title = "Descriptors of the simulated scenarios"
  )
  # saveWidget() inlines every dependency but leaves their copies behind.
  unlink(paste0(file_stem, "_files"), recursive = TRUE)
  printed <- values
  printed[numeric_cols] <- lapply(printed[numeric_cols], .signif_label)
  printed$Group <- NULL
  # Headers in LaTeX; the cells hold plain labels and numbers only.
  names(printed) <- c(
    "Scenario",
    unname(.scenario_latex_names[table$columns$id])
  )
  groups <- rle(as.character(values$Group))
  starts <- cumsum(c(1L, utils::head(groups$lengths, -1L)))
  column_groups <- split(
    seq_along(table$columns$id) + 1L,
    factor(
      .scenario_group_titles[table$columns$group],
      unique(.scenario_group_titles[table$columns$group])
    )
  )
  pdf <- tinytable::tt(
    printed,
    notes = list(
      .signature_note(data, latex = TRUE),
      paste(
        "Green: easiest value of a column; red: hardest. Easiest is the highest",
        "number of cells, Shannon entropy, Hoyer sparsity, Hellinger distance,",
        "covariance share and $\\lambda_{\\min}$; the lowest cell-size ratio,",
        "concentration, edge density, mean degree, overlap and condition number.",
        "Shading and size: Shannon entropy, mean overlap and condition number."
      )
    )
  ) |>
    tinytable::group_tt(
      i = stats::setNames(as.list(starts), groups$values),
      j = column_groups
    ) |>
    .style_scenario_cells(values, table$columns, starts) |>
    tinytable::style_tt(
      i = "groupi",
      bold = TRUE,
      fontsize = 1.25,
      align = "l",
      line = "t",
      line_color = "grey70",
      line_width = 0.05
    ) |>
    tinytable::style_tt(i = "groupj", bold = TRUE)
  # tinytable's own PDF wraps the table in a text-width minipage, which is
  # too narrow for twenty columns; the same preamble with a free width
  # keeps the table at its natural size.
  preamble <- sub(
    "\\\\minipage\\{\\\\textwidth\\}\n%s\n\\\\endminipage",
    "%s",
    get("latex_standalone", asNamespace("tinytable"))
  )
  tex <- paste0(file_stem, ".tex")
  writeLines(
    sprintf(
      sub(
        "\\\\documentclass\\{standalone\\}",
        "\\\\documentclass[border=8pt]{standalone}\n\\\\usepackage{amsmath}",
        preamble
      ),
      # the float of tinytable has no place in a standalone page
      gsub(
        "\\\\begin\\{table\\}|\\\\end\\{table\\}|\\\\centering",
        "",
        tinytable::save_tt(pdf, output = "latex")
      )
    ),
    tex
  )
  tinytex::latexmk(
    tex,
    engine = "xelatex",
    pdf_file = paste0(file_stem, ".pdf")
  )
  unlink(c(tex, paste0(file_stem, c(".log", ".aux"))))
  invisible(file_stem)
}

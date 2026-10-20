# Build DeCovarT function call network with pkgapi + visNetwork
# - Excludes R/zzz.R and R/data.R
# - Edges only between DeCovarT-defined functions (no base/stats/pkg:: calls)
# - One colour per source file; node tooltips include formal parameters
# see also https://pkgmap.app to list dependencies assoicated with a package
#
# Run from the package root:
#   Rscript scripts/auxiliary/generate_package_network_functions.R

stop_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Missing package '",
      pkg,
      "'. Install with e.g. pak::pak(\"",
      pkg,
      "\").",
      call. = FALSE
    )
  }
}

stop_if_missing("pkgapi")
stop_if_missing("visNetwork")

project_path <- getwd()
excluded_files <- c("zzz.R", "data.R", "utils-ui.R")
# Tiny I/O / argparse helpers clutter the call graph without clarifying
# the statistical API; keep them out of the visNetwork report.
excluded_functions <- c(
  ".match_arg_case_insensitive",
  ".write_artefact",
  ".ensure_file_suffix"
)

out_dir <- file.path(project_path, "output", "package_network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Two graphs. Core is estimation and intervals for the Gaussian
# convolution. Simulation is the generative scenarios and the benchmark.
network_groups <- list(
  list(
    files = c(
      "02_01_first_generation_deconvolution.R",
      "02_02_covariance_backends.R",
      "02_03_decovart_estimate_ratios.R",
      "02_04_decovart_fit.R",
      "03_01_decovart_inference.R"
    ),
    csv = "decovart_functions_core.csv",
    html = "decovart_function_network_core.html",
    legacy_html = "decovart_function_network.html"
  ),
  list(
    files = c(
      "01_01_toy_simulation.R",
      "01_02_generate_synthetic_networks.R",
      "01_03_scenario_descriptors.R",
      "01_04_infer_gaussian_overlap.R",
      "01_05_target_gaussian_overlap.R",
      "04_01_run_benchmark.R",
      "04_02_benchmark_visualisation.R",
      "04_03_run_simulation_benchmark.R",
      "04_04_statistical_metrics.R"
    ),
    csv = "decovart_functions_simulation.csv",
    html = "decovart_function_network_simulation.html",
    legacy_html = NULL
  )
)

# -----------------------------------------------------------------------------
# Formal parameter names from AST (no eval; avoids failures on complex defs)
# -----------------------------------------------------------------------------

function_formal_names <- function(fun_expr) {
  if (!is.call(fun_expr) || !identical(fun_expr[[1]], quote(`function`))) {
    return(character(0))
  }
  pl <- fun_expr[[2]]
  if (is.null(pl)) {
    return(character(0))
  }
  if (!is.pairlist(pl)) {
    return(character(0))
  }
  nm <- names(pl)
  if (is.null(nm)) {
    return(character(0))
  }
  nm[nm != ""]
}

# -----------------------------------------------------------------------------
# 1) List all internal functions and extract their parameters
# -----------------------------------------------------------------------------

r_files <- list.files(
  path = file.path(project_path, "R"),
  pattern = "\\.R$",
  full.names = TRUE
)
r_files <- r_files[!basename(r_files) %in% excluded_files]

extract_function_metadata <- function(file_path) {
  exprs <- parse(file_path, keep.source = FALSE)
  function_calls <- Filter(
    function(expr) {
      is.call(expr) &&
        identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) &&
        is.call(expr[[3]]) &&
        identical(expr[[3]][[1]], as.name("function"))
    },
    as.list(exprs)
  )

  if (length(function_calls) == 0) {
    return(
      data.frame(
        function_name = character(0),
        file = character(0),
        parameters = character(0),
        stringsAsFactors = FALSE
      )
    )
  }

  do.call(
    rbind,
    lapply(function_calls, function(expr) {
      fname <- as.character(expr[[2]])
      fparams <- function_formal_names(expr[[3]])
      params_label <- if (length(fparams) == 0) {
        "<none>"
      } else {
        paste(fparams, collapse = ", ")
      }
      data.frame(
        function_name = fname,
        file = basename(file_path),
        parameters = params_label,
        stringsAsFactors = FALSE
      )
    })
  )
}

internal_functions <- do.call(rbind, lapply(r_files, extract_function_metadata))
internal_functions <- internal_functions[
  !duplicated(internal_functions$function_name),
]
internal_functions <- internal_functions[
  !internal_functions$function_name %in% excluded_functions,
]

if (nrow(internal_functions) == 0L) {
  stop("No top-level functions found in R/ (after exclusions).", call. = FALSE)
}

message(
  "Internal functions (excluding ",
  paste(excluded_files, collapse = ", "),
  "): ",
  nrow(internal_functions)
)

# -----------------------------------------------------------------------------
# 2) Build call edges from pkgapi (see https://github.com/r-lib/pkgapi)
# -----------------------------------------------------------------------------

pkg_map <- pkgapi::map_package(path = project_path)
calls <- pkg_map$calls

calls$file_base <- basename(calls$file)
calls <- calls[!calls$file_base %in% excluded_files, ]
calls$target <- calls$str

internal_names <- internal_functions$function_name
# Only calls where both ends are functions defined in this package (R/ sources).
calls <- calls[
  calls$from %in% internal_names & calls$target %in% internal_names,
]

edges <- unique(data.frame(
  from = calls$from,
  to = calls$target,
  arrows = "to",
  color = "#8ea0b3",
  stringsAsFactors = FALSE
))

# S3 methods on `decovart_fit` are not *called* by fit_decovart(); they
# dispatch on objects it returns. Add synthetic edges so the network
# shows the constructor ↔ accessor family.
s3_fit_methods <- grep(
  "\\.decovart_fit$",
  internal_functions$function_name,
  value = TRUE
)
if (
  "fit_decovart" %in%
    internal_functions$function_name &&
    length(s3_fit_methods) > 0L
) {
  edges <- rbind(
    edges,
    data.frame(
      from = "fit_decovart",
      to = s3_fit_methods,
      arrows = "to",
      color = "#c4a35a",
      stringsAsFactors = FALSE
    )
  )
}

for (group in network_groups) {
  fns <- internal_functions[internal_functions$file %in% group$files, ]
  if (nrow(fns) == 0L) {
    stop("No functions in network group ", group$html, ".", call. = FALSE)
  }
  names_keep <- fns$function_name
  group_edges <- edges[
    edges$from %in% names_keep & edges$to %in% names_keep,
  ]
  message(
    group$html,
    ": ",
    nrow(fns),
    " functions, ",
    nrow(group_edges),
    " edges"
  )
  utils::write.csv(
    fns,
    file.path(out_dir, group$csv),
    row.names = FALSE
  )

  file_levels <- unique(fns$file)
  palette <- grDevices::hcl.colors(
    max(1L, length(file_levels)),
    palette = "Set 2"
  )
  file_colours <- stats::setNames(palette, file_levels)
  nodes <- data.frame(
    id = fns$function_name,
    label = fns$function_name,
    type = "DeCovarT",
    file = fns$file,
    parameters = fns$parameters,
    stringsAsFactors = FALSE
  )
  nodes$shape <- "ellipse"
  nodes$color.background <- unname(file_colours[nodes$file])
  nodes$color.border <- "#2f3e4f"
  nodes$title <- paste0(
    "<b>",
    nodes$id,
    "</b><br/>",
    "File: ",
    nodes$file,
    "<br/>",
    "Parameters: ",
    nodes$parameters
  )
  legend_nodes <- rbind(
    data.frame(
      label = "DeCovarT function",
      shape = "ellipse",
      color = "#b8d8b8",
      stringsAsFactors = FALSE
    ),
    data.frame(
      label = paste0("File: ", names(file_colours)),
      shape = "dot",
      color = unname(file_colours),
      stringsAsFactors = FALSE
    )
  )
  plot_edges <- group_edges
  if (nrow(plot_edges) == 0L) {
    plot_edges <- data.frame(
      from = character(0),
      to = character(0),
      arrows = character(0),
      color = character(0),
      stringsAsFactors = FALSE
    )
  }
  visnetwork_plot <- visNetwork::visNetwork(
    nodes = nodes,
    edges = plot_edges,
    width = "100%",
    height = "800px"
  ) |>
    visNetwork::visLegend(
      useGroups = FALSE,
      addNodes = legend_nodes,
      addEdges = data.frame(
        label = "Function call",
        color = "#8ea0b3",
        stringsAsFactors = FALSE
      )
    ) |>
    visNetwork::visInteraction(navigationButtons = TRUE) |>
    visNetwork::visOptions(highlightNearest = TRUE) |>
    visNetwork::visHierarchicalLayout(
      direction = "LR",
      levelSeparation = 220,
      nodeSpacing = 140
    ) |>
    visNetwork::visEdges(
      smooth = list(enabled = TRUE, type = "cubicBezier", roundness = 0.35)
    )
  html_path <- file.path(out_dir, group$html)
  visNetwork::visSave(visnetwork_plot, html_path)
  message("Wrote: ", html_path)
  if (!is.null(group$legacy_html)) {
    legacy_path <- file.path(out_dir, group$legacy_html)
    file.copy(html_path, legacy_path, overwrite = TRUE)
    message("Wrote: ", legacy_path)
  }
}

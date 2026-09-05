# Optional `cli` UI (Suggests). One cached install check; plain
# message() / warning() / stop() when `cli` is not installed.

.ui_cli_state <- local({
  cache <- new.env(parent = emptyenv())
  cache$value <- NULL
  list(
    has_cli = function() {
      if (is.null(cache$value)) {
        cache$value <- requireNamespace("cli", quietly = TRUE)
      }
      isTRUE(cache$value)
    },
    set = function(value) {
      cache$value <- value
      invisible(value)
    }
  )
})

.has_cli <- function() {
  .ui_cli_state$has_cli()
}

#' @keywords internal
#' @noRd
.ui_set_has_cli <- function(value) {
  .ui_cli_state$set(value)
}

.ui_eval_txt <- function(expr_txt, envir) {
  val <- tryCatch(
    eval(parse(text = expr_txt, keep.source = FALSE)[[1L]], envir = envir),
    error = function(e) expr_txt
  )
  if (is.null(val)) {
    return("NULL")
  }
  if (length(val) > 1L) {
    return(toString(val))
  }
  paste(as.character(val), collapse = "")
}

.ui_style_plain <- function(style, value) {
  switch(
    style,
    arg = ,
    code = ,
    fun = ,
    fn = ,
    field = ,
    pkg = paste0("`", value, "`"),
    val = paste0('"', value, '"'),
    path = ,
    file = ,
    url = value,
    value
  )
}

.ui_sub_pattern <- function(text, pattern, fun) {
  repeat {
    m <- regexec(pattern, text, perl = TRUE)
    loc <- m[[1L]]
    if (length(loc) < 1L || loc[[1L]] == -1L) {
      break
    }
    caps <- regmatches(text, m)[[1L]]
    replacement <- do.call(fun, as.list(caps[-1L]))
    start <- loc[[1L]]
    len <- attr(loc, "match.length")[[1L]]
    text <- paste0(
      substr(text, 1L, start - 1L),
      replacement,
      substr(text, start + len, nchar(text))
    )
  }
  text
}

.ui_plain_one <- function(text, envir) {
  if (!length(text) || !nzchar(text)) {
    return(as.character(text))
  }
  text <- .ui_sub_pattern(
    text,
    "\\{\\.([A-Za-z0-9_]+) +\\{([^{}]+)\\}\\}",
    function(style, expr) {
      .ui_style_plain(style, .ui_eval_txt(expr, envir))
    }
  )
  text <- .ui_sub_pattern(
    text,
    "\\{\\.([A-Za-z0-9_]+) +([^{}]+)\\}",
    function(style, lit) .ui_style_plain(style, lit)
  )
  .ui_sub_pattern(
    text,
    "\\{([^{}]+)\\}",
    function(expr) .ui_eval_txt(expr, envir)
  )
}

.ui_plain <- function(message, envir) {
  if (length(message) == 1L && is.null(names(message))) {
    return(.ui_plain_one(message, envir))
  }
  nms <- names(message)
  if (is.null(nms)) {
    nms <- rep("", length(message))
  }
  lines <- character(length(message))
  for (i in seq_along(message)) {
    body <- .ui_plain_one(message[[i]], envir)
    prefix <- switch(
      nms[[i]],
      i = "i ",
      x = "x ",
      v = "* ",
      `*` = "* ",
      `!` = "! ",
      ""
    )
    lines[[i]] <- paste0(prefix, body)
  }
  paste(lines, collapse = "\n")
}

#' @keywords internal
#' @noRd
.ui_h1 <- function(message, ..., .envir = parent.frame()) {
  if (.has_cli()) {
    cli::cli_h1(message, ..., .envir = .envir)
    return(invisible(NULL))
  }
  title <- .ui_plain_one(message, .envir)
  rule <- paste(rep("-", 64L), collapse = "")
  message(rule)
  message(title)
  message(rule)
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.ui_info <- function(message, ..., .envir = parent.frame()) {
  if (.has_cli()) {
    cli::cli_alert_info(message, ..., .envir = .envir)
    return(invisible(NULL))
  }
  message(.ui_plain(message, .envir))
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.ui_success <- function(message, ..., .envir = parent.frame()) {
  if (.has_cli()) {
    cli::cli_alert_success(message, ..., .envir = .envir)
    return(invisible(NULL))
  }
  message(.ui_plain(message, .envir))
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.ui_warn <- function(message, ..., .envir = parent.frame()) {
  if (.has_cli()) {
    cli::cli_alert_warning(message, ..., .envir = .envir)
    return(invisible(NULL))
  }
  message(.ui_plain(message, .envir))
  invisible(NULL)
}

#' Signal a warning condition (`cli_warn` or base `warning()`)
#'
#' @keywords internal
#' @noRd
.ui_warning <- function(
  message,
  ...,
  call = rlang::caller_env(),
  .envir = parent.frame()
) {
  if (.has_cli()) {
    cli::cli_warn(message, ..., call = call, .envir = .envir)
    return(invisible(NULL))
  }
  warning(.ui_plain(message, .envir), call. = FALSE)
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.ui_abort <- function(
  message,
  ...,
  call = rlang::caller_env(),
  .envir = parent.frame()
) {
  if (.has_cli()) {
    cli::cli_abort(message, ..., call = call, .envir = .envir)
  }
  stop(.ui_plain(message, .envir), call. = FALSE)
}

#' Bind UI helpers into a script's calling environment
#'
#' @keywords internal
#' @noRd
.ui_attach_script <- function(env = parent.frame()) {
  assign(".ui_h1", .ui_h1, envir = env)
  assign(".ui_info", .ui_info, envir = env)
  assign(".ui_success", .ui_success, envir = env)
  assign(".ui_warn", .ui_warn, envir = env)
  assign(".ui_abort", .ui_abort, envir = env)
  invisible(NULL)
}

.ui_scenario_label <- function(row) {
  if (is.null(row) || !ncol(row)) {
    return("unnamed scenario")
  }
  if ("ID" %in% names(row)) {
    id <- row[["ID"]][[1L]]
    if (length(id) == 1L && !is.na(id) && nzchar(as.character(id))) {
      return(as.character(id))
    }
  }
  skip <- c("true_theta", "n", "true_parameters")
  cols <- setdiff(names(row), skip)
  if (!length(cols)) {
    return("unnamed scenario")
  }
  bits <- character(0L)
  for (nm in cols) {
    x <- row[[nm]]
    if (is.list(x)) {
      x <- x[[1L]]
    } else {
      x <- x[[1L]]
    }
    if (length(x) != 1L || is.na(x)) {
      next
    }
    bits <- c(bits, paste0(nm, "=", as.character(x)))
  }
  if (!length(bits)) {
    return("unnamed scenario")
  }
  paste(bits, collapse = " · ")
}

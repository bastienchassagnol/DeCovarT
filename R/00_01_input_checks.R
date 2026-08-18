#' Case-insensitive match to a set of allowed strings
#'
#' Like [match.arg()], but `arg` is matched after `tolower()` so `"Sigma"`
#' and `"sigma"` are equivalent. If `arg` is the full `choices` vector (the
#' usual default when the caller omits the argument), the first choice is
#' returned.
#'
#' @param arg Character vector supplied by the user.
#' @param choices Character vector of allowed values (canonical case).
#'
#' @return A character vector of canonical `choices` entries.
#'
#' @srrstats {G2.3} Univariate character arguments are restricted to a
#'   documented set of values.
#' @srrstats {G2.3a} Validation uses this `match.arg()` equivalent.
#' @srrstats {G2.3b} Matching is case-insensitive via `tolower()`.
#' @srrstats {G3.0} Numeric "equality" elsewhere uses
#'   `100 * .Machine$double.eps`, not `==` on floats.
#'
#' @keywords internal
#' @noRd
.match_arg_ci <- function(arg, choices) {
  if (!is.character(arg) || length(arg) < 1L) {
    stop("Argument must be a non-empty character vector.", call. = FALSE)
  }
  if (anyNA(arg)) {
    stop("Argument must not contain missing values.", call. = FALSE)
  }
  # match.arg() default: omitted arg is the full choices vector.
  if (
    length(arg) == length(choices) &&
      setequal(tolower(arg), tolower(choices))
  ) {
    arg <- arg[[1L]]
  }
  vapply(
    arg,
    function(a) {
      hit <- choices[tolower(choices) == tolower(a)]
      if (length(hit) != 1L) {
        stop(
          "Argument must be one of ",
          paste(choices, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
      hit[[1L]]
    },
    character(1L),
    USE.NAMES = FALSE
  )
}

#' Error if an object contains missing values
#'
#' @param x Atomic vector, matrix, or array.
#' @param name Argument name used in the error message.
#'
#' @srrstats {G2.13} Missing values are rejected before analytic routines.
#' @srrstats {G2.15} Downstream `mean` / `var` / `cor` therefore never see
#'   incomplete data (they keep default `na.rm = FALSE`).
#'
#' @keywords internal
#' @noRd
.assert_no_missing <- function(x, name) {
  if (anyNA(x)) {
    stop("`", name, "` must not contain missing values.", call. = FALSE)
  }
  invisible(TRUE)
}

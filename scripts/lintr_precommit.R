#!/usr/bin/env Rscript
# Local pre-commit lintr entry: silence renv sync chatter from the hook
# environment (DeCovarT itself does not use renv).

Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  RENV_CONFIG_SYNCHRONIZED_CHECK = "FALSE"
)

files <- commandArgs(trailingOnly = TRUE)
files <- files[grepl("\\.[Rr]$", files)]
files <- files[!grepl("^vignettes/", files)]

if (!length(files)) {
  quit(save = "no", status = 0L)
}

suppressPackageStartupMessages({
  if (!requireNamespace("lintr", quietly = TRUE)) {
    stop(
      "Package 'lintr' is required for the lintr pre-commit hook.",
      call. = FALSE
    )
  }
})

lints <- unlist(lapply(files, lintr::lint), recursive = FALSE)
if (!length(lints)) {
  quit(save = "no", status = 0L)
}

print(lints)
# Match upstream precommit lintr --warn_only: print but do not fail.
quit(save = "no", status = 0L)

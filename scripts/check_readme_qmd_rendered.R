# Fail if README.qmd is newer than README.md (Quarto source must be rendered).
# Analogous to {precommit}'s readme-rmd-rendered hook, for README.qmd.
readme_qmd <- "README.qmd"
readme_md <- "README.md"

if (!file.exists(readme_qmd)) {
  quit(save = "no", status = 0L)
}

if (!file.exists(readme_md)) {
  stop(
    "README.md is missing. Render README.qmd with Quarto before committing.",
    call. = FALSE
  )
}

qmd_mtime <- file.info(readme_qmd)$mtime
md_mtime <- file.info(readme_md)$mtime

if (isTRUE(qmd_mtime > md_mtime)) {
  stop(
    "README.qmd is newer than README.md. Render it before committing, e.g.\n",
    "  quarto::quarto_render('README.qmd')",
    call. = FALSE
  )
}

message("README.qmd is not newer than README.md.")

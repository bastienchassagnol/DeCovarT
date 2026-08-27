#!/usr/bin/env Rscript
# Generate inst/packages.bib (+ article/packages.bib) for R / Bioconductor
# software cited in vignettes/ and article/ — not Zotero literature.
#
# Preferred generator: knitr::write_bib() on a curated package list.
# Discovery: CRANpkg macros, library()/require()/pkg:: usage, plus an
# explicit SOFTWARE_TOOLS allowlist (prose-only mentions such as huge).
#
# Usage (from repository root):
#   Rscript scripts/curate_bibliography/generate_packages_bib.R

find_repo_root <- function() {
  candidates <- c(
    ".",
    "..",
    file.path("..", ".."),
    file.path("..", "..", "..")
  )
  for (cand in candidates) {
    if (file.exists(file.path(cand, "DESCRIPTION"))) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    "Run from the DeCovarT repository root ",
    "(or scripts/curate_bibliography/).",
    call. = FALSE
  )
}

root <- find_repo_root()

vignette_dir <- file.path(root, "vignettes")
article_dir <- file.path(root, "article")
out_inst <- file.path(root, "inst", "packages.bib")
out_article <- file.path(root, "article", "packages.bib")

# Packages discussed as tools in prose / tables (not always library()-loaded).
SOFTWARE_TOOLS <- c(
  "huge",
  "BDgraph",
  "bnlearn",
  "Matrix",
  "igraph",
  "RSpectra",
  "numDeriv",
  "mvtnorm",
  "compositions",
  "marqLevAlg",
  "ComplexHeatmap",
  "MixSim",
  "clusterGeneration",
  "limSolve",
  "glmnet",
  "AlgDesign",
  "GenSA",
  "GA",
  "mco",
  "expm",
  "SimBu"
)

# GitHub / non-CRAN Manual entries appended after knitr::write_bib().
EXTRA_MANUALS <- c(
  paste(
    "@Manual{R-pipeML,",
    "  title = {pipeML: A flexible and modular machine learning framework",
    "designed to support leakage-free model training through custom",
    "cross-validation fold construction},",
    "  author = {Marcelo Hurtado and Vera Pancaldi},",
    "  year = {2026},",
    "  note = {R package version 0.0.1},",
    "  url = {https://github.com/VeraPancaldiLab/pipeML},",
    "}",
    sep = "\n"
  ),
  paste(
    "@Manual{R-skpr,",
    "  title = {skpr: Design of Experiments Suite: Generate and Evaluate",
    "Optimal Designs},",
    "  author = {Tyler Morgan-Wall and George Khoury},",
    "  year = {2025},",
    "  note = {R package version 1.9.2},",
    "  url = {https://CRAN.R-project.org/package=skpr},",
    "}",
    sep = "\n"
  ),
  paste(
    "@Manual{R-SimBu,",
    "  title = {SimBu: Bias-aware simulation of bulk RNA-seq data with",
    "variable cell-type composition},",
    "  author = {Alexander Dietrich},",
    "  year = {2024},",
    "  note = {Bioconductor package;",
    "  https://bioconductor.org/packages/SimBu},",
    "  url = {https://doi.org/10.18129/B9.bioc.SimBu},",
    "}",
    sep = "\n"
  )
)

# Never treat these as citable analysis tools for packages.bib.
OMIT <- c(
  "DeCovarT",
  "grateful",
  "knitr",
  "rmarkdown",
  "quarto",
  "testthat",
  "devtools",
  "usethis",
  "roxygen2",
  "pkgdown",
  "lintr",
  "precommit",
  "styler",
  "air",
  "renv",
  "stats",
  "utils",
  "methods",
  "graphics",
  "grDevices",
  "datasets",
  "tools",
  "grid",
  "parallel",
  "splines",
  "tcltk",
  "compiler"
)

read_text <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

extract_matches <- function(text, pattern, perl = FALSE) {
  m <- gregexpr(pattern, text, perl = perl)
  raw <- regmatches(text, m)[[1]]
  if (length(raw) == 0L) {
    return(character())
  }
  raw
}

discover_from_files <- function(paths) {
  found <- character()
  for (path in paths) {
    if (!file.exists(path)) {
      next
    }
    text <- read_text(path)

    cran <- extract_matches(text, "\\\\CRANpkg\\{[^}]+\\}")
    if (length(cran) > 0L) {
      found <- c(found, gsub("^\\\\CRANpkg\\{|\\}$", "", cran))
    }

    libs <- extract_matches(
      text,
      "(?:library|require|requireNamespace)\\(\\s*[\"']?[A-Za-z][A-Za-z0-9.]+"
    )
    if (length(libs) > 0L) {
      found <- c(
        found,
        sub(
          "^(?:library|require|requireNamespace)\\(\\s*[\"']?",
          "",
          libs
        )
      )
    }

    dbl <- extract_matches(
      text,
      "(?<![[:alnum:]._])([A-Za-z][A-Za-z0-9.]*)::",
      perl = TRUE
    )
    if (length(dbl) > 0L) {
      found <- c(found, sub("::$", "", dbl))
    }
  }
  unique(found)
}

vignette_files <- list.files(
  vignette_dir,
  pattern = "\\.(qmd|Rmd|R)$",
  full.names = TRUE
)
article_files <- list.files(
  article_dir,
  pattern = "\\.(tex|R|qmd)$",
  full.names = TRUE
)

discovered <- discover_from_files(c(vignette_files, article_files))
# Normalise known CRAN spelling variants from prose / TeX macros.
PKG_ALIASES <- c(mixSim = "MixSim")
discovered <- ifelse(
  discovered %in% names(PKG_ALIASES),
  PKG_ALIASES[discovered],
  discovered
)
packages <- sort(unique(c(SOFTWARE_TOOLS, discovered)))
packages <- setdiff(packages, OMIT)

installed <- rownames(installed.packages())
missing <- setdiff(packages, installed)
if (length(missing) > 0L) {
  message(
    "Installing missing software packages for citation metadata:\n  ",
    paste(missing, collapse = ", ")
  )
  tryCatch(
    {
      install.packages(
        missing,
        repos = "https://cloud.r-project.org",
        quiet = TRUE
      )
    },
    error = function(e) {
      message("install.packages error: ", conditionMessage(e))
    }
  )
  installed <- rownames(installed.packages())
  still_missing <- setdiff(packages, installed)
  if (length(still_missing) > 0L) {
    warning(
      "Could not install; omitting from packages.bib: ",
      paste(still_missing, collapse = ", "),
      call. = FALSE
    )
    packages <- setdiff(packages, still_missing)
  }
}

to_write <- c("base", sort(packages))
message("Writing citations for: ", paste(to_write, collapse = ", "))

dir.create(dirname(out_inst), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_article), recursive = TRUE, showWarnings = FALSE)

header <- paste(
  "% Auto-generated R / Bioconductor software bibliography.",
  "% Do NOT edit by hand — regenerate with:",
  "%   Rscript scripts/curate_bibliography/generate_packages_bib.R",
  "% Literature (Zotero) lives in decovart_library.bib / REFERENCES.bib.",
  "% Quarto YAML should list both:",
  "%   bibliography:",
  "%     - ../inst/packages.bib",
  "%     - ../inst/REFERENCES.bib",
  "%",
  sep = "\n"
)

tmp <- tempfile(fileext = ".bib")
knitr::write_bib(to_write, file = tmp)
body <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
if (length(EXTRA_MANUALS)) {
  body <- paste(c(body, EXTRA_MANUALS), collapse = "\n\n")
}
# Exactly one trailing newline (end-of-file-fixer / sync-safe).
text <- paste0(header, "\n\n", body)
text <- paste0(sub("\n+$", "", text), "\n")
write_bib_copy <- function(path) {
  raw <- charToRaw(enc2utf8(text))
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(raw, con)
}
write_bib_copy(out_inst)
write_bib_copy(out_article)

keys <- gsub(
  "^@[^{]+\\{([^,]+),.*$",
  "\\1",
  grep("^@", readLines(out_inst), value = TRUE)
)
message("Wrote ", out_inst, " and ", out_article)
message(
  "Citation keys (",
  length(keys),
  "):\n  ",
  paste(keys, collapse = ", ")
)

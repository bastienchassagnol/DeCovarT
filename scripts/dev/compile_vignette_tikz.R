# Pre-compile the `{.tikz}` blocks of vignettes/*.qmd to SVG.
#
# Workflow: each diagram lives in its vignette inside `::: {.content-hidden}`,
# with one `%%| filename: <name>` line; the vignette shows
# `figures/tikz/<name>.svg`. CI has no TeX, so the SVGs are committed.
# Run from the repository root after editing a block or the preamble:
#   Rscript scripts/dev/compile_vignette_tikz.R [--force]
# Needs LuaLaTeX with tikz-bayesnet, Ghostscript (gs) and pdftocairo.
# Blocks whose sha1 (body + preamble) matches <name>.sha1 are skipped.

force <- "--force" %in% commandArgs(trailingOnly = TRUE)
out_dir <- normalizePath("vignettes/figures/tikz", mustWork = TRUE)
preamble_file <- file.path(out_dir, "preamble.tex")
preamble <- readLines(preamble_file)

lualatex <- local({
  has_bayesnet <- function(bin) {
    kpse <- file.path(dirname(bin), "kpsewhich")
    found <- if (file.exists(kpse)) {
      suppressWarnings(system2(
        kpse,
        "tikzlibrarybayesnet.code.tex",
        stdout = TRUE,
        stderr = FALSE
      ))
    }
    length(found) > 0L && nzchar(found[1])
  }
  candidates <- c(Sys.which("lualatex"), "/usr/bin/lualatex")
  found <- Filter(has_bayesnet, unname(candidates))
  if (length(found) == 0L) {
    stop(
      "No lualatex with the TikZ bayesnet library found; ",
      "run `tlmgr install tikz-bayesnet` or install texlive-pictures.",
      call. = FALSE
    )
  }
  found[[1]]
})

# Extract {.tikz} blocks: list of (file, filename, body).
read_blocks <- function(qmd) {
  lines <- readLines(qmd)
  starts <- which(lines == "```{.tikz}")
  lapply(starts, function(s) {
    end <- s + which(startsWith(lines[(s + 1):length(lines)], "```"))[1]
    block <- lines[(s + 1):(end - 1)]
    opts <- startsWith(block, "%%|")
    name <- sub(
      "^%%\\| filename:\\s*",
      "",
      grep("^%%\\| filename:", block, value = TRUE)
    )
    if (length(name) != 1L) {
      stop(qmd, ":", s, " needs one `%%| filename:` line", call. = FALSE)
    }
    list(file = qmd, name = trimws(name), body = block[!opts])
  })
}
blocks <- unlist(
  lapply(list.files("vignettes", "\\.qmd$", full.names = TRUE), read_blocks),
  recursive = FALSE
)
names_all <- vapply(blocks, `[[`, "", "name")
if (anyDuplicated(names_all)) {
  stop(
    "Duplicated filename: ",
    paste(unique(names_all[duplicated(names_all)]), collapse = ", "),
    call. = FALSE
  )
}

hash <- function(x) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(x, "sha1"))
  }
  tmp <- withr::local_tempfile()
  writeLines(x, tmp)
  unname(tools::md5sum(tmp))
}

compile_block <- function(b) {
  svg <- file.path(out_dir, paste0(b$name, ".svg"))
  stamp <- file.path(out_dir, paste0(b$name, ".sha1"))
  key <- hash(c(b$body, preamble))
  if (
    !force &&
      file.exists(svg) &&
      file.exists(stamp) &&
      identical(readLines(stamp), key)
  ) {
    return("up to date")
  }
  dir <- withr::local_tempdir()
  writeLines(
    c(
      "\\documentclass[tikz,border=4pt]{standalone}",
      preamble,
      "\\begin{document}",
      b$body,
      "\\end{document}"
    ),
    file.path(dir, "d.tex")
  )
  log <- system2(
    lualatex,
    c(
      "-interaction=nonstopmode",
      "-halt-on-error",
      "-output-directory",
      dir,
      file.path(dir, "d.tex")
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  pdf <- file.path(dir, "d.pdf")
  if (!file.exists(pdf)) {
    return(paste("FAILED\n", paste(utils::tail(log, 15), collapse = "\n ")))
  }
  # Outline fonts so the SVG does not depend on installed fonts.
  system2(
    "gs",
    c(
      "-q",
      "-dNOPAUSE",
      "-dBATCH",
      "-dNoOutputFonts",
      "-sDEVICE=pdfwrite",
      paste0("-sOutputFile=", file.path(dir, "o.pdf")),
      pdf
    )
  )
  status <- system2("pdftocairo", c("-svg", file.path(dir, "o.pdf"), svg))
  if (status != 0L) {
    return("FAILED (pdftocairo)")
  }
  writeLines(key, stamp)
  "compiled"
}

results <- vapply(blocks, compile_block, "")
for (i in seq_along(blocks)) {
  cat(sprintf("%-40s %s\n", blocks[[i]]$name, results[[i]]))
}
if (any(startsWith(results, "FAILED"))) {
  quit(status = 1L)
}

if (identical(Sys.getenv("GITHUB_PAT"), "")) {
  pat <- tryCatch(
    system("gh auth token", intern = TRUE),
    error = function(e) character(0)
  )
  if (length(pat) == 1L && nzchar(pat)) {
    Sys.setenv(GITHUB_PAT = pat)
  }
}

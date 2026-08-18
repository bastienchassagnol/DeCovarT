# Hunspell extras live in inst/WORDLIST: proper nouns, method/package
# names, and domain terms missing from en-GB. Do not add LaTeX commands,
# math indices, or URL fragments; prune after update_wordlist().
if (requireNamespace("spelling", quietly = TRUE)) {
  spelling::spell_check_test(
    vignettes = TRUE,
    error = FALSE,
    skip_on_cran = TRUE
  )
}

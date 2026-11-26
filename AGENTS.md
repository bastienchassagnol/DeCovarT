# Agent instructions

Standing guidance is `.cursor/rules/` (Cursor) and `.claude/rules/` (Claude Code): copies of the ubuntu-genesis bundle rules, refreshed by `agent-sync install-project` (section below).

## Package conventions (DeCovarT)

Scope: no unrelated tasks (docs, commits, refactors) unless asked; new non-trivial logic still gets one small test (`ponytail`), and exported functions their smoke tests. British English in prose. Repository-specific, rarely needed material is in `dev/package-dev-notes.md` (CRAN incoming lessons, CRAN Cookbook, GitHub Actions, pkgdown, load hooks, Roxygen, testing, plotting, bibliography streams, Nature Methods limits, scientific review). Read the relevant section on demand. Generic R, Python, Quarto and writing guidance comes from the bundle rules.

- **Layout.** `R/` (package API, numeric prefixes `02_01_…`, keep the convention), `scripts/` (one-off analyses), `tests/testthat/` (`R/foo.R` maps to `test-foo.R`), `vignettes/` (Quarto), `inst/` (assets, `CITATION`, bib), `article/` (LaTeX and TikZ manuscript), `data/` (`raw/` immutable, DVC pointers only), `output/`, `docs/` (pkgdown build, gitignored).
- **Package code.** R only in `R/`; Python only in `scripts/`. `R/` and `NAMESPACE` stay ASCII. Document with roxygen2 (`devtools::document()`); a documented function carries `@export`, and `@noRd` only on undocumented helpers.
- **Checks before a push** touching package code: `devtools::check()` and `pkgcheck::pkgcheck()`, warnings included; a local `pkgdown::build_site()` when exports, `_pkgdown.yml`, vignettes or docs workflows change. Before a CRAN submission follow `dev/package-dev-notes.md` (CRAN incoming section).
- **pre-commit.** `pre-commit clean`, stage everything, `pre-commit run`, re-stage autofixes, commit. Autofix hooks fail once by design. No `--no-verify` unless asked. The project does not use `renv`: the R hooks are local `language: system` hooks running the scripts of the `precommit` R package (`scripts/dev/precommit_hook.sh`), so that package must be installed. `slides/_extensions/` is excluded from the whitespace hooks.
- **Diagrams (TikZ).** Sources stay in the vignette inside `::: {.content-hidden}`; the SVG is pre-compiled locally with `Rscript scripts/dev/compile_vignette_tikz.R` (needs LuaLaTeX with tikz-bayesnet, `gs`, `pdftocairo`) and committed under `vignettes/figures/tikz/` (preamble: `vignettes/figures/tikz/preamble.tex`). CI never compiles TikZ. If CI ever gets a TeX distribution (e.g. a cached TinyTeX with tikz-bayesnet), the pandoc-ext/diagram filter with `TIKZ_BIN` pointing at `scripts/dev/tikz-lualatex.sh` could compile the blocks at render time instead.
- **Plots.** Faceted `ggplot2` figures use `theme_decovart_facets()`.
- **Two citation streams, never mixed.** Literature: Zotero export into `inst/REFERENCES.bib` and `article/decovart_library.bib`, checked with `scripts/curate_bibliography/check_bibliography.py --clean`. Software: `packages.bib` from `scripts/curate_bibliography/generate_packages_bib.R`. Never invent citekeys, DOIs or BibTeX.
- **Manuscript (`article/`)**: Nature Methods Article: abstract at most 150 words and unreferenced, main text about 3,000 words, no footnotes, at most 50 references, at most 6 display items. Details in `dev/package-dev-notes.md`.

## Rules and skills: where they live

Source of truth: `Config/.cursor/` in
[ubuntu-genesis](https://github.com/bastienchassagnol/ubuntu-genesis).

```yaml
# read by agent-sync install-project
packs:
  - r-package-article
project_rules: copy
```

| Level | Rules | Skills | Written by |
| --- | --- | --- | --- |
| 1. Project | `.cursor/rules/*.mdc`, `.claude/rules/*.md` | none | `agent-sync install-project` (bundle copies, committed) |
| 2. User | `~/.cursor/rules`, `~/.claude/rules` | `~/.cursor/skills`, `~/.claude/skills` | `agent-sync install-user`; do not hand-edit |
| 3. Machine-local | none | `firecrawl*`, `~/.agents/skills` | not vendored |
| 4. Host built-in | none | `~/.cursor/skills-cursor`, plugins, `~/.claude/skills/synced/` | the tool itself |

- Project and user copies of a rule are the same bundle file. Claude Code loads
  both levels together (no winner), Cursor prefers the project copy; keeping
  them identical avoids any conflict. Refresh after a bundle change:
  `"$AGENT_BUNDLE_ROOT/agent-sync/agent-sync.sh" install-project .`
  (copies the pack into both folders, prunes rules removed from the bundle,
  keeps project-only `.cursor/rules` files and mirrors them to `.claude/rules`).
- Project-only rules go in `.cursor/rules/` under a name the bundle does not use.
- Claude Code reads `paths:` in a rule's frontmatter; Cursor reads `globs:` and
  `alwaysApply:`. Bundle rules carry both.
- Symbolic verification lives in `answer-cardamomOT/verification/`
  (`run_all.py`, `ledger.json`, `verification-ledger.tex`).
- Do not commit a symlink that leaves the repository (it breaks on clone, in CI
  and on Windows); this repository tracks none.

## Code graph (CodeGraph first, Graphify second)

CodeGraph is the default for code. Graphify is secondary and applies to mixed material. Routing follows the `code-intelligence` rule; this section adds the repository paths.

- **Code questions** (`R/`, `tests/testthat/`, `scripts/`, `answer-cardamomOT/verification/`: callers, callees, impact, how a call path reaches a function): where `.codegraph/` exists at the repository root, call `codegraph_explore` (MCP) or `codegraph explore "<question>"` before Grep, Glob or Read. Treat the returned source as already read. The index auto-syncs while the MCP server runs. `.codegraph/` is not committed; build it once per checkout with `codegraph init`. `codegraph.json` excludes `slides/_extensions/`, `vignettes/_extensions/` and `revdep/`.
- **Cross-document questions** (links across vignettes, `article/`, papers and design notes): use `graphify query "<q>"` only where `graphify-out/graph.json` exists. The folder is not committed and is rebuilt at milestones (annotated tag or version bump, via `Config/.cursor/code-knowledge/scripts/graphify-milestone.sh` in the ubuntu-genesis bundle), not on every commit. Do not run `graphify hook install` for code sync: CodeGraph covers it. `.graphifyignore` already excludes generated and vendored paths (rendered output, vendored extensions).
- **Neither index exists**: for a code question, tell the user to run `codegraph init` and wait; for a cross-document question, tell them to run the `/graphify` skill from the repository root. Do not substitute a long manual exploration.

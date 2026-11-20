# Agent instructions

Standing guidance is `.cursor/rules/` (Cursor) and `.claude/rules/` (Claude Code). The shared rules are also deployed from the ubuntu-genesis bundle at the user level; the levels below explain which copy applies.

## Package conventions (DeCovarT)

Scope: do not add tasks (tests, docs, commits, refactors) unless asked; ask first. British English in prose. The deprecated `.cursorrules` was removed; its repository-specific, rarely needed material is in `dev/package-dev-notes.md` (CRAN incoming lessons, CRAN Cookbook, GitHub Actions, pkgdown, load hooks, Roxygen, testing, plotting, bibliography streams, Nature Methods limits, scientific review). Read the relevant section on demand. Generic R, Python, Quarto and writing guidance comes from the bundle rules.

- **Layout.** `R/` (package API, numeric prefixes `02_01_…`, keep the convention), `scripts/` (one-off analyses), `tests/testthat/` (`R/foo.R` maps to `test-foo.R`), `vignettes/` (Quarto), `inst/` (assets, `CITATION`, bib), `article/` (LaTeX and TikZ manuscript), `data/` (`raw/` immutable, DVC pointers only), `output/`, `docs/` (pkgdown build, gitignored).
- **Package code.** R only in `R/`; Python only in `scripts/`. `R/` and `NAMESPACE` stay ASCII. Document with roxygen2 (`devtools::document()`); a documented function carries `@export`, and `@noRd` only on undocumented helpers.
- **Checks before a push** touching package code: `devtools::check()` and `pkgcheck::pkgcheck()`, warnings included; a local `pkgdown::build_site()` when exports, `_pkgdown.yml`, vignettes or docs workflows change. Before a CRAN submission follow `dev/package-dev-notes.md` (CRAN incoming section).
- **pre-commit.** `pre-commit clean`, stage everything, `pre-commit run`, re-stage autofixes, commit. Autofix hooks fail once by design. No `--no-verify` unless asked. The project does not use `renv`; lintr hooks set `RENV_CONFIG_AUTOLOADER_ENABLED=FALSE`. `slides/_extensions/` is excluded from the whitespace hooks.
- **Plots.** Faceted `ggplot2` figures use `theme_decovart_facets()`.
- **Two citation streams, never mixed.** Literature: Zotero export into `inst/REFERENCES.bib` and `article/decovart_library.bib`, checked with `scripts/curate_bibliography/check_bibliography.py --clean`. Software: `packages.bib` from `scripts/curate_bibliography/generate_packages_bib.R`. Never invent citekeys, DOIs or BibTeX.
- **Manuscript (`article/`)**: Nature Methods Article: abstract at most 150 words and unreferenced, main text about 3,000 words, no footnotes, at most 50 references, at most 6 display items. Details in `dev/package-dev-notes.md`.

## Skills and rules: where they live, which one wins

Source of truth: `Config/.cursor/` in
[ubuntu-genesis](https://github.com/bastienchassagnol/ubuntu-genesis), cloned to
`~/repos/ubuntu-genesis`. `agent-sync.sh` copies it to the user level of both tools.

| Level | Rules | Skills | Written by |
| --- | --- | --- | --- |
| 1. Project | `.cursor/rules/*.mdc`, `.claude/rules/*.md` | `.cursor/skills`, `.claude/skills` (none yet) | this repository, committed |
| 2. User | `~/.cursor/rules`, `~/.claude/rules` | `~/.cursor/skills`, `~/.claude/skills` | `agent-sync.sh install-user`; do not hand-edit |
| 3. Machine-local | none | `firecrawl*`, `~/.agents/skills` (`npx skills`) | not vendored, not in git |
| 4. Host built-in | none | `~/.cursor/skills-cursor`, plugins, `~/.claude/skills/synced/` | the tool itself |

The tools disagree on a name clash, so a name is defined at one level only:

- Cursor merges rules with the project over the user level (Team, Project, User).
- Claude Code loads user and project rules together and does not define a winner.
  For skills of the same name, a personal skill (`~/.claude/skills`) beats a project skill.
- Claude Code reads only `paths:` in a rule's frontmatter; Cursor reads `globs:` and
  `alwaysApply:`. Bundle rules carry both. A rule without either loads in every session.
- Claude Code does not read `.cursor/` or `.agents/`.

Project copies in `.claude/rules/` and `.cursor/rules/` have the same names as bundle
rules and may differ from them. Read the project copy as the authority for this
repository. Do not add a project skill whose name the bundle already uses.

### Keeping the levels in step

- Edit in the bundle, commit there, `git pull`; the `post-merge` hook runs `install-user`.
- Edit under `~/.cursor` or `~/.claude`: the `afterFileEdit` (Cursor) and `PostToolUse`
  (Claude Code) hook, `~/.cursor/hooks/agent-sync-after-edit.sh`, runs `publish`
  and redeploys. Review `git diff` in the bundle, then commit.
- Log: `~/.cache/agent-sync/hook.log`. Check: `$AGENT_BUNDLE_ROOT/agent-sync/agent-sync.sh verify`,
  or open a shell with `CLAUDE_AGENT_SYNC=1` (set in `~/.zshrc`; output in `verify.log`).
- `AGENT_BUNDLE_ROOT` is exported in `~/.zshrc` and in the `env` block of
  `~/.claude/settings.json`. Hook shells do not read `~/.zshrc`.

### Symlinks

- Do not commit a symlink that leaves the repository (for example `.cursor/rules`
  pointing at a bundle outside it): it breaks on clone, in CI and on Windows, and a
  symlinked `skills/` puts every bundle skill at project level. Use copies.
  A relative link inside the repository is acceptable.
- `~/.bashrc` is a separate file that hands interactive shells to zsh. Do not link it
  to `~/.zshrc`: that file uses zsh-only syntax (`setopt`, `typeset -g`).
- Audit of `dtoo_project/` on 2026-11-18: this repository tracks no symlink.

## Code graph (CodeGraph first, Graphify second)

CodeGraph is the default for code. Graphify is secondary and applies to mixed material. Where this section and the `graphify` rule (Cursor `.cursor/rules/graphify.mdc`, Claude Code `~/.claude/rules/graphify.md`) disagree about a code question, this section applies.

- **Code questions** (`R/`, `tests/testthat/`, `scripts/`, `answer-cardamomOT/verification/`: callers, callees, impact, how a call path reaches a function): where `.codegraph/` exists at the repository root, call `codegraph_explore` (MCP) or `codegraph explore "<question>"` before Grep, Glob or Read. Treat the returned source as already read. The index auto-syncs while the MCP server runs. `.codegraph/` is not committed; build it once per checkout with `codegraph init`. `codegraph.json` excludes `slides/_extensions/`, `vignettes/_extensions/` and `revdep/`.
- **Cross-document questions** (links across vignettes, `article/`, papers and design notes): use `/graphify query` only where `graphify-out/graph.json` exists. The folder is not committed and is rebuilt at milestones (annotated tag or version bump, via `Config/.cursor/code-knowledge/scripts/graphify-milestone.sh` in the ubuntu-genesis bundle), not on every commit. Do not run `graphify hook install` for code sync: CodeGraph covers it. `.graphifyignore` already excludes generated and vendored paths (rendered output, vendored extensions).
- **Neither index exists**: for a code question, tell the user to run `codegraph init` and wait; for a cross-document question, tell them to run the `/graphify` skill from the repository root. Do not substitute a long manual exploration.

---
description: Working scope, data layout and host settings for R/Python projects
alwaysApply: true
---
# R/Python projects: scope and data

Package work adds `r-package-development`; language details are in
`r-development` and `python-development`.

## Scope

- British English unless the repository or user uses another variant.
- Do what was asked. No unrelated refactors, features, dependencies, files,
  analyses, tests, docs or commits. New non-trivial logic still gets its one
  small test or assert (`ponytail`), and package exports their smoke tests.
- Inspect the repository first and follow its architecture, naming,
  formatting, testing, documentation and dependency conventions.
- Do not invent files, APIs, functions, results or package behaviour.
  Label claims as verified, inferred, recommended or uncertain.
- For unstable technical claims, check current primary documentation; drop a
  workaround once its reason is gone.
- Keep backwards compatibility unless a break is requested or justified.
- Never weaken tests, validation, statistical safeguards or error handling to
  pass a check. Never silently discard user data, citations, metadata or
  meaningful output.
- Before finishing: inspect the diff, run the narrowest relevant checks, and
  report any check that could not run.

## Repository, files and data

After [Data Science Team Training — Data Organisation](https://dstt.stephenturner.us/dataorg.html#additional-best-practices):

- One project = one working directory. Raw data is immutable; derive cleaned
  data with scripts or pipelines.
- Keep generated, cache, credential and machine-local files out of git; DVC or
  equivalent for large binaries that need history.
- Plain-text artefacts (CSV/TSV/Parquet) over `.xlsx`.
- No `final`, `final2`, `latest`: use `snake_case`, `YYYY-MM-DD-slug`, and
  zero-padded ordered scripts.
- Keep a codebook for important datasets (meaning, units, types, missing
  coding, provenance, transformations).
- Do not silently repair implausible data; report it and script deterministic
  rules.

## Laptop and cluster (Cursor / VS Code over SSH)

User settings are reused on SSH remotes, but interpreter and terminal paths
differ. Precedence: User → Remote (per SSH host) → workspace
`.vscode/settings.json`. Keep host paths (`r.rterm.linux`, `r.rpath.linux`,
`python.defaultInterpreterPath`) in User or Remote settings, never in
committed workspace settings. Never commit secrets or credentials.

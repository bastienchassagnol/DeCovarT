---
description: Zotero literature workflow, curation scripts, software citations and web sources
globs: "**/*.{bib,qmd,Rmd,tex},**/_quarto.yml,**/CITATION*"
alwaysApply: false
paths:
  - "**/*.{bib,qmd,Rmd,tex}"
  - "**/_quarto.yml"
  - "**/CITATION*"
---
# Bibliography and citations

Quarto labels, callouts and cross-references: `quarto-presentations`. Never
invent BibTeX, citekeys, DOIs or quotations (`scientific-writing`).

## Literature: Zotero → citation keys → Git

- **Source of truth:** the author's
  [Zotero library](https://www.zotero.org/bastienchassagnol/library)
  (library ID `7172053`). Add in Zotero, then export.
- **Citation keys:** Zotero 8+ native `citationKey` (`fields.fieldName =
  'citationKey'` in `zotero.sqlite`; formula under Settings → Cite → Citation
  Key). Better BibTeX is optional. Open the DB with `immutable=1` when Zotero
  holds a lock.
- **Collections:** one leaf collection per item (the most specific); never
  parent and child.
- **zotcurate (`zotc`)** bridges `.bib`/`.qmd`/CSV keys and collections:
  `~/.local/bin/zotc` wraps `-z ~/Zotero/zotero.sqlite`; credentials in
  `~/.zotc/` (chmod 600, never committed).

Workflow:

1. Zotero → export `.bib` → curate (below).
2. Git/docs → Zotero: `zotc keys extract`, then `collection add|create`
   (dry run first).
3. Software → `packages.bib` (or the repo equivalent), never the Zotero
   literature library.

## Curation scripts (`scripts/curate_bibliography/`, uv project)

They live in the package repo (DeCovarT); other repos vendor a copy or call
the sibling repo with `uv run --project <path>/scripts/curate_bibliography`.

```bash
uv run --project scripts/curate_bibliography \
  python scripts/curate_bibliography/bib_enrich.py PATH.bib --write --backup \
  --min-score 92 --report scripts/curate_bibliography/bib-enrichment.json
npx bibtex-tidy@1.14.0 PATH.bib --modify --curly --numeric --sort=key
```

- `bib_enrich.py` fills missing fields only; never overwrite verified fields
  without review. The report is gitignored.
- `check_bibliography.py --clean` (field whitelist, DOI dedupe, paired
  `article/` and `inst/` copies) only in package layouts; never on a book or
  website bibliography.
- `add_zotero_perspectives.py`: one-off Web API inserts, one leaf collection
  each.
- Required fields per type follow `check_bibliography.py` (`@article`: title,
  author, year, journal; `@book`: title, year, publisher, author or editor).
  Allowed fields when a whitelist applies: `title`, `author`, `editor`, `year`,
  `journal`, `booktitle`, `publisher`, `howpublished`, `doi`, `volume`,
  `number`, `pages`, `isbn`, `issn`, `edition`, `school`, `institution`, `url`.
- Without a repo copy: `uv venv ~/.venvs/bib-enrich` with `bibtexparser httpx
  rapidfuzz typer rich`.

## Software citations

- Generate from official metadata (`citation()`, `CITATION`,
  `knitr::write_bib()`); cite software that materially supports methods or
  reproducibility, not every transitive dependency.
- Split bibliographies in Quarto:

```yaml
bibliography:
  - path/to/packages.bib
  - path/to/REFERENCES.bib
link-citations: true
```

## Web resources

Prefer DOI, ISBN or official package pages. Try the live URL first; use the
[Wayback Machine](https://web.archive.org/) for dead links. Never fabricate
archive timestamps.

---
description: Model routing, answer caches and prefix stability to save tokens
alwaysApply: true
---
# Token budget: routing, caches, smaller models

Short, always-on. Code lookups: `code-intelligence`; documents:
`document-extraction`.

## Route each task to the cheapest model that can do it

| Task | Model |
| --- | --- |
| Search, excerpting, file lookup, renaming, formatting, boilerplate | Haiku: Explore subagent (or a local `doc-retriever` agent when installed), or `/model haiku` |
| Routine code and prose edits with a clear spec | Sonnet |
| Design, statistical review, debugging across files, final synthesis | Opus |

- Delegate reading to a Haiku subagent and keep its brief, not the files.
- `/effort low` for mechanical edits; raise effort only for reasoning.
- Cursor has no per-subagent model: switch to a small model for search turns.

## Reuse before recomputing

- Documents: `doc-search search` checks cached answers (notes) first.
- After an answer has been checked against the excerpts it cites, save it:
  `doc-search note add --source F --question Q --answer A --pages P`.
  Never save an answer that was not checked. Re-indexing a changed file
  drops its notes.
- Reuse `.firecrawl/extract/<slug>/content.md` and `figures.json` when they
  exist.

## Keep the cached prefix stable

Do not edit rules, `CLAUDE.md`, or the MCP server list during a session: the
prompt-cache prefix changes and the next turn is billed in full. Use `/clear`
between unrelated tasks; run `/context` before `/compact`.

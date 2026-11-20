---
description: Route code vs corpus questions to CodeGraph or Graphify
alwaysApply: true
---
# Code intelligence (solo repos)

Default stack for individual repositories: **CodeGraph** for source structure;
**Graphify** for mixed corpora at **milestones** only—not on every commit.

Do **not** default to pxpipe or Codebase Memory MCP unless the user opts in.

## Routing

| Question | Prefer | When |
| --- | --- | --- |
| Callers, callees, impact, “how does X reach Y” in source | `codegraph_explore` (MCP) or `codegraph explore` | `.codegraph/` exists (`codegraph init`) |
| Links across Quarto, papers, PDFs, design notes, multi-format repos | `/graphify query` or the graphify skill | `graphify-out/graph.json` exists |
| Neither index | Targeted search, LSP, small reads; suggest `codegraph init` or a milestone Graphify run | — |

## Graphify boundaries

- **Do not** install `graphify hook install` or `graphify claude install` for
  always-on rebuilds. CodeGraph auto-syncs code via MCP.
- Run a full or incremental Graphify build on **git tags** or **R `DESCRIPTION`
  version bumps** via [`Config/.cursor/code-knowledge/`](../code-knowledge/README.md)
  scripts when the user releases.
- For code-only exploration with a live CodeGraph index, do not rebuild Graphify.

## Claude and Cursor

User skills live under `~/.claude/skills/` and `~/.cursor/skills/` after
`agent-sync.sh install-user`. Deploy from git in `Config/.cursor/`, not from
project symlinks alone.
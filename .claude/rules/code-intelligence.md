---
description: Route code questions to CodeGraph and cross-document questions to Graphify
alwaysApply: true
---
# Code intelligence

| Question | Tool | Needs |
| --- | --- | --- |
| Callers, callees, impact, "how does X reach Y" in source | `codegraph_explore` (MCP) or `codegraph explore` (CLI), before Grep, Glob or Read | `.codegraph/` (`codegraph init`) |
| Links across Quarto chapters, vignettes, papers, PDFs, design notes | `graphify query "<q>"`, `graphify path "<A>" "<B>"`, `graphify explain "<concept>"` | `graphify-out/graph.json` |
| Neither index | Targeted search, LSP, small reads; suggest `codegraph init` | — |

- Never send code questions to Graphify, and never rebuild the graph for one.
- With a Graphify index, read `graphify-out/wiki/index.md` when it exists
  before raw files; `GRAPH_REPORT.md` only for broad architecture review.
  Switch to Read/Grep once oriented.
- Graphify is built at milestones only (annotated git tag, R `DESCRIPTION`
  version bump) with the `code-knowledge` scripts of the agent bundle. No
  `graphify hook install`, `graphify claude install` or `--watch`; run
  `graphify update .` (AST only) only when asked.
- No pxpipe or Codebase Memory MCP unless the user opts in.
- Repeat this routing in every subagent prompt that explores a repository.

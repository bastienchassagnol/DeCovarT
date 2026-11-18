# Agent instructions

Standing guidance is `.cursor/rules/` (Cursor) and `.claude/rules/` (Claude Code). The shared rules are also deployed from the ubuntu-genesis bundle at the user level; the levels below explain which copy applies.

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

## Code graph (graphify)

The `graphify` rule applies only where `graphify-out/` exists, and the folder is not committed. If `graphify-out/graph.json` is missing, tell the user to run the `/graphify` skill from the repository root before exploring the code, and wait for the build. Do not substitute a long manual exploration. `.graphifyignore` already excludes generated and vendored paths (`renv/`, rendered output, vendored extensions).

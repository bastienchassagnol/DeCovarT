# Symbolic verification of equations (SymPy + Lean)

Use this rule whenever a report, article, or derivation states an equation that
is new, rearranged, or copied from a source and then used. The model writes the
check; the engine decides. This is the neuro-symbolic loop of SymCode
([arXiv:2510.25975](https://arxiv.org/abs/2510.25975)): no equation enters the
text as "verified" on the model's say-so.

## Engines and versions

| Engine | Use for | Pin |
| --- | --- | --- |
| SymPy via `uv run --script` | every labelled equation: identities, limits, series, ranks, numeric evaluation | `sympy==1.14.0` in PEP 723 metadata |
| mpmath (ships with SymPy) | high-precision numerics: matrix exponentials, quadrature, pmf sums | 40 to 50 digits |
| Lean 4 + Mathlib | field/ring identities and elementary inequalities that carry an argument | `lean-toolchain` of the build project |

Each repository keeps its own `verification/` tree. The Lean build stays
**outside OneDrive** (`.lake` is about 7 GB). The tree
`~/lean-projects/dynamic_pdes_verification`, with canonical source
`verification/lean/TEOU/Lemmas.lean`, belongs to that dynamic-PDE project.
Do not reuse it for another repository. `verification/run_all.py` in the
repository that owns the proof copies the Lean sources into its build project
before `lake build`.

## Workflow

0. Before writing a check, compare the equation's symbols against the shared
   notation file (`vignettes/_notation.qmd` in DeCovarT). A check is written
   in that notation, so a symbol clash is caught before the ledger runs.
1. Give every display equation a `\label{eq:...}`. Each claim gets a ledger id
   (`T1a`, `O2e`, ...) and a check in `verification/sympy/check_<group>_*.py`.
2. Derive, do not copy. Recompute the quantity from first principles (a CGF,
   a generator, a Jacobian) and compare it with the printed formula. A check
   that only re-types the paper's formula proves nothing.
3. Test in this order, and record which one passed (`kind`):
   - `exact`: `sympy.simplify(lhs - rhs) == 0`, else `cancel(together(expand(...)))`;
   - `series`: for "approximately" claims, expand in a bookkeeping parameter
     (`u -> eps*u`) and compare the stated order only;
   - `rational-sampling`: random positive rationals, 25 points, fixed seed —
     a fallback, never reported as a proof;
   - `numeric`: mpmath at stated precision, with the tolerance in `detail`.
4. Status is `pass`, `fail`, or `flag`. Use `flag` when the engine confirms a
   **deviation from the source** (a dropped term, a typo). Report the deviation
   in the text with the source page; never silently "fix" the source.
5. Formalise in Lean only statements that are pure algebra or arithmetic
   (closed forms, Schur complements, polynomial identities after clearing
   denominators, `norm_num` facts). Prefer `linear_combination` for "these
   equations imply that one". Prefer `rw [div_eq_iff h]; ring` over `field_simp`
   when denominators are compound, because the current `field_simp` can leave
   unreduced inverses.
6. Lean acceptance: `lake build` succeeds; `rg -w sorry verification/lean` is
   empty; `#print axioms` (in `verification/lean/Axioms.lean`) lists only
   `propext`, `Classical.choice`, `Quot.sound`.
7. Run `python3 verification/run_all.py`. It writes `verification/ledger.json`
   and `report/verification-ledger.tex`; exit status 0 means no `fail` and a
   clean Lean build. Re-run after every change to an equation.

## Writing about results

- A statement is "verified" only if its ledger id has `pass` or `flag`.
- Probabilistic content (expectations, asymptotics) is verified only
  conditional on the stated inputs (for example Wishart moments, or a closure
  assumption). Say which.
- When a planned claim fails, report the rank, value, or residual the engine
  returned, and rewrite the claim. One precedent: a trajectory rank came out
  4, where the draft said 5.

## Optional tooling for later sessions

- `lean-lsp-mcp` (MIT): `claude mcp add lean-lsp -s project uvx lean-lsp-mcp`
  with `LEAN_PROJECT_PATH` set to the build project; gives goal states and
  Mathlib search to the agent.
- `cameronfreer/lean4-skills` (MIT): `/plugin marketplace add cameronfreer/lean4-skills`
  then `/plugin install lean4`; prove/review/golf workflows.
- Avoid `sympy-mcp`: unmaintained and evaluates input with `eval`.

# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""Run every SymPy check of the CardamomOT note and write the ledger.

Usage (from any directory):
    uv run --script verification/run_all.py

Writes verification/ledger.json and ../verification-ledger.tex (next to
cardamomot_decovart_scopes.tex, which inputs it in Appendix B) as one
longtable. Claims and details are plain text in which $...$ segments are
LaTeX (built with sympy.latex) and are kept verbatim; a raw '**' or an odd
number of '$' makes the record fail. Exits with status 1 if any check fails or
raises.
"""

from __future__ import annotations

import datetime as dt
import importlib.util
import json
import platform
import re
import sys
import traceback
from pathlib import Path

import numpy
import sympy

sys.dont_write_bytecode = True  # keep verification/sympy/ free of __pycache__
HERE = Path(__file__).resolve().parent
NOTE = HERE.parent
TEX = NOTE / "cardamomot_decovart_scopes.tex"
CHECK_DIR = HERE / "sympy"
VIGNETTES = NOTE.parent / "vignettes"  # package vignettes cite ids as <!-- chk: V13a -->
KEYS = ("id", "claim", "kind", "status", "detail", "tex_line")
# optional key: "detail_full", the full engine output kept in ledger.json; the
# short "detail" (plain text, with $...$ formulas from sympy.latex) goes into the
# table
KINDS = {"exact", "series", "rational-sampling", "numeric"}
STATUSES = {"pass", "fail", "flag"}


def _load(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _chk_lines(tex: Path) -> dict[str, list[int]]:
    """Map each ledger id cited in \\chk{...} to the lines that cite it."""
    found: dict[str, list[int]] = {}
    if not tex.exists():
        return found
    for i, line in enumerate(tex.read_text(encoding="utf-8").splitlines(), start=1):
        for arg in re.findall(r"\\chk\{([^}]*)\}", line):
            for tok in re.split(r"[,\s]+", arg):
                if tok:
                    found.setdefault(tok, []).append(i)
    return found


def _qmd_lines(folder: Path) -> dict[str, list[str]]:
    """Map each ledger id cited in a vignette comment <!-- chk: ... --> to
    'file-prefix:line' locations (file prefix = text before the first '-' pair,
    e.g. '02-1')."""
    found: dict[str, list[str]] = {}
    if not folder.is_dir():
        return found
    for qmd in sorted(folder.glob("*.qmd")):
        prefix = "-".join(qmd.stem.split("-")[:2])
        for i, line in enumerate(qmd.read_text(encoding="utf-8").splitlines(), start=1):
            for arg in re.findall(r"<!--\s*chk:\s*([^>]*?)\s*-->", line):
                for tok in re.split(r"[,\s]+", arg):
                    if tok:
                        found.setdefault(tok, []).append(f"{prefix}:{i}")
    return found


_TEX_ESC = {
    "\\": r"\textbackslash{}",
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}


def tex_escape(s: str) -> str:
    s = "".join(_TEX_ESC.get(ch, ch) for ch in str(s))
    return s.replace("\n", " ")


GROUPS = {
    "V1": "Linear mixed model: covariance, BLUP and Wiener filter",
    "V2": "Stein's loss minimised at the true covariance",
    "V3": "Poisson log-normal moment map",
    "V4": "Fisher information against the ANOVA ratio, $J=2$",
    "V5": "Two-gene Lyapunov covariance and precision; fork precision",
    "V6": "Single-cell (linear) against population (quadratic) covariance",
    "V7": "Maximum likelihood of the linear mixed model against DeCovarT's fit",
    "V8": "Hellinger and 2-Wasserstein distances against the Mahalanobis distance",
    "V9": "Depth and Poisson sampling in the within-type covariance",
    "V10": "Rank of the empirical covariance and existence of the Gaussian MLE",
    "V11": "Linear shrinkage is positive definite",
    "V12": "Quantile universal threshold against Bonferroni, FWER and false edges",
    "V13": "Raw or Poisson-denoised within-type covariance in the bulk",
    "V14": "Denoised within-type covariance: technical size factor and finite-sample bias",
    "V15": "Gaussian transport map (Bisque's gene-wise map as its diagonal case)",
    "V16": "Bulk sampling diagonal against the variance scale: identifiability and profile",
    "V17": "Presence gate: zero-variance and single-detection genes",
    "V18": "D-, E- and A-optimal gene panels from the Fisher information for $\\bm p$",
    "V19": "Explicit tail bound for RNA-Sieve's range filter",
    "V20": "Negative-binomial counts are sub-exponential: moments, Chernoff and Bernstein range bounds",
    "V21": "Unico: confounded variance components and spectral entropy",
    "V22": "CSNet's shared-gamma negative binomials and the factor covariance",
    "V23": "Multivariate Fisher information for $\\bm p$ and the Bhattacharyya overlap, $J=2$",
    "V24": "Lyapunov covariance, conditional precision and the triangular factor of a DAG",
    "V25": "Intrinsic noise, bulk scales and naive intervals",
    "V26": "Dynamic and mixed Gaussian Bayesian networks",
}


def render(s: str) -> str:
    """Escape text for LaTeX but keep $...$ segments (typeset with sympy.latex)
    verbatim."""
    parts = str(s).split("$")
    return "".join(
        ("$" + part + "$") if i % 2 else tex_escape(part) for i, part in enumerate(parts)
    )


def _layout_problem(r: dict) -> str | None:
    """Return a reason if the claim or detail cannot go into the table."""
    for key in ("claim", "detail"):
        txt = str(r[key])
        if "**" in txt:
            return f"raw '**' in {key}"
        if txt.count("$") % 2:
            return f"unbalanced '$' in {key}"
    return None


def write_tex(records: list[dict], path: Path, meta: dict) -> None:
    counts = {s: sum(r["status"] == s for r in records) for s in ("pass", "flag", "fail")}
    col = r">{\raggedright\arraybackslash}"
    spec = (
        "@{}"
        + col + "p{1cm}"
        + col + r"p{0.29\linewidth}"
        + col + r"p{\dimexpr0.71\linewidth-1cm-1.45cm-0.8cm-0.85cm-10\tabcolsep\relax}"
        + col + "p{1.45cm}"
        + col + "p{0.8cm}"
        + r">{\raggedleft\arraybackslash}p{0.85cm}"
        + "@{}"
    )
    header = r"ID & Claim & Detail & Kind & Status & Line \\"
    lines = [
        "% generated by verification/run_all.py; do not edit by hand",
        r"{\scriptsize",
        r"\noindent Generated "
        + tex_escape(meta["date"])
        + r" with SymPy "
        + tex_escape(meta["sympy"])
        + r", NumPy "
        + tex_escape(meta["numpy"])
        + r", Python "
        + tex_escape(meta["python"])
        + f". {len(records)} checks: {counts['pass']} pass, {counts['flag']} flag, "
        f"{counts['fail']} fail. \\enquote{{flag}} means the engine confirms a deviation "
        r"from the original wording, now corrected in the text. The text line is where "
        r"the identifier is cited; \texttt{ledger.json} keeps the full engine output of "
        r"each check (key \texttt{detail\_full}).\par\medskip",
        r"\setlength{\LTpre}{0pt}\setlength{\LTpost}{0pt}",
        r"\begin{longtable}{" + spec + "}",
        r"\toprule",
        header,
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        header,
        r"\midrule",
        r"\endhead",
        r"\midrule",
        r"\multicolumn{6}{r@{}}{\textit{continued on the next page}} \\",
        r"\endfoot",
        r"\bottomrule",
        r"\endlastfoot",
    ]
    current = None
    for r in records:
        g = re.match(r"V\d+", r["id"])
        g = g.group(0) if g else r["id"]
        if g != current:
            if current is not None:
                # space, not a rule: a rule here could meet the foot's rule at a page break
                lines.append(r"\addlinespace[4pt]")
            title = GROUPS.get(g, "")
            lines.append(
                r"\multicolumn{6}{@{}l@{}}{\textbf{"
                + tex_escape(g)
                + "}"
                + (r"\enspace\textit{" + render(title) + "}" if title else "")
                + r"} \\*"
            )
            current = g
        lines.append(
            " & ".join(
                [
                    r"\texttt{" + tex_escape(r["id"]) + "}",
                    render(r["claim"]),
                    render(r["detail"]),
                    tex_escape(r["kind"]),
                    (r"\textbf{" if r["status"] != "pass" else "{")
                    + tex_escape(r["status"])
                    + "}",
                    tex_escape(r["tex_line"]),
                ]
            )
            + r" \\"
        )
    lines += [r"\end{longtable}", "}"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    records: list[dict] = []
    # numeric order of the group number (a plain sort puts check_V10 before check_V2)
    def _group_no(p: Path) -> int:
        m = re.match(r"check_V(\d+)", p.name)
        return int(m.group(1)) if m else 10**6

    for path in sorted(CHECK_DIR.glob("check_V*.py"), key=lambda p: (_group_no(p), p.name)):
        try:
            res = _load(path).run()
        except Exception:  # a crashing check is a failed check
            res = [
                {
                    "id": path.stem,
                    "claim": "check raised an exception",
                    "kind": "exact",
                    "status": "fail",
                    "detail": traceback.format_exc(limit=3),
                    "tex_line": "",
                }
            ]
        for r in res:
            r.setdefault("detail_full", r.get("detail"))
            missing = [k for k in KEYS if k not in r]
            if missing or r["kind"] not in KINDS or r["status"] not in STATUSES:
                r["status"] = "fail"
                r["detail"] = f"malformed record (missing {missing})"
            problem = None if missing else _layout_problem(r)
            if problem:
                r["status"] = "fail"
                r["detail"] = f"malformed record ({problem})"
            r["script"] = f"verification/sympy/{path.name}"
        records += res

    cited = _chk_lines(TEX)
    cited_qmd = _qmd_lines(VIGNETTES)
    for key, locs in _qmd_lines(NOTE).items():  # the note's own Quarto files
        cited_qmd.setdefault(key, []).extend(locs)
    for r in records:
        if r["id"] in cited_qmd:
            r["qmd_lines"] = cited_qmd[r["id"]]
        if r["id"] in cited:
            r["tex_line"] = ", ".join(str(i) for i in cited[r["id"]])
        elif r["id"] in cited_qmd:
            r["tex_line"] = "qmd"  # cited in a vignette only; locations in ledger.json
        else:
            r["tex_line"] = f"{r['tex_line']} (not cited)" if r["tex_line"] else "not cited"

    meta = {
        "date": dt.date.today().isoformat(),
        "sympy": sympy.__version__,
        "numpy": numpy.__version__,
        "python": platform.python_version(),
    }
    (HERE / "ledger.json").write_text(
        json.dumps({"meta": meta, "records": records}, indent=2) + "\n", encoding="utf-8"
    )
    write_tex(records, NOTE / "verification-ledger.tex", meta)

    width = max(len(r["id"]) for r in records)
    for r in records:
        qmd = f"  qmd {', '.join(r['qmd_lines'])}" if r.get("qmd_lines") else ""
        print(f"{r['id']:<{width}}  {r['status']:<4}  {r['kind']:<17}  line {r['tex_line']}{qmd}")
    n_fail = sum(r["status"] == "fail" for r in records)
    print(f"{len(records)} checks, {n_fail} fail")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())

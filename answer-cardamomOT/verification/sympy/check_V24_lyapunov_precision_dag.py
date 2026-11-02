# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V24: what the diagonal of a precision matrix is (note, section on the
diagonal and self-regulation).

- V24a (exact): if A is symmetric (and stable) and Q = q I, then
  Sigma = -(q/2) A^{-1} solves the Lyapunov equation A Sigma + Sigma A' + Q = 0,
  so Omega = -(2/q) A and omega_gg = 2 a_g / q with A_gg = -a_g. G = 3,
  symbolic symmetric A.
- V24b (exact): linear SEM x = B x + e, e ~ N(0, D), B strictly lower
  triangular (G = 3, parents first): Omega = (I - B)' D^{-1} (I - B), obtained
  by inverting Cov x = (I - B)^{-1} D (I - B)^{-T}, and
  omega_gg = 1/d_g + sum over children h of b_hg^2 / d_h.
- V24c (exact): omega_gg = 1 / {sigma_gg (1 - R^2_g)}, with R^2_g the squared
  multiple correlation of gene g on the others,
  R^2_g = s_g,-g' S_-g,-g^{-1} s_g,-g / sigma_gg; G = 3 with a rational
  positive-definite covariance and a symbolic one-parameter family.
"""

from __future__ import annotations

import json

import sympy as sp


def _rec(id_, claim, kind, status, detail, tex_line, detail_full=None):
    return {
        "id": id_,
        "claim": claim,
        "kind": kind,
        "status": status,
        "detail": detail,
        "detail_full": detail if detail_full is None else detail_full,
        "tex_line": tex_line,
    }


def run() -> list[dict]:
    out = []

    # ---- V24a ----------------------------------------------------------
    a1, a2, a3, b12, b13, b23, q = sp.symbols("a1 a2 a3 b12 b13 b23 q", real=True)
    A = sp.Matrix([[-a1, b12, b13], [b12, -a2, b23], [b13, b23, -a3]])
    Sig = -(q / 2) * A.inv()
    lyap = (A * Sig + Sig * A.T + q * sp.eye(3)).applyfunc(sp.simplify)
    Om = Sig.inv().applyfunc(sp.simplify)
    ok_a = lyap == sp.zeros(3, 3) and (Om + (2 / q) * A).applyfunc(sp.simplify) == sp.zeros(3, 3)
    out.append(
        _rec(
            "V24a",
            r"Symmetric $\bm A$ and $\bm Q=q\mathbf I$: $\bm\Sigma=-\tfrac q2\bm A^{-1}$ solves the Lyapunov equation, "
            r"so $\bm\Omega=-\tfrac2q\bm A$ and $\omega_{gg}=2a_g/q$",
            "exact",
            "pass" if ok_a else "fail",
            r"$G=3$, symbolic symmetric $\bm A$; residual of the Lyapunov equation and $\bm\Omega+\tfrac2q\bm A$ "
            r"simplify to zero",
            "",
        )
    )

    # ---- V24b ----------------------------------------------------------
    b21, b31, b32 = sp.symbols("b21 b31 b32", real=True)
    d1, d2, d3 = sp.symbols("d1 d2 d3", positive=True)
    B = sp.Matrix([[0, 0, 0], [b21, 0, 0], [b31, b32, 0]])
    D = sp.diag(d1, d2, d3)
    I3 = sp.eye(3)
    Cov = (I3 - B).inv() * D * (I3 - B).inv().T
    Om_b = Cov.inv().applyfunc(sp.simplify)
    target = ((I3 - B).T * D.inv() * (I3 - B)).applyfunc(sp.simplify)
    diag_formula = [1 / d1 + b21**2 / d2 + b31**2 / d3, 1 / d2 + b32**2 / d3, 1 / d3]
    ok_b = (Om_b - target).applyfunc(sp.simplify) == sp.zeros(3, 3) and all(
        sp.simplify(Om_b[g, g] - diag_formula[g]) == 0 for g in range(3)
    )
    out.append(
        _rec(
            "V24b",
            r"Linear SEM on a directed acyclic graph, parents first: $\bm\Omega=(\mathbf I-\bm B)^\top\bm D^{-1}"
            r"(\mathbf I-\bm B)$ and $\omega_{gg}=1/d_g+\sum_{g\to h}b_{hg}^2/d_h$",
            "exact",
            "pass" if ok_b else "fail",
            r"$G=3$, symbolic strictly lower triangular $\bm B$ and diagonal $\bm D$; inverse of "
            r"$(\mathbf I-\bm B)^{-1}\bm D(\mathbf I-\bm B)^{-\top}$",
            "",
        )
    )

    # ---- V24c ----------------------------------------------------------
    t = sp.Symbol("t", real=True)
    Sg = sp.Matrix([[4, 1, t], [1, 3, sp.Rational(1, 2)], [t, sp.Rational(1, 2), 2]])
    Om_c = Sg.inv()
    ok_c = True
    for g in range(3):
        rest = [h for h in range(3) if h != g]
        s_gr = Sg.extract([g], rest)
        S_rr = Sg.extract(rest, rest)
        R2 = (s_gr * S_rr.inv() * s_gr.T)[0, 0] / Sg[g, g]
        if sp.simplify(Om_c[g, g] - 1 / (Sg[g, g] * (1 - R2))) != 0:
            ok_c = False
    out.append(
        _rec(
            "V24c",
            r"$\omega_{gg}=1/\operatorname{Var}(x_g\mid\bm x_{-g})=1/\{\sigma_{gg}(1-R^2_g)\}>0$",
            "exact",
            "pass" if ok_c else "fail",
            r"$G=3$, covariance with one symbolic entry $t$; Schur complement against the multiple correlation, "
            r"every gene",
            "",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

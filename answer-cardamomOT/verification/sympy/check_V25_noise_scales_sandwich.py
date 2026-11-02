# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V25: intrinsic noise, bulk scales and naive intervals
(summary_side_comparison.qmd).

- V25a (exact): linear-noise approximation of one gene. Birth-death
  (production k, degradation gamma): drift a = gamma, noise intensity
  q = k + gamma x* at the fixed point x* = k/gamma; the one-gene Lyapunov
  equation -2 a sigma + q = 0 gives sigma = x*, Fano factor 1. Bursty
  production (bursts at rate k, geometric size with mean b, so E B = b and
  E B^2 = b + 2 b^2): x* = k b / gamma, q = k (b + 2 b^2) + gamma x*, and
  sigma = x* (1 + b), Fano factor 1 + b.
- V25b (exact): with c_i = a_i n_i and phi_i = a_i^2 n_i, a_i = phi_i / c_i
  and n_i = c_i^2 / phi_i: depth and cell number are separately identified
  from the mean and variance scales.
- V25c (exact): one proportion, two genes with the same signature entry m
  and within-type correlation rho, working independence W = I: the true
  (sandwich) variance (M'M)^-1 M'VM (M'M)^-1 = (1 + rho)/(2 m^2) equals the
  GLS variance (M'V^-1 M)^-1, whereas the naive (M'M)^-1 = 1/(2 m^2): same
  point-estimate efficiency, naive interval too narrow by sqrt(1 + rho).
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


def zero(e) -> bool:
    return sp.simplify(e) == 0


def run() -> list[dict]:
    out = []

    # ---- V25a ----------------------------------------------------------
    k, g, b, s = sp.symbols("k gamma b sigma", positive=True)
    xs = k / g
    q = k + g * xs
    sol_bd = sp.solve(sp.Eq(-2 * g * s + q, 0), s)[0]
    ok_bd = zero(sol_bd - xs)
    xs_b = k * b / g
    EB2 = b + 2 * b**2  # geometric on {0,1,...} with mean b
    # check the geometric second moment from its pmf
    t = sp.Symbol("t", positive=True)
    pgf = 1 / (1 + b - b * t)  # E t^B
    EB = sp.diff(pgf, t).subs(t, 1)
    EBB = (sp.diff(pgf, t, 2) + sp.diff(pgf, t)).subs(t, 1)
    ok_geo = zero(EB - b) and zero(EBB - EB2)
    q_b = k * EB2 + g * xs_b
    sol_b = sp.solve(sp.Eq(-2 * g * s + q_b, 0), s)[0]
    ok_burst = zero(sol_b - xs_b * (1 + b))
    ok_a = ok_bd and ok_geo and ok_burst
    out.append(
        _rec(
            "V25a",
            r"Linear-noise approximation of one gene: Fano factor 1 for birth--death, $1+b$ for geometric bursts "
            r"of mean size $b$",
            "exact",
            "pass" if ok_a else "fail",
            r"Fixed point, noise intensity $q$ and the one-gene Lyapunov equation $-2a\sigma+q=0$; geometric moments "
            r"from the generating function",
            "",
            detail_full=f"birth-death {ok_bd}; geometric moments {ok_geo}; bursts {ok_burst}",
        )
    )

    # ---- V25b ----------------------------------------------------------
    a, n = sp.symbols("a n", positive=True)
    c, phi = a * n, a**2 * n
    ok_b = zero(phi / c - a) and zero(c**2 / phi - n)
    out.append(
        _rec(
            "V25b",
            r"With $c_i=a_in_i$ and $\phi_i=a_i^2n_i$: $a_i=\phi_i/c_i$ and $n_i=c_i^2/\phi_i$",
            "exact",
            "pass" if ok_b else "fail",
            r"Direct substitution",
            "",
        )
    )

    # ---- V25c ----------------------------------------------------------
    m, rho = sp.symbols("m rho", positive=True)
    M = sp.Matrix([m, m])
    V = sp.Matrix([[1, rho], [rho, 1]])
    naive = (M.T * M).inv()[0, 0]
    sandwich = ((M.T * M).inv() * M.T * V * M * (M.T * M).inv())[0, 0]
    gls = (M.T * V.inv() * M).inv()[0, 0]
    ok_c = zero(naive - 1 / (2 * m**2)) and zero(sandwich - (1 + rho) / (2 * m**2)) and zero(gls - sandwich)
    out.append(
        _rec(
            "V25c",
            r"Two genes, correlation $\rho$, same signature entry: OLS and GLS share the variance $(1+\rho)/(2m^2)$, "
            r"the naive independence variance is $1/(2m^2)$",
            "exact",
            "pass" if ok_c else "fail",
            r"Sandwich, GLS and naive variances computed symbolically",
            "",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

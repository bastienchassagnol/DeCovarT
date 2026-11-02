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
- V25d (exact): two genes, symmetric A and Q = q I: closed-form Sigma and
  correlation b/a; one-way regulation or Q = diag(q1, q2) break Omega ~ A.
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

    # ---- V25d ----------------------------------------------------------
    # Two genes, symmetric regulation, isotropic noise: closed form and
    # correlation b/a; then the two counter-examples (one-way regulation,
    # gene-specific noise) where Omega is not proportional to A.
    a_, b_, q_ = sp.symbols("a b q", positive=True)
    s11, s12, s22 = sp.symbols("s11 s12 s22", real=True)
    S = sp.Matrix([[s11, s12], [s12, s22]])

    def lyap_solve(A2, Q2):
        eqs = list(A2 * S + S * A2.T + Q2)
        sol = sp.solve(eqs, [s11, s12, s22], dict=True)[0]
        return S.subs(sol).applyfunc(sp.simplify)

    A_sym = sp.Matrix([[-a_, b_], [b_, -a_]])
    Sig2 = lyap_solve(A_sym, q_ * sp.eye(2))
    target2 = q_ / (2 * (a_**2 - b_**2)) * sp.Matrix([[a_, b_], [b_, a_]])
    corr = sp.simplify(Sig2[0, 1] / sp.sqrt(Sig2[0, 0] * Sig2[1, 1]))
    # Sig2[0,0] > 0 requires a > b; check corr on that branch
    ok_d1 = (Sig2 - target2).applyfunc(sp.simplify) == sp.zeros(2, 2) and zero(
        corr.subs(a_, a_ + b_) - b_ / (a_ + b_)
    )
    # one-way regulation 1 -> 2, isotropic noise: Omega symmetric, A not
    A_one = sp.Matrix([[-a_, 0], [b_, -a_]])
    Om_one = lyap_solve(A_one, q_ * sp.eye(2)).inv().applyfunc(sp.simplify)
    ok_d2 = not zero(Om_one[0, 1] - (-2 / q_) * A_one[0, 1]) or not zero(Om_one[1, 0] - (-2 / q_) * A_one[1, 0])
    # symmetric A, gene-specific noise diag(q1, q2): Omega != c A for any c
    q1, q2, cc = sp.symbols("q1 q2 c", positive=True)
    Om_an = lyap_solve(A_sym, sp.diag(q1, q2)).inv().applyfunc(sp.simplify)
    sol_c = sp.solve(list(Om_an - cc * A_sym), cc, dict=True)
    ok_d3 = len(sol_c) == 0 and zero(
        (Om_an.subs(q2, q1) + (2 / q1) * A_sym).applyfunc(sp.simplify).norm()
    )
    ok_d = ok_d1 and ok_d2 and ok_d3
    out.append(
        _rec(
            "V25d",
            r"Two genes, symmetric $\bm A$ (decay $a$, mutual effect $b$), $\bm Q=q\mathbf I$: "
            r"$\bm\Sigma=\tfrac{q}{2(a^2-b^2)}\bigl(\begin{smallmatrix}a&b\\b&a\end{smallmatrix}\bigr)$, "
            r"correlation $b/a$; one-way regulation or gene-specific noise break $\bm\Omega\propto\bm A$",
            "exact",
            "pass" if ok_d else "fail",
            r"Lyapunov equation solved symbolically for $G=2$; proportionality constant solved for and absent "
            r"when $q_1\neq q_2$",
            "",
            detail_full=f"closed form and corr {ok_d1}; one-way {ok_d2}; anisotropic {ok_d3}",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

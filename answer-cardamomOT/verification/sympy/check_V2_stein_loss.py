# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V2: the expected Gaussian negative log-score (Stein's loss up to constants)
L(Sigma) = 1/2 {log det Sigma + tr(Sigma^-1 Sigma_true)} is minimised at
Sigma = Sigma_true, whatever the law of the cells (finite second moments,
known mean).

Derivation route: the expected score is computed from the log-density of
N(0, Sigma) under a non-Gaussian law (exact expectation over a finite
support); stationarity and the Hessian are derived by differentiation of
L in the free entries of Sigma; a numeric check at random SPD matrices uses
mpmath (shipped with SymPy) and a fixed seed.
"""

from __future__ import annotations

import json
import random

import mpmath as mp
import sympy as sp


def _rec(id_, claim, kind, status, detail, tex_line, detail_full=None):
    """detail: short text for the ledger table, with $...$ LaTeX (sympy.latex);
    detail_full: the full engine output, kept in ledger.json."""
    return {
        "id": id_,
        "claim": claim,
        "kind": kind,
        "status": status,
        "detail": detail,
        "detail_full": detail if detail_full is None else detail_full,
        "tex_line": tex_line,
    }


def _loss(S: sp.Matrix, T: sp.Matrix):
    return sp.Rational(1, 2) * (sp.log(S.det()) + (S.inv() * T).trace())


def run() -> list[dict]:
    out: list[dict] = []

    # ------------------------------------------------ scalar case
    s, t, x = sp.symbols("s t x", positive=True)
    f = _loss(sp.Matrix([[s]]), sp.Matrix([[t]]))
    crit = sp.solve(sp.diff(f, s), s)
    f2 = sp.simplify(sp.diff(f, s, 2).subs(s, t))
    # global: f(s) - f(t) = (x - 1 - log x)/2 with x = t/s, convex in x, zero at 1
    h = sp.simplify((f - f.subs(s, t)).subs(s, t / x))
    h_target = (x - 1 - sp.log(x)) / 2
    glob = (
        sp.simplify(sp.expand_log(h - h_target, force=True)) == 0
        and sp.diff(h_target, x).subs(x, 1) == 0
        and sp.simplify(sp.diff(h_target, x, 2) - 1 / (2 * x**2)) == 0
    )
    ok = crit == [t] and sp.simplify(f2 - 1 / (2 * t**2)) == 0 and glob
    out.append(
        _rec(
            "V2a",
            r"Scalar case: $L(s)=\tfrac12(\log s+t/s)$ has its unique minimum at $s=t$",
            "exact",
            "pass" if ok else "fail",
            r"$L'(s)=0$ only at $s=" + ", ".join(sp.latex(c) for c in crit) + r"$; $L''(t)=" + sp.latex(f2) + r"$ $>0$; $L(s)-L(t)=\tfrac12(x-1-\log x)\ge0$ with $x=t/s$",
            453,
            detail_full=(
                f"L'(s) = 0 only at s = {crit}; L''(t) = {sp.sstr(f2)} > 0; "
                "L(s) - L(t) = (x - 1 - log x)/2 with x = t/s, convex with minimum 0 at x = 1, "
                "so the minimum is global"
            ),
        )
    )

    # ------------------------------------------------ 2 x 2 case
    s11, s22 = sp.symbols("s11 s22", positive=True)
    s12 = sp.Symbol("s12", real=True)
    l11, l22 = sp.symbols("l11 l22", positive=True)
    l21 = sp.Symbol("l21", real=True)
    Lt = sp.Matrix([[l11, 0], [l21, l22]])
    T = Lt * Lt.T  # Sigma_true, SPD by construction
    S = sp.Matrix([[s11, s12], [s12, s22]])
    Lf = _loss(S, T)
    vars_ = [s11, s12, s22]
    at_true = {s11: T[0, 0], s12: T[0, 1], s22: T[1, 1]}
    grad = [sp.simplify(sp.diff(Lf, v).subs(at_true)) for v in vars_]
    H = sp.hessian(Lf, vars_).subs(at_true).applyfunc(sp.simplify)
    minors = [sp.factor(H[:k, :k].det()) for k in (1, 2, 3)]
    # positivity from SymPy's assumption system (l11, l22 > 0, l21 real)
    pos = [m.is_positive is True for m in minors]
    ok = all(g == 0 for g in grad) and all(pos)
    out.append(
        _rec(
            "V2b",
            r"$2\times2$: the gradient of $L$ vanishes at $\bm\Sigma=\bm\Sigma_{\mathrm{true}}$ and the Hessian is positive definite there",
            "exact",
            "pass" if ok else "fail",
            r"gradient zero at $\bm\Sigma_{\mathrm{true}}=\bm L\bm L^{\mathsf T}$ (symbolic Cholesky factor); leading principal minors of the Hessian " + ", ".join("$" + sp.latex(m) + "$" for m in minors) + ", all positive",
            453,
            detail_full=(
                "Sigma_true = L L^T (symbolic Cholesky factor); gradient in (s11, s12, s22) at "
                f"Sigma_true = {grad}; leading principal minors of the Hessian: "
                f"{[sp.sstr(m) for m in minors]}, all > 0"
            ),
        )
    )

    # ------------------------------------------------ numeric, random SPD
    mp.mp.dps = 40
    rng = random.Random(20261030)
    d, n_true, n_pert = 4, 20, 50
    worst = mp.inf
    for _ in range(n_true):
        A = mp.matrix([[rng.gauss(0, 1) for _ in range(d)] for _ in range(d)])
        Tn = A * A.T + mp.eye(d) * mp.mpf("0.1")

        def Ln(Sm):
            return (mp.log(mp.det(Sm)) + sum((mp.inverse(Sm) * Tn)[i, i] for i in range(d))) / 2

        L0 = Ln(Tn)
        for _ in range(n_pert):
            Bm = mp.matrix([[rng.gauss(0, 1) for _ in range(d)] for _ in range(d)])
            Sm = Bm * Bm.T + mp.eye(d) * mp.mpf("0.05")
            worst = min(worst, Ln(Sm) - L0)
    ok = worst > 0
    out.append(
        _rec(
            "V2c",
            r"Numeric: $L(\bm\Sigma)>L(\bm\Sigma_{\mathrm{true}})$ at random SPD $\bm\Sigma$",
            "numeric",
            "pass" if ok else "fail",
            f"mpmath at 40 digits, seed 20261030, $d={d}$, {n_true} true matrices with {n_pert} random competitors each; smallest gap {mp.nstr(worst, 6)}",
            453,
            detail_full=(
                f"mpmath, 40 digits, seed 20261030, d = {d}, {n_true} true matrices x {n_pert} "
                f"random SPD competitors; smallest L(Sigma) - L(Sigma_true) = "
                f"{mp.nstr(worst, 6)} (> 0)"
            ),
        )
    )

    # ------------------------------------------------ "whatever the law of the cells"
    # Non-Gaussian law with finite support, mean zero: expected negative
    # Gaussian log-score E[-log N(x; 0, Sigma)] computed exactly.
    pts = [sp.Matrix(v) for v in ([3, 1], [-1, 2], [-2, -2], [0, -1])]
    w = [sp.Rational(1, 6), sp.Rational(1, 3), sp.Rational(1, 4), sp.Rational(1, 4)]
    mean = sum((wi * pi for wi, pi in zip(w, pts)), sp.zeros(2, 1))
    pts = [pi - mean for pi in pts]  # centre exactly
    T_disc = sum((wi * pi * pi.T for wi, pi in zip(w, pts)), sp.zeros(2, 2))
    Sinv = S.inv()
    exp_score = sum(
        wi * (2 * sp.log(2 * sp.pi) + sp.log(S.det()) + (pi.T * Sinv * pi)[0, 0]) / 2
        for wi, pi in zip(w, pts)
    )
    target = sp.log(2 * sp.pi) + _loss(S, T_disc)
    same = sp.simplify(exp_score - target) == 0
    grad_d = [sp.simplify(sp.diff(exp_score, v).subs({s11: T_disc[0, 0], s12: T_disc[0, 1], s22: T_disc[1, 1]})) for v in vars_]
    ok = same and all(g == 0 for g in grad_d)
    out.append(
        _rec(
            "V2d",
            r"The expected Gaussian log-score depends on the cell law only through its covariance, so the minimiser is $\bm\Sigma_{\mathrm{true}}$ for any law with finite second moments",
            "exact",
            "pass" if ok else "fail",
            r"centred 4-point non-Gaussian law: expected score $=\log 2\pi+L(\bm\Sigma)$; gradient zero at its covariance $" + sp.latex(T_disc) + r"$; the mean is taken as known",
            453,
            detail_full=(
                "4-point non-Gaussian law (rational weights, centred): expected score equals "
                "log(2 pi) + L(Sigma) with Sigma_true = its covariance "
                f"{sp.sstr(T_disc)}; gradient at that covariance = {grad_d}. "
                "Conditional on a known (or separately estimated) mean"
            ),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

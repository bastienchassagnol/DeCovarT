# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V8: Hellinger and 2-Wasserstein distances between two Gaussian type laws
with a common covariance, against the Mahalanobis distance (gene-selection
section, after the Fisher-information paragraph).

The note says: the Hellinger distance between the two type laws is monotone in
the Mahalanobis distance Delta' Sigma^-1 Delta; the 2-Wasserstein distance is
not, since for equal covariances it reduces to ||Delta||.

Derivation routes:
- V8a: the Bhattacharyya coefficient BC = int sqrt(f1 f2) is derived by
  completing the square, an exact polynomial identity in x for a general
  symmetric 2 x 2 Sigma, which leaves exp(-M/8) times a normal density that
  integrates to 1; one dimension is cross-checked with sympy.integrate.
  H^2 = 1 - BC, and dH^2/dM > 0.
- V8b: for any coupling (X, Y) with X ~ N(mu1, Sigma), Y ~ N(mu2, Sigma),
  E||X - Y||^2 = ||Delta||^2 + tr Cov(X - Y) (exact identity with a symbolic
  cross-covariance C); the trace is >= 0 and is 0 for the translation
  coupling Y = X - Delta, so W2 = ||Delta||, free of Sigma.
- V8c: exact counter-example with one Sigma and two mean differences, in
  which W2 and the Mahalanobis distance (hence Hellinger) rank the pairs in
  opposite orders.
"""

from __future__ import annotations

import json

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


def run() -> list[dict]:
    out: list[dict] = []

    # ---------------- V8a: Hellinger, general 2 x 2 Sigma ----------------
    xv = sp.Matrix(sp.symbols("x1 x2", real=True))
    m1 = sp.Matrix(sp.symbols("m1_1 m1_2", real=True))
    m2 = sp.Matrix(sp.symbols("m2_1 m2_2", real=True))
    s11, s22 = sp.symbols("s11 s22", positive=True)
    s12 = sp.Symbol("s12", real=True)
    S = sp.Matrix([[s11, s12], [s12, s22]])
    P = S.inv()
    q = lambda v: (v.T * P * v)[0]  # noqa: E731
    Dl = m1 - m2
    M = q(Dl)  # squared Mahalanobis distance
    mbar = (m1 + m2) / 2
    # log sqrt(f1 f2) = -log(2 pi) - log det S / 2 - {q(x - m1) + q(x - m2)} / 4
    # claim: {q(x - m1) + q(x - m2)} / 4 = q(x - mbar) / 2 + M / 8
    ident = sp.simplify(sp.expand((q(xv - m1) + q(xv - m2)) / 4 - q(xv - mbar) / 2 - M / 8)) == 0
    # so sqrt(f1 f2) = exp(-M/8) * N(x; mbar, S), and BC = exp(-M/8)
    # 1-D cross-check by direct integration
    x = sp.Symbol("x", real=True)
    a, b = sp.symbols("a b", real=True)
    s = sp.Symbol("s", positive=True)
    f = lambda m: sp.exp(-((x - m) ** 2) / (2 * s)) / sp.sqrt(2 * sp.pi * s)  # noqa: E731
    integrand = sp.exp(-((x - a) ** 2 + (x - b) ** 2) / (4 * s)) / sp.sqrt(2 * sp.pi * s)
    assert sp.simplify(sp.sqrt(f(a) * f(b)) - integrand) == 0
    bc1 = sp.simplify(sp.integrate(integrand, (x, -sp.oo, sp.oo), conds="none"))
    ok_1d = sp.simplify(bc1 - sp.exp(-((a - b) ** 2) / (8 * s))) == 0
    Ms = sp.Symbol("M", nonnegative=True)
    H2 = 1 - sp.exp(-Ms / 8)
    dH2 = sp.diff(H2, Ms)
    mono = sp.simplify(dH2 - sp.exp(-Ms / 8) / 8) == 0  # > 0 for every M
    out.append(
        _rec(
            "V8a",
            r"Common $\bm\Sigma$: $H^2=1-\exp(-\tfrac18\bm\Delta^{\mathsf T}\bm\Sigma^{-1}\bm\Delta)$, "
            r"increasing in the Mahalanobis distance",
            "exact",
            "pass" if ident and ok_1d and mono else "fail",
            r"completing the square: $\sqrt{f_1f_2}=e^{-M/8}\,\mathcal N(\bm x;\bar{\bm\mu},\bm\Sigma)$ "
            r"(polynomial identity, general $2\times2$ $\bm\Sigma$); one-dimensional integral $="
            + sp.latex(sp.exp(sp.factor(sp.log(bc1).expand(force=True))))
            + r"$; $\mathrm dH^2/\mathrm dM="
            + sp.latex(dH2)
            + r">0$",
            155,
            detail_full=(
                f"square-completion identity: {ident}; 1-D Bhattacharyya integral = {sp.sstr(bc1)} "
                f"(matches exp(-(a-b)^2/(8 s)): {ok_1d}); dH2/dM = {sp.sstr(dH2)}"
            ),
        )
    )

    # ---------------- V8b: W2 = ||Delta|| for a common Sigma -----------------
    # coupling with cross-covariance C = Cov(X, Y)
    C = sp.Matrix(2, 2, sp.symbols("c11 c12 c21 c22", real=True))
    EXX = S + m1 * m1.T
    EYY = S + m2 * m2.T
    EXY = C + m1 * m2.T
    cost = (EXX - EXY - EXY.T + EYY).trace()  # E||X - Y||^2 = tr E(X-Y)(X-Y)'
    cov_diff = 2 * S - C - C.T  # Cov(X - Y)
    ident_w = sp.expand(cost - (Dl.T * Dl)[0] - cov_diff.trace()) == 0
    # translation coupling Y = X - Delta: C = Cov(X, X) = S, Cov(X - Y) = 0
    trans = sp.simplify(cost.subs({C[i, j]: S[i, j] for i in range(2) for j in range(2)}) - (Dl.T * Dl)[0]) == 0
    out.append(
        _rec(
            "V8b",
            r"Common $\bm\Sigma$: $W_2=\lVert\bm\Delta\rVert$, which does not depend on $\bm\Sigma$",
            "exact",
            "pass" if ident_w and trans else "fail",
            r"for any coupling, $\mathbb E\lVert\bm X-\bm Y\rVert^2=\lVert\bm\Delta\rVert^2+"
            r"\operatorname{tr}\operatorname{Cov}(\bm X-\bm Y)\ge\lVert\bm\Delta\rVert^2$ (exact identity, "
            r"symbolic cross-covariance); the translation $\bm Y=\bm X-\bm\Delta$ attains it",
            156,
            detail_full=(
                f"E||X-Y||^2 = {sp.sstr(sp.expand(cost))}; identity with ||Delta||^2 + tr Cov(X-Y): "
                f"{ident_w}; translation coupling gives ||Delta||^2: {trans}"
            ),
        )
    )

    # ---------------- V8c: rankings disagree ---------------------------------
    S0 = sp.diag(1, 100)
    pairs = {"A": sp.Matrix([1, 0]), "B": sp.Matrix([0, 3])}
    w2 = {k: sp.sqrt((v.T * v)[0]) for k, v in pairs.items()}
    mah = {k: (v.T * S0.inv() * v)[0] for k, v in pairs.items()}
    hel = {k: 1 - sp.exp(-mah[k] / 8) for k in pairs}
    opposite = (w2["A"] < w2["B"]) and (mah["A"] > mah["B"]) and bool(sp.N(hel["A"] - hel["B"], 30) > 0)
    out.append(
        _rec(
            "V8c",
            r"Counter-example: $W_2$ and the Mahalanobis distance can rank two pairs of type laws in "
            r"opposite orders",
            "exact",
            "pass" if opposite else "fail",
            r"$\bm\Sigma=" + sp.latex(S0) + r"$; $\bm\Delta_A=(1,0)$: $W_2="
            + sp.latex(w2["A"])
            + r"$, $M="
            + sp.latex(mah["A"])
            + r"$; $\bm\Delta_B=(0,3)$: $W_2="
            + sp.latex(w2["B"])
            + r"$, $M="
            + sp.latex(mah["B"])
            + r"$. $W_2$ prefers $B$; $M$ and $H$ prefer $A$",
            156,
            detail_full=(
                f"W2: {w2}; Mahalanobis^2: {mah}; H^2: "
                f"{ {k: sp.N(v, 8) for k, v in hel.items()} }"
            ),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

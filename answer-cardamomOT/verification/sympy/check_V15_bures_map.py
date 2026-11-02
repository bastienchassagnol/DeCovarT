# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V15: the Gaussian transport map of vignette 05-1 (eq-bures-map,
sec-bisque-ot), promoted from agent C's scratch check.

Claims: A = S1^{-1/2} (S1^{1/2} S2 S1^{1/2})^{1/2} S1^{-1/2} is symmetric
positive definite and pushes N(m1, S1) onto N(m2, S2) (A S1 A = S2); the
expected squared displacement of T(x) = m2 + A(x - m1) is
||m1 - m2||^2 + tr{S1 + S2 - 2 (S1^{1/2} S2 S1^{1/2})^{1/2}}; for diagonal
covariances A is diagonal with entries s2g/s1g (Bisque's gene-wise map) and,
in one dimension, the cost is (m1 - m2)^2 + (s1 - s2)^2.
Optimality of a symmetric positive definite linear map (Brenier's theorem)
is quoted, not checked.

Derivation route:
- V15a: exact 2 x 2 instance with rational S1 = diag(1, 4) and
  S2 = [[2, 1], [1, 2]]; the square root of a 2 x 2 SPD matrix M is
  (M + sqrt(det M) I) / sqrt(tr M + 2 sqrt(det M)), verified by squaring.
  The displacement cost is computed as E||(A - I)(x - m1)||^2 + ||m2 - m1||^2
  = tr((A - I) S1 (A - I)) + ||m2 - m1||^2 and compared with the formula.
- V15b: G = 4, random SPD pairs at 40-digit precision (mpmath), seed fixed:
  pushforward, symmetry, positivity, cost identity.
- V15c: diagonal and one-dimensional cases, symbolic.
"""

from __future__ import annotations

import json
import random

import sympy as sp
import mpmath as mp  # ships with SymPy


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


def _sqrt2(M: sp.Matrix) -> sp.Matrix:
    d = sp.sqrt(M.det())
    return (M + d * sp.eye(2)) / sp.sqrt(M.trace() + 2 * d)


def _sci(x) -> str:
    mant, ex = f"{float(x):.1e}".split("e")
    return f"${mant}\\times10^{{{int(ex)}}}$"


def run() -> list[dict]:
    out: list[dict] = []
    # ---------------- V15a: exact 2 x 2 ----------------------------------
    S1 = sp.diag(1, 4)
    S2 = sp.Matrix([[2, 1], [1, 2]])
    r1 = sp.diag(1, 2)
    ok_r1 = r1 * r1 == S1
    M = r1 * S2 * r1
    rM = _sqrt2(M)
    ok_rM = all(sp.simplify(e) == 0 for e in rM * rM - M)
    A = (r1.inv() * rM * r1.inv()).applyfunc(sp.simplify)
    ok_push = all(sp.simplify(e) == 0 for e in A * S1 * A - S2)
    ok_sym = all(sp.simplify(e) == 0 for e in A - A.T)
    ok_pd = sp.simplify(A[0, 0]) > 0 and sp.simplify(A.det()) > 0
    m1v = sp.Matrix([1, 0])
    m2v = sp.Matrix([0, 3])
    I2 = sp.eye(2)
    cost_map = sp.simplify(((A - I2) * S1 * (A - I2)).trace() + (m2v - m1v).dot(m2v - m1v))
    cost_formula = sp.simplify((m1v - m2v).dot(m1v - m2v) + (S1 + S2 - 2 * rM).trace())
    ok_cost = sp.simplify(cost_map - cost_formula) == 0
    ok_a = ok_r1 and ok_rM and ok_push and ok_sym and ok_pd and ok_cost
    out.append(
        _rec(
            "V15a",
            r"eq-bures-map: $\bm A$ symmetric positive definite, $\bm A\bm\Sigma_1\bm A=\bm\Sigma_2$, and the "
            r"displacement cost equals $\lVert\bm m_1-\bm m_2\rVert^2+\operatorname{tr}\{\bm\Sigma_1+\bm\Sigma_2"
            r"-2(\bm\Sigma_1^{1/2}\bm\Sigma_2\bm\Sigma_1^{1/2})^{1/2}\}$",
            "exact",
            "pass" if ok_a else "fail",
            r"$\bm\Sigma_1=\operatorname{diag}(1,4)$, $\bm\Sigma_2=\begin{pmatrix}2&1\\1&2\end{pmatrix}$; closed-form "
            r"$2\times2$ square root checked by squaring; both costs equal exactly ($\approx" + f"{float(cost_formula):.4f}" + r"$)",
            "",
            detail_full=(
                f"sqrt(S1): {ok_r1}; sqrt(M): {ok_rM}; A = {sp.sstr(A)}; A S1 A = S2: {ok_push}; symmetric: {ok_sym}; "
                f"PD: {ok_pd}; cost by the map = {cost_map}; formula = {cost_formula}; equal: {ok_cost}"
            ),
        )
    )

    # ---------------- V15b: G = 4, 40 digits -------------------------------
    mp.mp.dps = 40
    rnd = random.Random(20261031)
    G = 4
    worst = mp.mpf(0)
    ok_b = True
    for _ in range(20):
        mats = []
        for _k in range(2):
            X = mp.matrix([[mp.mpf(rnd.gauss(0, 1)) for _ in range(G)] for _ in range(G)])
            mats.append(X * X.T + G * mp.eye(G))
        Sa, Sb = mats
        ra = mp.sqrtm(Sa)
        rai = mp.inverse(ra)
        rM_ = mp.sqrtm(ra * Sb * ra)
        A_ = rai * rM_ * rai
        push = mp.mnorm(A_ * Sa * A_ - Sb, 1)
        sym = mp.mnorm(A_ - A_.T, 1)
        eig = min(mp.re(e) for e in mp.eig((A_ + A_.T) / 2)[0])
        Ig = mp.eye(G)
        P = (A_ - Ig) * Sa * (A_ - Ig)
        cost_m = sum(P[i, i] for i in range(G))
        cost_f = sum((Sa + Sb - 2 * rM_)[i, i] for i in range(G))
        gap = abs(cost_m - cost_f)
        worst = max(worst, push, sym, gap)
        ok_b &= eig > 0
    ok_b = ok_b and worst < mp.mpf("1e-30")
    out.append(
        _rec(
            "V15b",
            r"Same identities for $G=4$",
            "numeric",
            "pass" if ok_b else "fail",
            r"mpmath at 40 digits, 20 random SPD pairs (seed 20261031); largest residual "
            + _sci(worst) + r", tolerance $10^{-30}$; all eigenvalues of $\bm A$ positive",
            "",
            detail_full=f"max residual (pushforward, symmetry, cost) = {mp.nstr(worst, 5)}",
        )
    )

    # ---------------- V15c: diagonal and 1-D cases ------------------------
    s1, s2 = sp.symbols("sigma1 sigma2", positive=True)
    m1, m2 = sp.symbols("m1 m2", real=True)
    Ad = sp.sqrt(1 / s1**2) * sp.sqrt(s1**2 * s2**2) * sp.sqrt(1 / s1**2)
    ok_diag = sp.simplify(Ad - s2 / s1) == 0
    w2 = (m1 - m2) ** 2 + s1**2 + s2**2 - 2 * sp.sqrt(s1**2 * s2**2)
    ok_1d = sp.simplify(w2 - ((m1 - m2) ** 2 + (s1 - s2) ** 2)) == 0
    ok_c = ok_diag and ok_1d
    out.append(
        _rec(
            "V15c",
            r"Diagonal covariances: $\bm A=\operatorname{diag}(\sigma_{2g}/\sigma_{1g})$, Bisque's gene-wise map; "
            r"one dimension: $W_2^2=(m_1-m_2)^2+(\sigma_1-\sigma_2)^2$",
            "exact",
            "pass" if ok_c else "fail",
            r"entrywise square roots of diagonal matrices, symbolic",
            "",
            detail_full=f"diagonal entry = s2/s1: {ok_diag}; 1-D cost: {ok_1d}",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

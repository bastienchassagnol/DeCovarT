# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V11: the linear shrinkage estimate (1 - rho) V + rho T is positive definite
for any C >= 2 when the target T is positive definite and rho is in (0, 1]
(identifiability box, part (a)).

Derivation routes:
- V11a (2 x 2, exact): write V = L L' and T = M M' (Cholesky, any positive
  semi-definite V, positive-definite T). det(A + B) = det A + det B +
  tr(adj(A) B) for 2 x 2 matrices, and tr(adj(L L') M M') = ||adj(L) M||_F^2,
  so det{(1 - rho) V + rho T} = (1 - rho)^2 det V + rho^2 det T
  + rho (1 - rho) ||adj(L) M||_F^2 >= rho^2 det T > 0; the trace is
  (1 - rho) tr V + rho tr T > 0. Positive trace and determinant give a
  positive-definite 2 x 2 matrix. The identity is checked as a polynomial.
- V11b (any G, exact): x' V x = ||X_c x||^2 / C >= 0 (symbolic data, G = 3,
  C = 2, so V has rank 1), hence x'{(1 - rho) V + rho T} x
  >= rho lambda_min(T) ||x||^2 (Weyl's inequality in its simplest form).
- V11c (numeric): G = 40, C = 2 (rank-1 V), mpmath at 40 digits, seed
  20261030: the smallest eigenvalue of the shrinkage estimate with the
  Schaefer-Strimmer diagonal target T = diag(V) and with the Ledoit-Wolf
  target (tr V / G) I is >= rho lambda_min(T) > 0 for rho in {0.01, 0.5, 1}.
  The diagonal target is positive definite only if no gene is constant
  within the type, which is the role of the presence gate.
"""

from __future__ import annotations

import json
import random

import mpmath as mp
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
    out: list[dict] = []
    rho = sp.Symbol("rho", positive=True)

    # ---------------- V11a: 2 x 2, exact ----------------
    l1, l2, l3 = sp.symbols("l1 l2 l3", real=True)
    m1, m3 = sp.symbols("m1 m3", positive=True)
    m2 = sp.Symbol("m2", real=True)
    L = sp.Matrix([[l1, 0], [l2, l3]])
    M = sp.Matrix([[m1, 0], [m2, m3]])
    V, T = L * L.T, M * M.T
    S = (1 - rho) * V + rho * T
    adjL = L.adjugate()
    mixed = sum(e**2 for e in (adjL * M))  # ||adj(L) M||_F^2
    ident = sp.expand(S.det() - ((1 - rho) ** 2 * V.det() + rho**2 * T.det() + rho * (1 - rho) * mixed)) == 0
    adj_ok = sp.expand((V.adjugate() * T).trace() - mixed) == 0
    trace = sp.expand(S.trace() - ((1 - rho) * V.trace() + rho * T.trace())) == 0
    out.append(
        _rec(
            "V11a",
            r"$2\times2$: $(1-\rho)\bm V+\rho\bm T\succ0$ for $\bm V\succeq0$, $\bm T\succ0$, $\rho\in(0,1]$",
            "exact",
            "pass" if ident and adj_ok and trace else "fail",
            r"$\bm V=\bm L\bm L^{\mathsf T}$, $\bm T=\bm M\bm M^{\mathsf T}$: $\det=(1-\rho)^2\det\bm V+\rho^2\det\bm T"
            r"+\rho(1-\rho)\lVert\operatorname{adj}(\bm L)\bm M\rVert_F^2\ge\rho^2\det\bm T>0$ "
            r"(polynomial identity); trace $(1-\rho)\operatorname{tr}\bm V+\rho\operatorname{tr}\bm T>0$",
            0,
            detail_full=(
                f"determinant identity: {ident}; tr(adj(V) T) = ||adj(L) M||_F^2: {adj_ok}; trace identity: {trace}; "
                f"||adj(L) M||_F^2 = {sp.sstr(sp.expand(mixed))}"
            ),
        )
    )

    # ---------------- V11b: any G, quadratic form ----------------
    G, C = 3, 2
    X = sp.Matrix(C, G, sp.symbols("x1:7", real=True))
    H = sp.eye(C) - sp.ones(C, C) / C
    Xc = H * X
    Vs = Xc.T * Xc / C
    xv = sp.Matrix(sp.symbols("z1:4", real=True))
    quad = sp.expand((xv.T * Vs * xv)[0] - sum(e**2 for e in (Xc * xv)) / C) == 0
    rank1 = Vs.subs({s: i + 1 for i, s in enumerate(sorted(X.free_symbols, key=str))}).rank() == 1
    out.append(
        _rec(
            "V11b",
            r"Any $G$: $\lambda_{\min}\{(1-\rho)\bm V+\rho\bm T\}\ge\rho\,\lambda_{\min}(\bm T)$",
            "exact",
            "pass" if quad and rank1 else "fail",
            r"$\bm x^{\mathsf T}\bm V\bm x=\lVert\bm X_c\bm x\rVert^2/C\ge0$ (identity with symbolic data, $G=3$, "
            r"$C=2$, $\operatorname{rank}\bm V=1$), so $\bm x^{\mathsf T}\hat{\bm\Sigma}\bm x\ge\rho\,\bm x^{\mathsf T}"
            r"\bm T\bm x\ge\rho\,\lambda_{\min}(\bm T)\lVert\bm x\rVert^2$",
            0,
            detail_full=f"x'Vx = ||X_c x||^2 / C: {quad}; rank V = 1 at an integer point: {rank1}",
        )
    )

    # ---------------- V11c: numeric, G = 40, C = 2 ----------------
    mp.mp.dps = 40
    rng = random.Random(20261030)
    Gn, Cn = 40, 2
    Xn = mp.matrix(Cn, Gn)
    for i in range(Cn):
        for j in range(Gn):
            Xn[i, j] = mp.mpf(rng.gauss(0, 1)) * (1 + j / 10)
    mean = [sum(Xn[i, j] for i in range(Cn)) / Cn for j in range(Gn)]
    Xc_n = mp.matrix(Cn, Gn)
    for i in range(Cn):
        for j in range(Gn):
            Xc_n[i, j] = Xn[i, j] - mean[j]
    Vn = (Xc_n.T * Xc_n) / Cn
    ev_V = mp.eigsy(Vn)[0]
    rank_num = sum(1 for e in ev_V if abs(e) > mp.mpf(10) ** -30)
    targets = {
        "diag": mp.diag([Vn[j, j] for j in range(Gn)]),
        "scaled identity": mp.eye(Gn) * (sum(Vn[j, j] for j in range(Gn)) / Gn),
    }
    worst = None
    ok_num = rank_num == 1
    rows = []
    for name, Tn in targets.items():
        tmin = min(Tn[j, j] for j in range(Gn))
        for r in (mp.mpf("0.01"), mp.mpf("0.5"), mp.mpf(1)):
            Sn = (1 - r) * Vn + r * Tn
            lmin = min(mp.eigsy(Sn)[0])
            gap = lmin - r * tmin
            ok_num &= lmin > 0 and gap > -mp.mpf(10) ** -30
            rows.append((name, float(r), mp.nstr(lmin, 6), mp.nstr(r * tmin, 6)))
            worst = gap if worst is None else min(worst, gap)
    out.append(
        _rec(
            "V11c",
            r"Numeric check of V11b at $G=40$, $C=2$, with the diagonal and the scaled-identity targets",
            "numeric",
            "pass" if ok_num else "fail",
            r"mpmath, 40 digits, seed 20261030; $\operatorname{rank}\bm V=" + str(rank_num)
            + r"$; for $\rho\in\{0.01,0.5,1\}$ and both targets, $\lambda_{\min}(\hat{\bm\Sigma})\ge\rho\,"
            r"\lambda_{\min}(\bm T)>0$ (smallest gap $"
            + (mp.nstr(worst, 3) if abs(worst) > mp.mpf(10) ** -30 else "0")
            + r"$ to 30 digits, attained at $\rho=1$)",
            0,
            detail_full=f"numerical rank of V: {rank_num}; (target, rho, lambda_min, rho*lambda_min(T)): {rows}",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

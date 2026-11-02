# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V21: Unico's moment model and the spectral (von Neumann) entropy
(vignettes 05-1 #sec-unico and 03-1 #sec-fs-vn-entropy).

Unico (Chen et al. 2025, Genome Biol., Eqs 3-6, p. 11) writes, per gene,
X_i = w_i' Z_i + e_i with Var Z_i = Sigma (k x k, across cell types) and
Var e_i = tau^2, so Var X_i = w_i' Sigma w_i + tau^2, with sum_h w_ih = 1.

- V21a (exact, our remark): w'(Sigma + c 1 1')w = w'Sigma w + c whenever
  1'w = 1, so the constant part c 1 1' of Sigma and tau^2 enter the bulk
  variance identically: the moment conditions cannot separate them.
- V21b (numeric): for a k x k correlation matrix R (trace k), the entropy
  H = -sum_h l_h log l_h of the eigenvalues l_h of R/k satisfies
  0 <= H <= log k, with H = log k at R = I and H = 0 at the rank-one
  R = 1 1'. Random correlation matrices (seed 20261101) and the two
  extremes are checked; the general bound is the entropy bound for a
  probability vector, quoted, not proved here. The Renyi-2 lower bound
  H >= log[k / (1 + (k - 1) mean rho^2)] (mean over off-diagonal pairs;
  Shannon >= Renyi-2 and sum l_h^2 = tr(R^2)/k^2) is checked on the same
  matrices.
"""

from __future__ import annotations

import json

import numpy as np
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


def _entropy(R: np.ndarray) -> float:
    lam = np.linalg.eigvalsh(R) / R.shape[0]
    lam = lam[lam > 1e-15]
    return float(-(lam * np.log(lam)).sum())


def run() -> list[dict]:
    out = []

    # ---- V21a ----------------------------------------------------------
    k = 4
    w = sp.Matrix(sp.symbols("w1:5", real=True))
    S = sp.Matrix(k, k, lambda i, j: sp.Symbol(f"s{min(i, j) + 1}{max(i, j) + 1}"))
    c = sp.Symbol("c", real=True)
    one = sp.ones(k, 1)
    lhs = (w.T * (S + c * one * one.T) * w)[0, 0]
    rhs = (w.T * S * w)[0, 0] + c
    on_simplex = {w[k - 1]: 1 - sum(w[:k - 1])}
    ok_a = sp.simplify(sp.expand((lhs - rhs).subs(on_simplex))) == 0
    out.append(
        _rec(
            "V21a",
            r"If $\mathbf{1}^\top\bm w=1$, then $\bm w^\top(\bm\Sigma_g+c\mathbf{1}\mathbf{1}^\top)\bm w"
            r"=\bm w^\top\bm\Sigma_g\bm w+c$: a constant shift of Unico's cross-type covariance "
            r"is confounded with $\tau_g^2$",
            "exact",
            "pass" if ok_a else "fail",
            r"Symbolic $4\times4$ symmetric $\bm\Sigma$, $w_4=1-w_1-w_2-w_3$; difference expands to 0",
            "",
        )
    )

    # ---- V21b ----------------------------------------------------------
    rng = np.random.default_rng(20261101)
    ok_b = True
    worst_lo, worst_hi = np.inf, -np.inf
    for kk in (2, 3, 5, 10, 50):
        for _ in range(200):
            A = rng.standard_normal((kk, kk + rng.integers(0, 5)))
            C = A @ A.T
            d = np.sqrt(np.diag(C))
            R = C / np.outer(d, d)
            H = _entropy(R)
            worst_lo = min(worst_lo, H)
            worst_hi = max(worst_hi, H / np.log(kk))
            if H < -1e-12 or H > np.log(kk) + 1e-12:
                ok_b = False
            rho2 = (R[~np.eye(kk, dtype=bool)] ** 2).mean()
            if H < np.log(kk / (1 + (kk - 1) * rho2)) - 1e-12:
                ok_b = False
        if abs(_entropy(np.eye(kk)) - np.log(kk)) > 1e-12:
            ok_b = False
        if abs(_entropy(np.ones((kk, kk)))) > 1e-10:
            ok_b = False
    out.append(
        _rec(
            "V21b",
            r"Spectral entropy of a $k\times k$ correlation matrix: $0\le H\le\log k$, "
            r"$H=\log k$ at $\bm R=\mathbf I$, $H=0$ at $\bm R=\mathbf{1}\mathbf{1}^\top$, and the Rényi-2 bound "
            r"$H\ge\log[k/\{1+(k-1)\overline{\rho^2}\}]$",
            "numeric",
            "pass" if ok_b else "fail",
            rf"$k\in\{{2,3,5,10,50\}}$, 200 random correlation matrices each, seed 20261101; "
            rf"smallest $H$ {worst_lo:.3f}, largest $H/\log k$ {worst_hi:.3f}; tolerance $10^{{-12}}$",
            "",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

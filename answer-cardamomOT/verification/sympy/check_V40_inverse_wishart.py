# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V40: mean of the inverse empirical covariance (vignette 03-2,
eq-inverse-wishart-bias).

With V_j = X_c' X_c / C_j and C_j V_j ~ W_G(C_j - 1, Sigma), the inverse
Wishart mean E[W^{-1}] = Sigma^{-1} / (n - G - 1), n = C_j - 1, gives
E[V_j^{-1}] = C_j / (C_j - G - 2) Omega.
Checks: (i) G = 1 exactly, from the chi-square density (W = sigma^2 chi^2_n);
(ii) the algebra from the Wishart mean to the stated factor; (iii) a Monte
Carlo check for G = 3, C_j = 15 (fixed seed, numeric, 3 % tolerance);
(iv) the factor at C_j = 2G tends to 2.
"""

from __future__ import annotations

import json

import numpy as np
import sympy as sp


def run():
    n, x = sp.symbols("n x", positive=True)
    chi2 = x ** (n / 2 - 1) * sp.exp(-x / 2) / (2 ** (n / 2) * sp.gamma(n / 2))
    e_inv = sp.simplify(
        sp.integrate(chi2 / x, (x, 0, sp.oo), conds="none").subs(n, 7)
    )
    ok_i = e_inv == sp.Rational(1, 5)  # 1/(n - 2) at n = 7, i.e. 1/(n - G - 1), G = 1

    C, G = sp.symbols("C G", positive=True)
    factor = C * 1 / ((C - 1) - G - 1)
    ok_ii = sp.simplify(factor - C / (C - G - 2)) == 0
    lim = sp.limit((C / (C - G - 2)).subs(C, 2 * G), G, sp.oo)

    rng = np.random.default_rng(20261110)
    g, c, reps = 3, 15, 40000
    L = np.array([[1.0, 0, 0], [0.5, 1.0, 0], [-0.3, 0.2, 1.0]])
    sigma = L @ L.T
    acc = np.zeros((g, g))
    for _ in range(reps):
        X = rng.standard_normal((c, g)) @ L.T
        Xc = X - X.mean(axis=0)
        acc += np.linalg.inv(Xc.T @ Xc / c)
    mc = acc / reps
    target = c / (c - g - 2) * np.linalg.inv(sigma)
    rel = float(np.max(np.abs(mc - target)) / np.max(np.abs(target)))
    ok_iii = rel < 0.03

    ok = ok_i and ok_ii and ok_iii and lim == 2
    return [
        {
            "id": "V40",
            "claim": r"$\mathbb E[\boldsymbol V_j^{-1}]=\tfrac{C_j}{C_j-G-2}\boldsymbol\Omega_j$ for $C_j\boldsymbol V_j\sim W_G(C_j-1,\boldsymbol\Sigma)$; factor $\to2$ at $C_j=2G$",
            "kind": "numeric",
            "status": "pass" if ok else "fail",
            "detail": r"$G=1$ exact from the $\chi^2$ density; $G=3$, $C_j=15$, 40000 draws: max relative error "
            + f"{rel:.4f}"
            + r" (tolerance 0.03)",
            "tex_line": "",
        }
    ]


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

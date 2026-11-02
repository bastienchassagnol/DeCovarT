# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V23: multivariate Fisher information for p and the Bhattacharyya overlap
of two Gaussian type laws, J = 2 (note, appendix A.4, eqs fisher-multi and
bhattacharyya).

- V23a (exact): y ~ N(c m(p), phi S(p)), m(p) = p mu1 + (1-p) mu2,
  S(p) = p Sigma1 + (1-p) Sigma2. The Fisher information is derived as
  -E[d^2 l / dp^2] from the Gaussian log-likelihood itself, the expectation
  taken by replacing r r' with V = phi S(p) (E r = 0). Compared with
  n (c^2/phi = n) Delta' S^-1 Delta + 1/2 tr[(S^-1 (Sigma1 - Sigma2))^2], for
  G = 2 with rational covariances and symbolic p, c, phi.
- V23b (numeric): D_B = 1/8 Delta' Sbar^-1 Delta + 1/2 log(det Sbar /
  sqrt(det S1 det S2)) against -log of the numerically integrated
  Bhattacharyya coefficient int sqrt(f1 f2), G = 1 (mpmath quadrature) and
  G = 2 (tensor grid); equal covariances give Delta' S^-1 Delta / 8 (V8a).
"""

from __future__ import annotations

import json

import mpmath as mp
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


def _db(d, S1, S2):
    Sb = (S1 + S2) / 2
    return d @ np.linalg.solve(Sb, d) / 8 + 0.5 * np.log(np.linalg.det(Sb) / np.sqrt(np.linalg.det(S1) * np.linalg.det(S2)))


def run() -> list[dict]:
    out = []

    # ---- V23a ----------------------------------------------------------
    p, c, phi = sp.symbols("p c phi", positive=True)
    mu1 = sp.Matrix([3, 1])
    mu2 = sp.Matrix([1, 2])
    S1 = sp.Matrix([[2, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    S2 = sp.Matrix([[1, -sp.Rational(1, 3)], [-sp.Rational(1, 3), 3]])
    m = c * (p * mu1 + (1 - p) * mu2)
    S = p * S1 + (1 - p) * S2
    V = phi * S
    y = sp.Matrix(sp.symbols("y1 y2", real=True))
    r = y - m
    ell = -sp.log(V.det()) / 2 - (r.T * V.inv() * r)[0, 0] / 2
    d2 = sp.expand(sp.diff(ell, p, 2))
    # expectation: y_a y_b -> V_ab + m_a m_b, y_a -> m_a
    poly = sp.Poly(d2, *y)
    E = 0
    for mon, coef in poly.terms():
        e1, e2 = mon
        if e1 + e2 == 0:
            E += coef
        elif e1 + e2 == 1:
            E += coef * (m[0] if e1 else m[1])
        elif (e1, e2) == (2, 0):
            E += coef * (V[0, 0] + m[0] ** 2)
        elif (e1, e2) == (0, 2):
            E += coef * (V[1, 1] + m[1] ** 2)
        else:
            E += coef * (V[0, 1] + m[0] * m[1])
    I_derived = sp.simplify(-E)
    n = c**2 / phi
    Delta = mu1 - mu2
    Si = S.inv()
    Mtr = Si * (S1 - S2)
    I_formula = n * (Delta.T * Si * Delta)[0, 0] + (Mtr * Mtr).trace() / 2
    ok_a = sp.simplify(I_derived - I_formula) == 0
    out.append(
        _rec(
            "V23a",
            r"Multivariate Fisher information for $p$, $J=2$: $I(p)=n_i\bm\Delta^\top\bm S(p)^{-1}\bm\Delta+"
            r"\tfrac12\operatorname{tr}[\{\bm S(p)^{-1}(\bm\Sigma_{W,1}-\bm\Sigma_{W,2})\}^2]$",
            "exact",
            "pass" if ok_a else "fail",
            r"$-\mathbb E\,\partial^2\ell/\partial p^2$ from the Gaussian log-likelihood, $G=2$, rational covariances, "
            r"symbolic $p$, $c$, $\phi$ with $n_i=c^2/\phi$",
            "",
        )
    )

    # ---- V23b ----------------------------------------------------------
    mp.mp.dps = 30
    worst = 0.0
    # G = 1
    for (m1, m2, s1, s2) in [(0.0, 1.5, 1.0, 2.5), (2.0, 2.0, 0.5, 3.0), (-1.0, 1.0, 1.0, 1.0)]:
        f = lambda x: mp.sqrt(mp.npdf(x, m1, mp.sqrt(s1)) * mp.npdf(x, m2, mp.sqrt(s2)))
        bc = mp.quad(f, [-mp.inf, 0, mp.inf])
        db_num = float(-mp.log(bc))
        db_form = float(_db(np.array([m1 - m2]), np.array([[s1]]), np.array([[s2]])))
        worst = max(worst, abs(db_num - db_form))
    # G = 2 on a grid
    m1v, m2v = np.array([0.0, 0.0]), np.array([1.2, -0.6])
    A1 = np.array([[1.0, 0.6], [0.6, 1.5]])
    A2 = np.array([[2.0, -0.4], [-0.4, 0.8]])
    xs = np.linspace(-9, 9, 1201)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    P = np.stack([X, Y], axis=-1)

    def dens(mu, S):
        Si = np.linalg.inv(S)
        d = P - mu
        q = np.einsum("...i,ij,...j->...", d, Si, d)
        return np.exp(-q / 2) / (2 * np.pi * np.sqrt(np.linalg.det(S)))

    h = xs[1] - xs[0]
    bc2 = np.sqrt(dens(m1v, A1) * dens(m2v, A2)).sum() * h * h
    err2 = abs(-np.log(bc2) - _db(m1v - m2v, A1, A2))
    worst = max(worst, err2)
    # equal covariances: reduces to Mahalanobis / 8
    eq = abs(_db(m1v - m2v, A1, A1) - (m1v - m2v) @ np.linalg.solve(A1, m1v - m2v) / 8)
    ok_b = worst < 1e-6 and eq < 1e-12
    out.append(
        _rec(
            "V23b",
            r"Bhattacharyya distance of two Gaussian laws: $D_B=\tfrac18\bm\Delta^\top\bar{\bm\Sigma}^{-1}\bm\Delta+"
            r"\tfrac12\log\{\det\bar{\bm\Sigma}/(\det\bm\Sigma_1\det\bm\Sigma_2)^{1/2}\}$; equal covariances give "
            r"$\bm\Delta^\top\bm\Sigma^{-1}\bm\Delta/8$",
            "numeric",
            "pass" if ok_b else "fail",
            rf"$-\log\int\sqrt{{f_1f_2}}$ by quadrature ($G=1$, three cases, 30 digits) and on a $1201^2$ grid "
            rf"($G=2$); largest difference {worst:.1e}",
            "",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

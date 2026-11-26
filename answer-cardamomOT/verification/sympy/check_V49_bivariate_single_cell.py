# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V49: the bivariate single-cell scenario (vignette R3-bivariate-single-cell,
script scripts/fig05_bivariate_single_cell.R).

Model at known scales (a = 1, so c = phi = n): y ~ N(n mu p, n S_W(p)),
S_W(p) = p1 Sigma_1 + (1 - p1) Sigma_2, G = J = 2.

Derivation routes (nothing is copied from the vignette):
* V49a: S_W(p)/n at p = (1/2, 1/2), n = 1000 and Sigma_j = 500 R_j is
  compared with the population covariance sum_j p_j^2 R_j of fig02.
* V49b: E[(y - m)' A (y - m)] = tr(A E[(y - m)(y - m)']) with
  E[(y - m)(y - m)'] = V* + (m* - m)(m* - m)' (first and second moments of
  y), so the expected Gaussian log-density is built from moments and then
  compared with the closed form coded in .expected_loglik_single_cell_known().
* V49c: the Fisher information of p1 is minus the second derivative of that
  expected log-likelihood at p1 = p1*, compared with the split
  n d' S_W^-1 d + tr[(S_W^-1 (Sigma_1 - Sigma_2))^2] / 2.
* V49d: its first derivative vanishes at p1 = p1* (the maximiser is p*).
"""

from __future__ import annotations

import json
import random

import sympy as sp


def _sym2(prefix):
    a, b, c = sp.symbols(f"{prefix}11 {prefix}12 {prefix}22", real=True)
    return sp.Matrix([[a, b], [b, c]])


def run():
    records = []

    # V49a: calibration of the per-cell variance.
    r1, r2 = sp.symbols("rho1 rho2", real=True)
    R1 = sp.Matrix([[1, r1], [r1, 1]])
    R2 = sp.Matrix([[1, r2], [r2, 1]])
    half = sp.Rational(1, 2)
    single_cell = (half * 500 * R1 + half * 500 * R2) / 1000
    population = half**2 * R1 + half**2 * R2
    ok_a = sp.simplify(single_cell - population) == sp.zeros(2, 2)
    records.append(
        {
            "id": "V49a",
            "claim": r"With $\boldsymbol\Sigma_{W,j}=500\,\boldsymbol R_j$, $n=1000$ and $\boldsymbol p=(1/2,1/2)$: $\boldsymbol S_W(\boldsymbol p)/n=\sum_j p_j^2\boldsymbol R_j$",
            "kind": "exact",
            "status": "pass" if ok_a else "fail",
            "detail": r"symbolic in $\varrho_1,\varrho_2$; residual $\boldsymbol 0$",
            "tex_line": "",
        }
    )

    # Generic 2 x 2 model for V49b-d.
    S1, S2 = _sym2("s"), _sym2("t")
    m1, m2, m3, m4 = sp.symbols("mu11 mu21 mu12 mu22", real=True)
    mu = sp.Matrix([[m1, m3], [m2, m4]])
    p1, q1, n = sp.symbols("p1 q1 n", positive=True)

    def moments(x):
        p = sp.Matrix([x, 1 - x])
        return n * mu * p, n * (x * S1 + (1 - x) * S2)

    m_star, V_star = moments(q1)
    m, V = moments(p1)
    A = V.inv()
    gap = m_star - m
    second = V_star + gap * gap.T
    expected_from_moments = -sp.Rational(1, 2) * (
        sp.log(V.det()) + (A * second).trace()
    )
    S = p1 * S1 + (1 - p1) * S2
    S_star = q1 * S1 + (1 - q1) * S2
    g = mu * sp.Matrix([q1 - p1, p1 - q1])
    closed_form = -sp.Rational(1, 2) * (
        2 * sp.log(n)
        + sp.log(S.det())
        + (S.inv() * S_star).trace()
        + n**2 / n * (g.T * S.inv() * g)[0, 0]
    )
    diff_b = expected_from_moments - closed_form

    # Exact for the rational part; log det V = 2 log n + log det S for G = 2.
    rng = random.Random(20261125)
    worst_b = 0
    worst_c = 0
    worst_d = 0
    for _ in range(25):
        vals = {}
        for M in (S1, S2):
            a = sp.Rational(rng.randint(2, 9), 1)
            c = sp.Rational(rng.randint(2, 9), 1)
            b = sp.Rational(rng.randint(-9, 9), 10)
            vals.update({M[0, 0]: a, M[0, 1]: b, M[1, 1]: c})
        vals.update(
            {
                m1: rng.randint(10, 30),
                m2: rng.randint(10, 30),
                m3: rng.randint(10, 30),
                m4: rng.randint(10, 30),
                n: rng.choice([100, 1000, 10000]),
                q1: sp.Rational(rng.randint(1, 9), 10),
            }
        )
        at_p = dict(vals)
        at_p[p1] = sp.Rational(rng.randint(1, 9), 10)
        worst_b = max(worst_b, abs(sp.N(diff_b.subs(at_p), 40)))
        # V49c and V49d at p1 = q1.
        loglik_q = closed_form.subs(vals)
        first = sp.diff(loglik_q, p1).subs(p1, vals[q1])
        info = -sp.diff(loglik_q, p1, 2).subs(p1, vals[q1])
        Sq = S_star.subs(vals)
        d = mu.subs(vals) * sp.Matrix([1, -1])
        Delta = (S1 - S2).subs(vals)
        split = vals[n] * (d.T * Sq.inv() * d)[0, 0] + sp.Rational(1, 2) * (
            (Sq.inv() * Delta) ** 2
        ).trace()
        worst_c = max(worst_c, abs(sp.nsimplify(sp.simplify(info - split))))
        worst_d = max(worst_d, abs(sp.simplify(first)))
    records.append(
        {
            "id": "V49b",
            "claim": r"$\mathrm E_{\boldsymbol p^\star}\ell(\boldsymbol p)=-\tfrac12\{G\log\phi+\log|\boldsymbol S_W(\boldsymbol p)|+\mathrm{tr}(\boldsymbol S_W(\boldsymbol p)^{-1}\boldsymbol S_W(\boldsymbol p^\star))+(c^2/\phi)\,\boldsymbol g^\top\boldsymbol S_W(\boldsymbol p)^{-1}\boldsymbol g\}$, $\boldsymbol g=\tilde{\boldsymbol\mu}(\boldsymbol p^\star-\boldsymbol p)$, up to $-\tfrac G2\log2\pi$",
            "kind": "rational-sampling",
            "status": "pass" if worst_b < sp.Float("1e-30") else "fail",
            "detail": rf"25 random rational $2\times2$ models, $c=\phi=n$; largest residual ${sp.latex(sp.Float(worst_b, 3))}$ (mpmath, 40 digits)",
            "tex_line": "",
        }
    )
    records.append(
        {
            "id": "V49c",
            "claim": r"$I(p_1^\star)=n\,\boldsymbol d^\top\boldsymbol S_W^{-1}\boldsymbol d+\tfrac12\mathrm{tr}[\{\boldsymbol S_W^{-1}(\boldsymbol\Sigma_{W,1}-\boldsymbol\Sigma_{W,2})\}^2]$ at known scales $c=\phi=n$",
            "kind": "rational-sampling",
            "status": "pass" if worst_c == 0 else "fail",
            "detail": rf"minus the second derivative of V49b at $p_1=p_1^\star$; 25 random rational models; largest residual ${sp.latex(worst_c)}$",
            "tex_line": "",
        }
    )
    records.append(
        {
            "id": "V49d",
            "claim": r"The expected log-likelihood of V49b is stationary at $\boldsymbol p=\boldsymbol p^\star$",
            "kind": "rational-sampling",
            "status": "pass" if worst_d == 0 else "fail",
            "detail": rf"first derivative in $p_1$ at $p_1^\star$; 25 random rational models; largest value ${sp.latex(worst_d)}$",
            "tex_line": "",
        }
    )
    return records


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

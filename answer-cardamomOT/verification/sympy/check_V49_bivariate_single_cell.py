# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V49: the bivariate single-cell scenario (vignette R3-bivariate-single-cell,
script scripts/fig05_bivariate_single_cell.R).

Model at known scales (a = 1, so c = phi = n): y ~ N(n mu p, n S_W(p)),
S_W(p) = p1 Sigma_1 + (1 - p1) Sigma_2, G = J = 2.

Derivation routes (nothing is copied from the vignette):
* V49a: with Sigma_W,j = n R_j, the covariance of y / n is
  n S_W(p) / n^2 = p1 R1 + (1 - p1) R2 for every n (n cancels).
* V49b: E[(y - m)' A (y - m)] = tr(A E[(y - m)(y - m)']) with
  E[(y - m)(y - m)'] = V* + (m* - m)(m* - m)' (first and second moments of
  y), so the expected Gaussian log-density is built from moments and then
  compared with the closed form coded in .expected_loglik_single_cell_known().
* V49c: the Fisher information of p1 is minus the second derivative of that
  expected log-likelihood at p1 = p1*, compared with the split
  n d' S_W^-1 d + tr[(S_W^-1 (Sigma_1 - Sigma_2))^2] / 2.
* V49d: its first derivative vanishes at p1 = p1* (the maximiser is p*).
* V51a: per-cell scale, y / n ~ N(mu p, V(p)) with the population law
  V = p1^2 R1 + (1 - p1)^2 R2 (fig02) or the single-cell law
  V = p1 R1 + (1 - p1) R2: minus the second derivative of the expected
  log-likelihood (built from moments, as V49b) at p1* equals
  d' V^-1 d + tr[(V^-1 V')^2] / 2 for both laws.
* V51b: with R1 = R2 = R (G = 2), the covariance part is 0 under the
  single-cell law and 4 (2 p1 - 1)^2 / (p1^2 + (1 - p1)^2)^2 under the
  population law.
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

    # V49a: component matching Sigma_W,j = n R_j.
    r1, r2 = sp.symbols("rho1 rho2", real=True)
    R1 = sp.Matrix([[1, r1], [r1, 1]])
    R2 = sp.Matrix([[1, r2], [r2, 1]])
    x, nn = sp.symbols("x n", positive=True)
    per_cell = nn * (x * nn * R1 + (1 - x) * nn * R2) / nn**2
    ok_a = sp.simplify(per_cell - (x * R1 + (1 - x) * R2)) == sp.zeros(2, 2)
    records.append(
        {
            "id": "V49a",
            "claim": r"With $\boldsymbol\Sigma_{W,j}=n\boldsymbol R_j$ and $c=\phi=n$: $\mathrm{Var}(\boldsymbol y/n)=p_1\boldsymbol R_1+p_2\boldsymbol R_2$ for every $n$",
            "kind": "exact",
            "status": "pass" if ok_a else "fail",
            "detail": r"symbolic in $p_1,n,\varrho_1,\varrho_2$; residual $\boldsymbol 0$",
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
    # V51a: Fisher split on the per-cell scale for both laws.
    rr1, rr2 = sp.symbols("r1 r2", real=True)
    Q1 = sp.Matrix([[1, rr1], [rr1, 1]])
    Q2 = sp.Matrix([[1, rr2], [rr2, 1]])
    laws = {
        "linear": (lambda x: x * Q1 + (1 - x) * Q2),
        "quadratic": (lambda x: x**2 * Q1 + (1 - x) ** 2 * Q2),
    }
    worst_51a = 0
    for law, Vf in laws.items():
        Vp = Vf(p1)
        Vq = Vf(q1)
        g51 = mu * sp.Matrix([q1 - p1, p1 - q1])
        Ap = Vp.inv()
        expected = -sp.Rational(1, 2) * (
            sp.log(Vp.det()) + (Ap * (Vq + g51 * g51.T)).trace()
        )
        for _ in range(10):
            vals = {
                rr1: sp.Rational(rng.randint(-9, 9), 10),
                rr2: sp.Rational(rng.randint(-9, 9), 10),
                m1: rng.randint(10, 40),
                m2: rng.randint(10, 40),
                m3: rng.randint(10, 40),
                m4: rng.randint(10, 40),
                q1: sp.Rational(rng.randint(1, 9), 10),
            }
            e = expected.subs({k: v for k, v in vals.items() if k != q1})
            e = e.subs(q1, vals[q1])
            info = -sp.diff(e, p1, 2).subs(p1, vals[q1])
            Vs = Vf(vals[q1]).subs(vals)
            dV = sp.diff(Vf(p1), p1).subs(p1, vals[q1]).subs(vals)
            dvec = mu.subs(vals) * sp.Matrix([1, -1])
            split = (dvec.T * Vs.inv() * dvec)[0, 0] + sp.Rational(1, 2) * (
                (Vs.inv() * dV) ** 2
            ).trace()
            worst_51a = max(worst_51a, abs(sp.simplify(info - split)))
    records.append(
        {
            "id": "V51a",
            "claim": r"For $\boldsymbol y/n\sim\mathcal N(\tilde{\boldsymbol\mu}\boldsymbol p,\boldsymbol V(p_1))$ with $\boldsymbol V=\sum_jp_j\boldsymbol R_j$ or $\sum_jp_j^2\boldsymbol R_j$: $I(p_1^\star)=\boldsymbol d^\top\boldsymbol V^{-1}\boldsymbol d+\tfrac12\mathrm{tr}[(\boldsymbol V^{-1}\boldsymbol V')^2]$",
            "kind": "rational-sampling",
            "status": "pass" if worst_51a == 0 else "fail",
            "detail": rf"minus the second derivative of the expected log-likelihood at $p_1^\star$; 10 random rational models per law; largest residual ${sp.latex(worst_51a)}$",
            "tex_line": "",
        }
    )

    # V51b: equal covariances.
    rho = sp.symbols("rho", real=True)
    R = sp.Matrix([[1, rho], [rho, 1]])
    cov_part = {}
    for law, w in {"linear": (p1, 1 - p1), "quadratic": (p1**2, (1 - p1) ** 2)}.items():
        Vp = w[0] * R + w[1] * R
        dV = sp.diff(Vp, p1)
        cov_part[law] = sp.simplify(sp.Rational(1, 2) * ((Vp.inv() * dV) ** 2).trace())
    target = 4 * (2 * p1 - 1) ** 2 / (p1**2 + (1 - p1) ** 2) ** 2
    ok_b = cov_part["linear"] == 0 and sp.simplify(cov_part["quadratic"] - target) == 0
    records.append(
        {
            "id": "V51b",
            "claim": r"With $\boldsymbol R_1=\boldsymbol R_2$ ($G=2$): covariance part $0$ under $\sum_jp_j\boldsymbol R_j$ and $4(2p_1-1)^2/\{p_1^2+(1-p_1)^2\}^2$ under $\sum_jp_j^2\boldsymbol R_j$",
            "kind": "exact",
            "status": "pass" if ok_b else "fail",
            "detail": r"symbolic in $p_1,\varrho$; residuals $0$",
            "tex_line": "",
        }
    )
    return records


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

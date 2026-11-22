# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V47: claims of the reply to Elias's review
(review-decovart-bastien-elias-marielle-anais/answer_elias.qmd).

- V47a (exact): for y ~ N(m(p), V(p)) the score splits into the GLS
  quasi-score m'^T V^-1 r and a covariance score
  1/2 [r^T V^-1 V' V^-1 r - tr(V^-1 V')]; the covariance score has mean zero,
  the two scores are uncorrelated, Var(GLS score) = I_mean = m'^T V^-1 m' and
  Var(cov score) = I_cov = 1/2 tr[(V^-1 V')^2]. Expectations are taken over
  z ~ N(0, I) with r = L z (L the Cholesky factor of V), at a rational point,
  with exact radicals. So the GLS estimating equation has asymptotic variance
  I_mean^-1 and the ML one (I_mean + I_cov)^-1.
- V47b (exact): single-cell law y ~ N(c mu p, phi S_W(p)), J = 2, G = 2,
  symbolic n = c^2/phi: I_mean is proportional to n, I_cov does not depend on
  n, hence I_cov / I_mean = O(1/n) and n * I_cov / I_mean has a finite limit.
- V47c (exact): two-level sum. Each type j contributes m_j = n p_j cells
  x = mu_j + b_j + e_c with a shared type-level offset b_j (variance Sigma_B,j)
  and cell noise e_c (variance Sigma_W,j); y = a sum x gives
  Cov(y) = a^2 n^2 sum p_j^2 Sigma_B,j + a^2 n sum p_j Sigma_W,j
  = c^2 S_B(p) + phi S_W(p) (Elias's formula, eq-av-four-layer without D_i);
  G = 1: f_cov = I_cov / (I_mean + I_cov) has a positive limit as n -> inf, and
  the effective number n_hat = c^2 / phi_hat (phi_hat = Var(y) / S_W) tends to
  S_W / S_B. Relative noise sd/mean scales as n^(-1/2) without Sigma_B.
- V47d (exact): MEAD counter-example with perfect markers:
  gamma~ Lambda~ U p~ = gamma Lambda U p for any delta > 0.
- V47e (exact): a gene with zero variance in type 1 makes det Sigma_W,1 = 0 but
  S_W(p) = p Sigma_W,1 + (1-p) Sigma_W,2 stays positive definite on (0, 1).
- V47f (exact): Gamma(alpha, beta)-Poisson posterior mean after zero counts
  over exposure L is alpha / (beta + L): 1/L for the flat prior (alpha = 1,
  beta -> 0), 1/(2L) for Jeffreys (alpha = 1/2).
- V47g (exact): two genes, correlation rho, loadings x. The variance reported
  under independence is (x^T x)^-1 sigma^2 while the actual variance of that
  estimator is (x^T x)^-1 x^T Sigma x (x^T x)^-1: ratio 1 + rho for x = (1, 1),
  1 - rho for x = (1, -1). Positive correlation makes the univariate interval
  too narrow in the first case and too wide in the second.
"""

from __future__ import annotations

import sympy as sp


def _rec(id_, claim, kind, status, detail, tex_line="", detail_full=None):
    return {
        "id": id_,
        "claim": claim,
        "kind": kind,
        "status": status,
        "detail": detail,
        "detail_full": detail if detail_full is None else detail_full,
        "tex_line": tex_line,
    }


def _gauss_expect(expr, z):
    """E over independent standard normals z of a polynomial expression."""
    poly = sp.Poly(sp.expand(expr), *z)
    mom = {0: 1, 1: 0, 2: 1, 3: 0, 4: 3, 5: 0, 6: 15}
    tot = 0
    for mon, coef in poly.terms():
        w = coef
        for k in mon:
            w *= mom[k]
        tot += w
    return sp.nsimplify(sp.simplify(tot)) if tot != 0 else sp.Integer(0)


def run() -> list[dict]:
    out = []

    # ---- V47a: score split, expectations exact ---------------------------
    p = sp.symbols("p", positive=True)
    c, phi = sp.Integer(2), sp.Rational(3, 2)
    mu1, mu2 = sp.Matrix([3, 1]), sp.Matrix([1, 2])
    S1 = sp.Matrix([[2, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    S2 = sp.Matrix([[1, -sp.Rational(1, 3)], [-sp.Rational(1, 3), 3]])
    m = c * (p * mu1 + (1 - p) * mu2)
    V = phi * (p * S1 + (1 - p) * S2)
    dm, dV = m.diff(p), V.diff(p)
    p0 = sp.Rational(1, 3)
    V0, dm0, dV0 = V.subs(p, p0), dm.subs(p, p0), dV.subs(p, p0)
    W0 = V0.inv()
    L0 = V0.cholesky(hermitian=False)
    z = sp.symbols("z1 z2", real=True)
    r = L0 * sp.Matrix(z)
    # full score from the log-likelihood, evaluated at p0 with y = m(p0) + r
    yv = sp.Matrix(sp.symbols("y1 y2", real=True))
    ell = -sp.log(V.det()) / 2 - ((yv - m).T * V.inv() * (yv - m))[0, 0] / 2
    score = sp.diff(ell, p).subs(p, p0).subs(dict(zip(yv, m.subs(p, p0) + r)))
    u_gls = (dm0.T * W0 * r)[0, 0]
    u_cov = ((r.T * W0 * dV0 * W0 * r)[0, 0] - (W0 * dV0).trace()) / 2
    split_ok = sp.simplify(sp.expand(score - u_gls - u_cov)) == 0
    I_mean = sp.nsimplify(sp.simplify((dm0.T * W0 * dm0)[0, 0]))
    I_cov = sp.nsimplify(sp.simplify(((W0 * dV0) ** 2).trace() / 2))
    e_cov = _gauss_expect(u_cov, z)
    e_cross = _gauss_expect(u_gls * u_cov, z)
    v_gls = _gauss_expect(u_gls**2, z)
    v_cov = _gauss_expect(u_cov**2, z)
    ok_a = (
        split_ok
        and e_cov == 0
        and e_cross == 0
        and sp.simplify(v_gls - I_mean) == 0
        and sp.simplify(v_cov - I_cov) == 0
    )
    out.append(
        _rec(
            "V47a",
            r"Gaussian score = GLS quasi-score $\bm m'^\top\bm V^{-1}\bm r$ + covariance score "
            r"$\tfrac12\{\bm r^\top\bm V^{-1}\bm V'\bm V^{-1}\bm r-\operatorname{tr}(\bm V^{-1}\bm V')\}$; "
            r"uncorrelated, with variances $I_{\mathrm{mean}}$ and $I_{\mathrm{cov}}=\tfrac12\operatorname{tr}\{(\bm V^{-1}\bm V')^2\}$",
            "exact",
            "pass" if ok_a else "fail",
            rf"$G=2$, $J=2$, $p=1/3$, $c=2$, $\phi=3/2$; Gaussian moments of $\bm z$ with $\bm r=\bm L\bm z$: "
            rf"$I_{{\mathrm{{mean}}}}={sp.latex(I_mean)}$, $I_{{\mathrm{{cov}}}}={sp.latex(I_cov)}$",
        )
    )

    # ---- V47b: I_cov / I_mean = O(1/n) -----------------------------------
    n, cc, pp = sp.symbols("n c p", positive=True)
    phis = cc**2 / n
    S = pp * S1 + (1 - pp) * S2
    Si = S.inv()
    Delta = mu1 - mu2
    Im = sp.simplify(cc**2 / phis * (Delta.T * Si * Delta)[0, 0])
    Ic = sp.simplify(((Si * (S1 - S2)) ** 2).trace() / 2)
    lim = sp.simplify(sp.limit(n * Ic / Im, n, sp.oo))
    ok_b = (
        sp.simplify(sp.diff(Im / n, n)) == 0
        and sp.diff(Ic, n) == 0
        and lim.is_finite is not False
        and sp.limit(Ic / Im, n, sp.oo) == 0
    )
    out.append(
        _rec(
            "V47b",
            r"Single-cell law with $\phi=c^2/n$: $I_{\mathrm{mean}}\propto n$, $I_{\mathrm{cov}}$ free of $n$, so "
            r"$I_{\mathrm{cov}}/I_{\mathrm{mean}}=O(1/n)$ and ML and GLS agree to first order as $n\to\infty$",
            "exact",
            "pass" if ok_b else "fail",
            r"$G=2$, $J=2$, symbolic $p$, $c$, $n$; $n\,I_{\mathrm{cov}}/I_{\mathrm{mean}}$ has a finite limit",
        )
    )

    # ---- V47c: two-level covariance, f_cov limit, effective n ------------
    a, nn = sp.symbols("a n", positive=True)
    p1 = sp.symbols("p1", positive=True)
    B1, B2, W1, W2 = sp.symbols("B1 B2 W1 W2", positive=True)
    m1, m2 = nn * p1, nn * (1 - p1)
    # var of a * (sum over m_j cells of (b_j + e_c)) = a^2 (m_j^2 B_j + m_j W_j)
    var_sum = a**2 * (m1**2 * B1 + m1 * W1 + m2**2 * B2 + m2 * W2)
    c2, phi2 = a * nn, a**2 * nn
    SB = p1**2 * B1 + (1 - p1) ** 2 * B2
    SW = p1 * W1 + (1 - p1) * W2
    elias_ok = sp.expand(var_sum - (c2**2 * SB + phi2 * SW)) == 0
    # G = 1 information in p1 (mu1 != mu2 scalars)
    u1, u2 = sp.symbols("u1 u2", real=True)
    Vt = c2**2 * SB + phi2 * SW
    I_m = (c2 * (u1 - u2)) ** 2 / Vt
    I_c = (sp.diff(Vt, p1) / Vt) ** 2 / 2
    f_lim = sp.simplify(sp.limit(I_c / (I_m + I_c), nn, sp.oo))
    vals = {p1: sp.Rational(1, 3), B1: 1, B2: 2, W1: 3, W2: 5, u1: 4, u2: 1}
    f_num = f_lim.subs(vals)
    n_hat = sp.simplify(c2**2 / (Vt / SW))
    n_lim = sp.simplify(sp.limit(n_hat, nn, sp.oo))
    # relative noise without Sigma_B: sd / mean = sqrt(phi S_W)/(c m) ~ n^(-1/2)
    rel = sp.sqrt(phi2 * SW) / (c2 * (p1 * u1 + (1 - p1) * u2))
    ok_rel = sp.simplify(rel * sp.sqrt(nn) - sp.sqrt(SW) / (p1 * u1 + (1 - p1) * u2)) == 0
    ratio_500 = sp.sqrt(sp.Rational(10**4, 500))
    ok_c = elias_ok and f_num > 0 and f_num < 1 and sp.simplify(n_lim - SW / SB) == 0 and ok_rel
    out.append(
        _rec(
            "V47c",
            r"Two-level sum: $\operatorname{Cov}\bm y=a^2n^2\sum_jp_j^2\bm\Sigma_{B,j}+a^2n\sum_jp_j\bm\Sigma_{W,j}"
            r"=c^2\bm S_B(\bm p)+\phi\bm S_W(\bm p)$; with $\bm\Sigma_B\neq0$, $f_{\mathrm{cov}}$ keeps a positive limit "
            r"and $\hat n=c^2/\hat\phi\to S_W/S_B$ ($G=1$); without it the relative noise is $O(n^{-1/2})$",
            "exact",
            "pass" if ok_c else "fail",
            rf"$J=2$, symbolic; at $p_1=1/3$, $B=(1,2)$, $W=(3,5)$, $\mu=(4,1)$: $\lim f_{{\mathrm{{cov}}}}={sp.latex(f_num)}"
            rf"\approx{float(f_num):.3f}$; noise ratio $\sqrt{{10^4/500}}\approx{float(ratio_500):.2f}$",
        )
    )

    # ---- V47d: MEAD perfect-marker counter-example -----------------------
    lam = sp.symbols("l1:5", positive=True)  # genes 1,2 -> type 1; genes 3,4 -> type 2
    uu = sp.symbols("w1:5", positive=True)
    d1, d2, g, q = sp.symbols("delta1 delta2 gamma q", positive=True)
    U = sp.Matrix([[uu[0], 0], [uu[1], 0], [0, uu[2]], [0, uu[3]]])
    P = sp.Matrix([q, 1 - q])
    Lam = sp.diag(*lam)
    dk = [d1, d1, d2, d2]
    Lam_t = sp.diag(*[lam[i] / dk[i] for i in range(4)])
    tot = d1 * q + d2 * (1 - q)
    P_t = sp.Matrix([d1 * q / tot, d2 * (1 - q) / tot])
    g_t = g * tot
    ok_d = sp.simplify(g_t * Lam_t * U * P_t - g * Lam * U * P) == sp.zeros(4, 1)
    ok_d = ok_d and sp.simplify(P_t[0] - q) != 0
    out.append(
        _rec(
            "V47d",
            r"MEAD counter-example: with perfect markers, $\tilde\gamma\tilde{\bm\Lambda}\bm U\tilde{\bm p}="
            r"\gamma\bm\Lambda\bm U\bm p$ for $\tilde\lambda_g=\lambda_g/\delta_k$, "
            r"$\tilde p_k=\delta_kp_k/\sum_l\delta_lp_l$, $\tilde\gamma=\gamma\sum_l\delta_lp_l$",
            "exact",
            "pass" if ok_d else "fail",
            r"$G=4$ genes, $K=2$ types, symbolic $\lambda$, $u$, $\delta$, $\gamma$, $p$; $\tilde{\bm p}\neq\bm p$ when $\delta_1\neq\delta_2$",
        )
    )

    # ---- V47e: zero-variance gene in one type ----------------------------
    s, b, dd, rho = sp.symbols("s b d r", positive=True)
    t = sp.symbols("t", positive=True)
    Z1 = sp.Matrix([[0, 0], [0, s]])
    Z2 = sp.Matrix([[b, rho], [rho, dd]])
    Sp = t * Z1 + (1 - t) * Z2
    det_p = sp.factor(Sp.det())
    # det = (1-t) [ t b s + (1-t)(b d - r^2) ] > 0 when b d > r^2 and 0 < t < 1
    target = (1 - t) * (t * b * s + (1 - t) * (b * dd - rho**2))
    ok_e = Z1.det() == 0 and sp.simplify(det_p - target) == 0 and sp.simplify(Sp[0, 0] - (1 - t) * b) == 0
    out.append(
        _rec(
            "V47e",
            r"A gene with zero variance in type 1 gives $\det\bm\Sigma_{W,1}=0$, but "
            r"$\det\bm S_W(p)=(1-p)\{pbs+(1-p)(bd-r^2)\}>0$ on $(0,1)$: singular only on the face $p=1$",
            "exact",
            "pass" if ok_e else "fail",
            r"$G=2$, $\bm\Sigma_{W,1}=\operatorname{diag}(0,s)$, $\bm\Sigma_{W,2}=\begin{pmatrix}b&r\\r&d\end{pmatrix}$ positive definite",
        )
    )

    # ---- V47f: pseudo-count as a posterior mean --------------------------
    al, be, L, lmb = sp.symbols("alpha beta L lambda", positive=True)
    prior = lmb ** (al - 1) * sp.exp(-be * lmb)
    like0 = sp.exp(-lmb * L)  # Poisson probability of zero counts over exposure L
    post = prior * like0
    mean = sp.simplify(sp.integrate(lmb * post, (lmb, 0, sp.oo), conds="none") / sp.integrate(post, (lmb, 0, sp.oo), conds="none"))
    ok_f = (
        sp.simplify(mean - al / (be + L)) == 0
        and sp.simplify(mean.subs({al: 1, be: 0}) - 1 / L) == 0
        and sp.simplify(mean.subs({al: sp.Rational(1, 2), be: 0}) - 1 / (2 * L)) == 0
    )
    out.append(
        _rec(
            "V47f",
            r"Zero counts over exposure $L$ with a $\mathrm{Gamma}(\alpha,\beta)$ prior: posterior mean $\alpha/(\beta+L)$, "
            r"that is $1/L$ for the flat prior and $1/(2L)$ for Jeffreys",
            "exact",
            "pass" if ok_f else "fail",
            r"conjugate integral computed by SymPy, symbolic $\alpha$, $\beta$, $L$",
        )
    )

    # ---- V47g: univariate interval, sign depends on geometry -------------
    sg, rh = sp.symbols("sigma rho", real=True)
    Sig = sg**2 * sp.Matrix([[1, rh], [rh, 1]])
    res = {}
    for name, x in {"same": sp.Matrix([1, 1]), "opposite": sp.Matrix([1, -1])}.items():
        xtx = (x.T * x)[0, 0]
        reported = sg**2 / xtx
        actual = (x.T * Sig * x)[0, 0] / xtx**2
        res[name] = sp.simplify(actual / reported)
    ok_g = sp.simplify(res["same"] - (1 + rh)) == 0 and sp.simplify(res["opposite"] - (1 - rh)) == 0
    out.append(
        _rec(
            "V47g",
            r"Actual over reported variance of the independence (univariate) estimator: $1+\rho$ for loadings $(1,1)$, "
            r"$1-\rho$ for $(1,-1)$; positive correlation makes the interval too narrow or too wide depending on the loadings",
            "exact",
            "pass" if ok_g else "fail",
            r"two genes, $\bm\Sigma=\sigma^2\begin{pmatrix}1&\rho\\\rho&1\end{pmatrix}$, one coefficient",
        )
    )
    return out


if __name__ == "__main__":
    import json

    for rec in run():
        print(json.dumps({k: rec[k] for k in ("id", "status", "detail")}, ensure_ascii=False))

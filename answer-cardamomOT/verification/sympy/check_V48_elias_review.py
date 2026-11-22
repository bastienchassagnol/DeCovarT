# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V48: second batch of claims of the reply to Elias
(review-decovart-bastien-elias-marielle-anais/answer_elias.qmd).

Single-cell law y ~ N(c mu p, phi S(p)), S(p) = sum_j p_j Sigma_j, J = 2,
p = (p, 1 - p), contrast d = (1, -1).

- V48a (exact): with phi = c^2/n and y ~ N(c mu p0, phi S(p0)), the expected
  log-likelihood E l(p) is stationary at p0 with second derivative -I(p0),
  I = n I_mean_unit + I_cov (Fisher consistency of the small-dispersion limit).
- V48b (exact): non-concavity with known c, phi and distinct means:
  mu1 = (2,1), mu2 = (1,2), Sigma1 = I, Sigma2 = 4 I, c = phi = 1,
  y = mu (1/2, 1/2): l''(1/2) = 16/25 > 0 at the true composition.
- V48c (exact): two-layer bulk Cov y = c^2 S_B(p) + phi S_W(p), p known:
  E(phi_hat)/c^2 = tr(S_W^-1 S_B)/G + 1/n with phi_hat = Q/G, hence
  E(1/n_hat) = tr(S_W^-1 S_B)/G + 1/n (G = 2, symbolic).
- V48d (exact): cross information I_{d phi} = sum_j d_j tr(S^-1 Sigma_j)/(2 phi)
  and I_{dc} = c (mu d)' S^-1 mu p / phi, derived from -E of second
  derivatives; I_{d phi} = 0 when Sigma1 = Sigma2, non-zero otherwise.
- V48e (numeric): at the GLS fixed point mu' W r = 0, so m' W r = 0; at the
  maximum-likelihood estimate mu' W r = -(covariance score) != 0.
- V48f (exact): independence (univariate) estimator, actual over reported
  variance = x' D^-1 Sigma D^-1 x / x' D^-1 x; equals 1 + rho and 1 - rho for
  x = (1, 1) and (1, -1) with D = diag(Sigma).
- V48g (exact): one perfect marker gene (mean c m p, variance phi s p):
  information n m^2/(p s) + 1/(2 p^2); with a free gene-wise factor b the
  profiled information is 1/(2 p^2) (the O(n) mean part vanishes).
- V48h (exact): Sigma_j = S + eps_j E with p0 . eps = 0 keeps S(p0) = S, hence
  I_mean(p0) and I_{dc}(p0); I_cov = (d . eps)^2 tr((S^-1 E)^2)/2; if
  tr(S^-1 E) = 0 then I_{d phi} = 0, so the profiled information gains I_cov.
"""

from __future__ import annotations

import numpy as np
import sympy as sp


def _rec(id_, claim, kind, status, detail, tex_line=""):
    return {"id": id_, "claim": claim, "kind": kind, "status": status, "detail": detail,
            "detail_full": detail, "tex_line": tex_line}


def _expect_quadratic(expr, y, mean, cov):
    """E over y ~ N(mean, cov) of a polynomial of degree <= 2 in y."""
    poly = sp.Poly(sp.expand(expr), *y)
    tot = 0
    for mon, coef in poly.terms():
        idx = [k for k, e in enumerate(mon) for _ in range(e)]
        if len(idx) == 0:
            tot += coef
        elif len(idx) == 1:
            tot += coef * mean[idx[0]]
        else:
            a, b = idx
            tot += coef * (cov[a, b] + mean[a] * mean[b])
    return tot


def run() -> list[dict]:
    out = []
    p, p0, c, n = sp.symbols("p p0 c n", positive=True)
    y = sp.Matrix(sp.symbols("y1 y2", real=True))
    mu1, mu2 = sp.Matrix([3, 1]), sp.Matrix([1, 2])
    S1 = sp.Matrix([[2, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    S2 = sp.Matrix([[1, -sp.Rational(1, 3)], [-sp.Rational(1, 3), 3]])
    M = sp.Matrix.hstack(mu1, mu2)

    def loglik(pp, cc, phi, Sa, Sb, yy):
        S = phi * (pp * Sa + (1 - pp) * Sb)
        r = yy - cc * (pp * mu1 + (1 - pp) * mu2)
        return -sp.log(S.det()) / 2 - (r.T * S.inv() * r)[0, 0] / 2

    # ---- V48a ----------------------------------------------------------
    phi = c**2 / n
    ell = loglik(p, c, phi, S1, S2, y)
    mean0 = c * (p0 * mu1 + (1 - p0) * mu2)
    cov0 = phi * (p0 * S1 + (1 - p0) * S2)
    El = _expect_quadratic(ell, y, mean0, cov0)
    g1 = sp.simplify(sp.diff(El, p).subs(p, p0))
    h1 = sp.simplify(sp.diff(El, p, 2).subs(p, p0))
    S0 = p0 * S1 + (1 - p0) * S2
    Delta = mu1 - mu2
    I_th = n * (Delta.T * S0.inv() * Delta)[0, 0] + ((S0.inv() * (S1 - S2)) ** 2).trace() / 2
    ok_a = g1 == 0 and sp.simplify(h1 + I_th) == 0
    out.append(_rec(
        "V48a",
        r"Small dispersion: with $\phi_i=c_i^2/n_i$ the expected log-likelihood is stationary at $p_0$ with curvature "
        r"$-I(p_0)=-\{n_i\bm\Delta^\top\bm S_W^{-1}\bm\Delta+\tfrac12\operatorname{tr}[(\bm S_W^{-1}(\bm\Sigma_{W,1}-\bm\Sigma_{W,2}))^2]\}$",
        "exact", "pass" if ok_a else "fail",
        r"$G=2$, $J=2$, symbolic $p_0$, $c$, $n$; expectation of a quadratic form under the true law",
    ))

    # ---- V48b ----------------------------------------------------------
    m1, m2 = sp.Matrix([2, 1]), sp.Matrix([1, 2])
    yb = (m1 + m2) / 2
    v = p * 1 + (1 - p) * 4
    rb = yb - (p * m1 + (1 - p) * m2)
    lb = -sp.log((v * sp.eye(2)).det()) / 2 - (rb.T * rb)[0, 0] / (2 * v)
    l2 = sp.simplify(sp.diff(lb, p, 2).subs(p, sp.Rational(1, 2)))
    # numeric scan: local maxima of lb on (0, 1)
    f = sp.lambdify(p, lb, "numpy")
    grid = np.linspace(1e-4, 1 - 1e-4, 20001)
    vals = f(grid)
    argmax = float(grid[np.argmax(vals)])
    ok_b = l2 == sp.Rational(16, 25)
    out.append(_rec(
        "V48b",
        r"Single-cell log-likelihood not concave: $\bm\mu_{\cdot1}=(2,1)$, $\bm\mu_{\cdot2}=(1,2)$, "
        r"$\bm\Sigma_{W,1}=\mathbf I$, $\bm\Sigma_{W,2}=4\mathbf I$, $c=\phi=1$, $\bm y=\tilde{\bm\mu}(1/2,1/2)$: "
        r"$\ell''(1/2)=16/25>0$ at the true composition",
        "exact", "pass" if ok_b else "fail",
        rf"exact second derivative; grid maximiser on $(0,1)$ at $p\approx{argmax:.3f}$ (boundary $p=1$ attracts the fit)",
    ))

    # ---- V48c ----------------------------------------------------------
    pp = sp.symbols("pp", positive=True)
    B1 = sp.Matrix([[1, sp.Rational(1, 4)], [sp.Rational(1, 4), 2]])
    B2 = sp.Matrix([[3, 0], [0, 1]])
    SW = pp * S1 + (1 - pp) * S2
    SB = pp**2 * B1 + (1 - pp) ** 2 * B2
    covr = c**2 * SB + (c**2 / n) * SW
    EQ = (SW.inv() * covr).trace()  # E r' S_W^-1 r with E r = 0
    lhs = sp.simplify(EQ / 2 / c**2)
    rhs = sp.simplify((SW.inv() * SB).trace() / 2 + 1 / n)
    ok_c = sp.simplify(lhs - rhs) == 0
    out.append(_rec(
        "V48c",
        r"Effective number of cells: with $\hat\phi_i=Q/G$ and $\operatorname{Cov}\bm y=c_i^2\bm S_B+\phi_i\bm S_W$, "
        r"$\mathbb E(1/\hat n_i)=\mathbb E(\hat\phi_i)/c_i^2=\operatorname{tr}(\bm S_W^{-1}\bm S_B)/G+1/n_i$",
        "exact", "pass" if ok_c else "fail",
        r"$G=2$, $J=2$, symbolic $p$, $c$, $n$, composition known",
    ))

    # ---- V48d ----------------------------------------------------------
    ph = sp.symbols("phi", positive=True)
    ell_d = loglik(p, c, ph, S1, S2, y)
    mean_p = c * (p * mu1 + (1 - p) * mu2)
    cov_p = ph * (p * S1 + (1 - p) * S2)
    I_pphi = sp.simplify(-_expect_quadratic(sp.diff(ell_d, p, ph), y, mean_p, cov_p))
    I_pc = sp.simplify(-_expect_quadratic(sp.diff(ell_d, p, c), y, mean_p, cov_p))
    Sp = p * S1 + (1 - p) * S2
    I_pphi_th = ((Sp.inv() * S1).trace() - (Sp.inv() * S2).trace()) / (2 * ph)
    I_pc_th = c * (Delta.T * Sp.inv() * (p * mu1 + (1 - p) * mu2))[0, 0] / ph
    ell_eq = loglik(p, c, ph, S1, S1, y)
    I_eq = sp.simplify(-_expect_quadratic(sp.diff(ell_eq, p, ph), y, mean_p, ph * S1))
    ok_d = (sp.simplify(I_pphi - I_pphi_th) == 0 and sp.simplify(I_pc - I_pc_th) == 0
            and I_eq == 0 and sp.simplify(I_pphi.subs(p, sp.Rational(1, 2))) != 0)
    out.append(_rec(
        "V48d",
        r"Scales and composition are not orthogonal: $I_{d\phi}=\sum_jd_j\operatorname{tr}(\bm S_W^{-1}\bm\Sigma_{W,j})/(2\phi_i)$, "
        r"$I_{dc}=c_i(\tilde{\bm\mu}\bm d)^\top\bm S_W^{-1}\tilde{\bm\mu}\bm p/\phi_i$; $I_{d\phi}=0$ for identical $\bm\Sigma_{W,j}$ only",
        "exact", "pass" if ok_d else "fail",
        rf"$-\mathbb E$ of mixed second derivatives, $G=2$; at $p=1/2$: $I_{{d\phi}}={sp.latex(sp.simplify(I_pphi.subs(p, sp.Rational(1, 2))))}$",
    ))

    # ---- V48e ----------------------------------------------------------
    Mn = np.array([[3.0, 1.0], [1.0, 2.0]])
    A1 = np.array([[2.0, 0.5], [0.5, 1.0]])
    A2 = np.array([[1.0, -1 / 3], [-1 / 3, 3.0]])
    yv = np.array([2.6, 1.9])
    phin = 0.05

    def nll(pv, cv):
        Sn = phin * (pv * A1 + (1 - pv) * A2)
        r = yv - cv * (Mn @ np.array([pv, 1 - pv]))
        return 0.5 * np.log(np.linalg.det(Sn)) + 0.5 * r @ np.linalg.solve(Sn, r)

    best = min(((nll(a, b), a, b) for a in np.linspace(0.01, 0.99, 981) for b in np.linspace(0.5, 1.5, 1001)))
    _, ph_, ch_ = best
    # Newton polish with central finite differences, so the stationarity error is negligible
    xk, h = np.array([ph_, ch_]), 1e-5
    for _ in range(20):
        f0 = lambda z: nll(z[0], z[1])
        gr = np.array([(f0(xk + h * e) - f0(xk - h * e)) / (2 * h) for e in np.eye(2)])
        He = np.array([[(f0(xk + h * (e + u)) - f0(xk + h * (e - u)) - f0(xk - h * (e - u)) + f0(xk - h * (e + u))) / (4 * h * h)
                        for u in np.eye(2)] for e in np.eye(2)])
        xk = xk - np.linalg.solve(He, gr)
    ph_, ch_ = xk
    grad_norm = float(np.max(np.abs(gr)))
    W = np.linalg.inv(phin * (ph_ * A1 + (1 - ph_) * A2))
    r = yv - ch_ * (Mn @ np.array([ph_, 1 - ph_]))
    gls_part = Mn.T @ W @ r
    # GLS fixed point with W frozen at the ML point: beta = (M'WM)^-1 M'Wy
    beta = np.linalg.solve(Mn.T @ W @ Mn, Mn.T @ W @ yv)
    r_g = yv - Mn @ beta
    orth = (Mn @ beta) @ W @ r_g
    ok_e = abs(orth) < 1e-10 and np.max(np.abs(gls_part)) > 1e-3 and grad_norm < 1e-6
    out.append(_rec(
        "V48e",
        r"Whitened residuals are orthogonal to fitted values at the GLS fixed point ($\bm m^\top\bm W\bm r=0$), "
        r"not at the maximum-likelihood estimate, where $\tilde{\bm\mu}^\top\bm W\bm r\neq0$ balances the covariance score",
        "numeric", "pass" if ok_e else "fail",
        rf"$G=2$, $J=2$, ML, grid then Newton ($\hat p={ph_:.3f}$, $\hat c={ch_:.3f}$): $\max|\tilde{{\bm\mu}}^\top\bm W\bm r|={np.max(np.abs(gls_part)):.3g}$; "
        rf"GLS $|\bm m^\top\bm W\bm r|={abs(orth):.1e}$",
    ))

    # ---- V48f ----------------------------------------------------------
    sg, rho = sp.symbols("sigma rho", real=True)
    Sig = sg**2 * sp.Matrix([[1, rho], [rho, 1]])
    D = sp.diag(Sig[0, 0], Sig[1, 1])
    res = []
    for x in (sp.Matrix([1, 1]), sp.Matrix([1, -1])):
        a = (x.T * D.inv() * x)[0, 0]
        actual = (x.T * D.inv() * Sig * D.inv() * x)[0, 0] / a**2
        res.append(sp.simplify(actual / (1 / a)))
    ok_f = sp.simplify(res[0] - (1 + rho)) == 0 and sp.simplify(res[1] - (1 - rho)) == 0
    out.append(_rec(
        "V48f",
        r"Independence estimator: actual over reported variance $=\bm x^\top\bm D^{-1}\bm\Sigma\bm D^{-1}\bm x/\bm x^\top\bm D^{-1}\bm x$, "
        r"$1+\rho$ for $\bm x=(1,1)$ and $1-\rho$ for $\bm x=(1,-1)$",
        "exact", "pass" if ok_f else "fail", r"two genes, $\bm D=\operatorname{diag}(\bm\Sigma)$, symbolic $\sigma$, $\rho$",
    ))

    # ---- V48g ----------------------------------------------------------
    m, s, b = sp.symbols("m s b", positive=True)
    yg = sp.symbols("yg", real=True)
    phg = c**2 / n
    lg = -sp.log(phg * s * p) / 2 - (yg - c * b * m * p) ** 2 / (2 * phg * s * p)
    mg, vg = c * b * m * p, phg * s * p

    def info(u, w):
        return sp.simplify(-_expect_quadratic(sp.diff(lg, u, w), [yg], [mg], sp.Matrix([[vg]])))

    Ipp, Ipb, Ibb = info(p, p), info(p, b), info(b, b)
    full = sp.simplify(Ipp.subs(b, 1))
    prof = sp.simplify((Ipp - Ipb**2 / Ibb).subs(b, 1))
    ok_g = sp.simplify(full - (n * m**2 / (p * s) + 1 / (2 * p**2))) == 0 and sp.simplify(prof - 1 / (2 * p**2)) == 0
    out.append(_rec(
        "V48g",
        r"One perfect marker gene: information $n_i m^2/(p s)+1/(2p^2)$ on $p$; with a free gene-wise factor the profiled "
        r"information is $1/(2p^2)$, the $O(n_i)$ mean part vanishes",
        "exact", "pass" if ok_g else "fail", r"mean $c_i b m p$, variance $\phi_i s p$, $\phi_i=c_i^2/n_i$, Schur complement in $b$",
    ))

    # ---- V48h ----------------------------------------------------------
    e1, e2, q0 = sp.symbols("epsilon1 epsilon2 q0", real=True)
    Sbase = sp.Matrix([[2, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    E = sp.Matrix([[1, 0], [0, -2]])  # tr(Sbase^-1 E) computed below
    trSE = sp.simplify((Sbase.inv() * E).trace())
    Ez = E - trSE / (Sbase.inv()).trace() * sp.eye(2) * 0  # keep E; test both cases via trSE
    Sa, Sb_ = Sbase + e1 * E, Sbase + e2 * E
    cond = {e2: -q0 * e1 / (1 - q0)}
    Sp0 = (q0 * Sa + (1 - q0) * Sb_).subs(cond)
    keep = sp.simplify(Sp0 - Sbase) == sp.zeros(2, 2)
    Icov = sp.simplify((((Sbase.inv() * (Sa - Sb_)) ** 2).trace() / 2).subs(cond))
    Icov_th = sp.simplify((((e1 - e2) ** 2) * ((Sbase.inv() * E) ** 2).trace() / 2).subs(cond))
    Idphi = sp.simplify((((Sbase.inv() * Sa).trace() - (Sbase.inv() * Sb_).trace()) / 2).subs(cond))
    E0 = sp.Matrix([[1, 0], [0, 0]])
    E0 = E0 - (Sbase.inv() * E0).trace() / (Sbase.inv() * Sbase).trace() * Sbase  # tr(S^-1 E0) = 0
    Idphi0 = sp.simplify(((Sbase.inv() * (e1 - e2) * E0).trace() / 2).subs(cond))
    ok_h = keep and sp.simplify(Icov - Icov_th) == 0 and Idphi0 == 0 and sp.simplify(Idphi - (e1 - e2).subs(cond) * trSE / 2) == 0
    out.append(_rec(
        "V48h",
        r"$\bm\Sigma_{W,j}=\bm S+\varepsilon_j\bm E$ with $\bm p_0^\top\bm\varepsilon=0$ keeps $\bm S_W(\bm p_0)=\bm S$ (Mahalanobis term fixed) "
        r"while $I^{\mathrm{cov}}=\tfrac12(\bm d^\top\bm\varepsilon)^2\operatorname{tr}\{(\bm S^{-1}\bm E)^2\}$; "
        r"$\operatorname{tr}(\bm S^{-1}\bm E)=0$ gives $I_{d\phi}=0$",
        "exact", "pass" if ok_h else "fail", r"$G=2$, $J=2$, symbolic $\varepsilon$ and $p_0$; $\phi_i$ factor omitted (common to both cases)",
    ))

    # ---- V48i: calibration of t for a target f_cov --------------------------
    Mm, Kk, fs = sp.symbols("M K f", positive=True)
    tstar = sp.sqrt(fs * Mm / ((1 - fs) * Kk))
    fcov = Kk * tstar**2 / (Mm + Kk * tstar**2)
    # whitened identity Sigma_j(t) = S^{1/2}(I + t a_j H)S^{1/2}, checked as S + t a_j E with E = S^{1/2} H S^{1/2}
    tt, aj = sp.symbols("t a_j", real=True)
    Hs = sp.Matrix([[1, sp.Rational(1, 3)], [sp.Rational(1, 3), -1]])
    Sh = sp.Matrix([[2, 0], [0, 1]])  # diagonal S so that S^{1/2} is exact
    Shalf = sp.diag(sp.sqrt(2), 1)
    Ew = Shalf * Hs * Shalf
    ok_whiten = sp.simplify(Shalf * (sp.eye(2) + tt * aj * Hs) * Shalf - (Sh + tt * aj * Ew)) == sp.zeros(2, 2)
    ok_tr = sp.simplify(((Sh.inv() * Ew) ** 2).trace() - (Hs**2).trace()) == 0
    ok_i = sp.simplify(fcov - fs) == 0 and ok_whiten and ok_tr
    out.append(_rec(
        "V48i",
        r"Calibration: $t^\star=\sqrt{f^\star M/\{(1-f^\star)K\}}$ gives $f_{\mathrm{cov}}=f^\star$; "
        r"$\bm\Sigma_{W,j}(t)=\bm S^{1/2}(\mathbf I+ta_j\bm H)\bm S^{1/2}$ and $\operatorname{tr}\{(\bm S^{-1}\bm E)^2\}=\lVert\bm H\rVert_F^2$",
        "exact", "pass" if ok_i else "fail", r"symbolic $M$, $K$, $f^\star$; whitening identity at $G=2$",
    ))

    # ---- V48j: multi-mode tangent covariance information ---------------------
    J3 = 3
    Bh = sp.Matrix([[1, 1], [-1, 1], [0, -2]])
    Bh = sp.Matrix.hstack(*[Bh[:, k] / Bh[:, k].norm() for k in range(2)])
    q1, q2 = sp.symbols("t1 t2", real=True)
    Ab = sp.Matrix([[1, 2], [-1, 0], [0, -3]])  # columns a_q (p0 orthogonality irrelevant to the identity)
    Sm = sp.Matrix([[2, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    E1 = sp.Matrix([[1, 0], [0, -1]])
    E2 = sp.Matrix([[0, 1], [1, 0]])
    Es = [E1, E2]
    Tm = sp.diag(q1, q2)
    Sig = [Sm + q1 * Ab[j, 0] * E1 + q2 * Ab[j, 1] * E2 for j in range(J3)]
    Smi = Sm.inv()
    D = [sum((Bh[j, r] * Sig[j] for j in range(J3)), sp.zeros(2, 2)) for r in range(2)]
    Icov_direct = sp.Matrix(2, 2, lambda r, s: (Smi * D[r] * Smi * D[s]).trace() / 2)
    GE = sp.Matrix(2, 2, lambda a, b: (Smi * Es[a] * Smi * Es[b]).trace())
    Icov_formula = (Bh.T * Ab * Tm) * GE * (Bh.T * Ab * Tm).T / 2
    ok_j = sp.simplify(Icov_direct - Icov_formula) == sp.zeros(2, 2)
    out.append(_rec(
        "V48j",
        r"Several modes: $\mathbf I_T^{\mathrm{cov}}=\tfrac12(\mathbf B^\top\mathbf A\mathbf T)\mathbf G_E(\mathbf B^\top\mathbf A\mathbf T)^\top$ "
        r"with $\mathbf G_E(q,r)=\operatorname{tr}(\bm S^{-1}\bm E_q\bm S^{-1}\bm E_r)$",
        "exact", "pass" if ok_j else "fail", r"$J=3$, $Q=2$, $G=2$, normalised Helmert basis, symbolic $t_1$, $t_2$",
    ))
    return out


if __name__ == "__main__":
    for r_ in run():
        print(r_["id"], r_["status"], r_["detail"])

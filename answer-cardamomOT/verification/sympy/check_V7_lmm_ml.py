# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V7: "maximum likelihood is DeCovarT's fit" (side axis, Wiener-filter item).

The note writes the population model as a linear mixed model
    y = mu p + sum_j p_j eps_j + e,  eps_j ~ N(0, Sigma_j),  e ~ N(0, D),
and says that its maximum likelihood is DeCovarT's fit. DeCovarT's likelihood,
as stated in article/main_layout.tex (eq:convolution and
eq-loglikelihood-multivariate-gaussian, assumptions A1-A5), is
    y ~ N(mu p, Sigma(p)),  Sigma(p) = sum_j p_j^2 Sigma_j,
    l(p) = C - 1/2 log det Sigma(p) - 1/2 (y - mu p)' Sigma(p)^-1 (y - mu p),
with no extra measurement noise (A5).

Derivation routes (none re-types the mixed model's marginal):
- V7a, G = 1, J = 2: the latent profiles x_1, x_2 are integrated out of the
  joint density in log space with the one-dimensional Gaussian integral
  int exp(-A x^2 + B x) dx = sqrt(pi/A) exp(B^2/(4A)), itself obtained from
  sympy.integrate; the result is compared with DeCovarT's density at D = 0.
- V7a, G = 2, J = 2, general symmetric 2 x 2 Sigma_j and diagonal D: mean and
  covariance from the cumulant generating function of y (sum of independent
  Gaussian terms); the log-likelihood built from them is compared with l(p)
  at D = 0.
- V7b: one-gene counter-example with a nugget d = 1: the stationary points of
  both log-likelihoods in p are the real roots of polynomials, isolated
  exactly; the maximisers in (0, 1) differ.
- V7c: restricted maximum likelihood. With p the fixed effect, the error
  contrasts K'y (K'Delta = 0) have a law free of p only if K' V(p) K is; its
  derivative in p is computed for G = 3, J = 2.
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


def _gauss_lemma():
    """Return (A, B) -> log int exp(-A x^2 + B x) dx, derived by sympy."""
    x = sp.Symbol("x", real=True)
    A = sp.Symbol("A", positive=True)
    B = sp.Symbol("B", real=True)
    val = sp.simplify(sp.integrate(sp.exp(-A * x**2 + B * x), (x, -sp.oo, sp.oo), conds="none"))
    ok = sp.simplify(val - sp.sqrt(sp.pi / A) * sp.exp(B**2 / (4 * A))) == 0
    return ok, (lambda a_, b_: sp.log(val.subs({A: a_, B: b_}, simultaneous=True)))


def _integrate_out(logf, x, lemma):
    """log int exp(logf) dx for logf quadratic in x with negative leading term."""
    poly = sp.Poly(sp.expand(logf), x)
    A = -poly.coeff_monomial(x**2)
    B = poly.coeff_monomial(x)
    C = poly.coeff_monomial(1)
    return C + lemma(A, B)


def _lognorm(z, m, v):
    return -((z - m) ** 2) / (2 * v) - sp.log(2 * sp.pi * v) / 2


def _stationary_max(ell, p):
    """Maximiser of ell on (0, 1) among the real roots of the score numerator."""
    num = sp.numer(sp.together(sp.diff(ell, p)))
    roots = [r for r in sp.Poly(num, p).real_roots() if 0 < r < 1]
    best = max(roots, key=lambda r: sp.N(ell.subs(p, r), 30))
    return best, roots


def run() -> list[dict]:
    out: list[dict] = []

    # ---------------- V7a, G = 1, J = 2: explicit integration ----------------
    ok_lemma, lemma = _gauss_lemma()
    x1, x2, y, m1, m2 = sp.symbols("x1 x2 y m1 m2", real=True)
    s1, s2, d = sp.symbols("s1 s2 d", positive=True)
    p = sp.Symbol("p", positive=True)  # p_1 = p, p_2 = 1 - p
    joint = _lognorm(y, p * x1 + (1 - p) * x2, d) + _lognorm(x1, m1, s1) + _lognorm(x2, m2, s2)
    marg = _integrate_out(_integrate_out(joint, x1, lemma), x2, lemma)
    v_lmm = p**2 * s1 + (1 - p) ** 2 * s2 + d
    ok_g1_lmm = sp.simplify(sp.expand_log(marg - _lognorm(y, p * m1 + (1 - p) * m2, v_lmm), force=True)) == 0
    # DeCovarT, G = 1: l(p) = C - 1/2 log v - (y - m)^2 / (2 v), C = -1/2 log 2 pi
    ell_dec_1 = (
        -sp.log(2 * sp.pi) / 2
        - sp.log(p**2 * s1 + (1 - p) ** 2 * s2) / 2
        - (y - p * m1 - (1 - p) * m2) ** 2 / (2 * (p**2 * s1 + (1 - p) ** 2 * s2))
    )
    # marg has d in denominators; it equals the closed form for every d > 0
    # (ok_g1_lmm), which is continuous at d = 0, so the d -> 0 limit of marg is
    # the closed form at d = 0
    ell_lmm_1 = _lognorm(y, p * m1 + (1 - p) * m2, v_lmm)
    ok_g1_dec = ok_g1_lmm and sp.simplify(sp.expand_log(ell_lmm_1.subs(d, 0) - ell_dec_1, force=True)) == 0

    # ---------------- V7a, G = 2, J = 2: CGF route, general 2 x 2 ------------
    ts = sp.symbols("t1 t2", real=True)
    tv = sp.Matrix(ts)
    mus, Sigs = [], []
    for j in (1, 2):
        mus.append(sp.Matrix(sp.symbols(f"mu{j}_1 mu{j}_2", real=True)))
        a11, a22 = sp.symbols(f"S{j}_11 S{j}_22", positive=True)
        a12 = sp.Symbol(f"S{j}_12", real=True)
        Sigs.append(sp.Matrix([[a11, a12], [a12, a22]]))
    dd = sp.symbols("d1 d2", positive=True)
    Dm = sp.diag(*dd)
    w = [p, 1 - p]

    def K(mu, S, tt):  # Gaussian cumulant generating function
        return (tt.T * mu)[0] + (tt.T * S * tt)[0] / 2

    K_y = sum(K(mus[j], Sigs[j], w[j] * tv) for j in range(2)) + K(sp.zeros(2, 1), Dm, tv)
    at0 = {t: 0 for t in ts}
    mean_y = sp.Matrix([sp.diff(K_y, t).subs(at0) for t in ts])
    cov_y = sp.hessian(K_y, ts).subs(at0)
    third = [sp.diff(K_y, ts[i], ts[j], ts[k]) for i in range(2) for j in range(2) for k in range(2)]
    gaussian = all(sp.simplify(c) == 0 for c in third)  # y is Gaussian: cumulants > 2 vanish
    yv = sp.Matrix(sp.symbols("y1 y2", real=True))
    mu_mat = sp.Matrix.hstack(mus[0], mus[1])

    def loglik(V):
        r = yv - mean_y
        return -sp.log(2 * sp.pi) - sp.log(V.det()) / 2 - (r.T * V.inv() * r)[0] / 2

    ell_lmm0 = loglik(cov_y.subs({di: 0 for di in dd}))
    Sig_p = p**2 * Sigs[0] + (1 - p) ** 2 * Sigs[1]  # main_layout.tex, Sigma(p)
    r_dec = yv - mu_mat * sp.Matrix(w)
    ell_dec = -sp.log(2 * sp.pi) - sp.log(Sig_p.det()) / 2 - (r_dec.T * Sig_p.inv() * r_dec)[0] / 2
    ok_g2 = sp.simplify(sp.together(sp.expand_log(ell_lmm0 - ell_dec, force=True))) == 0
    ok_mean = sp.simplify(mean_y - mu_mat * sp.Matrix(w)) == sp.zeros(2, 1)
    ok_cov = sp.simplify(cov_y - Sig_p - Dm) == sp.zeros(2, 2)
    ok_a = ok_lemma and ok_g1_lmm and ok_g1_dec and gaussian and ok_g2 and ok_mean and ok_cov
    out.append(
        _rec(
            "V7a",
            r"The mixed model's marginal log-likelihood in $\bm p$ equals DeCovarT's Gaussian "
            r"log-likelihood, $\bm y\sim\mathcal N(\bm\mu\bm p,\sum_jp_j^2\bm\Sigma_j)$, when $\bm D=\bm 0$",
            "exact",
            "pass" if ok_a else "fail",
            r"$G=1$, $J=2$: integrating $x_1,x_2$ out of the joint density gives "
            r"$\mathcal N(y;\,m(p),\,p^2s_1+(1-p)^2s_2+d)$, whose limit at $d=0$ is DeCovarT's density; "
            r"$G=J=2$, general $2\times2$ $\bm\Sigma_j$: the cumulant generating function gives "
            r"$\bm V=\sum_jp_j^2\bm\Sigma_j+\bm D$ and the two log-likelihoods agree at $\bm D=\bm 0$",
            400,
            detail_full=(
                f"Gaussian lemma from sympy.integrate: {ok_lemma}; G = 1 marginal after two "
                f"integrations = N(y; p m1 + (1-p) m2, {sp.sstr(v_lmm)}): {ok_g1_lmm}; limit d -> 0 "
                f"equals DeCovarT l(p): {ok_g1_dec}. G = 2: mean {sp.sstr(mean_y.T)}, covariance "
                f"{sp.sstr(cov_y)}; third cumulants zero: {gaussian}; mean = mu p: {ok_mean}; "
                f"cov = sum p_j^2 Sigma_j + D: {ok_cov}; l_LMM(D = 0) - l_DeCovarT = 0: {ok_g2}"
            ),
        )
    )

    # ---------------- V7b: a nugget D != 0 moves the maximiser ---------------
    num = {m1: 0, m2: 1, s1: 1, s2: 4, y: sp.Rational(1, 2)}
    ell_d0 = _lognorm(y, p * m1 + (1 - p) * m2, p**2 * s1 + (1 - p) ** 2 * s2).subs(num)
    ell_d1 = _lognorm(y, p * m1 + (1 - p) * m2, p**2 * s1 + (1 - p) ** 2 * s2 + 1).subs(num)
    p0, roots0 = _stationary_max(ell_d0, p)
    p1, roots1 = _stationary_max(ell_d1, p)
    differ = sp.N(p0 - p1, 30) != 0
    out.append(
        _rec(
            "V7b",
            r"With a nugget $\bm D\neq\bm 0$ the maximum-likelihood estimate of $\bm p$ differs from "
            r"DeCovarT's in general",
            "exact",
            "flag" if differ else "fail",
            r"one gene, $\mu_1=0$, $\mu_2=1$, $\sigma_1^2=1$, $\sigma_2^2=4$, $y=1/2$: $\hat p="
            + f"{sp.N(p0, 6)}"
            + r"$ at $d=0$ (DeCovarT) and $\hat p="
            + f"{sp.N(p1, 6)}"
            + r"$ at $d=1$ (exact real roots of the score). Deviation: the note's claim holds only at "
            r"$\bm D=\bm 0$, which is DeCovarT's assumption (A5)",
            400,
            detail_full=(
                f"stationary points in (0,1): d = 0 -> {[sp.N(r, 12) for r in roots0]}, "
                f"d = 1 -> {[sp.N(r, 12) for r in roots1]}; maximisers {sp.N(p0, 15)} and "
                f"{sp.N(p1, 15)}. The LMM with D != 0 is N(mu p, sum p_j^2 Sigma_j + D), not "
                "DeCovarT's eq:convolution; its ML fit coincides with DeCovarT's only at D = 0"
            ),
        )
    )

    # ---------------- V7c: REML error contrasts still depend on p ------------
    G3 = 3
    dl = sp.Matrix(sp.symbols("delta1:4", real=True))  # Delta = mu_1 - mu_2, fixed-effect design
    sA = sp.diag(*sp.symbols("a1:4", positive=True))  # Sigma_1 (diagonal suffices)
    sB = sp.diag(*sp.symbols("b1:4", positive=True))  # Sigma_2
    dg = sp.diag(*sp.symbols("g1:4", positive=True))  # D
    Vp = p**2 * sA + (1 - p) ** 2 * sB + dg
    Kc = sp.Matrix.hstack(*dl.T.nullspace())  # error contrasts: K' Delta = 0
    assert Kc.shape == (G3, G3 - 1) and sp.simplify(Kc.T * dl) == sp.zeros(G3 - 1, 1)
    dV = (Kc.T * Vp * Kc).diff(p).applyfunc(sp.simplify)
    dV_e1 = dV.subs({dl[1]: 0, dl[2]: 0}).applyfunc(sp.factor)  # Delta along e_1
    depends = any(sp.simplify(e) != 0 for e in dV)
    out.append(
        _rec(
            "V7c",
            r"Restricted maximum likelihood does not separate $\bm p$ from the variance: the error "
            r"contrasts' covariance $\bm K^{\mathsf T}\bm V(\bm p)\bm K$ depends on $\bm p$",
            "exact",
            "flag" if depends else "fail",
            r"$G=3$, $J=2$, $\bm K^{\mathsf T}\bm\Delta=\bm 0$, diagonal $\bm\Sigma_1=\operatorname{diag}(a_g)$, "
            r"$\bm\Sigma_2=\operatorname{diag}(b_g)$; for $\bm\Delta\propto\bm e_1$, $\partial_p"
            r"(\bm K^{\mathsf T}\bm V\bm K)=" + sp.latex(dV_e1) + r"$, not identically $\bm 0$. "
            r"Deviation: the note's `REML is then an option' assumed a variance free of the fixed effect",
            401,
            detail_full=(
                f"K = {sp.sstr(Kc)}; d/dp (K' V(p) K) = {sp.sstr(dV)}; standard REML needs the law of "
                "K'y to be free of the fixed effects, which fails here because the random-effect design "
                "Z = (p_1 I, ..., p_J I) depends on p"
            ),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

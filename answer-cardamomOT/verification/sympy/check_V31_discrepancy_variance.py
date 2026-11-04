# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V31: moments of the discrepancy D and of the module sum S_M (vignette
03-1, #sec-fs-five-bounds, @eq-fs-var-d and @eq-fs-module-var).

Derivation route:
- V31a: D is written as an explicit linear form a'X over individual cells
  (toy: n = 5 bulk cells with m = (2, 3), reference C = (2, 3)), with the
  per-cell mean and variance taken from the NB MGF (M'(0), M''(0) - M'(0)^2);
  E D and Var D = sum_c a_c^2 sigma_c^2 (independent cells) are computed from
  that form and compared with the printed sums, which use sigma_gj^2 =
  mu + mu^2/psi (the vignette does not define sigma_gj^2 in this section; this
  is the per-cell NB variance of the spec). The same with symbolic J = 3 and
  m_j = p_j n via group sums. Numeric: sd D / mu_bar at the script parameters
  (low-expression gene) against the printed 0.0600.
- V31b: two genes per cell with within-cell covariance c_j (type j),
  independent cells: Cov(a'X_g, a'X_h) = sum_c a_c^2 c_{j(c)} for the same
  coefficients a, compared with the printed sums (toy and symbolic J = 3).
  Var S_M = 1' Corr 1 for a 5-gene graph with symbolic edge correlations and
  zero correlation off the edges equals |M| + 2 sum_edges Corr. Cauchy-Schwarz:
  Lagrange's identity (sum w a^2)(sum w b^2) - (sum w a b)^2 =
  sum_{j<l} w_j w_l (a_j b_l - a_l b_j)^2 (J = 3, exact), so
  Corr(D_g, D_h) = rho C / sqrt(A B) <= rho, with equality iff the
  vectors sigma_g. and sigma_h. are proportional.
"""

from __future__ import annotations

import json

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


def zero(e) -> bool:
    if sp.cancel(sp.together(sp.expand(e))) == 0:
        return True
    return sp.simplify(e) == 0


psi = sp.Symbol("psi", positive=True)
v = sp.Symbol("v", real=True)


def nb_moments(mu):
    """Mean and variance of one NB cell from its MGF (pi = mu/(mu + psi))."""
    pi_ = mu / (mu + psi)
    M = ((1 - pi_) / (1 - pi_ * sp.exp(v))) ** psi
    m1 = sp.simplify(sp.diff(M, v).subs(v, 0))
    m2 = sp.simplify(sp.diff(M, v, 2).subs(v, 0))
    return m1, sp.simplify(m2 - m1**2)


def cells(m_bulk, C_ref):
    """Coefficients a_c and types of the cells of D = (1/n) sum_bulk - sum_j (p_j/C_j) sum_ref_j."""
    n = sum(m_bulk)
    p = [sp.Rational(mj, n) for mj in m_bulk]
    coef, typ = [], []
    for j, mj in enumerate(m_bulk):
        coef += [sp.Rational(1, n)] * mj
        typ += [j] * mj
    for j, cj in enumerate(C_ref):
        coef += [-p[j] / cj] * cj
        typ += [j] * cj
    return coef, typ, p, n


def run() -> list[dict]:
    out = []
    J = 3

    # ---- V31a ---------------------------------------------------------------
    mu2 = sp.symbols("mu1:3", positive=True)
    mom = [nb_moments(m) for m in mu2]
    ok_mom = all(zero(mom[j][0] - mu2[j]) for j in range(2))
    sig2_def = [mom[j][1] for j in range(2)]
    coef, typ, p_t, n_t = cells((2, 3), (2, 3))
    C_t = (2, 3)
    ED = sum(a * mom[j][0] for a, j in zip(coef, typ))
    VD = sum(a**2 * sig2_def[j] for a, j in zip(coef, typ))
    sig2 = [mu2[j] + mu2[j] ** 2 / psi for j in range(2)]
    printed_toy = sum(p_t[j] * sig2[j] / n_t for j in range(2)) + sum(p_t[j] ** 2 * sig2[j] / C_t[j] for j in range(2))
    ok_toy = zero(ED) and zero(VD - printed_toy)
    # symbolic J = 3 with group sums (m_j = p_j n cells of coefficient 1/n, C_j cells of -p_j/C_j)
    mu = sp.symbols("mu1:4", positive=True)
    p = sp.symbols("p1:4", positive=True)
    Cs = sp.symbols("C1:4", positive=True)
    n = sp.Symbol("n", positive=True)
    mom3 = [nb_moments(m) for m in mu]
    ED3 = sum(p[j] * n * mom3[j][0] / n for j in range(J)) - sum(Cs[j] * p[j] / Cs[j] * mom3[j][0] for j in range(J))
    VD3 = sum(p[j] * n * mom3[j][1] / n**2 for j in range(J)) + sum(
        Cs[j] * (p[j] / Cs[j]) ** 2 * mom3[j][1] for j in range(J)
    )
    sig3 = [mu[j] + mu[j] ** 2 / psi for j in range(J)]
    printed = sum(p[j] * sig3[j] / n for j in range(J)) + sum(p[j] ** 2 * sig3[j] / Cs[j] for j in range(J))
    ok_sym = zero(ED3) and zero(VD3 - printed)
    # numeric: sd D / mu_bar for the low-expression gene of fig_03_1_tail_bounds.R
    sub = {psi: 2, n: 10_000}
    sub.update({p[j]: sp.nsimplify(x) for j, x in enumerate((0.6, 0.3, 0.1))})
    sub.update({Cs[j]: x for j, x in enumerate((2000, 300, 30))})
    sub.update({mu[j]: sp.nsimplify(x) for j, x in enumerate((0.2, 0.6, 0.1))})
    mubar = sum(p[j] * mu[j] for j in range(J)).subs(sub)
    ratio = float(sp.sqrt(printed.subs(sub)) / mubar)
    ok_num = round(ratio, 4) == 0.0600
    ok_a = ok_mom and ok_toy and ok_sym and ok_num
    out.append(
        _rec(
            "V31a",
            r"$\operatorname{E}D=0$ and $\operatorname{Var}D=\sum_jp_j\sigma_{gj}^2/n+\sum_jp_j^2\sigma_{gj}^2/C_j$ "
            r"with the per-cell NB variance $\sigma_{gj}^2=\tilde\mu_{gj}+\tilde\mu_{gj}^2/\psi_g$",
            "exact",
            "pass" if ok_a else "fail",
            rf"Per-cell moments from the NB MGF; $D$ as a linear form over cells (toy $m=(2,3)$, $C=(2,3)$) and with "
            rf"symbolic $J=3$, $m_j=p_jn$: both match the printed sums; $\operatorname{{sd}}D/\bar\mu_g={ratio:.5f}$ "
            rf"for the low-expression gene (printed 0.0600)",
            "",
            detail_full=f"moments {ok_mom}; toy {ok_toy}; symbolic {ok_sym}; sd ratio {ratio} ({ok_num})",
        )
    )

    # ---- V31b ---------------------------------------------------------------
    cj = sp.symbols("c1:4", real=True)  # within-cell covariance of genes g, h in type j
    cov_toy = sum(a**2 * cj[j] for a, j in zip(coef, typ))  # same coefficients for both genes
    printed_cov_toy = sum(p_t[j] * cj[j] / n_t for j in range(2)) + sum(p_t[j] ** 2 * cj[j] / C_t[j] for j in range(2))
    cov3 = sum(p[j] * n * cj[j] / n**2 for j in range(J)) + sum(Cs[j] * (p[j] / Cs[j]) ** 2 * cj[j] for j in range(J))
    printed_cov = sum(p[j] * cj[j] / n for j in range(J)) + sum(p[j] ** 2 * cj[j] / Cs[j] for j in range(J))
    ok_cov = zero(cov_toy - printed_cov_toy) and zero(cov3 - printed_cov)
    # Var S_M on a 5-gene graph
    Mn = 5
    E = [(0, 1), (1, 2), (2, 3), (1, 4)]
    r = {e: sp.Symbol(f"r_{e[0]}{e[1]}", real=True) for e in E}
    R = sp.eye(Mn)
    for (a_, b_), rv in r.items():
        R[a_, b_] = R[b_, a_] = rv
    ones = sp.ones(Mn, 1)
    varS = (ones.T * R * ones)[0, 0]
    ok_varS = zero(varS - (Mn + 2 * sum(r.values())))
    # Cauchy-Schwarz via Lagrange's identity, weights w_j = p_j/n + p_j^2/C_j > 0
    w = [p[j] / n + p[j] ** 2 / Cs[j] for j in range(J)]
    sg = sp.symbols("s_g1:4", positive=True)
    sh = sp.symbols("s_h1:4", positive=True)
    A = sum(w[j] * sg[j] ** 2 for j in range(J))
    B = sum(w[j] * sh[j] ** 2 for j in range(J))
    Cx = sum(w[j] * sg[j] * sh[j] for j in range(J))
    lag = sum(w[j] * w[l] * (sg[j] * sh[l] - sg[l] * sh[j]) ** 2 for j in range(J) for l in range(j + 1, J))
    ok_lag = zero(A * B - Cx**2 - lag)
    lam = sp.Symbol("lambda", positive=True)
    ok_eq = zero((A * B - Cx**2).subs({sh[j]: lam * sg[j] for j in range(J)}))
    # Corr(D_g, D_h) = rho * Cx / sqrt(A B), using c_ghj = rho s_gj s_hj in printed_cov and Var D of V31a
    rho = sp.Symbol("rho", positive=True)
    var_g = sum(p[j] * sg[j] ** 2 / n for j in range(J)) + sum(p[j] ** 2 * sg[j] ** 2 / Cs[j] for j in range(J))
    var_h = sum(p[j] * sh[j] ** 2 / n for j in range(J)) + sum(p[j] ** 2 * sh[j] ** 2 / Cs[j] for j in range(J))
    corr = printed_cov.subs({cj[j]: rho * sg[j] * sh[j] for j in range(J)}) / sp.sqrt(var_g * var_h)
    ok_corr = zero(corr - rho * Cx / sp.sqrt(A * B))
    ok_b = ok_cov and ok_varS and ok_lag and ok_eq and ok_corr
    out.append(
        _rec(
            "V31b",
            r"$\operatorname{Cov}(D_g,D_h)=\sum_jp_jc_{ghj}/n+\sum_jp_j^2c_{ghj}/C_j$, "
            r"$\operatorname{Var}S_M=\lvert M\rvert+2\sum_{(g,h)\in\mathcal E}\operatorname{Corr}(D_g,D_h)$, and "
            r"$\operatorname{Corr}(D_g,D_h)\le\varrho$",
            "exact",
            "pass" if ok_b else "fail",
            r"Covariance of two linear forms with equal cell coefficients (toy and symbolic $J=3$); "
            r"$\bm 1^\top\bm R\bm 1$ on a 5-gene graph with symbolic edge correlations; Lagrange's identity "
            r"$(\sum w a^2)(\sum w b^2)-(\sum wab)^2=\sum_{j<l}w_jw_l(a_jb_l-a_lb_j)^2$ with "
            r"$w_j=p_j/n+p_j^2/C_j$, zero iff $\sigma_{h\cdot}\propto\sigma_{g\cdot}$",
            "",
            detail_full=f"cov {ok_cov}; varS {ok_varS}; lagrange {ok_lag}; equality {ok_eq}; corr form {ok_corr}",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

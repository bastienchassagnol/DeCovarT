# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V14: the denoised within-type covariance of vignette 05-1
(eq-denoised-within), under the measurement model of eq-count-covariance.

Model (inside one type): x_g | (s, lambda) ~ Poisson(b_g s lambda_g),
independent over genes; s = s_bio s_tech with s_bio, s_tech and lambda
mutually independent, E s_bio = E s_tech = 1, nu_bio = Var s_bio,
nu_tech = Var s_tech; E lambda = mu, Cov lambda = Sigma^lambda; B = diag(b).
Target: the s_tech-free, Poisson-free covariance
T = Cov(B s_bio lambda) = B[(1 + nu_bio) Sigma^lambda + nu_bio mu mu']B.
Printed estimator, with kappa = nu_tech / (1 + nu_tech) taken as known:
    E1 = V - diag(xbar) - kappa [V - diag(xbar) + xbar xbar'],
V the unbiased sample covariance and xbar the mean of C i.i.d. cells.

Derivation route (G = 2):
- V14a: E s^2 = E s_bio^2 E s_tech^2 by independence, so
  1 + nu = (1 + nu_bio)(1 + nu_tech), with nu = Var s and E s = 1.
- V14b: Cov x by the tower rule from Poisson moments summed from the pmf
  (as V9), equal to eq-count-covariance; the population version of E1
  (V -> Cov x, xbar -> E x) equals T exactly.
- V14c: finite C. E xbar xbar' = E x E x' + Cov x / C and E V = Cov x,
  both derived by expanding the double sum over cells with symbolic C;
  then E E1 - T = -kappa Cov x / C (a deviation from 'unbiased'), and the
  corrected E2, with xbar xbar' - V/C in place of xbar xbar', has E E2 = T.
- V14d: plugging the total extrinsic noise nu (the CardamomOT /
  Fang-Pachter size factor, which carries s_bio and s_tech) in place of
  nu_tech gives B Sigma^lambda B: the biological size-factor share
  nu_bio (Sigma^lambda + mu mu') is removed with the technical one.
- V14e: Monte Carlo (numpy, fixed seed): gamma size factors, log-normal
  lambda with correlated genes, C = 10 cells per replicate; the mean of E1
  matches T - kappa Cov x / C and the mean of E2 matches T.
"""

from __future__ import annotations

import json

import numpy as np
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


def _poisson_moments():
    k = sp.Symbol("k", integer=True, nonnegative=True)
    m = sp.Symbol("m", positive=True)
    pmf = sp.exp(-m) * m**k / sp.factorial(k)
    e1 = sp.simplify(sp.summation(k * pmf, (k, 0, sp.oo)))
    f2 = sp.simplify(sp.summation(k * (k - 1) * pmf, (k, 0, sp.oo)))
    return m, e1, f2


def run() -> list[dict]:
    out: list[dict] = []
    nb, nt = sp.symbols("nu_b nu_t", nonnegative=True)
    sb, st = sp.symbols("s_b s_t", positive=True)
    # ---------------- V14a: 1 + nu = (1 + nu_b)(1 + nu_t) --------------------
    mom_b = {0: 1, 1: 1, 2: 1 + nb}
    mom_t = {0: 1, 1: 1, 2: 1 + nt}

    def E_s(expr):
        poly = sp.Poly(sp.expand(expr), sb, st)
        return sum(c * mom_b[i] * mom_t[j] for (i, j), c in poly.terms())

    Es1 = E_s(sb * st)
    Es2 = E_s((sb * st) ** 2)
    nu = sp.expand(Es2 - Es1**2)
    ok_a = Es1 == 1 and sp.expand(1 + nu - (1 + nb) * (1 + nt)) == 0
    out.append(
        _rec(
            "V14a",
            r"$1+\nu=(1+\nu^{\mathrm{bio}})(1+\nu^{\mathrm{tech}})$ for $s=s^{\mathrm{bio}}s^{\mathrm{tech}}$, "
            r"independent factors of mean $1$",
            "exact",
            "pass" if ok_a else "fail",
            r"$\mathbb E s=1$ and $\operatorname{Var}s=" + sp.latex(sp.factor(nu), symbol_names={nb: r"\nu^{\mathrm{bio}}", nt: r"\nu^{\mathrm{tech}}"})
            + r"$ from the product of moments",
            "",
            detail_full=f"E s = {Es1}; Var s = {nu}",
        )
    )

    # ---------------- V14b: Cov x and the population identity --------------
    m, e1, f2 = _poisson_moments()
    s = sp.Symbol("s", positive=True)
    l1, l2 = sp.symbols("l1 l2", positive=True)
    lam = [l1, l2]
    b = sp.symbols("b1 b2", positive=True)
    mu = sp.symbols("mu1 mu2", positive=True)
    S11, S22 = sp.symbols("S11 S22", nonnegative=True)
    S12 = sp.Symbol("S12", real=True)
    Sig = sp.Matrix([[S11, S12], [S12, S22]])
    B = sp.diag(*b)
    mu_v = sp.Matrix(mu)

    def lam_mom(b1, b2):
        table = {
            (0, 0): 1, (1, 0): mu[0], (0, 1): mu[1],
            (2, 0): S11 + mu[0] ** 2, (0, 2): S22 + mu[1] ** 2, (1, 1): S12 + mu[0] * mu[1],
        }
        return table[(b1, b2)]

    def E_full(expr, smom):
        """expectation over (s, lambda), s independent of lambda, degree <= 2."""
        poly = sp.Poly(sp.expand(expr), s, l1, l2)
        return sp.expand(sum(c * smom[a_] * lam_mom(b1_, b2_) for (a_, b1_, b2_), c in poly.terms()))

    def cond_E(g):
        return e1.subs(m, b[g] * s * lam[g])

    def cond_E2(g, h):
        if g == h:
            return f2.subs(m, b[g] * s * lam[g]) + e1.subs(m, b[g] * s * lam[g])
        return cond_E(g) * cond_E(h)

    s_mom_total = {0: 1, 1: 1, 2: 1 + nu}  # s = s_b s_t
    Ex = sp.Matrix([E_full(cond_E(g), s_mom_total) for g in range(2)])
    Exx = sp.Matrix(2, 2, lambda g, h: E_full(cond_E2(g, h), s_mom_total))
    Cov = (Exx - Ex * Ex.T).applyfunc(sp.expand)
    printed_cov = sp.diag(*(B * mu_v)) + B * ((1 + nu) * Sig + nu * mu_v * mu_v.T) * B
    ok_cov = all(sp.expand(e) == 0 for e in Cov - printed_cov) and all(sp.expand(e) == 0 for e in Ex - B * mu_v)
    # target: Cov(B s_b lambda), s_b of mean 1 and variance nu_b
    s_mom_bio = {0: 1, 1: 1, 2: 1 + nb}
    m_b = sp.Matrix([b[g] * s * lam[g] for g in range(2)])
    E_mb = sp.Matrix([E_full(m_b[g], s_mom_bio) for g in range(2)])
    T = (sp.Matrix(2, 2, lambda g, h: E_full(m_b[g] * m_b[h], s_mom_bio)) - E_mb * E_mb.T).applyfunc(sp.expand)
    ok_T = all(sp.expand(e) == 0 for e in T - B * ((1 + nb) * Sig + nb * mu_v * mu_v.T) * B)
    kappa = nt / (1 + nt)
    A_pop = Cov - sp.diag(*Ex)
    est_pop = A_pop - kappa * (A_pop + Ex * Ex.T)
    ok_pop = all(sp.simplify(e) == 0 for e in est_pop - T)
    ok_b = ok_cov and ok_T and ok_pop
    out.append(
        _rec(
            "V14b",
            r"eq-count-covariance, and the population version of eq-denoised-within equals the target "
            r"$\operatorname{Cov}(\bm B s^{\mathrm{bio}}\bm\lambda)=\bm B[(1+\nu^{\mathrm{bio}})\bm\Sigma^\lambda"
            r"+\nu^{\mathrm{bio}}\bm\mu\bm\mu^{\mathsf T}]\bm B$",
            "exact",
            "pass" if ok_b else "fail",
            r"Poisson moments from the pmf; tower rule over $(s,\bm\lambda)$, $G=2$, gene-wise capture "
            r"$\bm B$; $\bm V\to\operatorname{Cov}\bm x$, $\bar{\bm x}\to\mathbb E\bm x$ gives the target exactly",
            "",
            detail_full=(
                f"Cov x = eq-count-covariance and E x = B mu: {ok_cov}; target T = B[(1+nu_b)Sigma + nu_b mu mu']B: "
                f"{ok_T}; population estimator = T: {ok_pop}"
            ),
        )
    )

    # ---------------- V14c: finite C, bias and correction -------------------
    C = sp.Symbol("C", integer=True, positive=True)
    # E sum_c x_c x_c' = C Exx; E sum_{c != c'} x_c x_c'' = C(C-1) Ex Ex'
    E_xbar_xbar = (C * Exx + C * (C - 1) * Ex * Ex.T) / C**2
    E_V = (C * Exx - C * E_xbar_xbar) / (C - 1)
    ok_moms = all(sp.simplify(e) == 0 for e in E_V - Cov) and all(
        sp.simplify(e) == 0 for e in E_xbar_xbar - (Ex * Ex.T + Cov / C)
    )
    A_fin = E_V - sp.diag(*Ex)  # E diag(xbar) = diag(E x)
    E_est1 = A_fin - kappa * (A_fin + E_xbar_xbar)
    bias1 = (E_est1 - T).applyfunc(sp.simplify)
    ok_bias_form = all(sp.simplify(e) == 0 for e in bias1 + kappa * Cov / C)
    nonzero = sp.simplify(bias1[0, 1].subs({nt: 1, C: 10, b[0]: 1, b[1]: 1, mu[0]: 1, mu[1]: 1, S12: 0, nb: 0})) != 0
    E_est2 = A_fin - kappa * (A_fin + E_xbar_xbar - E_V / C)
    ok_corr = all(sp.simplify(e) == 0 for e in E_est2 - T)
    ok_c = ok_moms and ok_bias_form and nonzero and ok_corr
    out.append(
        _rec(
            "V14c",
            r"eq-denoised-within as first printed is biased: $\mathbb E\hat{\bm\Sigma}-\text{target}"
            r"=-\kappa\operatorname{Cov}\bm x/C_j$, $\kappa=\nu^{\mathrm{tech}}/(1+\nu^{\mathrm{tech}})$; "
            r"replacing $\bar{\bm x}\bar{\bm x}^{\mathsf T}$ by $\bar{\bm x}\bar{\bm x}^{\mathsf T}-\bm V/C_j$ "
            r"makes it unbiased",
            "exact",
            "flag" if ok_c else "fail",
            r"$\mathbb E\bar{\bm x}\bar{\bm x}^{\mathsf T}=\mathbb E\bm x\,\mathbb E\bm x^{\mathsf T}+\operatorname{Cov}\bm x/C$ "
            r"and $\mathbb E\bm V=\operatorname{Cov}\bm x$ from the double sum over cells, symbolic $C$. "
            r"Deviation: the printed form is consistent, biased at order $1/C_j$; corrected in 05-1",
            "",
            detail_full=(
                f"E xbar xbar' and E V: {ok_moms}; bias of printed estimator = -kappa Cov x / C: {ok_bias_form}; "
                f"bias non-zero at a test point: {nonzero}; corrected estimator unbiased: {ok_corr}; "
                f"bias(1,2) = {sp.factor(bias1[0, 1])}"
            ),
        )
    )

    # ---------------- V14d: plugging the total nu --------------------------
    est_total = A_pop - (nu / (1 + nu)) * (A_pop + Ex * Ex.T)
    diff_d = (est_total - B * Sig * B).applyfunc(sp.simplify)
    lost = (T - est_total).applyfunc(sp.simplify)
    ok_d = all(e == 0 for e in diff_d) and all(
        sp.simplify(e) == 0 for e in lost - nb * B * (Sig + mu_v * mu_v.T) * B
    )
    out.append(
        _rec(
            "V14d",
            r"Plugging the total extrinsic noise $\nu$ (a size factor carrying $s^{\mathrm{bio}}$ and "
            r"$s^{\mathrm{tech}}$) for $\nu^{\mathrm{tech}}$ returns $\bm B\bm\Sigma^\lambda\bm B$ and drops "
            r"$\nu^{\mathrm{bio}}\bm B(\bm\Sigma^\lambda+\bm\mu\bm\mu^{\mathsf T})\bm B$",
            "exact",
            "pass" if ok_d else "fail",
            r"population estimator with $\kappa=\nu/(1+\nu)$, $\nu$ from V14a; difference from the target "
            r"factorised",
            "",
            detail_full=f"estimate with total nu = B Sigma B: {all(e == 0 for e in diff_d)}; lost share: {sp.sstr(lost)}",
        )
    )

    # ---------------- V14e: Monte Carlo ------------------------------------
    rng = np.random.default_rng(20261031)
    reps, Cn = 200_000, 10
    nb_v, nt_v = 0.3, 0.5
    b_v = np.array([0.8, 0.5])
    mz = np.array([np.log(20.0), np.log(5.0)])
    Sz = np.array([[0.25, 0.10], [0.10, 0.36]])
    # log-normal lambda: exact moments
    mu_l = np.exp(mz + np.diag(Sz) / 2)
    Sig_l = np.outer(mu_l, mu_l) * (np.exp(Sz) - 1)
    Bm = np.diag(b_v)
    T_v = Bm @ ((1 + nb_v) * Sig_l + nb_v * np.outer(mu_l, mu_l)) @ Bm
    nu_v = (1 + nb_v) * (1 + nt_v) - 1
    Cov_v = np.diag(b_v * mu_l) + Bm @ ((1 + nu_v) * Sig_l + nu_v * np.outer(mu_l, mu_l)) @ Bm
    kap_v = nt_v / (1 + nt_v)
    z = rng.multivariate_normal(mz, Sz, size=(reps, Cn))
    lam_s = np.exp(z)
    s_b = rng.gamma(1 / nb_v, nb_v, size=(reps, Cn, 1))
    s_t = rng.gamma(1 / nt_v, nt_v, size=(reps, Cn, 1))
    x = rng.poisson(b_v * s_b * s_t * lam_s)
    xbar = x.mean(axis=1)
    xc = x - xbar[:, None, :]
    V = np.einsum("rcg,rch->rgh", xc, xc) / (Cn - 1)
    xx = np.einsum("rg,rh->rgh", xbar, xbar)
    Dg = np.einsum("rg,gh->rgh", xbar, np.eye(2))
    A = V - Dg
    E1 = A - kap_v * (A + xx)
    E2 = A - kap_v * (A + xx - V / Cn)
    m1, se1 = E1.mean(0), E1.std(0) / np.sqrt(reps)
    m2, se2 = E2.mean(0), E2.std(0) / np.sqrt(reps)
    z1 = np.abs(m1 - (T_v - kap_v * Cov_v / Cn)) / se1
    z2 = np.abs(m2 - T_v) / se2
    zb = np.abs(m1 - T_v) / se1
    ok_e = bool(z1.max() < 4.5 and z2.max() < 4.5 and zb.max() > 4.5)
    out.append(
        _rec(
            "V14e",
            r"Monte Carlo of V14c: gamma size factors ($\nu^{\mathrm{bio}}=0.3$, $\nu^{\mathrm{tech}}=0.5$), "
            r"log-normal $\bm\lambda$ with correlated genes, $C_j=10$ cells",
            "numeric",
            "pass" if ok_e else "fail",
            r"numpy seed 20261031, $2\times10^5$ replicates; largest $|z|$: printed form against "
            r"target minus bias " + f"{z1.max():.2f}" + r", corrected form against target " + f"{z2.max():.2f}"
            + r" (tolerance 4.5); printed form against target " + f"{zb.max():.1f}",
            "",
            detail_full=(
                f"T = {T_v.round(3).tolist()}; mean printed = {m1.round(3).tolist()}; mean corrected = "
                f"{m2.round(3).tolist()}; |z| printed vs T - kappa Cov/C = {z1.round(2).tolist()}; |z| corrected "
                f"vs T = {z2.round(2).tolist()}; |z| printed vs T = {zb.round(1).tolist()}"
            ),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

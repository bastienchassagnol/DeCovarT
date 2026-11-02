# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V4: Fisher information for the fraction p (two types, one gene) against the
ANOVA F ratio, and the exact condition under which they rank genes alike.

Derivation route:
- the information is -E d^2/dp^2 log f(y; p), with y = m(p) + L(p) z and the
  expectation taken over standard normal z by their moments (nothing typed
  from the note); the univariate and bivariate cases use the same operator;
- the ANOVA mean squares are built from symbolic data (C = 3 cells per type),
  and their expectations from the sampling laws of the group means and of
  the pooled variance (inverse chi-square moment integrated from the density);
- agreement conditions are solved as polynomial identities in p.
"""

from __future__ import annotations

import json
from math import prod

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


def _normal_moment(k: int):
    return 0 if k % 2 else sp.factorial2(k - 1) if k > 0 else 1


def expect_normal(expr, zs):
    """E over independent standard normal zs of a polynomial in zs."""
    P = sp.Poly(sp.expand(expr), *zs)
    return sp.Add(*[c * prod(_normal_moment(k) for k in mon) for mon, c in P.terms()])


def fisher_1d(m, v, p):
    """-E d2/dp2 log N(y; m(p), v(p)) at y = m + sqrt(v) z."""
    y, z = sp.symbols("y z", real=True)
    ll = -sp.log(v) / 2 - (y - m) ** 2 / (2 * v)
    d2 = sp.diff(ll, p, 2)
    d2 = sp.together(d2.subs(y, m + sp.sqrt(v) * z))
    num, den = sp.fraction(d2)
    return sp.simplify(-expect_normal(num, [z]) / den)


def run() -> list[dict]:
    out: list[dict] = []
    p = sp.Symbol("p", positive=True)
    mu1, mu2 = sp.symbols("mu1 mu2", real=True)
    s1, s2, sig, a, n = sp.symbols("sigma1 sigma2 sigma a n", positive=True)
    Delta = mu1 - mu2
    c = a * n  # c_i = a_i n_i
    phi = a**2 * n  # phi_i = a_i^2 n_i

    # ------------------------------------------------ single-cell model
    m = c * (p * mu1 + (1 - p) * mu2)
    v = phi * (p * s1**2 + (1 - p) * s2**2)
    I_sc = fisher_1d(m, v, p)
    note_I = c**2 * Delta**2 / v + sp.Rational(1, 2) * (
        (s1**2 - s2**2) / (p * s1**2 + (1 - p) * s2**2)
    ) ** 2
    ok_gen = sp.simplify(I_sc - note_I) == 0
    I_eq = sp.simplify(I_sc.subs({s1: sig, s2: sig}))
    ok_eq = sp.simplify(I_eq - n * Delta**2 / sig**2) == 0
    out.append(
        _rec(
            "V4a",
            r"Single-cell model, $J=2$, one gene: $I(p)=c^2\Delta^2/v+\tfrac12\{(\sigma_1^2-\sigma_2^2)/v\}^2$, and $I(p)=n_i\Delta^2/\sigma^2$ when $\sigma_1=\sigma_2$",
            "exact",
            "pass" if ok_gen and ok_eq else "fail",
            r"$-\mathbb E\,\partial_p^2\ell$ with $y=m+\sqrt v\,z$ equals the printed form; with equal variances, $c=an$ and $\phi=a^2n$: $I=" + sp.latex(sp.factor(I_eq)) + "$",
            168,
            detail_full=(
                "information derived as -E d2/dp2 log-likelihood with y = m + sqrt(v) z; equals "
                f"the note's form; with equal variances and c = a n, phi = a^2 n: I = {sp.sstr(I_eq)}"
            ),
        )
    )

    # ------------------------------------------------ ANOVA ratio, J = 2
    C = 3
    y1 = sp.symbols(f"u1:{C + 1}", real=True)
    y2 = sp.symbols(f"w1:{C + 1}", real=True)
    m1, m2 = sum(y1) / C, sum(y2) / C
    grand = (sum(y1) + sum(y2)) / (2 * C)
    msb = C * ((m1 - grand) ** 2 + (m2 - grand) ** 2) / (2 - 1)
    ident = sp.expand(msb - C * (m1 - m2) ** 2 / 2) == 0
    # expectations: ybar1 - ybar2 ~ N(Delta, 2 sigma^2/C); pooled variance
    # s^2 ~ sigma^2 chi2_nu / nu, nu = 2C - 2, independent of the means
    Cs, nu = sp.symbols("C nu", positive=True)
    z = sp.Symbol("z", real=True)
    E_msb = expect_normal(Cs * (Delta + sp.sqrt(2 * sig**2 / Cs) * z) ** 2 / 2, [z])
    xx = sp.Symbol("x", positive=True)
    k = sp.Integer(2 * C - 2)
    chi2 = xx ** (k / 2 - 1) * sp.exp(-xx / 2) / (2 ** (k / 2) * sp.gamma(k / 2))
    E_inv = sp.simplify(sp.integrate(chi2 / xx, (xx, 0, sp.oo)))  # E 1/chi2_k
    E_F = sp.simplify(E_msb.subs(Cs, C) * k * E_inv / sig**2)
    lam = sp.Rational(C, 2) * Delta**2 / sig**2
    target = k / (k - 2) * (1 + lam)
    ok_EF = sp.simplify(E_F - target) == 0 and sp.simplify(E_msb - (sig**2 + Cs * Delta**2 / 2)) == 0
    out.append(
        _rec(
            "V4b",
            r"ANOVA ratio, $C$ cells per type: $F$ is the plug-in of $(C/2)\Delta^2/\sigma^2$, but $\mathbb E F=\frac{\nu}{\nu-2}\{1+(C/2)\Delta^2/\sigma^2\}$",
            "exact",
            "flag" if ident and ok_EF else "fail",
            r"$\mathrm{MS}_B=C(\bar y_1-\bar y_2)^2/2$ on symbolic data; $\mathbb E(1/\chi^2_4)=" + sp.latex(E_inv) + r"$; at $C=3$, $\mathbb E F=" + sp.latex(sp.factor(E_F)) + r"$. Deviation: $I(p)\approx(2n_i/C)(F-1)$, not $(2n_i/C)F$",
            168,
            detail_full=(
                "MS_between = C (ybar1 - ybar2)^2 / 2 exactly (symbolic data, C = 3); "
                f"E MS_between = {sp.sstr(E_msb)}; E 1/chi2_4 = {E_inv}; E F = {sp.sstr(E_F)} at C = 3. "
                "Deviation from the note: F estimates 1 + (C/2) Delta^2/sigma^2 up to nu/(nu-2), so "
                "I(p) ~ (2 n_i/C)(F - 1), not (2 n_i/C) F; rankings by F are unaffected"
            ),
        )
    )

    # ------------------------------------------------ agreement condition
    # ANOVA with equal C targets Delta^2 / sbar^2, sbar^2 = (s1^2 + s2^2)/2.
    # Agreement for every p: I(p) = n Delta^2 / sbar^2 identically in p.
    S1, S2 = sp.symbols("S1 S2", positive=True)  # variances
    I_v = I_sc.subs({s1: sp.sqrt(S1), s2: sp.sqrt(S2)})
    gap = sp.together(I_v - n * Delta**2 / ((S1 + S2) / 2))
    num = sp.Poly(sp.numer(gap), p)
    sols = sp.solve(num.coeffs(), [S2], dict=True)
    ok_agree = sols == [{S2: S1}]
    out.append(
        _rec(
            "V4c",
            r"Fisher information and ANOVA ratio agree for every $p$ ($J=2$) only with equal within-type variances",
            "exact",
            "pass" if ok_agree else "fail",
            r"$I(p)-n\Delta^2/\bar s^2$ is a polynomial in $p$ whose coefficients all vanish only at " + ", ".join("$" + sp.latex(k) + "=" + sp.latex(v) + "$" for sol in sols for k, v in sol.items()) + r" ($\Delta\neq0$)",
            223,
            detail_full=(
                "I(p) - n Delta^2/sbar^2 written as a polynomial in p; its coefficients vanish "
                f"together only for {sols} (Delta != 0); at p = 1/2 with S1 != S2 the gap is the "
                "positive term (1/2){(S1 - S2)/sbar^2}^2"
            ),
        )
    )

    # ------------------------------------------------ multivariate, common Sigma_W
    G = 2
    zs = sp.symbols("z1 z2", real=True)
    l11, l22 = sp.symbols("l11 l22", positive=True)
    l21 = sp.Symbol("l21", real=True)
    L = sp.Matrix([[l11, 0], [l21, l22]])
    SigW = L * L.T
    D1, D2 = sp.symbols("D1 D2", real=True)
    Dv = sp.Matrix([D1, D2])
    mu2v = sp.Matrix(sp.symbols("nu1 nu2", real=True))
    mv = c * (p * (mu2v + Dv) + (1 - p) * mu2v)
    V = phi * SigW  # p * SigW + (1-p) * SigW
    yv = sp.Matrix(sp.symbols("y1 y2", real=True))
    ll = -sp.log(V.det()) / 2 - ((yv - mv).T * V.inv() * (yv - mv))[0] / 2
    d2 = sp.diff(ll, p, 2).subs(dict(zip(yv, mv + sp.sqrt(phi) * L * sp.Matrix(zs))))
    I_mv = sp.simplify(-expect_normal(sp.expand(sp.cancel(d2)), list(zs)))
    ok_mv = sp.simplify(I_mv - n * (Dv.T * SigW.inv() * Dv)[0]) == 0
    out.append(
        _rec(
            "V4d",
            r"Multivariate information with a common $\bm\Sigma_W$ is $n_i\bm\Delta^{\mathsf T}\bm\Sigma_W^{-1}\bm\Delta$",
            "exact",
            "pass" if ok_mv else "fail",
            r"bivariate $-\mathbb E\,\partial_p^2\ell$, expectation over $\bm y=\bm m+\sqrt\phi\,\bm L\bm z$, equals $n\bm\Delta^{\mathsf T}\bm\Sigma_W^{-1}\bm\Delta$",
            170,
            detail_full=(
                "bivariate Gaussian log-likelihood differentiated twice in p, expectation over "
                "y = m + sqrt(phi) L z; equals n Delta' Sigma_W^-1 Delta"
            ),
        )
    )

    # ------------------------------------------------ additivity over genes
    rho = sp.Symbol("rho", real=True)
    d1, d2_ = sp.symbols("d1 d2", real=True)  # standardised differences Delta_g/sigma_g
    R = sp.Matrix([[1, rho], [rho, 1]])
    dv = sp.Matrix([d1, d2_])
    gap_add = sp.factor(sp.simplify((dv.T * R.inv() * dv)[0] - (d1**2 + d2_**2)))
    roots = sp.solve(sp.numer(gap_add), rho)
    zero_at_0 = gap_add.subs(rho, 0) == 0
    extra = [r for r in roots if r != 0]
    st_add = "fail" if not zero_at_0 else ("flag" if extra else "pass")
    out.append(
        _rec(
            "V4e",
            r"Per-gene informations add for independent genes; with correlated genes, only at isolated values of $\rho$",
            "exact",
            st_add,
            r"$n^{-1}$(multivariate $-$ sum of per-gene) $=" + sp.latex(gap_add) + r"$, zero at $\rho\in\{" + ", ".join(sp.latex(r) for r in roots) + r"\}$. Deviation: additivity also holds at the second root, a measure-zero set",
            504,
            detail_full=(
                f"n^-1 (multivariate - sum of per-gene) = {sp.sstr(gap_add)}; zero for rho in "
                f"{[sp.sstr(r) for r in roots]}. Deviation: additivity also holds at "
                "rho = 2 d1 d2/(d1^2 + d2^2) (d_g = Delta_g/sigma_g), a measure-zero set; "
                "the note now says 'in general'"
            ),
        )
    )

    # ------------------------------------------------ population model
    mp_ = p * mu1 + (1 - p) * mu2
    vp = sig**2 * (p**2 + (1 - p) ** 2)
    I_pop = fisher_1d(mp_, vp, p)
    offset = sp.simplify(I_pop.subs(mu1, mu2))  # information left when Delta = 0
    slope = sp.simplify((I_pop - offset) * sig**2 / Delta**2)
    not_prop = offset != 0 and sp.simplify(offset.subs(p, sp.Rational(1, 2))) == 0
    sigma_free = sp.diff(offset, sig) == 0
    ok_pop = not_prop and sigma_free and sp.diff(slope, p) != 0
    out.append(
        _rec(
            "V4f",
            r"Population model, equal variances: $v(p)$ depends on $p$, and $I(p)$ is not proportional to the ANOVA ratio",
            "exact",
            "pass" if ok_pop else "fail",
            r"$I(p)=\frac{\Delta^2}{\sigma^2}\cdot" + sp.latex(slope) + "+" + sp.latex(sp.factor(offset)) + r"$; the offset is free of $\Delta$ and $\sigma$, so the ranking at fixed $p$ is kept",
            505,
            detail_full=(
                f"I(p) = Delta^2/sigma^2 * {sp.sstr(slope)} + {sp.sstr(sp.factor(offset))}; the "
                "offset is non-zero for p != 1/2, so I is affine, not proportional, in "
                "Delta^2/sigma^2; the offset is free of sigma and Delta (common to all genes), so "
                "the ranking of genes at fixed p is unchanged"
            ),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

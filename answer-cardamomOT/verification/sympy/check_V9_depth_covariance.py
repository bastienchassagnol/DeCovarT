# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V9: eq:depth (noise-layer section and Appendix A.5): the within-type
covariance of UMI counts with a per-cell depth s and Poisson sampling,

    Cov(x) = E(s^2) Sigma^lam + Var(s) mu^lam mu^lam' + E(s) diag(mu^lam),

for x_g | (lam, s) independent Poisson(s lam_g), with s independent of lam,
E lam = mu^lam, Cov lam = Sigma^lam. Then Fang & Pachter's normalised
covariance (their eq. 5, p. 5).

Derivation route (G = 2, no distributional form for s or lam beyond their
first two moments):
- Poisson moments E X = m and E X(X - 1) = m^2 summed from the pmf;
- conditional moments of (x_1, x_2) given (lam, s) from those and the
  conditional independence of the genes;
- the tower rule applied to these polynomials in (s, lam_1, lam_2), with
  E s^a lam^b = E s^a E lam^b (independence) and the moments of lam
  expressed through mu^lam and Sigma^lam;
- the two terms of the law of total covariance are also computed separately
  and compared with the decomposition printed in Appendix A.5.
V9b: normalised covariance Cov(x_a, x_b)/(E x_a E x_b) for Sigma^lam_ab = 0,
and normalised variance (Var x_g - E x_g)/(E x_g)^2 for Sigma^lam_gg = 0 (a
"Poisson gene"), both against Var(s)/(E s)^2, with E s left free.
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


def _poisson_moments():
    k = sp.Symbol("k", integer=True, nonnegative=True)
    m = sp.Symbol("m", positive=True)
    pmf = sp.exp(-m) * m**k / sp.factorial(k)
    e1 = sp.simplify(sp.summation(k * pmf, (k, 0, sp.oo)))
    f2 = sp.simplify(sp.summation(k * (k - 1) * pmf, (k, 0, sp.oo)))
    return m, e1, f2


def run() -> list[dict]:
    out: list[dict] = []
    m, e1, f2 = _poisson_moments()
    ok_pois = sp.simplify(e1 - m) == 0 and sp.simplify(f2 - m**2) == 0

    s, l1, l2 = sp.symbols("s l1 l2", positive=True)
    lam = [l1, l2]
    Es1, Es2 = sp.symbols("Es1 Es2", positive=True)  # E s, E s^2
    mu = sp.symbols("mu1 mu2", positive=True)
    S11, S22 = sp.symbols("S11 S22", nonnegative=True)
    S12 = sp.Symbol("S12", real=True)
    Sig = sp.Matrix([[S11, S12], [S12, S22]])

    # conditional moments given (lam, s): x_g ~ Poisson(s lam_g), independent
    def cond_E(g):
        return e1.subs(m, s * lam[g])

    def cond_E2(g, h):
        if g == h:  # E X^2 = E X(X-1) + E X
            return f2.subs(m, s * lam[g]) + e1.subs(m, s * lam[g])
        return cond_E(g) * cond_E(h)  # conditional independence

    # expectation over (s, lam), s independent of lam, polynomial of degree <= 2
    s_mom = {0: 1, 1: Es1, 2: Es2}

    def lam_mom(b1, b2):
        if (b1, b2) == (0, 0):
            return 1
        if (b1, b2) == (1, 0):
            return mu[0]
        if (b1, b2) == (0, 1):
            return mu[1]
        if (b1, b2) == (2, 0):
            return S11 + mu[0] ** 2
        if (b1, b2) == (0, 2):
            return S22 + mu[1] ** 2
        if (b1, b2) == (1, 1):
            return S12 + mu[0] * mu[1]
        raise ValueError((b1, b2))

    def E(expr):
        poly = sp.Poly(sp.expand(expr), s, l1, l2)
        return sum(c * s_mom[a] * lam_mom(b1, b2) for (a, b1, b2), c in poly.terms())

    Ex = sp.Matrix([E(cond_E(g)) for g in range(2)])
    Exx = sp.Matrix(2, 2, lambda g, h: E(cond_E2(g, h)))
    Cov = (Exx - Ex * Ex.T).applyfunc(sp.expand)
    mu_v = sp.Matrix(mu)
    Vs = Es2 - Es1**2
    printed = Es2 * Sig + Vs * mu_v * mu_v.T + Es1 * sp.diag(*mu)
    ok_cov = all(sp.expand(e) == 0 for e in Cov - printed)
    ok_mean = all(sp.expand(e) == 0 for e in Ex - Es1 * mu_v)
    # the two terms of the law of total covariance, as written in A.5
    within = sp.Matrix(2, 2, lambda g, h: E(cond_E2(g, h) - cond_E(g) * cond_E(h)))
    m_cond = sp.Matrix([cond_E(g) for g in range(2)])  # = s lam
    between = sp.Matrix(2, 2, lambda g, h: E(m_cond[g] * m_cond[h])) - Ex * Ex.T
    ok_within = all(sp.expand(e) == 0 for e in within - Es1 * sp.diag(*mu))
    ok_between = all(sp.expand(e) == 0 for e in between - (Es2 * Sig + Vs * mu_v * mu_v.T))
    ok_a = ok_pois and ok_cov and ok_mean and ok_within and ok_between
    out.append(
        _rec(
            "V9a",
            r"eq:depth: $\operatorname{Cov}\bm x=\mathbb E(s^2)\bm\Sigma^\lambda+\operatorname{Var}(s)\,"
            r"\bm\mu^\lambda\bm\mu^{\lambda\mathsf T}+\mathbb E(s)\operatorname{diag}(\bm\mu^\lambda)$, "
            r"for Poisson counts given $(\bm\lambda,s)$ and $s\perp\bm\lambda$",
            "exact",
            "pass" if ok_a else "fail",
            r"Poisson moments from the pmf; conditional moments given $(\bm\lambda,s)$; tower rule with "
            r"only the first two moments of $s$ and $\bm\lambda$ ($G=2$); the within term is "
            r"$\mathbb E(s)\operatorname{diag}(\bm\mu^\lambda)$ and the between term "
            r"$\operatorname{Cov}(s\bm\lambda)$, as in A.5",
            286,
            detail_full=(
                f"Poisson pmf moments E X = {e1}, E X(X-1) = {f2}: {ok_pois}; E x = {sp.sstr(Ex.T)}; "
                f"Cov x = {sp.sstr(Cov)}; equal to eq:depth: {ok_cov}; within term = E(s) diag(mu): "
                f"{ok_within}; between term = E(s^2) Sigma + Var(s) mu mu': {ok_between}"
            ),
        )
    )

    # ---------------- V9b: Fang & Pachter's normalised (co)variance ---------
    target = Vs / Es1**2
    ncov = sp.simplify((Cov[0, 1] / (Ex[0] * Ex[1])).subs(S12, 0))
    nvar = sp.simplify(((Cov[0, 0] - Ex[0]) / Ex[0] ** 2).subs(S11, 0))
    ok_ncov = sp.simplify(ncov - target) == 0
    ok_nvar = sp.simplify(nvar - target) == 0
    # scale invariance: s -> c s, lam -> lam / c leaves the ratio unchanged
    c = sp.Symbol("c", positive=True)
    scaled = target.subs({Es1: c * Es1, Es2: c**2 * Es2}, simultaneous=True)
    ok_scale = sp.simplify(scaled - target) == 0
    ok_b = ok_ncov and ok_nvar and ok_scale
    names = {Es1: r"(\mathbb E s)", Es2: r"\mathbb E(s^2)"}
    out.append(
        _rec(
            "V9b",
            r"Fang & Pachter's eq. 5: for $\Sigma^\lambda_{ab}=0$, $\operatorname{Cov}(x_a,x_b)/"
            r"(\mathbb E x_a\,\mathbb E x_b)=\operatorname{Var}(s)/(\mathbb E s)^2$, which also equals the "
            r"normalised variance of a gene with $\Sigma^\lambda_{gg}=0$",
            "exact",
            "flag" if ok_b else "fail",
            r"normalised covariance $=" + sp.latex(sp.together(ncov), symbol_names=names)
            + r"$ and normalised variance $="
            + sp.latex(sp.together(nvar), symbol_names=names)
            + r"$, with $\mathbb E s$ free; both are invariant to $s\mapsto cs$. Deviation: the note's "
            r"`up to the normalisation $\mathbb E s=1$' is not needed",
            287,
            detail_full=(
                f"normalised covariance at S12 = 0: {sp.sstr(ncov)}; normalised variance at S11 = 0: "
                f"{sp.sstr(nvar)}; target Var(s)/(E s)^2 = {sp.sstr(target)}; scale invariant: {ok_scale}"
            ),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

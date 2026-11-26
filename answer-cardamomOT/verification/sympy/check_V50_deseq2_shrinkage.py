# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V50: DESeq2 box of the note on figures 1 and 4
(review-decovart-bastien-elias-marielle-anais/retour_figure_1_and_4.qmd).

- V50a (exact): ashr posterior mean for a point mass plus centred normals:
  E(beta | betahat) = sum_k w_k sigma_k^2/(sigma_k^2 + s^2) betahat, weights
  w_k proportional to pi_k N(betahat; 0, sigma_k^2 + s^2) normalised over
  k = 0..K (K = 2); the shrinkage factor tends to 0 when s -> infinity.
- V50b (exact): NB with mean mu and dispersion alpha: log P(K = 0) =
  -log(1 + alpha mu)/alpha, which tends to 0 when mu -> 0.
- V50c (exact): one-way design with k groups: full model k coefficients,
  reduced model 1, so the LRT has k - 1 degrees of freedom (k = 2..6).
- V50d (exact): the zero-inflated and hurdle (zero-truncated) laws sum to 1
  (Poisson and NB); for a group of n zero counts the zero-inflated
  likelihood is 1 at pi = 1 or mu = 0 (its maximum), and the hurdle
  likelihood pi^n is maximal at pi = 1.
"""

from __future__ import annotations

import sympy as sp


def _rec(id_, claim, kind, status, detail, tex_line=""):
    return {"id": id_, "claim": claim, "kind": kind, "status": status, "detail": detail,
            "detail_full": detail, "tex_line": tex_line}


def _norm(x, var):
    return sp.exp(-x**2 / (2 * var)) / sp.sqrt(2 * sp.pi * var)


def run() -> list[dict]:
    out = []

    # ---- V50a: ashr posterior mean ------------------------------------------
    b, bh = sp.symbols("beta betahat", real=True)
    s, s1, s2 = sp.symbols("s sigma1 sigma2", positive=True)
    p0, p1, p2 = sp.symbols("pi0 pi1 pi2", positive=True)
    lik = _norm(bh - b, s**2)
    # Each normal component: marginal N(betahat; 0, sigma^2 + s^2) and
    # posterior mean sigma^2/(sigma^2 + s^2) betahat (conjugacy).
    comp_ok = True
    parts = []
    for pk, sk in [(p1, s1), (p2, s2)]:
        i1 = sp.simplify(sp.integrate(lik * _norm(b, sk**2), (b, -sp.oo, sp.oo)))
        i2 = sp.simplify(sp.integrate(b * lik * _norm(b, sk**2), (b, -sp.oo, sp.oo)))
        comp_ok = comp_ok and sp.simplify(i1 / _norm(bh, sk**2 + s**2)) == 1
        comp_ok = comp_ok and sp.simplify(i2 / i1 - sk**2 / (sk**2 + s**2) * bh) == 0
        parts.append((pk, i1, i2))
    num = sum(pk * i2 for pk, _, i2 in parts)
    den = p0 * lik.subs(b, 0) + sum(pk * i1 for pk, i1, _ in parts)
    marg = [p0 * _norm(bh, s**2), p1 * _norm(bh, s1**2 + s**2), p2 * _norm(bh, s2**2 + s**2)]
    tot = sum(marg)
    formula = (marg[1] / tot) * s1**2 / (s1**2 + s**2) * bh + (marg[2] / tot) * s2**2 / (s2**2 + s**2) * bh
    # Mixture: check the assembled ratio at fixed rational points (fallback).
    pts = [{p0: sp.Rational(1, 2), p1: sp.Rational(1, 3), p2: sp.Rational(1, 6), s: sp.Rational(a, 4), s1: sp.Rational(1, 2), s2: sp.Integer(3), bh: sp.Rational(c, 3)}
           for a, c in [(1, 2), (3, -5), (7, 1), (2, 9)]]
    ok1 = comp_ok and all(abs(sp.N((num / den - formula).subs(pt), 30)) < 1e-25 for pt in pts)
    ok2 = sp.limit(s1**2 / (s1**2 + s**2), s, sp.oo) == 0
    out.append(_rec("V50a", r"ashr: $\mathrm E(\beta\mid\hat\beta)=\sum_k w_k\sigma_k^2/(\sigma_k^2+\hat s^2)\,\hat\beta$, factor $\to0$ as $\hat s\to\infty$",
                    "exact", "pass" if (ok1 and ok2) else "fail", "normal components exact (conjugacy); mixture of a point mass and two normals checked at 4 rational points, 30 digits"))

    # ---- V50b: NB zero probability ------------------------------------------
    mu, al = sp.symbols("mu alpha", positive=True)
    r = 1 / al  # size
    p0nb = (r / (r + mu)) ** r
    ok = sp.simplify(sp.log(p0nb) - (-sp.log(1 + al * mu) / al)) == 0 and sp.limit(-sp.log(1 + al * mu) / al, mu, 0) == 0
    out.append(_rec("V50b", r"NB: $\log\Pr(K=0)=-\alpha^{-1}\log(1+\alpha\mu)\to0$ as $\mu\to0$",
                    "exact", "pass" if ok else "fail", "size $1/\\alpha$, symbolic"))

    # ---- V50c: LRT degrees of freedom ----------------------------------------
    ok = True
    for k in range(2, 7):
        X = sp.zeros(k, k)  # one row per group: intercept + k-1 treatment dummies
        for i in range(k):
            X[i, 0] = 1
            if i > 0:
                X[i, i] = 1
        ok = ok and (X.rank() - sp.ones(k, 1).rank()) == k - 1
    out.append(_rec("V50c", r"One-way design with $k$ groups: LRT on $k-1$ degrees of freedom",
                    "exact", "pass" if ok else "fail", "rank of the design, $k=2,\\ldots,6$"))

    # ---- V50d: zero-inflated and hurdle laws ---------------------------------
    pi_, n = sp.symbols("pi n", positive=True)
    kk = sp.symbols("k", integer=True, nonnegative=True)
    pois = sp.exp(-mu) * mu**kk / sp.factorial(kk)
    zip_tot = pi_ + (1 - pi_) * sp.exp(-mu) + (1 - pi_) * sp.summation(pois, (kk, 1, sp.oo))
    hur_tot = pi_ + (1 - pi_) * sp.summation(pois, (kk, 1, sp.oo)) / (1 - sp.exp(-mu))
    ok1 = sp.simplify(zip_tot - 1) == 0 and sp.simplify(hur_tot - 1) == 0
    # NB version, numeric size and mean (exact rationals in the closed form)
    rr, mm = sp.Rational(2), sp.Rational(3)
    nb = sp.binomial(kk + rr - 1, kk) * (rr / (rr + mm)) ** rr * (mm / (rr + mm)) ** kk
    nb_pos = sp.summation(nb, (kk, 1, sp.oo))
    ok2 = sp.simplify(nb_pos - (1 - (rr / (rr + mm)) ** rr)) == 0
    zi_zero = (pi_ + (1 - pi_) * sp.exp(-mu)) ** n
    ok3 = zi_zero.subs(pi_, 1) == 1 and sp.limit(zi_zero, mu, 0) == 1
    ok4 = sp.diff(pi_**n, pi_).subs(pi_, 1) == n  # increasing up to pi = 1
    out.append(_rec("V50d", r"Zero-inflated and hurdle laws sum to 1; an all-zero group puts the zero-inflated optimum at the boundary and the hurdle optimum at $\pi=1$",
                    "exact", "pass" if (ok1 and ok2 and ok3 and ok4) else "fail",
                    "Poisson symbolic, NB with size 2 and mean 3; likelihoods of $n$ zeros"))
    return out


if __name__ == "__main__":
    for r_ in run():
        print(r_["id"], r_["status"], r_["detail"])

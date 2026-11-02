# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V20: negative-binomial counts are sub-exponential; Chernoff and Bernstein
versions of RNA-Sieve's range bound (vignette 03-1, #sec-fs-subexp).

Parametrisation: X ~ NB(r, pi), failures before the r-th success,
P(X = x) = C(x + r - 1, x) pi^x (1 - pi)^r, r > 0, 0 < pi < 1;
mu = r pi / (1 - pi), Var X = mu + mu^2 / r.

- V20a (exact + series): the generating-function identity
  sum_x C(x+r-1, x) z^x = (1 - z)^(-r) (coefficients up to z^10, symbolic r)
  gives M(lambda) = ((1 - pi) / (1 - pi e^lambda))^r for pi e^lambda < 1,
  i.e. lambda < -log pi (ratio test: the series diverges beyond); M'(0) = mu
  and M''(0) - M'(0)^2 = mu + mu^2/r.
- V20b (numeric): the moment growth of Vershynin's Prop. 2.8.1(ii):
  (E|X - mu|^p)^(1/p) / p stays bounded in p (p = 1..60) for a grid of
  (r, pi); log-space summation in float64 to x where the summand is
  negligible.
- V20c (exact + numeric): the Cramér-Chernoff rate
  I(a) = a log[a (mu + r) / (mu (a + r))] + r log[(mu + r) / (a + r)]
  equals sup_lambda {lambda a - log M(lambda)} (stationary point
  e^lambda* = a (mu + r) / (mu (a + r)) checked symbolically), with I(mu) = 0,
  I(0) = r log(1 + mu/r) = -log P(X = 0) and the Poisson limit as r -> oo.
  Monte Carlo (seed 20261101): the tails of the mean of n NB counts, and the
  range distance delta of RNA-Sieve's filter with NB counts, stay below the
  Chernoff bounds e^{-n I}; the mixture rate I_0 is maximised on a grid of
  lambda, which under-estimates the supremum and so over-estimates the bound
  (conservative).
"""

from __future__ import annotations

import json
import math

import numpy as np
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
    return sp.simplify(sp.expand(e)) == 0


def _logpmf(x: np.ndarray, r: float, pi: float) -> np.ndarray:
    lg = np.vectorize(math.lgamma)
    return lg(x + r) - math.lgamma(r) - lg(x + 1.0) + x * math.log(pi) + r * math.log1p(-pi)


def _rate(a, mu, r):
    """Closed-form NB Chernoff rate I(a) (vectorised, a >= 0)."""
    a = np.asarray(a, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        t1 = np.where(a > 0, a * np.log(a * (mu + r) / (mu * (a + r))), 0.0)
    return t1 + r * np.log((mu + r) / (a + r))


def run() -> list[dict]:
    out = []
    r, pi, lam, z = sp.symbols("r pi lambda z", positive=True)

    # ---- V20a ----------------------------------------------------------
    order = 10
    ser = sp.series((1 - z) ** (-r), z, 0, order).removeO()
    ok_series = all(
        zero(ser.coeff(z, x) - sp.rf(r, x) / sp.factorial(x)) for x in range(order)
    )
    M = ((1 - pi) / (1 - pi * sp.exp(lam))) ** r
    mu = r * pi / (1 - pi)
    m1 = sp.diff(M, lam).subs(lam, 0)
    m2 = sp.diff(M, lam, 2).subs(lam, 0)
    ok_mean = zero(sp.simplify(m1 - mu))
    ok_var = zero(sp.simplify(m2 - m1**2 - (mu + mu**2 / r)))
    # domain: ratio of consecutive terms of the MGF series -> pi e^lambda
    x = sp.Symbol("x", positive=True, integer=True)
    term = sp.rf(r, x) / sp.factorial(x) * (pi * sp.exp(lam)) ** x
    ratio = sp.limit(sp.combsimp(term.subs(x, x + 1) / term), x, sp.oo)
    ok_dom = zero(ratio - pi * sp.exp(lam))
    ok_a = ok_series and ok_mean and ok_var and ok_dom
    out.append(
        _rec(
            "V20a",
            r"For $X\sim\mathrm{NB}(r,\pi)$, $M(\lambda)=\{(1-\pi)/(1-\pi e^\lambda)\}^r$, finite exactly for "
            r"$\lambda<-\log\pi$, with mean $\mu=r\pi/(1-\pi)$ and variance $\mu+\mu^2/r$: the NB is sub-exponential",
            "exact",
            "pass" if ok_a else "fail",
            r"Coefficients of $(1-z)^{-r}$ to $z^9$ equal $(r)_x/x!$ (symbolic $r$); $M'(0)$, $M''(0)-M'(0)^2$ "
            r"simplify; ratio test limit $\pi e^\lambda$",
            "",
            detail_full=f"series {ok_series}; mean {ok_mean}; var {ok_var}; domain {ok_dom}",
        )
    )

    # ---- V20b ----------------------------------------------------------
    ps = np.arange(1, 61)
    worst = []
    ok_b = True
    for rr in (0.5, 2.0, 10.0):
        for pp in (0.3, 0.6, 0.9):
            muv = rr * pp / (1 - pp)
            xmax = int(muv + 200 / -math.log(pp) + 60 / -math.log(pp) * 10)
            xs = np.arange(0, xmax, dtype=float)
            lp = _logpmf(xs, rr, pp)
            dev = np.abs(xs - muv)
            ratios = []
            for p in ps:
                with np.errstate(divide="ignore"):
                    w = lp + p * np.log(dev)
                mx = w.max()
                logm = mx + math.log(np.exp(w - mx).sum())
                ratios.append(math.exp(logm / p) / p)
            ratios = np.array(ratios)
            # bounded: the ratio does not grow; its tail values stay below its early maximum
            if not (ratios[30:].max() <= ratios[:30].max() + 1e-9 and np.isfinite(ratios).all()):
                ok_b = False
            # tail check: the last summand is negligible at the highest p
            last = lp[-1] + ps[-1] * math.log(dev[-1])
            if last > logm - 40:
                ok_b = False
            worst.append((rr, pp, ratios.max(), ratios[-1]))
    big = max(worst, key=lambda t: t[2])
    out.append(
        _rec(
            "V20b",
            r"Moment growth of Prop. 2.8.1(ii) for the centred NB: $(\mathbb E|X-\mu|^p)^{1/p}/p$ stays bounded in $p$",
            "numeric",
            "pass" if ok_b else "fail",
            rf"$r\in\{{0.5,2,10\}}$, $\pi\in\{{0.3,0.6,0.9\}}$, $p=1,\ldots,60$, log-space sums in float64; largest "
            rf"ratio {big[2]:.2f} at $(r,\pi)=({big[0]},{big[1]})$; ratios for $p>30$ never exceed those for $p\le30$",
            "",
            detail_full=json.dumps([[a, b, round(c, 4), round(d, 4)] for a, b, c, d in worst]),
        )
    )

    # ---- V20c ----------------------------------------------------------
    a, m_, rr_ = sp.symbols("a mu r", positive=True)
    pis = m_ / (m_ + rr_)  # pi in terms of mu and r
    Msym = ((1 - pis) / (1 - pis * sp.exp(lam))) ** rr_
    lstar = sp.log(a * (m_ + rr_) / (m_ * (a + rr_)))
    stationary = sp.diff(lam * a - sp.log(Msym), lam).subs(lam, lstar)
    ok_stat = zero(sp.simplify(stationary))
    I_closed = a * sp.log(a * (m_ + rr_) / (m_ * (a + rr_))) + rr_ * sp.log((m_ + rr_) / (a + rr_))
    I_at = sp.expand_log(sp.simplify((lam * a - sp.log(Msym)).subs(lam, lstar)), force=True)
    ok_val = zero(sp.expand_log(I_at - I_closed, force=True))
    ok_mu = zero(sp.simplify(I_closed.subs(a, m_)))
    I0 = sp.limit(I_closed, a, 0)
    ok_zero = zero(sp.expand_log(I0 - rr_ * sp.log(1 + m_ / rr_), force=True))
    ok_pois = zero(sp.limit(I_closed, rr_, sp.oo) - (a * sp.log(a / m_) - a + m_))
    ok_sym = ok_stat and ok_val and ok_mu and ok_zero and ok_pois

    rng = np.random.default_rng(20261101)
    # (1) mean of n NB counts
    muv, rv, n, reps = 3.0, 2.0, 50, 200_000
    piv = muv / (muv + rv)
    xbar = rng.negative_binomial(rv, 1 - piv, size=(reps, n)).mean(axis=1)
    ok_mc1 = True
    rows1 = []
    for u in np.linspace(0.25, 2.0, 8):
        up, lo = (xbar - muv >= u).mean(), (muv - xbar >= u).mean()
        bu, bl = math.exp(-n * _rate(muv + u, muv, rv)), math.exp(-n * _rate(max(muv - u, 0), muv, rv))
        for emp, bd in ((up, bu), (lo, bl)):
            se = math.sqrt(max(emp * (1 - emp), 1e-12) / reps)
            if emp > min(bd, 1.0) + 3 * se:
                ok_mc1 = False
        rows1.append((u, up, bu))
    # (2) RNA-Sieve range distance with NB counts
    K, nb, ck, reps2 = 3, 100, 100, 20_000
    mus = np.array([1.0, 3.0, 6.0])
    rs = np.array([2.0, 2.0, 2.0])
    alpha = np.array([0.98, 0.01, 0.01])
    pis_k = mus / (mus + rs)
    b = float(alpha @ mus)
    types = rng.choice(K, size=(reps2, nb), p=alpha)
    bulk = rng.negative_binomial(rs[types], 1 - pis_k[types]).mean(axis=1)
    ref = np.stack(
        [rng.negative_binomial(rs[k], 1 - pis_k[k], size=(reps2, ck)).mean(axis=1) for k in range(K)], axis=1
    )
    dist = np.maximum.reduce([ref.min(axis=1) - bulk, bulk - ref.max(axis=1), np.zeros(reps2)])
    lam_hi = float(np.min(-np.log(pis_k)))
    grid = np.concatenate([-np.geomspace(50, 1e-6, 4000), np.linspace(1e-6, lam_hi * (1 - 1e-9), 4000)])
    logmix = np.log(
        np.sum(alpha[:, None] * ((1 - pis_k[:, None]) / (1 - pis_k[:, None] * np.exp(grid[None, :]))) ** rs[:, None], axis=0)
    )

    def I0(av):
        return float(np.max(grid * av - logmix)) if av >= 0 else math.inf

    ok_mc2 = True
    rows2 = []
    for tt in np.linspace(0.1, 2.5, 25):
        emp = (dist >= tt).mean()
        bound = math.exp(-nb * I0(b + tt / 2)) + math.exp(-nb * I0(b - tt / 2))
        for k in range(K):
            bound += math.exp(-ck * _rate(mus[k] + tt / 2, mus[k], rs[k]))
            lo = mus[k] - tt / 2
            bound += math.exp(-ck * _rate(lo, mus[k], rs[k])) if lo >= 0 else 0.0
        se = math.sqrt(max(emp * (1 - emp), 1e-12) / reps2)
        if emp > min(bound, 1.0) + 3 * se:
            ok_mc2 = False
        rows2.append((tt, emp, bound))
    t_inf = next((r_ for r_ in rows2 if r_[2] < 1), rows2[-1])
    ok_c = ok_sym and ok_mc1 and ok_mc2
    out.append(
        _rec(
            "V20c",
            r"Cramér--Chernoff rate of the NB, $I(a)=a\log\frac{a(\mu+r)}{\mu(a+r)}+r\log\frac{\mu+r}{a+r}$, with "
            r"$I(\mu)=0$, $I(0)=r\log(1+\mu/r)$ and the Poisson limit; Monte Carlo tails of the mean and of the range "
            r"distance $\delta$ stay below the Chernoff bounds",
            "numeric",
            "pass" if ok_c else "fail",
            rf"Stationary point and closed form symbolic; MC seed 20261101: mean of $n=50$ NB($\mu=3,r=2$), "
            rf"$2\times10^5$ replicates, at $u=1$ tail {rows1[3][1]:.4f} vs bound {rows1[3][2]:.4f}; range filter "
            rf"$K=3$, $\mu_k=(1,3,6)$, $r=2$, $\bm\alpha=(0.98,0.01,0.01)$, $n=c_k=100$, $2\times10^4$ replicates; "
            rf"first bound below 1 at $t={t_inf[0]:.1f}$: {t_inf[2]:.3f} vs {t_inf[1]:.4f}",
            "",
            detail_full=(
                f"stationary {ok_stat}; value {ok_val}; I(mu)=0 {ok_mu}; I(0) {ok_zero}; Poisson {ok_pois}; "
                f"mc mean {ok_mc1}; mc range {ok_mc2}; rows1 {rows1}; rows2 {rows2}"
            ),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

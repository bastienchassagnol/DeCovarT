# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V19: an explicit tail bound for RNA-Sieve's range filter (vignette 03-1,
#sec-fs-cross-platform).

RNA-Sieve (Erdmann-Pham et al. 2021, Genome Res., Eq. 12, PDF p. 11) states
that P[ delta( b_hat_g, [min_k m~_gk, max_k m~_gk] ) >= t ] "decays in t with
sub-Gaussian tails" and prints no constant; delta(p, A) = inf_{a in A}|p - a|.
The bound below is ours, not the paper's.

Setting: per-cell expression of gene g bounded in [0, B]; b_hat = bulk mean
over n cells drawn i.i.d. from the mixture sum_k alpha_k nu_gk, so
E b_hat = b = sum_k alpha_k m_k; m~_k = mean of c_k cells of type k.

Derivation route:
- V19a (exact, linear certificates): (i) b - min_k m_k = sum_k alpha_k
  (m_k - min m) >= 0 when sum alpha = 1, and likewise max m - b >= 0;
  (ii) if |b_hat - b| < t/2 and |m~_k - m_k| < t/2 for all k, then
  min m~ - b_hat < t and b_hat - max m~ < t, so delta < t. Each implication is
  an identity t - (L - b_hat) = sum of terms that are positive under the
  hypotheses, checked by expansion. Hence
  {delta >= t} is contained in {|b_hat - b| >= t/2} or {max_k |m~_k - m_k| >= t/2}.
- V19b (exact): Hoeffding's inequality for a mean of N variables in [0, B],
  P(|X_bar - E X_bar| >= u) <= 2 exp(-2 N u^2 / B^2), at u = t/2 gives
  2 exp(-N t^2 / (2 B^2)); with the union bound,
  P(delta >= t) <= 2 exp(-n t^2/(2B^2)) + 2 sum_k exp(-c_k t^2/(2B^2)).
  Hoeffding's inequality itself is quoted, not checked.
- V19c (numeric): Monte Carlo with bounded (binomial) counts, K = 3 types,
  n = c_k = 100, B = 2, alpha = (0.98, 0.01, 0.01), numpy seed 20261101,
  20,000 replicates; the empirical tail never exceeds the bound on a grid of t.
"""

from __future__ import annotations

import json

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


def run() -> list[dict]:
    out = []

    # ---- V19a ----------------------------------------------------------
    K = 3
    al = sp.symbols("alpha1:4", nonnegative=True)
    m = sp.symbols("m1:4", real=True)
    mmin, mmax = sp.symbols("m_min m_max", real=True)
    b = sum(a * mk for a, mk in zip(al, m))
    # (i) with sum alpha = 1: b - m_min = sum alpha_k (m_k - m_min)
    sub = {al[2]: 1 - al[0] - al[1]}
    ok_i = zero((b - mmin).subs(sub) - sum(a * (mk - mmin) for a, mk in zip(al, m)).subs(sub)) and zero(
        (mmax - b).subs(sub) - sum(a * (mmax - mk) for a, mk in zip(al, m)).subs(sub)
    )
    # (ii) lower side: L = min m~ <= m~_* (k* = argmin m, so m_* = m_min).
    t, bh, bb, mt, ms, L, U = sp.symbols("t bhat b mtilde_star m_star L U", real=True)
    # t - (L - bhat) = (m_* + t/2 - m~_*) + (b - m_min)|_{m_*=m_min} + (bhat - b + t/2) + (m~_* - L)
    lower = (ms + t / 2 - mt) + (bb - ms) + (bh - bb + t / 2) + (mt - L)
    ok_low = zero(lower - (t - (L - bh)))
    # upper side: U = max m~ >= m~^* (k^* = argmax m)
    upper = (mt - ms + t / 2) + (ms - bb) + (bb - bh + t / 2) + (U - mt)
    ok_up = zero(upper - (t - (bh - U)))
    ok_a = ok_i and ok_low and ok_up
    out.append(
        _rec(
            "V19a",
            r"The bulk mean $b_g=\sum_k\alpha_km_{gk}$ lies in $[\min_km_{gk},\max_km_{gk}]$, and "
            r"$\delta(\hat b_g,[\min_k\tilde m_{gk},\max_k\tilde m_{gk}])\ge t$ implies $|\hat b_g-b_g|\ge t/2$ or "
            r"$\max_k|\tilde m_{gk}-m_{gk}|\ge t/2$",
            "exact",
            "pass" if ok_a else "fail",
            r"Linear certificates: $t-(L-\hat b)$ and $t-(\hat b-U)$ written as sums of terms positive under the "
            r"hypotheses, checked by expansion",
            "",
            detail_full=f"convexity {ok_i}; lower {ok_low}; upper {ok_up}",
        )
    )

    # ---- V19b ----------------------------------------------------------
    N, B, u = sp.symbols("N B u", positive=True)
    hoeff = 2 * sp.exp(-2 * N * u**2 / B**2)
    at_half = sp.simplify(hoeff.subs(u, t / 2))
    ok_b = zero(at_half - 2 * sp.exp(-N * t**2 / (2 * B**2)))
    out.append(
        _rec(
            "V19b",
            r"With expression in $[0,B]$, Hoeffding and the union bound give "
            r"$P(\delta\ge t)\le2e^{-nt^2/(2B^2)}+2\sum_ke^{-c_kt^2/(2B^2)}$",
            "exact",
            "pass" if ok_b else "fail",
            r"Hoeffding $2e^{-2Nu^2/B^2}$ at $u=t/2$ equals $2e^{-Nt^2/(2B^2)}$; Hoeffding itself quoted",
            "",
        )
    )

    # ---- V19c ----------------------------------------------------------
    rng = np.random.default_rng(20261101)
    Bn, n, ck, reps = 2, 100, 100, 20000
    alpha = np.array([0.98, 0.01, 0.01])  # bulk mean close to the lower end of the range
    pis = np.array([0.2, 0.3, 0.5])  # binomial(B, pi_k): means 0.4, 0.6, 1.0
    types = rng.choice(3, size=(reps, n), p=alpha)
    bulk = rng.binomial(Bn, pis[types]).mean(axis=1)
    ref = np.stack([rng.binomial(Bn, pis[k], size=(reps, ck)).mean(axis=1) for k in range(3)], axis=1)
    lo, hi = ref.min(axis=1), ref.max(axis=1)
    dist = np.maximum.reduce([lo - bulk, bulk - hi, np.zeros(reps)])
    ok_c = True
    rows = []
    for tt in np.linspace(0.02, 0.8, 40):
        emp = (dist >= tt).mean()
        bound = 2 * np.exp(-n * tt**2 / (2 * Bn**2)) + 3 * 2 * np.exp(-ck * tt**2 / (2 * Bn**2))
        se = np.sqrt(max(emp * (1 - emp), 1e-12) / reps)
        if emp > min(bound, 1.0) + 3 * se:
            ok_c = False
        rows.append((tt, emp, bound))
    t_first = next(r for r in rows if r[2] < 1)  # first t at which the bound is informative
    emp_small = rows[2]
    out.append(
        _rec(
            "V19c",
            r"Monte Carlo: the empirical tail of $\delta$ stays below the bound of V19b",
            "numeric",
            "pass" if ok_c else "fail",
            rf"Binomial counts, $B=2$, $K=3$, $n=c_k=100$, $\bm\alpha=(0.98,0.01,0.01)$, seed 20261101, "
            rf"$2\times10^4$ replicates, 40 values of $t$; tail {emp_small[1]:.3f} at $t={emp_small[0]:.2f}$; first "
            rf"bound below 1 at $t={t_first[0]:.2f}$: {t_first[2]:.3f} against {t_first[1]:.4f} (conservative)",
            "",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

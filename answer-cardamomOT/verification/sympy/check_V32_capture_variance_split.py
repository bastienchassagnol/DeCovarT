# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V32: binomial capture of a mixed-Poisson count and the three-term variance
split (article app:dynamic, eq:dyn-variance-split).

Model. The molecule count of gene g in cell c is Poisson given a random mean
Lambda (upstream biology); each molecule is captured with probability b
(b = b_c b_g), b independent of Lambda; the observed count is
x | x_mol, b ~ Binomial(x_mol, b).

Derivation route:
- V32a: the binomial MGF (1 - b + b e^t)^n, summed from the pmf for small n,
  turns a Poisson(Lambda) count into exp(Lambda((1 - b + b e^t) - 1)) =
  exp(b Lambda (e^t - 1)), the MGF of Poisson(b Lambda): binomial capture of a
  Poisson count is Poisson with mean b Lambda.
- V32b: with x | (Lambda, b) ~ Poisson(b Lambda), the mean and variance are
  the first two derivatives of the CGF of x, obtained from E exp(b Lambda
  (e^t - 1)) with b and Lambda independent and generic moments
  (E b, E b^2, E Lambda, E Lambda^2), expanded to second order in t; compared
  with E b E Lambda + E b^2 Var Lambda + Var b (E Lambda)^2.
- V32c: the law of total variance of Lambda over a two-state label zeta with
  probabilities (w, 1 - w), state means m1, m2 and variances v1, v2.
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


TEX = "article app:dynamic"


def run():
    out = []

    # V32a: binomial thinning of a Poisson count, through MGFs summed from pmfs
    b, t, lam = sp.symbols("b t Lambda", positive=True)
    ok_bin = True
    for n in range(0, 7):
        k = sp.Symbol("k", integer=True, nonnegative=True)
        mgf = sp.summation(sp.binomial(n, k) * b**k * (1 - b) ** (n - k) * sp.exp(t * k), (k, 0, n))
        ok_bin = ok_bin and sp.expand(mgf - (1 - b + b * sp.exp(t)) ** n) == 0
    # Poisson MGF in the variable u: E u^X = exp(Lambda (u - 1)); replace u by the binomial factor
    u = sp.Symbol("u", positive=True)
    kk = sp.Symbol("kk", integer=True, nonnegative=True)
    pgf = sp.simplify(sp.summation(sp.exp(-lam) * (lam * u) ** kk / sp.factorial(kk), (kk, 0, sp.oo)))
    ok_pgf = sp.simplify(pgf - sp.exp(lam * (u - 1))) == 0
    thinned = pgf.subs(u, 1 - b + b * sp.exp(t))
    ok_thin = sp.simplify(sp.log(thinned) - b * lam * (sp.exp(t) - 1)) == 0
    status = "pass" if (ok_bin and ok_pgf and ok_thin) else "fail"
    out.append(
        _rec(
            "V32a",
            r"Binomial capture with probability $b$ of a $\mathrm{Poisson}(\Lambda)$ count is $\mathrm{Poisson}(b\Lambda)$",
            "exact",
            status,
            r"binomial MGF $(1-b+be^{t})^{n}$ summed for $n\le6$; Poisson PGF $e^{\Lambda(u-1)}$ summed; "
            r"composite CGF $=b\Lambda(e^{t}-1)$",
            TEX,
            detail_full=f"binomial MGF ok={ok_bin}; Poisson PGF ok={ok_pgf}; log composite MGF - b*Lambda*(e^t-1) == 0: {ok_thin}",
        )
    )

    # V32b: three-term variance split, from the CGF with generic moments of independent b and Lambda
    Eb, Eb2, EL, EL2 = sp.symbols("E_b E_b2 E_L E_L2", positive=True)
    s = sp.Symbol("s")  # s = e^t - 1, s = t + t^2/2 + O(t^3)
    # E exp(b Lambda s) to second order in s: 1 + s E(b)E(L) + s^2/2 E(b^2)E(L^2)
    mgf_x = 1 + Eb * EL * s + sp.Rational(1, 2) * Eb2 * EL2 * s**2
    mgf_t = sp.series(mgf_x.subs(s, sp.exp(t) - 1), t, 0, 3).removeO()
    cgf = sp.series(sp.log(mgf_t), t, 0, 3).removeO()
    mean = sp.expand(cgf.coeff(t, 1))
    var = sp.expand(2 * cgf.coeff(t, 2))
    VarL = EL2 - EL**2
    Varb = Eb2 - Eb**2
    split = Eb * EL + Eb2 * VarL + Varb * EL**2
    ok_mean = sp.simplify(mean - Eb * EL) == 0
    ok_var = sp.simplify(var - split) == 0
    out.append(
        _rec(
            "V32b",
            r"$\operatorname{Var}x=\mathbb{E}b\,\mathbb{E}\Lambda+\mathbb{E}b^{2}\operatorname{Var}\Lambda"
            r"+\operatorname{Var}b\,(\mathbb{E}\Lambda)^{2}$ for $x\mid(\Lambda,b)\sim\mathrm{Poisson}(b\Lambda)$, $b$ independent of $\Lambda$",
            "exact",
            "pass" if (ok_mean and ok_var) else "fail",
            r"second-order CGF of $x$ from generic moments: mean $" + sp.latex(mean)
            + r"$, variance $" + sp.latex(sp.factor_terms(var)) + r"$",
            TEX,
            detail_full=f"mean={mean}; var={var}; split={sp.expand(split)}; equal={ok_var}",
        )
    )

    # V32c: law of total variance over a two-state label
    w = sp.Symbol("w", positive=True)
    m1, m2, v1, v2 = sp.symbols("m1 m2 v1 v2", positive=True)
    EL_mix = w * m1 + (1 - w) * m2
    EL2_mix = w * (v1 + m1**2) + (1 - w) * (v2 + m2**2)
    var_direct = sp.expand(EL2_mix - EL_mix**2)
    within = w * v1 + (1 - w) * v2
    between = w * (m1 - EL_mix) ** 2 + (1 - w) * (m2 - EL_mix) ** 2
    ok_c = sp.simplify(var_direct - (within + between)) == 0
    out.append(
        _rec(
            "V32c",
            r"$\operatorname{Var}\Lambda=\mathbb{E}_{\zeta}\operatorname{Var}(\Lambda\mid\zeta)+\operatorname{Var}_{\zeta}\mathbb{E}(\Lambda\mid\zeta)$, two states",
            "exact",
            "pass" if ok_c else "fail",
            r"between-state term $" + sp.latex(sp.factor(sp.expand(between))) + r"$",
            TEX,
            detail_full=f"direct={var_direct}; within+between={sp.expand(within + between)}",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V33: binomial capture of a negative binomial and the limiting capture
(article app:dynamic, eq:dyn-limiting-capture).

Model. The molecule count is negative binomial with mean mu and size psi,
variance mu + mu^2/psi (the article's NB convention), built as a Poisson count
whose mean is Gamma(shape psi, rate psi/mu). Each molecule is captured with
probability b.

Derivation route:
- V33a: the PGF of the NB is the Gamma integral of the Poisson PGF
  exp(lambda(u - 1)); capture replaces u by 1 - b + b u; the result is the PGF
  of an NB with mean b mu and the same size psi; its variance, from the PGF,
  is b mu + b^2 mu^2 / psi.
- V33b: the Poisson term b mu exceeds the extra term b^2 mu^2/psi iff
  b < psi/mu; with the bursty NB of CardamomOT, psi = k_on/k_deg and
  mu = k_on k_syn/(k_deg k_off), so psi/mu = k_off/k_syn, the inverse mean
  burst size; at b = 1/10 the threshold burst size is 10.
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


def _var_from_pgf(G, u):
    g1 = sp.diff(G, u).subs(u, 1)
    g2 = sp.diff(G, u, 2).subs(u, 1)
    return sp.simplify(g1), sp.simplify(g2 + g1 - g1**2)


def run():
    out = []
    mu, psi, b, lam = sp.symbols("mu psi b lambda", positive=True)
    u = sp.Symbol("u", positive=True)

    # V33a: NB PGF from the Gamma-Poisson mixture, then capture
    rate = psi / mu
    gamma_pdf = rate**psi * lam ** (psi - 1) * sp.exp(-rate * lam) / sp.gamma(psi)
    # integrate for u in (0, 1), written u = 1 - v with v > 0, so the exponent -lambda v is negative
    v = sp.Symbol("v", positive=True)
    G_v = sp.simplify(sp.integrate(gamma_pdf * sp.exp(-lam * v), (lam, 0, sp.oo), conds="none"))
    G = G_v.subs(v, 1 - u)
    G_ref = (1 + mu * (1 - u) / psi) ** (-psi)
    ok_pgf = sp.simplify(sp.powsimp(G_v / (1 + mu * v / psi) ** (-psi), force=True)) == 1
    G_cap = G_ref.subs(u, 1 - b + b * u)
    G_cap_ref = (1 + b * mu * (1 - u) / psi) ** (-psi)
    ok_cap = sp.simplify(G_cap - G_cap_ref) == 0
    mean_cap, var_cap = _var_from_pgf(G_cap, u)
    ok_mom = sp.simplify(mean_cap - b * mu) == 0 and sp.simplify(var_cap - (b * mu + b**2 * mu**2 / psi)) == 0
    out.append(
        _rec(
            "V33a",
            r"Capture with probability $b$ of $\mathrm{NB}(\mu,\psi)$ gives $\mathrm{NB}(b\mu,\psi)$, of variance $b\mu+b^{2}\mu^{2}/\psi$",
            "exact",
            "pass" if (ok_pgf and ok_cap and ok_mom) else "fail",
            r"Gamma--Poisson PGF $" + sp.latex(G_ref) + r"$; after capture mean $" + sp.latex(mean_cap)
            + r"$, variance $" + sp.latex(sp.expand(var_cap)) + r"$",
            TEX,
            detail_full=f"G={G}; pgf ok={ok_pgf}; capture ok={ok_cap}; mean={mean_cap}; var={var_cap}",
        )
    )

    # V33b: limiting capture and burst size
    diff_terms = sp.factor(b * mu - b**2 * mu**2 / psi)  # = b mu (psi - b mu)/psi
    root = sp.solve(sp.Eq(b * mu, b**2 * mu**2 / psi), b)
    threshold = [r for r in root if r != 0]
    ok_thr = threshold == [psi / mu]
    kon, koff, ksyn, kdeg = sp.symbols("k_on k_off k_syn k_deg", positive=True)
    psi_c = kon / kdeg
    mu_c = kon * ksyn / (kdeg * koff)
    ratio = sp.simplify(psi_c / mu_c)
    ok_ratio = sp.simplify(ratio - koff / ksyn) == 0
    burst_at_tenth = sp.solve(sp.Eq(sp.Rational(1, 10), 1 / sp.Symbol("B", positive=True)), sp.Symbol("B", positive=True))
    ok_num = burst_at_tenth == [10]
    out.append(
        _rec(
            "V33b",
            r"Poisson term larger iff $b<\psi_{g}/\mu_{g}=k^{\mathrm{off}}_{g}/k^{\mathrm{syn}}_{g}$, the inverse mean burst size; "
            r"threshold burst size $10$ at $b=0.1$",
            "exact",
            "pass" if (ok_thr and ok_ratio and ok_num) else "fail",
            r"$b\mu-b^{2}\mu^{2}/\psi=" + sp.latex(diff_terms) + r"$; $\psi/\mu=" + sp.latex(ratio) + r"$",
            TEX,
            detail_full=f"difference={diff_terms}; nonzero root={threshold}; psi/mu={ratio}; burst at b=1/10: {burst_at_tenth}",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

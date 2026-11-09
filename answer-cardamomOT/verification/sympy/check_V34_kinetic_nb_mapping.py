# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V34: kinetic parameters of CardamomOT and of the Gorin et al. mixture limit
mapped to the negative-binomial size psi_g (article app:dynamic).

Derivation route:
- V34a: CardamomOT, Eq. 3: mRNA given proteins ~ Gamma(shape A, rate B)
  mixed into a Poisson, A = k_on/k_deg, B = k_off/k_syn. Mean and variance
  from the law of total variance (Poisson given lambda, Gamma moments
  integrated symbolically): mean A/B = k_on k_syn/(k_deg k_off), Fano factor
  1 + k_syn/k_off, and the article's size psi (variance mu + mu^2/psi) equals A.
- V34b: Gorin et al., mixture limit: stationary rate K ~ Gamma(shape
  k_arr/k_rev, scale varsigma); mature RNA ~ Poisson(K/k_deg) given K; the
  rescaled rate K/k_deg is Gamma(shape k_arr/k_rev, scale varsigma/k_deg), so
  mature counts are NB with size psi = k_arr/k_rev and mean
  k_arr varsigma/(k_rev k_deg).
- V34c: consistency with Gorin et al. Table 1 (their notation kappa = k_rev,
  theta = varsigma, beta = k_spl): sigma_N^2 - mu_N = mu_N theta/(kappa + beta)
  tends, as kappa -> 0, to the NB excess mu_N^2/psi of the nascent marginal
  of the mixture model, psi = a/kappa, mu_N = a theta/(kappa beta).
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


def _gamma_moments(shape, rate):
    lam = sp.Symbol("lambda", positive=True)
    pdf = rate**shape * lam ** (shape - 1) * sp.exp(-rate * lam) / sp.gamma(shape)
    m1 = sp.simplify(sp.integrate(lam * pdf, (lam, 0, sp.oo), conds="none"))
    m2 = sp.simplify(sp.integrate(lam**2 * pdf, (lam, 0, sp.oo), conds="none"))
    return m1, m2


def run():
    out = []
    kon, koff, ksyn, kdeg = sp.symbols("k_on k_off k_syn k_deg", positive=True)

    # V34a: CardamomOT NB moments
    A, B = kon / kdeg, koff / ksyn
    m1, m2 = _gamma_moments(A, B)
    mean = sp.simplify(m1)  # E x = E lambda
    var = sp.simplify(m1 + (m2 - m1**2))  # Var x = E lambda + Var lambda
    fano = sp.simplify(var / mean)
    psi = sp.simplify(mean**2 / (var - mean))
    ok_a = (
        sp.simplify(mean - kon * ksyn / (kdeg * koff)) == 0
        and sp.simplify(fano - (1 + ksyn / koff)) == 0
        and sp.simplify(psi - kon / kdeg) == 0
    )
    out.append(
        _rec(
            "V34a",
            r"CardamomOT negative binomial: mean $k^{\mathrm{on}}k^{\mathrm{syn}}/(k^{\mathrm{deg}}k^{\mathrm{off}})$, "
            r"Fano $1+k^{\mathrm{syn}}/k^{\mathrm{off}}$, size $\psi=k^{\mathrm{on}}/k^{\mathrm{deg}}$",
            "exact",
            "pass" if ok_a else "fail",
            r"mean $" + sp.latex(mean) + r"$, Fano $" + sp.latex(fano) + r"$, $\psi=" + sp.latex(psi) + r"$",
            TEX,
            detail_full=f"mean={mean}; var={var}; fano={fano}; psi={psi}",
        )
    )

    # V34b: Gorin mixture limit for mature RNA
    karr, krev, vs = sp.symbols("k_arr k_rev varsigma", positive=True)
    # K ~ Gamma(shape karr/krev, scale vs) -> K/kdeg ~ Gamma(shape karr/krev, rate kdeg/vs)
    m1b, m2b = _gamma_moments(karr / krev, kdeg / vs)
    mean_b = sp.simplify(m1b)
    var_b = sp.simplify(m1b + (m2b - m1b**2))
    psi_b = sp.simplify(mean_b**2 / (var_b - mean_b))
    ok_b = sp.simplify(psi_b - karr / krev) == 0 and sp.simplify(mean_b - karr * vs / (krev * kdeg)) == 0
    out.append(
        _rec(
            "V34b",
            r"Gorin et al., slow reversion: mature counts negative binomial with $\psi=k^{\mathrm{arr}}/k^{\mathrm{rev}}$, "
            r"mean $k^{\mathrm{arr}}\varsigma/(k^{\mathrm{rev}}k^{\mathrm{deg}})$",
            "exact",
            "pass" if ok_b else "fail",
            r"mean $" + sp.latex(mean_b) + r"$, $\psi=" + sp.latex(psi_b) + r"$",
            TEX,
            detail_full=f"mean={mean_b}; var={var_b}; psi={psi_b}",
        )
    )

    # V34c: Table 1 of Gorin et al. against the mixture NB, kappa -> 0
    a, theta, kappa, beta = sp.symbols("a theta kappa beta", positive=True)
    muN = a * theta / (kappa * beta)
    excess_table = muN * theta / (kappa + beta)  # sigma_N^2 - mu_N, their Table 1
    excess_nb = muN**2 / (a / kappa)  # NB excess with psi = a/kappa
    ratio = sp.simplify(excess_table / excess_nb)
    lim = sp.limit(ratio, kappa, 0)
    first = sp.series(ratio, kappa, 0, 2).removeO()
    ok_c = lim == 1
    out.append(
        _rec(
            "V34c",
            r"Gorin et al.\ Table~1, $\sigma_{N}^{2}-\mu_{N}=\mu_{N}\theta/(\kappa+\beta)$, tends to the mixture excess "
            r"$\mu_{N}^{2}/\psi$, $\psi=a/\kappa$, as $\kappa\to0$",
            "series",
            "pass" if ok_c else "fail",
            r"ratio $" + sp.latex(ratio) + r"=" + sp.latex(first) + r"+O(\kappa^{2})$, limit $" + sp.latex(lim) + r"$",
            TEX,
            detail_full=f"ratio={ratio}; series={first}; limit={lim}",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

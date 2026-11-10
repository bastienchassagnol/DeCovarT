# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V37: size factors of the measurement model (vignette 03-2,
sec-cov-size-factors).

Model. x_gc | s_c, lambda_c ~ Poisson(s_c b_g lambda_gc), s_c = s^bio_c s^tech_c,
with s^bio, s^tech and lambda independent.

V37a: the law depends on (s_c, lambda_c) only through s_c lambda_c, so
(k s_c, lambda_c / k) gives the same Poisson mean (and hence the same law);
the scale of s_c is fixed by E[s^bio | j] = E[s^tech | j] = 1.
V37b: under that constraint, derived on explicit two-point laws (no moment
retyped), E[s_c | j] = 1, E[x_c | j] = B E[lambda_c | j] = mu_j and
E[x_{+c} | s_c, j] = s_c S_j with S_j = 1' mu_j; the plug-in
s_hat_c = x_{+c} / mean(x_{+c'}) averages exactly one over the cells.
V37c: with E s = 1, nu = E s^2 - 1 and 1 + nu = (1 + nu_bio)(1 + nu_tech)
under independence; a dependent pair (s^tech = s^bio) breaks it.
"""

from __future__ import annotations

import json

import sympy as sp


def _rec(id_, claim, kind, status, detail, tex_line="", detail_full=None):
    return {
        "id": id_,
        "claim": claim,
        "kind": kind,
        "status": status,
        "detail": detail,
        "detail_full": detail if detail_full is None else detail_full,
        "tex_line": tex_line,
    }


def _two_point(d):
    """Two-point law with mean 1: values 1 -/+ d, probability 1/2 each."""
    return [(1 - d, sp.Rational(1, 2)), (1 + d, sp.Rational(1, 2))]


def run():
    out = []
    # V37a: scale invariance of the Poisson mean
    s, lam, b, k = sp.symbols("s lambda b k", positive=True)
    ok_a = sp.simplify((k * s) * b * (lam / k) - s * b * lam) == 0
    out.append(
        _rec(
            "V37a",
            r"The Poisson mean $s_cb_g\lambda_{gc}$ is invariant under $(s_c,\boldsymbol\lambda_{\cdot c})\mapsto(ks_c,\boldsymbol\lambda_{\cdot c}/k)$",
            "exact",
            "pass" if ok_a else "fail",
            r"$(ks)\,b\,(\lambda/k)-s\,b\,\lambda=0$",
        )
    )

    # V37b: expectations on explicit laws, G = 2
    d1, d2 = sp.symbols("d_1 d_2", positive=True)
    b1, b2 = sp.symbols("b_1 b_2", positive=True)
    q = sp.Symbol("q", positive=True)
    l = sp.symbols("l11 l12 l21 l22", positive=True)
    lam_law = [(sp.Matrix([l[0], l[1]]), q), (sp.Matrix([l[2], l[3]]), 1 - q)]
    B = sp.diag(b1, b2)
    E_s = sp.expand(
        sum(pb * pt * vb * vt for vb, pb in _two_point(d1) for vt, pt in _two_point(d2))
    )
    E_lam = sum((v * p for v, p in lam_law), sp.zeros(2, 1))
    mu = B * E_lam
    # E[x | s, lambda] = s B lambda (Poisson mean); average over lambda at fixed s
    E_x = sp.zeros(2, 1)
    for vb, pb in _two_point(d1):
        for vt, pt in _two_point(d2):
            for v, p in lam_law:
                E_x += pb * pt * p * vb * vt * (B * v)
    sv = sp.Symbol("s_c", positive=True)
    E_lib_given_s = sum(p * sv * (sp.ones(1, 2) * B * v)[0] for v, p in lam_law)
    S_j = (sp.ones(1, 2) * mu)[0]
    xs = sp.symbols("x_1:5", positive=True)
    shat = [xi / (sum(xs) / 4) for xi in xs]
    ok_b = (
        sp.simplify(E_s - 1) == 0
        and sp.simplify(E_x - mu) == sp.zeros(2, 1)
        and sp.simplify(E_lib_given_s - sv * S_j) == 0
        and sp.simplify(sum(shat) / 4 - 1) == 0
    )
    out.append(
        _rec(
            "V37b",
            r"With $\mathbb E[s^{\mathrm{bio}}]=\mathbb E[s^{\mathrm{tech}}]=1$: $\mathbb E[s_c]=1$, $\mathbb E\boldsymbol x_{\cdot c}=\boldsymbol{\mathcal B}\mathbb E\boldsymbol\lambda_{\cdot c}$, $\mathbb E[x_{+c}\mid s_c]=s_cS_j$, and $\hat s_c$ averages one",
            "exact",
            "pass" if ok_b else "fail",
            r"two-point laws for $s^{\mathrm{bio}}$, $s^{\mathrm{tech}}$ and $\boldsymbol\lambda$ ($G=2$); four cells for $\hat s_c$",
        )
    )

    # V37c: product rule for nu, and its failure under dependence
    E_s2 = sp.expand(
        sum(pb * pt * (vb * vt) ** 2 for vb, pb in _two_point(d1) for vt, pt in _two_point(d2))
    )
    nu = E_s2 - E_s**2
    nu_b = sp.expand(sum(p * v**2 for v, p in _two_point(d1)) - 1)
    nu_t = sp.expand(sum(p * v**2 for v, p in _two_point(d2)) - 1)
    ok_c = (
        sp.simplify(1 + nu - (1 + nu_b) * (1 + nu_t)) == 0
        and sp.simplify(nu - (nu_b + nu_t + nu_b * nu_t)) == 0
    )
    # dependent case: s_tech = s_bio
    E_s_dep = sum(p * v * v for v, p in _two_point(d1))
    E_s2_dep = sum(p * v**4 for v, p in _two_point(d1))
    nu_dep = sp.simplify(E_s2_dep / E_s_dep**2 - 1)
    gap = sp.simplify(1 + nu_dep - (1 + nu_b) ** 2)
    ok_c = ok_c and gap != 0
    out.append(
        _rec(
            "V37c",
            r"$1+\nu=(1+\nu^{\mathrm{bio}})(1+\nu^{\mathrm{tech}})$ for independent factors with mean one; not for $s^{\mathrm{tech}}=s^{\mathrm{bio}}$",
            "exact",
            "pass" if ok_c else "fail",
            r"independent: residual $0$; dependent: $1+\nu-(1+\nu^{\mathrm{bio}})^2=" + sp.latex(gap) + r"$",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

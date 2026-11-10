# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V38: what the raw within-type covariance contains (vignette 03-2,
sec-cov-raw-decomposition and sec-cov-cv).

Derivation route: explicit two-point laws (mean one) for s^bio and s^tech, a
two-point law for lambda (G = 2), the law of total covariance with the
Poisson conditional covariance diag(eta), eta = s^tech s^bio B lambda.
V38a: Cov x = diag(mu) + (1 + nu_tech) Sigma_bio + nu_tech mu mu', with
Sigma_bio = Cov(s^bio B lambda) = (1 + nu_bio) B Sigma_lambda B + nu_bio mu mu'.
V38b: solving for Sigma_bio with kappa = nu_tech / (1 + nu_tech) gives
Sigma_raw - diag(mu) - kappa (Sigma_raw - diag(mu) + mu mu').
V38c: the scoring reconstruction Sigma_bio + diag(mu) + kappa/(1-kappa)
(Sigma_bio + mu mu') returns Sigma_raw.
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
    return [(1 - d, sp.Rational(1, 2)), (1 + d, sp.Rational(1, 2))]


def _moments(atoms):
    """Mean and covariance of a finite law given as [(vector, prob)]."""
    m = sum((p * v for v, p in atoms), sp.zeros(2, 1))
    c = sum((p * (v - m) * (v - m).T for v, p in atoms), sp.zeros(2, 2))
    return m, c


def _zero(M):
    return all(sp.simplify(e) == 0 for e in M)


def run():
    d1, d2, q = sp.symbols("d_1 d_2 q", positive=True)
    b1, b2 = sp.symbols("b_1 b_2", positive=True)
    l = sp.symbols("l11 l12 l21 l22", positive=True)
    B = sp.diag(b1, b2)
    lam_law = [(sp.Matrix([l[0], l[1]]), q), (sp.Matrix([l[2], l[3]]), 1 - q)]
    nu_b, nu_t = d1**2, d2**2  # CV^2 of the two-point laws

    eta = [
        (vt * vb * B * v, pt * pb * p)
        for vb, pb in _two_point(d1)
        for vt, pt in _two_point(d2)
        for v, p in lam_law
    ]
    mu, cov_eta = _moments(eta)
    sigma_raw = sp.diag(*mu) + cov_eta  # law of total covariance, Poisson
    _, sigma_lam = _moments(lam_law)
    mu_b, sigma_bio = _moments(
        [(vb * B * v, pb * p) for vb, pb in _two_point(d1) for v, p in lam_law]
    )
    claim_bio = (1 + nu_b) * B * sigma_lam * B + nu_b * mu * mu.T
    claim_raw = sp.diag(*mu) + (1 + nu_t) * sigma_bio + nu_t * mu * mu.T
    ok_a = _zero(sigma_bio - claim_bio) and _zero(sigma_raw - claim_raw) and _zero(mu_b - mu)

    kappa = nu_t / (1 + nu_t)
    D = sp.diag(*mu)
    inv = sigma_raw - D - kappa * (sigma_raw - D + mu * mu.T)
    ok_b = _zero(inv - sigma_bio)
    rec = sigma_bio + D + kappa / (1 - kappa) * (sigma_bio + mu * mu.T)
    ok_c = _zero(rec - sigma_raw)

    out = [
        _rec(
            "V38a",
            r"$\boldsymbol\Sigma^{\mathrm{raw}}_j=\operatorname{diag}(\tilde{\boldsymbol\mu}_{\cdot j})+(1+\nu^{\mathrm{tech}})\boldsymbol\Sigma^{\mathrm{bio}}_j+\nu^{\mathrm{tech}}\tilde{\boldsymbol\mu}_{\cdot j}\tilde{\boldsymbol\mu}_{\cdot j}^{\top}$, with $\boldsymbol\Sigma^{\mathrm{bio}}_j=(1+\nu^{\mathrm{bio}})\boldsymbol{\mathcal B}\boldsymbol\Sigma^{\lambda}_j\boldsymbol{\mathcal B}+\nu^{\mathrm{bio}}\tilde{\boldsymbol\mu}_{\cdot j}\tilde{\boldsymbol\mu}_{\cdot j}^{\top}$",
            "exact",
            "pass" if ok_a else "fail",
            r"two-point laws, $G=2$; residual matrices $\bm 0$",
        ),
        _rec(
            "V38b",
            r"Inversion with $\kappa=\nu^{\mathrm{tech}}/(1+\nu^{\mathrm{tech}})$ returns $\boldsymbol\Sigma^{\mathrm{bio}}_j$",
            "exact",
            "pass" if ok_b else "fail",
            r"residual $\bm 0$",
        ),
        _rec(
            "V38c",
            r"Scoring reconstruction $\boldsymbol\Sigma^{\mathrm{bio}}+\operatorname{diag}(\tilde{\boldsymbol\mu})+\tfrac{\kappa}{1-\kappa}(\boldsymbol\Sigma^{\mathrm{bio}}+\tilde{\boldsymbol\mu}\tilde{\boldsymbol\mu}^{\top})=\boldsymbol\Sigma^{\mathrm{raw}}$",
            "exact",
            "pass" if ok_c else "fail",
            r"residual $\bm 0$",
        ),
    ]
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

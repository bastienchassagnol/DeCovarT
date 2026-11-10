# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V41-V45: worked numbers and small identities of vignette 03-2.

V41 (sec-grn-attenuation): two genes, log-normal latent with mu_z = log 0.5,
unit variances, correlation 0.6, Poisson counts at s = 1. Moments derived from
the bivariate normal MGF exp(t'mu + t'Sigma t / 2): mean count ~0.82, raw
correlation ~0.28, correlation after removing the Poisson diagonal
(e^0.6 - 1)/(e - 1) ~0.48.
V42 (sec-cov-regimes, layer 3; book table): d^2 log G and d log G with d = 10,
natural log, for G = 1417, 3420, 5821.
V43 (sec-cov-preprocessing): the PLN averaged moments without the Poisson
term, with s-bar = 1 and mean(s^2) = 1 + nu, equal (1 + nu) C + nu m m',
C = (m m') o (e^Sigma_z - 1 1'), which is Sigma_raw - diag(mu) at
Sigma_lambda = C (the kappa = 0 end); with s fixed at 1 they equal C (the
other end).
V44 (tbl-cov-normalisation): dividing each cell by its library turns the
signatures into theta_j = mu_j / S_j; the bulk mean sum_j p_j mu_j equals
T sum_j q_j theta_j with q_j = p_j S_j / T, T = sum_j p_j S_j.
V45 (eq-bulk-depth-per-cell): E y_+ = c 1' mu p = a n sum_j p_j S_j; the
illustration 2e7 / (1e5 * 5e3) = 0.04 and the Fano threshold 1/a = 25.
"""

from __future__ import annotations

import json
import math

import sympy as sp


def _rec(id_, claim, kind, status, detail, tex_line=""):
    return {
        "id": id_,
        "claim": claim,
        "kind": kind,
        "status": status,
        "detail": detail,
        "tex_line": tex_line,
    }


def _v41():
    rho = sp.Rational(3, 5)
    mu0 = sp.log(sp.Rational(1, 2))
    Sig = sp.Matrix([[1, rho], [rho, 1]])
    mu = sp.Matrix([mu0, mu0])

    def mgf(t):
        t = sp.Matrix(t)
        return sp.exp((t.T * mu)[0] + (t.T * Sig * t)[0] / 2)

    m = mgf([1, 0])
    var_lat = mgf([2, 0]) - m**2
    cov_lat = mgf([1, 1]) - m**2
    raw = cov_lat / (m + var_lat)  # Poisson adds the mean to each variance
    den = cov_lat / var_lat
    vals = [float(m), float(raw), float(den)]
    ok = (
        round(vals[0], 2) == 0.82
        and round(vals[1], 2) == 0.28
        and round(vals[2], 2) == 0.48
        and sp.simplify(den - (sp.exp(rho) - 1) / (sp.E - 1)) == 0
    )
    return _rec(
        "V41",
        r"Attenuation example: mean count $0.82$, raw correlation $0.28$, denoised $(e^{0.6}-1)/(e-1)\approx0.48$",
        "numeric",
        "pass" if ok else "fail",
        "values " + ", ".join(f"{v:.4f}" for v in vals) + r" (rounded to 2 decimals)",
    )


def _v42():
    vals = {g: (100 * math.log(g), 10 * math.log(g)) for g in (1417, 3420, 5821)}
    expected = {1417: (726, 73), 3420: (814, 81), 5821: (867, 87)}
    ok = all(round(vals[g][0]) == e[0] and round(vals[g][1]) == e[1] for g, e in expected.items())
    return _rec(
        "V42",
        r"Support-recovery orders with $d=10$: $d^2\log G$ and $d\log G$ at $G=1417,3420,5821$",
        "numeric",
        "pass" if ok else "fail",
        "; ".join(f"G={g}: {a:.1f}, {b:.2f}" for g, (a, b) in vals.items()),
    )


def _v43():
    m1, m2, s11, s12, s22, nu = sp.symbols("m1 m2 s11 s12 s22 nu", positive=True)
    m = sp.Matrix([m1, m2])
    E = sp.Matrix([[sp.exp(s11), sp.exp(s12)], [sp.exp(s12), sp.exp(s22)]])
    one = sp.ones(2, 2)
    C = sp.matrix_multiply_elementwise(m * m.T, E - one)
    avg_no_poisson = sp.matrix_multiply_elementwise(m * m.T, (1 + nu) * E - one)
    # Sigma_raw - diag(mu) with Sigma_lambda = C, total extrinsic noise nu
    target_k0 = (1 + nu) * C + nu * m * m.T
    fixed_s = sp.matrix_multiply_elementwise(m * m.T, E - one)
    ok = all(sp.simplify(e) == 0 for e in avg_no_poisson - target_k0) and fixed_s == C
    return _rec(
        "V43",
        r"PLN export brackets: averaged map without the Poisson term $=(1+\nu)\boldsymbol C+\nu\boldsymbol m\boldsymbol m^{\top}$ ($\kappa=0$); $s$ fixed gives $\boldsymbol C$",
        "exact",
        "pass" if ok else "fail",
        r"$G=2$; residual $\bm 0$",
    )


def _v44():
    p1, S1, S2 = sp.symbols("p1 S1 S2", positive=True)
    p = [p1, 1 - p1]
    S = [S1, S2]
    th = [sp.Matrix(sp.symbols(f"t{j}_1:4", positive=True)) for j in (1, 2)]
    th = [t / sum(t) for t in th]  # unit-sum signatures
    mu = [S[j] * th[j] for j in range(2)]
    T = sum(p[j] * S[j] for j in range(2))
    q = [p[j] * S[j] / T for j in range(2)]
    lhs = sum((p[j] * mu[j] for j in range(2)), sp.zeros(3, 1))
    rhs = T * sum((q[j] * th[j] for j in range(2)), sp.zeros(3, 1))
    ok = all(sp.simplify(e) == 0 for e in lhs - rhs)
    return _rec(
        "V44",
        r"Library-normalised signatures turn cell proportions into $q_j=p_jS_j/\sum_\ell p_\ell S_\ell$",
        "exact",
        "pass" if ok else "fail",
        r"$G=3$, $J=2$; residual $\bm 0$",
    )


def _v45():
    a, n, p1, S1, S2 = sp.symbols("a n p1 S1 S2", positive=True)
    c = a * n
    Ey = c * (p1 * S1 + (1 - p1) * S2)  # c 1' mu p with 1' mu_j = S_j
    ok = sp.simplify(Ey / (n * (p1 * S1 + (1 - p1) * S2)) - a) == 0
    a_ill = sp.Rational(2 * 10**7, 10**5 * 5 * 10**3)
    ok = ok and a_ill == sp.Rational(1, 25) and 1 / a_ill == 25
    return _rec(
        "V45",
        r"$a_i=\mathbb E\,y_{+i}/(n_i\sum_jp_{ji}S_j)$; illustration $a_i=0.04$, Fano threshold $25$",
        "exact",
        "pass" if ok else "fail",
        r"$2\times10^7/(10^5\cdot5\times10^3)=1/25$",
    )


def run():
    return [_v41(), _v42(), _v43(), _v44(), _v45()]


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

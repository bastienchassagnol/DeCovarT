# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V22: CSNet's multivariate negative-binomial simulator and the factor
covariance (vignette S5, #sec-s5-review-csnet and #sec-s5-review-hadaca).

CSNet (Su et al., bioRxiv 10.1101/2021.12.21.473558, supplement A1,
Algorithm 1, PDF pp. 33-35): theta_g = phi * (sum of T Exp(1) draws), of
which T_gg' are shared with gene g'; x_g | theta_g ~ Poisson(s theta_g).

- V22a (exact, our derivation): with theta_g = phi (U + A_g),
  theta_g' = phi (U + A_g'), U ~ Gamma(T_gg', 1), A_g ~ Gamma(T - T_gg', 1)
  independent, the law of total (co)variance gives E x = s T phi,
  Var x = s^2 T phi^2 + s T phi, Cov(x_g, x_g') = s^2 phi^2 T_gg', hence
  Cor = (T_gg'/T) b with b = s phi / (1 + s phi); Var x = sigma^2 solves to
  s phi = sqrt(sigma^2/T + 1/4) - 1/2.
- V22b (numeric): b = 0.77 at s = 6e7, sigma^2 = e^8, T = 200 (paper p. 35)
  and E x ~ 678; Monte Carlo of Algorithm 1 (MA(1), rho = 0.5, p = 6,
  seed 20261101) matches mean, variance and lag-one correlation 0.5 b;
  the row sums of an AR-type block 0.9 * 0.9^(|d|-1) over 33 genes exceed 1,
  so Step 2.2 would need a negative number of private draws.
- V22c (exact): Woodbury, (Psi + B B')^{-1} = Psi^{-1} - Psi^{-1} B
  (I_q + B' Psi^{-1} B)^{-1} B' Psi^{-1}, for G = 4, q = 2 symbolic.
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


def run() -> list[dict]:
    out = []

    # ---- V22a ----------------------------------------------------------
    T, Tgg, phi, s, sig2 = sp.symbols("T T_gg phi s sigma2", positive=True)
    # gamma(k, 1): mean k, variance k; Poisson(s theta): E = s theta, Var = s theta
    E_theta = phi * T
    Var_theta = phi**2 * T
    Cov_theta = phi**2 * Tgg  # only the shared Gamma(T_gg, 1) part is common
    E_x = s * E_theta
    Var_x = s * E_theta + s**2 * Var_theta  # E Var(x|theta) + Var E(x|theta)
    Cov_x = s**2 * Cov_theta  # conditional independence given the rates
    ok_mom = zero(E_x - s * T * phi) and zero(Var_x - (s**2 * T * phi**2 + s * T * phi)) and zero(
        Cov_x - s**2 * phi**2 * Tgg
    )
    cor = Cov_x / Var_x
    b = s * phi / (1 + s * phi)
    ok_cor = zero(sp.simplify(cor - Tgg / T * b))
    u = sp.Symbol("u", positive=True)  # u = s phi
    sol = sp.solve(sp.Eq(T * u**2 + T * u, sig2), u)
    target = sp.sqrt(sig2 / T + sp.Rational(1, 4)) - sp.Rational(1, 2)
    ok_phi = any(zero(sp.radsimp(sv - target)) for sv in sol)
    ok_a = ok_mom and ok_cor and ok_phi
    out.append(
        _rec(
            "V22a",
            r"CSNet's Algorithm 1: $\mathbb{E}x=sT\phi$, $\operatorname{Var}x=s^2T\phi^2+sT\phi$, "
            r"$\operatorname{Cov}(x_g,x_{g'})=s^2\phi^2T_{gg'}$, $\operatorname{Cor}=(T_{gg'}/T)\,s\phi/(1+s\phi)$; "
            r"$\operatorname{Var}x=\sigma^2$ gives $s\phi=\sqrt{\sigma^2/T+1/4}-1/2$",
            "exact",
            "pass" if ok_a else "fail",
            r"Law of total (co)variance with a shared $\mathrm{Gamma}(T_{gg'},1)$ component; positive root of "
            r"$Tu^2+Tu=\sigma^2$",
            "",
            detail_full=f"moments {ok_mom}; correlation {ok_cor}; phi {ok_phi}",
        )
    )

    # ---- V22b ----------------------------------------------------------
    Tn, sn, s2n = 200, 6e7, math.exp(8)
    uphi = math.sqrt(s2n / Tn + 0.25) - 0.5
    bn = uphi / (1 + uphi)
    mean_n = Tn * uphi
    ok_const = abs(bn - 0.77) < 0.005 and abs(mean_n - 678) < 1
    rng = np.random.default_rng(20261101)
    G, rho, reps = 6, 0.5, 200_000
    s_mc, phi_mc = 1.0, 2.0  # small depth so the Poisson part matters
    Tsh = int(Tn * rho)
    shared = rng.gamma(Tsh, 1.0, size=(reps, G - 1))  # one shared block per neighbour pair
    priv = np.empty((reps, G))
    for g in range(G):
        n_nb = (g > 0) + (g < G - 1)
        priv[:, g] = rng.gamma(Tn - n_nb * Tsh, 1.0, size=reps) if Tn - n_nb * Tsh > 0 else 0.0
    theta = priv.copy()
    theta[:, 1:] += shared
    theta[:, :-1] += shared
    theta *= phi_mc
    x = rng.poisson(s_mc * theta)
    e_mean = s_mc * Tn * phi_mc
    e_var = s_mc**2 * Tn * phi_mc**2 + s_mc * Tn * phi_mc
    b_mc = s_mc * phi_mc / (1 + s_mc * phi_mc)
    e_cor = rho * b_mc
    m_hat, v_hat = x.mean(), x.var(axis=0).mean()
    c_hat = np.mean([np.corrcoef(x[:, g], x[:, g + 1])[0, 1] for g in range(G - 1)])
    c_far = np.mean([np.corrcoef(x[:, g], x[:, g + 2])[0, 1] for g in range(G - 2)])
    ok_mc = (
        abs(m_hat - e_mean) / e_mean < 0.01
        and abs(v_hat - e_var) / e_var < 0.02
        and abs(c_hat - e_cor) < 0.01
        and abs(c_far) < 0.01
    )
    d = np.arange(1, 33)
    row_ar = 2 * np.sum(0.9 * 0.9 ** (d[:16] - 1))  # middle gene of a 33-gene block
    ok_feas = row_ar > 1
    ok_b = ok_const and ok_mc and ok_feas
    out.append(
        _rec(
            "V22b",
            r"CSNet constants and simulation: $b=0.77$ and $\mathbb{E}x\approx678$ at $s=6\times10^7$, "
            r"$\sigma^2=e^8$, $T=200$; Monte Carlo of Algorithm 1 matches the moments; an AR-type block "
            r"violates $\sum_{g'}T_{gg'}\le T$",
            "numeric",
            "pass" if ok_b else "fail",
            rf"$b={bn:.4f}$, $\mathbb{{E}}x={mean_n:.1f}$; MC MA(1), $\rho=0.5$, $G=6$, $s=1$, $\phi=2$, "
            rf"$2\times10^5$ draws, seed 20261101: mean {m_hat:.1f} (exact {e_mean:.0f}), variance {v_hat:.0f} "
            rf"(exact {e_var:.0f}), lag-one correlation {c_hat:.3f} (exact {e_cor:.3f}), lag two {c_far:.3f}; "
            rf"AR-type row sum {row_ar:.1f}",
            "",
        )
    )

    # ---- V22c ----------------------------------------------------------
    Gs, q = 4, 2
    B = sp.Matrix(Gs, q, lambda i, j: sp.Symbol(f"b{i}{j}", real=True))
    Psi = sp.diag(*sp.symbols("psi1:5", positive=True))
    lhs = (Psi + B * B.T).inv()
    Pi = Psi.inv()
    rhs = Pi - Pi * B * (sp.eye(q) + B.T * Pi * B).inv() * B.T * Pi
    diff = (lhs - rhs).applyfunc(lambda e: sp.cancel(sp.together(e)))
    ok_c = all(e == 0 for e in diff)
    out.append(
        _rec(
            "V22c",
            r"Woodbury: $(\bm\Psi+\bm B\bm B^\top)^{-1}=\bm\Psi^{-1}-\bm\Psi^{-1}\bm B(\mathbf I_q+\bm B^\top"
            r"\bm\Psi^{-1}\bm B)^{-1}\bm B^\top\bm\Psi^{-1}$",
            "exact",
            "pass" if ok_c else "fail",
            r"$G=4$, $q=2$, symbolic loadings and diagonal $\bm\Psi$; every entry of the difference cancels",
            "",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

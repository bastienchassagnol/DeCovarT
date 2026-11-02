# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V26: dynamic and mixed Gaussian Bayesian networks (note, appendix on purely
statistical approaches; summary_side_comparison.qmd, perspectives).

- V26a (exact): linear Gaussian DBN x_{t+1} = A x_t + e, e ~ N(0, Q). The
  stationary covariance solves the discrete Lyapunov (Stein) equation
  Sigma = A Sigma A' + Q, i.e. (I - A (x) A) vec Sigma = vec Q. Two genes with
  persistence a1, a2 and a feedback loop b12, b21 (edges 2_t -> 1_{t+1} and
  1_t -> 2_{t+1}); the solution is checked symbolically. A positive loop gives
  sigma_12 != 0; a symmetric negative loop (b12 = -b21, equal rates and noise)
  gives sigma_12 = 0 while the lag-one cross-covariance (A Sigma)_12 != 0: the
  loop is invisible in one snapshot and visible across time.
- V26b (numeric): snapshots without couplings. With M = Sigma_t^{1/2} and U
  orthogonal fixing M^{-1} mu_t, A' = A M U M^{-1} gives the same mu_{t+1} and
  Sigma_{t+1} as A (same Q), with A' != A: first and second moments of
  consecutive snapshots do not identify A (our remark). G = 3, seed 20261102.
- V26c (exact): mixture covariance, law of total covariance:
  Cov = sum_k pi_k Sigma_k + sum_k pi_k (mu_k - mubar)(mu_k - mubar)'; two
  classes with diagonal Sigma_k give a non-zero off-diagonal
  pi_1 pi_2 (mu_1 - mu_2)_1 (mu_1 - mu_2)_2.
- V26d (exact): X -> Y and Y -> X with matched parameters share the
  observational Gaussian (I-equivalence) but differ under do(X = x0).
- V26e (exact): v-structure X -> Z <- Y: Sigma_xy = 0 while J_xy != 0
  (Koller & Friedman, Ex. 7.5).
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
    return sp.simplify(e) == 0


def run() -> list[dict]:
    out = []

    # ---- V26a ----------------------------------------------------------
    a1, a2, b12, b21 = sp.symbols("a1 a2 b12 b21", real=True)
    q1, q2 = sp.symbols("q1 q2", positive=True)
    A = sp.Matrix([[a1, b12], [b21, a2]])
    Q = sp.diag(q1, q2)
    s11, s12, s22 = sp.symbols("s11 s12 s22", real=True)
    S = sp.Matrix([[s11, s12], [s12, s22]])
    eqs = list((S - A * S * A.T - Q).applyfunc(sp.expand))
    sol = sp.solve([eqs[0], eqs[1], eqs[3]], [s11, s12, s22], dict=True)[0]
    Ssol = S.subs(sol)
    resid = (Ssol - A * Ssol * A.T - Q).applyfunc(sp.simplify)
    # vec form: (I - A (x) A) vec S = vec Q
    kron = sp.kronecker_product(A, A)
    vecS = sp.Matrix([Ssol[0, 0], Ssol[1, 0], Ssol[0, 1], Ssol[1, 1]])
    vecQ = sp.Matrix([Q[0, 0], Q[1, 0], Q[0, 1], Q[1, 1]])
    resid_vec = ((sp.eye(4) - kron) * vecS - vecQ).applyfunc(sp.simplify)
    s12_expr = sp.simplify(Ssol[0, 1])
    # positive loop (1 -> 2 and 2 -> 1 both activating): contemporaneous covariance
    pos = {a1: sp.Rational(1, 2), a2: sp.Rational(1, 2), b12: sp.Rational(1, 4), b21: sp.Rational(1, 4), q1: 1, q2: 1}
    # negative loop (1 activates 2, 2 represses 1), symmetric rates: sigma_12 = 0,
    # yet the lag-one cross-covariance Cov(x_{t+1}, x_t) = A Sigma is not zero
    neg = {a1: sp.Rational(1, 2), a2: sp.Rational(1, 2), b12: -sp.Rational(1, 4), b21: sp.Rational(1, 4), q1: 1, q2: 1}
    rho_pos = max(abs(complex(v)) for v in A.subs(pos).eigenvals())
    rho_neg = max(abs(complex(v)) for v in A.subs(neg).eigenvals())
    s12_pos = sp.nsimplify(s12_expr.subs(pos))
    s12_neg = sp.nsimplify(s12_expr.subs(neg))
    lag_neg = sp.nsimplify((A * Ssol).subs(neg)[0, 1])
    ok_a = (
        resid == sp.zeros(2, 2)
        and resid_vec == sp.zeros(4, 1)
        and s12_pos != 0
        and s12_neg == 0
        and lag_neg != 0
        and rho_pos < 1
        and rho_neg < 1
    )
    out.append(
        _rec(
            "V26a",
            r"Linear Gaussian DBN: the stationary covariance solves $\bm\Sigma=\bm A\bm\Sigma\bm A^\top+\bm Q$, "
            r"$(\mathbf I-\bm A\otimes\bm A)\operatorname{vec}\bm\Sigma=\operatorname{vec}\bm Q$; a positive two-gene loop "
            r"between slices gives $\sigma_{12}\neq0$, a symmetric negative loop gives $\sigma_{12}=0$ but a non-zero "
            r"lag-one cross-covariance $(\bm A\bm\Sigma)_{12}$",
            "exact",
            "pass" if ok_a else "fail",
            rf"Symbolic $2\times2$ $\bm A$ with persistence $1/2$ and cross edges $\pm1/4$, $\bm Q=\mathbf I$; residuals "
            rf"vanish; positive loop $\sigma_{{12}}={s12_pos}$; negative loop $\sigma_{{12}}={s12_neg}$, "
            rf"$(\bm A\bm\Sigma)_{{12}}={lag_neg}$; spectral radii {rho_pos:.3f}, {rho_neg:.3f}",
            "",
        )
    )

    # ---- V26b ----------------------------------------------------------
    rng = np.random.default_rng(20261102)
    G = 3
    A0 = rng.normal(size=(G, G))
    A0 *= 0.8 / max(abs(np.linalg.eigvals(A0)))
    L = rng.normal(size=(G, G))
    St = L @ L.T + G * np.eye(G)
    mut = rng.normal(size=G)
    w, V = np.linalg.eigh(St)
    M = V @ np.diag(np.sqrt(w)) @ V.T
    v = np.linalg.solve(M, mut)
    v /= np.linalg.norm(v)
    # orthonormal basis with v first, rotate the complement by an angle
    B = np.linalg.qr(np.column_stack([v, rng.normal(size=(G, G - 1))]))[0]
    th = 0.9
    R = np.eye(G)
    R[1:, 1:] = [[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]]
    U = B @ R @ B.T
    A1 = A0 @ M @ U @ np.linalg.inv(M)
    err_mu = np.abs(A1 @ mut - A0 @ mut).max()
    err_S = np.abs(A1 @ St @ A1.T - A0 @ St @ A0.T).max()
    gap = np.abs(A1 - A0).max()
    ok_b = err_mu < 1e-10 and err_S < 1e-10 and gap > 1e-2
    out.append(
        _rec(
            "V26b",
            r"Without couplings, $\bm A'=\bm A\bm M\bm U\bm M^{-1}$ ($\bm M=\bm\Sigma_t^{1/2}$, $\bm U$ orthogonal fixing "
            r"$\bm M^{-1}\bm\mu_t$) gives the same $\bm\mu_{t+1}$ and $\bm\Sigma_{t+1}$ as $\bm A$",
            "numeric",
            "pass" if ok_b else "fail",
            rf"$G=3$, seed 20261102: mean error {err_mu:.1e}, covariance error {err_S:.1e}, "
            rf"$\max|\bm A'-\bm A|={gap:.2f}$",
            "",
        )
    )

    # ---- V26c ----------------------------------------------------------
    p1 = sp.Symbol("pi1", positive=True)
    p2 = 1 - p1
    m11, m12, m21, m22, d11, d12, d21, d22 = sp.symbols("m11 m12 m21 m22 d11 d12 d21 d22", real=True)
    mu1, mu2 = sp.Matrix([m11, m12]), sp.Matrix([m21, m22])
    S1, S2 = sp.diag(d11, d12), sp.diag(d21, d22)
    mbar = p1 * mu1 + p2 * mu2
    second = p1 * (S1 + mu1 * mu1.T) + p2 * (S2 + mu2 * mu2.T)
    cov = (second - mbar * mbar.T).applyfunc(sp.expand)
    formula = (p1 * S1 + p2 * S2 + p1 * (mu1 - mbar) * (mu1 - mbar).T + p2 * (mu2 - mbar) * (mu2 - mbar).T).applyfunc(
        sp.expand
    )
    off = sp.factor(cov[0, 1])
    ok_c = (cov - formula).applyfunc(sp.simplify) == sp.zeros(2, 2) and zero(
        off - p1 * p2 * (m11 - m21) * (m12 - m22)
    )
    out.append(
        _rec(
            "V26c",
            r"Mixture covariance $\sum_k\pi_k\bm\Sigma_k+\sum_k\pi_k(\bm\mu_k-\bar{\bm\mu})(\bm\mu_k-\bar{\bm\mu})^\top$: "
            r"diagonal classes give the off-diagonal $\pi_1\pi_2(\mu_{11}-\mu_{21})(\mu_{12}-\mu_{22})$",
            "exact",
            "pass" if ok_c else "fail",
            r"Second moments of a two-class mixture, $G=2$, symbolic",
            "",
        )
    )

    # ---- V26d ----------------------------------------------------------
    b, sx, se, x0 = sp.symbols("b s_x s_e x0", positive=True)
    # model 1: X ~ N(0, sx^2), Y = b X + e
    C1 = sp.Matrix([[sx**2, b * sx**2], [b * sx**2, b**2 * sx**2 + se**2]])
    # model 2: Y ~ N(0, sy^2), X = c Y + f, matched
    sy2 = b**2 * sx**2 + se**2
    c = b * sx**2 / sy2
    sf2 = sx**2 - c**2 * sy2
    C2 = sp.Matrix([[c**2 * sy2 + sf2, c * sy2], [c * sy2, sy2]])
    same_obs = (C1 - C2).applyfunc(sp.simplify) == sp.zeros(2, 2) and sp.simplify(sf2) != 0
    # do(X = x0): model 1 gives Y ~ N(b x0, se^2); model 2 leaves Y ~ N(0, sy2)
    differ = sp.simplify(b * x0) != 0 and sp.simplify(se**2 - sy2) != 0
    ok_d = same_obs and differ
    out.append(
        _rec(
            "V26d",
            r"$X\to Y$ and $Y\to X$ with matched parameters have the same observational law but differ under "
            r"$\mathrm{do}(X=x_0)$: $Y\sim\mathcal N(bx_0,s_e^2)$ against $\mathcal N(0,b^2s_x^2+s_e^2)$",
            "exact",
            "pass" if ok_d else "fail",
            r"Covariances of the two SEMs equal symbolically; interventional means and variances differ for $b\neq0$",
            "",
        )
    )

    # ---- V26e ----------------------------------------------------------
    d = sp.Symbol("d", positive=True)
    Bm = sp.Matrix([[0, 0, 0], [0, 0, 0], [1, 1, 0]])  # order X, Y, Z; Z = X + Y + e
    D = sp.diag(1, 1, d)
    I3 = sp.eye(3)
    Sig = (I3 - Bm).inv() * D * (I3 - Bm).inv().T
    J = Sig.inv().applyfunc(sp.simplify)
    ok_e = Sig[0, 1] == 0 and zero(J[0, 1] - 1 / d)
    out.append(
        _rec(
            "V26e",
            r"v-structure $X\to Z\leftarrow Y$: $\Sigma_{xy}=0$ but $J_{xy}=1/d\neq0$ (moralisation in the precision)",
            "exact",
            "pass" if ok_e else "fail",
            r"$Z=X+Y+e$, $e\sim\mathcal N(0,d)$; covariance and precision computed symbolically",
            "",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

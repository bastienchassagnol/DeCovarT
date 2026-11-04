# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V28: leverage and Cook's distance of cells and genes (vignette 03-1,
#sec-fs-leverage, @eq-fs-cook-cells and @eq-fs-hat-genes).

Derivation route:
- V28a: the GLM M_1 (log link, intercept + J - 1 type dummies, offset
  log s_c). The IRLS weight is derived as (d mu / d eta)^2 / V(mu) with
  mu = exp(eta) and the NB variance V(m) = m + m^2/psi (giving m^2/V(m)), then
  the hat matrix H = W^(1/2) X (X'WX)^(-1) X' W^(1/2) (handled as the similar
  matrix X (X'WX)^(-1) X' W, same diagonal, trace and idempotence) is built from the design
  with symbolic weights (J = 3 types of 3, 2 and 2 cells). Its diagonal is
  compared with w_c / sum_{c' in type} w_c'; trace = J, H^2 = H. The weight
  m^2/V(m) is increasing in m (derivative sign), so within a type h grows with
  the fitted mean; within a type the h sum to 1.
- V28b: Cook's distance as the one-step deletion distance
  (b - b_(c))' X'WX (b - b_(c)) / (J phi) at fixed IRLS weights. b_(c) is the
  weighted least-squares fit of the working response with row c removed,
  computed exactly (symbolic weights and working responses, the same design).
  With the working residual e = (x - m) d eta/d mu and the Pearson residual
  r = (x - m)/sqrt(V(m)), w e^2 = r^2. Compared with both printed forms,
  with r^SP = r / sqrt(phi (1 - h)).
- V28c: the weighted bulk-on-reference regression (no intercept), G x J
  design X, W = diag(1/tau_g^2). The hat diagonal is compared with
  mu_g'(X'WX)^(-1) mu_g / tau_g^2 (symbolic G = 4, J = 2), trace = J and
  idempotence (also for rational G = 7, J = 3), invariance of h_g when
  (mu_g, tau_g) -> (c mu_g, c tau_g), and Cook_g = r_g^2 h_g / (J (1 - h_g))
  against an exact refit without gene g (rational data, s^2 from the full fit).
"""

from __future__ import annotations

import json
import random

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
    """cancel first (rational functions of the weights), simplify as fallback."""
    if sp.cancel(sp.together(e)) == 0:
        return True
    return sp.simplify(e) == 0


def one_way_design(sizes):
    """Intercept + dummies for types 2..J (beta_g1 = 0), rows grouped by type."""
    J = len(sizes)
    rows, types = [], []
    for j, n in enumerate(sizes):
        for _ in range(n):
            rows.append([1] + [1 if j == l else 0 for l in range(1, J)])
            types.append(j)
    return sp.Matrix(rows), types


def wls(X, w, z):
    W = sp.diag(*w)
    A = X.T * W * X
    b = (A.inv(method="ADJ") * (X.T * W * z)).applyfunc(sp.cancel)
    return b, A


def run() -> list[dict]:
    out = []

    # ---- V28a ---------------------------------------------------------------
    eta, psi, m = sp.symbols("eta psi m", positive=True)
    mu_eta = sp.exp(eta)
    V = lambda mm: mm + mm**2 / psi  # NB variance function
    w_def = sp.diff(mu_eta, eta) ** 2 / V(mu_eta)
    ok_w = zero(w_def.subs(eta, sp.log(m)) - m**2 / V(m))
    dw = sp.simplify(sp.diff(m**2 / V(m), m))
    ok_incr = zero(dw - psi**2 / (m + psi) ** 2)  # > 0
    sizes = (3, 2, 2)
    J = len(sizes)
    X, types = one_way_design(sizes)
    C = X.shape[0]
    w = sp.symbols(f"w1:{C + 1}", positive=True)
    W = sp.diag(*w)
    A = X.T * W * X
    # P = X (X'WX)^(-1) X' W = W^(-1/2) H W^(1/2): same diagonal, trace and
    # idempotence as H, without square roots
    H = (X * A.inv(method="ADJ") * X.T * W).applyfunc(sp.cancel)
    ok_diag = all(
        zero(H[c, c] - w[c] / sum(w[k] for k in range(C) if types[k] == types[c])) for c in range(C)
    )
    ok_tr = zero(H.trace() - J)
    ok_idem = all(zero(e) for e in (H * H - H))
    ok_sum1 = all(
        zero(sum(H[c, c] for c in range(C) if types[c] == j) - 1) for j in range(J)
    )
    ok_a = ok_w and ok_incr and ok_diag and ok_tr and ok_idem and ok_sum1
    out.append(
        _rec(
            "V28a",
            r"One-way log-link GLM $M_1$: $h_{gc}=w_{gc}/\sum_{c':j(c')=j(c)}w_{gc'}$ with "
            r"$w_{gc}=\hat m_{gc}^2/V(\hat m_{gc})$; the hat values of a type sum to 1 and grow with the fitted mean",
            "exact",
            "pass" if ok_a else "fail",
            rf"IRLS weight $(d\mu/d\eta)^2/V$ at $\mu=e^\eta$ equals $m^2/V(m)$, "
            rf"$\partial_m\{{m^2/V(m)\}}={sp.latex(dw)}>0$; $\bm H$ from the design (types of 3, 2, 2 cells, "
            rf"symbolic weights): diagonal as printed, trace $J=3$, $\bm H^2=\bm H$",
            "",
            detail_full=(
                f"weight {ok_w}; increasing {ok_incr}; diag {ok_diag}; trace {ok_tr}; idempotent {ok_idem}; "
                f"type sums {ok_sum1}"
            ),
        )
    )

    # ---- V28b ---------------------------------------------------------------
    z = sp.Matrix(sp.symbols(f"z1:{C + 1}", real=True))  # working responses
    phi = sp.Symbol("phi", positive=True)
    b, A = wls(X, list(w), z)
    resid = z - X * b  # working residuals e_c
    ok_b_cells = True
    for c in (0, 3):  # one cell of a 3-cell type, one of a 2-cell type
        keep = [k for k in range(C) if k != c]
        Xd = X.extract(keep, list(range(J)))
        bd, _ = wls(Xd, [w[k] for k in keep], z.extract(keep, [0]))
        diff = (b - bd).applyfunc(sp.cancel)
        cook_def = sp.cancel((diff.T * A * diff)[0, 0]) / (J * phi)
        h = H[c, c]
        e = sp.cancel(resid[c])
        r2 = w[c] * e**2  # w e^2 = r^2 (Pearson residual squared)
        cook_pearson = r2 * h / (J * phi * (1 - h) ** 2)
        rsp2 = r2 / (phi * (1 - h))
        cook_sp = rsp2 / J * h / (1 - h)
        if not (zero(cook_def - cook_pearson) and zero(cook_def - cook_sp)):
            ok_b_cells = False
    # working residual and Pearson residual: w e^2 = r^2
    x_ = sp.Symbol("x", positive=True)
    e_work = (x_ - m) * sp.diff(sp.log(m), m)  # (x - m) d eta / d mu
    ok_wr = zero((m**2 / V(m)) * e_work**2 - (x_ - m) ** 2 / V(m))
    ok_bb = ok_b_cells and ok_wr
    out.append(
        _rec(
            "V28b",
            r"Cook's distance of a cell: $\mathrm{Cook}_{gc}=(r^{\mathrm{SP}}_{gc})^2h_{gc}/\{J(1-h_{gc})\}"
            r"=r_{gc}^2h_{gc}/\{J\hat\phi_g(1-h_{gc})^2\}$",
            "exact",
            "pass" if ok_bb else "fail",
            r"Exact one-step deletion refit at fixed IRLS weights (symbolic weights and working responses, "
            r"cells of a 3-cell and a 2-cell type): $(\bm b-\bm b_{(c)})^\top\bm X^\top\bm W\bm X(\bm b-\bm b_{(c)})/(J\phi)$ "
            r"equals both printed forms; $w\,e^2=r^2$ for the log link",
            "",
            detail_full=f"deletion {ok_b_cells}; working vs Pearson residual {ok_wr}",
        )
    )

    # ---- V28c ---------------------------------------------------------------
    G, Jg = 4, 2
    Mu = sp.Matrix(G, Jg, lambda i, j: sp.Symbol(f"mu{i + 1}{j + 1}", positive=True))
    tau = sp.symbols(f"tau1:{G + 1}", positive=True)
    Wg = sp.diag(*[1 / t**2 for t in tau])
    Ag = Mu.T * Wg * Mu
    Ainv = Ag.inv(method="ADJ").applyfunc(sp.cancel)
    Hg = (Mu * Ainv * Mu.T * Wg).applyfunc(sp.cancel)  # similar to H, see V28a
    ok_hdiag = all(
        zero(Hg[g, g] - (Mu[g, :] * Ainv * Mu[g, :].T)[0, 0] / tau[g] ** 2) for g in range(G)
    )
    ok_htr = zero(Hg.trace() - Jg)
    ok_hidem = all(zero(e) for e in (Hg * Hg - Hg))
    cst = sp.Symbol("c", positive=True)
    sub = {tau[0]: cst * tau[0], **{Mu[0, j]: cst * Mu[0, j] for j in range(Jg)}}
    ok_scale = all(zero(Hg[g, g].subs(sub, simultaneous=True) - Hg[g, g]) for g in range(G))
    # rational G = 7, J = 3: idempotence, trace, Cook against exact refit
    rng = random.Random(20261101)
    G2, J2 = 7, 3
    Xr = sp.Matrix(G2, J2, lambda i, j: sp.Rational(rng.randint(1, 40), rng.randint(1, 9)))
    t2 = [sp.Rational(rng.randint(1, 20), rng.randint(1, 5)) for _ in range(G2)]
    y = sp.Matrix([sp.Rational(rng.randint(1, 200), rng.randint(1, 9)) for _ in range(G2)])
    w2 = [1 / t**2 for t in t2]
    b2, A2 = wls(Xr, w2, y)
    H2 = Xr * A2.inv() * Xr.T * sp.diag(*w2)  # similar to H, see V28a
    ok_r = H2.trace() == J2 and (H2 * H2 - H2).applyfunc(sp.simplify) == sp.zeros(G2, G2)
    e2 = y - Xr * b2
    s2 = sum(w2[g] * e2[g] ** 2 for g in range(G2)) / (G2 - J2)
    ok_cook_g = True
    for g in range(G2):
        keep = [k for k in range(G2) if k != g]
        bd, _ = wls(Xr.extract(keep, list(range(J2))), [w2[k] for k in keep], y.extract(keep, [0]))
        d = b2 - bd
        cook_def = (d.T * A2 * d)[0, 0] / (J2 * s2)
        h = H2[g, g]
        rstd2 = w2[g] * e2[g] ** 2 / (s2 * (1 - h))
        if sp.simplify(cook_def - rstd2 / J2 * h / (1 - h)) != 0:
            ok_cook_g = False
    ok_c = ok_hdiag and ok_htr and ok_hidem and ok_scale and ok_r and ok_cook_g
    out.append(
        _rec(
            "V28c",
            r"Weighted bulk-on-reference regression: $h_g=\hat{\tilde{\bm\mu}}_{g\cdot}^\top(\bm X^\top\bm W\bm X)^{-1}"
            r"\hat{\tilde{\bm\mu}}_{g\cdot}/\tau_g^2$, $\sum_gh_g=J$, $\bm H^2=\bm H$, $h_g$ unchanged when "
            r"$\tau_g$ scales with the level; $\mathrm{Cook}_g=r_g^2h_g/\{J(1-h_g)\}$",
            "exact",
            "pass" if ok_c else "fail",
            r"$\bm H=\bm W^{1/2}\bm X(\bm X^\top\bm W\bm X)^{-1}\bm X^\top\bm W^{1/2}$ with symbolic $G=4$, $J=2$: "
            r"diagonal as printed, trace 2, idempotent, scale-invariant; rational $G=7$, $J=3$ (seed 20261101): "
            r"trace 3, idempotent, Cook equals the exact deletion distance for every gene",
            "",
            detail_full=(
                f"diag {ok_hdiag}; trace {ok_htr}; idempotent {ok_hidem}; scale {ok_scale}; rational {ok_r}; "
                f"cook {ok_cook_g}"
            ),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

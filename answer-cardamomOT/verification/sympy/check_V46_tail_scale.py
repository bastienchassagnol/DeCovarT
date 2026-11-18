# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V46a-V46g: effective cell count and multivariate discrepancy of vignette 03-1
(sec-fs-effective-cells, sec-fs-multivariate-discrepancy).

V46a (eq-fs-neff): with Var D = v / n, z^2 is linear in n, so the median
equation has the root n_eff = a q / median z^2(a).
V46b (eq-fs-n-from-scales): c = a n and phi = a^2 n give n = c^2 / phi.
V46c, V46d (eq-fs-mv-gene): Var(W r) = W when Var r = W^{-1}; the conditional
residual and variance through the Schur complement.
V46e: W_gg Sigma_gg >= 1 (power gain for a one-gene shift).
V46f: exact variance of the module sum of conditional statistics.
V46g: shared shift on two correlated genes; per-gene loss, equal module sums.
Same derivations as verification/sympy/check_tail_scale.py of the
reproducibility book (TS1-TS7).
"""

from __future__ import annotations

import random

import sympy as sp

SEED = 20261111


def _rec(id_, claim, kind, status, detail, tex_line=""):
    return {
        "id": id_,
        "claim": claim,
        "kind": kind,
        "status": status,
        "detail": detail,
        "tex_line": tex_line,
    }


def is_zero(expr):
    expr = sp.sympify(expr)
    if isinstance(expr, sp.MatrixBase):
        return all(is_zero(e) for e in expr)
    if sp.simplify(expr) == 0:
        return True
    return sp.cancel(sp.together(sp.expand(expr))) == 0


def spd_from(rows):
    a = sp.Matrix(rows)
    return a.T * a + sp.eye(a.shape[1])


def generic_symmetric(n, name):
    s = {}
    for i in range(n):
        for j in range(i, n):
            s[(i, j)] = s[(j, i)] = sp.Symbol(f"{name}{i + 1}{j + 1}", real=True)
    return sp.Matrix(n, n, lambda i, j: s[(i, j)])


def run():
    records = []
    # V46 (TS1). With no reference term, sd^2 D = v / n, so z^2(n) = z^2(a) n / a and
    # the median equation median z^2(n) = q has the root n = a q / median z^2(a).
    d, v, n, a, q, m_a = sp.symbols("d v n a q m_a", positive=True)
    z2 = lambda cells: d**2 / (v / cells)  # noqa: E731
    ok_scale = is_zero(z2(n) - z2(a) * n / a)
    root = sp.solve(sp.Eq(m_a * n / a, q), n)
    ok_root = len(root) == 1 and is_zero(root[0] - a * q / m_a)
    records.append(
        _rec(
            "V46a",
            "With Var D = v / n, z^2 scales linearly in n, so n_eff = a * qchisq(0.5, 1) / median z^2(a)",
            "exact",
            "pass" if ok_scale and ok_root else "fail",
            "z^2(n) - z^2(a) n / a simplified to 0; the median scales with every z^2, root from sympy.solve",
        )
    )

    # V46 (TS2). c = a n and phi = a^2 n give n = c^2 / phi.
    n_cells, scale = sp.symbols("n_cells a_platform", positive=True)
    c, phi = scale * n_cells, scale**2 * n_cells
    records.append(
        _rec(
            "V46b",
            "Under c = a n and phi = a^2 n, the cell count is c^2 / phi",
            "exact",
            "pass" if is_zero(c**2 / phi - n_cells) else "fail",
            "direct substitution",
        )
    )

    # V46 (TS3). If Var r = W^{-1}, then Var(W r) = W, so z_g = (W r)_g / sqrt(W_gg)
    # has variance 1 (conditional on the Gaussian model for the law N(0, 1)).
    w = generic_symmetric(3, "w")
    var_wr = w * w.inv() * w.T
    ok3 = is_zero(var_wr - w) and all(is_zero(var_wr[g, g] / w[g, g] - 1) for g in range(3))
    records.append(
        _rec(
            "V46c",
            "Var r = W^{-1} implies Var(W r) = W and Var z_g = 1 for z_g = (W r)_g / sqrt(W_gg)",
            "exact",
            "pass" if ok3 else "fail",
            "generic symmetric 3 x 3 W; Gaussianity of z_g follows from linearity, not checked here",
        )
    )

    # V46 (TS4). Conditional residual and variance: with W = Sigma^{-1},
    # r_1 - Sigma_{1,-1} Sigma_{-1}^{-1} r_{-1} = (W r)_1 / W_11 and
    # Sigma_11 - Sigma_{1,-1} Sigma_{-1}^{-1} Sigma_{-1,1} = 1 / W_11.
    s = generic_symmetric(3, "s")
    r = sp.Matrix(sp.symbols("r1 r2 r3", real=True))
    w = s.inv()
    s_12 = s[0, 1:]
    s_22 = s[1:, 1:]
    cond_resid = r[0] - (s_12 * s_22.inv() * r[1:, :])[0]
    cond_var = s[0, 0] - (s_12 * s_22.inv() * s_12.T)[0]
    ok4 = is_zero(cond_resid - (w * r)[0] / w[0, 0]) and is_zero(cond_var - 1 / w[0, 0])
    records.append(
        _rec(
            "V46d",
            "(W r)_g / W_gg = y_g - E(y_g | y_-g) and 1 / W_gg = Var(y_g | y_-g) (Schur complement)",
            "exact",
            "pass" if ok4 else "fail",
            "generic symmetric 3 x 3 Sigma, gene 1 against genes 2 and 3",
        )
    )

    # V46 (TS5). W_gg Sigma_gg >= 1, equal iff Sigma_{g,-g} = 0; a shift of k
    # marginal SDs in gene g moves z_cond by k sqrt(Sigma_gg W_gg) >= k.
    s11, s22 = sp.symbols("s11 s22", positive=True)
    s12 = sp.Symbol("s12", real=True)
    s2 = sp.Matrix([[s11, s12], [s12, s22]])
    w2 = s2.inv()
    gap = sp.simplify(w2[0, 0] - 1 / s11)
    ok_closed = is_zero(gap - s12**2 / (s11 * (s11 * s22 - s12**2)))
    k = sp.Symbol("k", positive=True)
    delta = sp.Matrix([k * sp.sqrt(s11), 0])
    mean_cond = (w2 * delta)[0] / sp.sqrt(w2[0, 0])
    ok_mean = is_zero(mean_cond**2 - k**2 * s11 * w2[0, 0])
    rng = random.Random(SEED)
    ok_draws, ok_equal = True, True
    for _ in range(25):
        sig = spd_from([[sp.Rational(rng.randint(-9, 9), rng.randint(1, 9)) for _ in range(5)] for _ in range(5)])
        wi = sig.inv()
        ok_draws &= all(wi[g, g] * sig[g, g] >= 1 for g in range(5))
        blk = sig.copy()
        blk[0, 1:] = sp.zeros(1, 4)
        blk[1:, 0] = sp.zeros(4, 1)
        ok_equal &= blk.inv()[0, 0] * blk[0, 0] == 1
    records.append(
        _rec(
            "V46e",
            "W_gg >= 1 / Sigma_gg (equality iff Sigma_{g,-g} = 0); a one-gene shift of k SDs moves z_cond by k sqrt(Sigma_gg W_gg) >= k",
            "exact",
            "pass" if ok_closed and ok_mean and ok_draws and ok_equal else "fail",
            "2 x 2: W_11 - 1/Sigma_11 = Sigma_12^2 / (Sigma_11 det Sigma) >= 0 (exact); "
            "5 x 5: 25 random rational SPD matrices and their block-diagonal versions (rational-sampling)",
        )
    )

    # V46 (TS6). Var(sum_{g in M} z_g) = 1' C_M 1, C = D^{-1/2} W D^{-1/2}.
    w = generic_symmetric(3, "w")
    dinv = sp.diag(*[1 / sp.sqrt(w[g, g]) for g in range(3)])
    sel = sp.Matrix([1, 1, 0])
    var_sum = (sel.T * dinv * w * dinv * sel)[0]
    c_corr = dinv * w * dinv
    records.append(
        _rec(
            "V46f",
            "Var(sum_{g in M} z_g) = 1' C_M 1 with C the correlation matrix of W",
            "exact",
            "pass" if is_zero(var_sum - sum(c_corr[i, j] for i in range(2) for j in range(2))) else "fail",
            "generic symmetric 3 x 3 W, module {1, 2}; uses Var(W r) = W (V46c)",
        )
    )

    # V46 (TS7). A shift of k SDs shared by two genes of correlation rho (unit
    # variances): each z_cond moves by k sqrt((1 - rho) / (1 + rho)) <= k, the
    # module sums move by the same sqrt(2) k / sqrt(1 + rho) under both tests.
    rho = sp.Symbol("rho", positive=True)
    s2 = sp.Matrix([[1, rho], [rho, 1]])
    w2 = s2.inv()
    delta = sp.Matrix([k, k])
    mean_g = (w2 * delta)[0] / sp.sqrt(w2[0, 0])
    c_w = sp.diag(*[1 / sp.sqrt(w2[g, g]) for g in range(2)]) * w2 * sp.diag(*[1 / sp.sqrt(w2[g, g]) for g in range(2)])
    sum_cond = 2 * mean_g / sp.sqrt(sum(c_w))
    sum_diag = 2 * k / sp.sqrt(sum(s2))
    ok7 = (
        is_zero(sp.simplify(mean_g.subs(rho, sp.Rational(1, 2))) - k * sp.sqrt(sp.Rational(1, 3)))
        and is_zero((mean_g**2 - k**2 * (1 - rho) / (1 + rho)).subs(rho, sp.Rational(3, 10)))
        and is_zero(sp.simplify(sum_cond**2 - 2 * k**2 / (1 + rho)).subs(rho, sp.Rational(3, 10)))
        and is_zero(sp.simplify(sum_diag**2 - 2 * k**2 / (1 + rho)))
    )
    rng = random.Random(SEED + 7)
    for _ in range(25):
        rr = sp.Rational(rng.randint(1, 19), 20)
        ok7 &= is_zero((mean_g**2 - k**2 * (1 - rho) / (1 + rho)).subs(rho, rr))
        ok7 &= is_zero((sum_cond**2 - 2 * k**2 / (1 + rho)).subs(rho, rr))
    records.append(
        _rec(
            "V46g",
            "Shared shift of k SDs on two genes of correlation rho > 0: per-gene z_cond moves by "
            "k sqrt((1 - rho)/(1 + rho)) < k; both module sums move by sqrt(2) k / sqrt(1 + rho)",
            "rational-sampling",
            "pass" if ok7 else "fail",
            "diagonal sum exact in rho; conditional forms at 25 rational rho in (0, 1) plus rho = 1/2 and 3/10",
        )
    )
    return records


if __name__ == "__main__":
    for r in run():
        print(r["id"], r["status"])

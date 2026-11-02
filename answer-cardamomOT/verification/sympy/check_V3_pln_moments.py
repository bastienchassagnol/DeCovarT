# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V3: Poisson log-normal moments on the linear scale, for a fixed offset s.
x_g | z ~ Poisson(s e^{z_g}) independently, z ~ N(mu, Sigma).

Derivation route:
- Poisson conditional moments E(X | lam) = lam and E{X(X-1) | lam} = lam^2
  are summed from the probability mass function;
- the normal moment generating function E e^{tZ} = exp(t m + t^2 v / 2) is
  integrated from the density, then applied to the scalar t'z ~ N(t'mu, t'Sigma t);
- the first two moments of x follow by the tower rule and are compared with
  the note's E x = s exp(mu + diag(Sigma)/2) and
  Cov x = diag(E x) + (E x E x^T) o (e^Sigma - 1).
A Monte Carlo check (numpy, fixed seed) confirms the result numerically.
"""

from __future__ import annotations

import json

import numpy as np
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


def _poisson_moments():
    k = sp.Symbol("k", integer=True, nonnegative=True)
    lam = sp.Symbol("lam", positive=True)
    pmf = sp.exp(-lam) * lam**k / sp.factorial(k)
    m1 = sp.simplify(sp.summation(k * pmf, (k, 0, sp.oo)))
    f2 = sp.simplify(sp.summation(k * (k - 1) * pmf, (k, 0, sp.oo)))
    return lam, m1, f2


def _normal_mgf():
    zz, t, m = sp.symbols("zz t m", real=True)
    v = sp.Symbol("v", positive=True)
    dens = sp.exp(-((zz - m) ** 2) / (2 * v)) / sp.sqrt(2 * sp.pi * v)
    mgf = sp.simplify(sp.integrate(sp.exp(t * zz) * dens, (zz, -sp.oo, sp.oo), conds="none"))
    return (t, m, v), mgf


def run() -> list[dict]:
    out: list[dict] = []
    lam, m1, f2 = _poisson_moments()
    (t_, m_, v_), mgf = _normal_mgf()

    G = 2
    s = sp.Symbol("s", positive=True)
    mu = sp.Matrix(sp.symbols("mu1 mu2", real=True))
    sg11, sg22 = sp.symbols("S11 S22", positive=True)
    sg12 = sp.Symbol("S12", real=True)
    Sig = sp.Matrix([[sg11, sg12], [sg12, sg22]])

    def M(tvec: sp.Matrix):
        # E exp(t'z) with t'z ~ N(t'mu, t'Sigma t)
        return mgf.subs({t_: 1, m_: (tvec.T * mu)[0], v_: (tvec.T * Sig * tvec)[0]})

    e = [sp.Matrix([1 if i == g else 0 for i in range(G)]) for g in range(G)]
    # tower rule with lam_g = s exp(z_g); the pmf sums gave E X = lam and
    # E X(X-1) = lam^2, which the lines below use
    assert m1 == lam and f2 == lam**2
    Ex = [s * M(e[g]) for g in range(G)]
    Exx = sp.zeros(G, G)
    for g in range(G):
        for h in range(G):
            if g == h:
                # E X^2 = E{X(X-1)} + E X = E lam^2 + E lam
                Exx[g, h] = s**2 * M(2 * e[g]) + s * M(e[g])
            else:
                # conditional independence: E X_g X_h = E lam_g lam_h
                Exx[g, h] = s**2 * M(e[g] + e[h])
    Ex_v = sp.Matrix(Ex)
    Cov = Exx - Ex_v * Ex_v.T

    note_mean = sp.Matrix([s * sp.exp(mu[g] + Sig[g, g] / 2) for g in range(G)])
    hadamard = sp.Matrix(G, G, lambda g, h: note_mean[g] * note_mean[h] * (sp.exp(Sig[g, h]) - 1))
    note_cov = sp.diag(*note_mean) + hadamard
    ok_mean = all(sp.simplify(a - b) == 0 for a, b in zip(Ex_v, note_mean))
    ok_cov = all(sp.simplify(sp.expand(a - b)) == 0 for a, b in zip(Cov, note_cov))
    out.append(
        _rec(
            "V3a",
            r"Poisson log-normal moments for fixed $s$: $\mathbb E\bm x=s\,e^{\bm\mu+\operatorname{diag}(\bm\Sigma)/2}$, $\operatorname{Cov}\bm x=\operatorname{diag}(\mathbb E\bm x)+(\mathbb E\bm x\,\mathbb E\bm x^{\mathsf T})\circ(e^{\bm\Sigma}-1)$, exponential entrywise",
            "exact",
            "pass" if ok_mean and ok_cov else "fail",
            r"Poisson moments summed from the pmf ($\mathbb E X=" + sp.latex(m1, symbol_names={lam: r"\lambda"}) + r"$, $\mathbb E X(X-1)=" + sp.latex(f2, symbol_names={lam: r"\lambda"}) + r"$); normal MGF $" + sp.latex(mgf) + r"$ integrated from the density; tower rule, $G=2$",
            323,
            detail_full=(
                f"Poisson moments summed from the pmf (E X = {m1}, E X(X-1) = {f2}); normal MGF "
                f"integrated from the density: {sp.sstr(mgf)}; tower rule gives mean and "
                "covariance equal to the note's forms (G = 2, symbolic mu, Sigma)"
            ),
        )
    )

    # matrix-exponential reading of e^{Sigma} is not the result
    Sig_num = sp.Matrix([[sp.Rational(1, 2), sp.Rational(1, 5)], [sp.Rational(1, 5), sp.Rational(3, 10)]])
    entry = Sig_num.applyfunc(sp.exp)
    matexp = Sig_num.exp()
    diff = (entry - matexp).applyfunc(lambda q: sp.N(q, 8))
    differs = any(abs(q) > 1e-6 for q in diff)
    out.append(
        _rec(
            "V3b",
            r"$e^{\bm\Sigma_z}$ in the covariance is the entrywise exponential, not the matrix exponential",
            "numeric",
            "flag" if differs else "pass",
            r"for $\bm\Sigma=\left(\begin{smallmatrix}0.5&0.2\\0.2&0.3\end{smallmatrix}\right)$, entrywise minus matrix exponential $=" + sp.latex(diff.applyfunc(lambda q: sp.Float(round(float(q), 3), 3))) + r"$; the derivation uses $e^{\Sigma_{gh}}$ entry by entry",
            323,
            detail_full=(
                "for Sigma = [[0.5, 0.2], [0.2, 0.3]], entrywise exp minus matrix exp = "
                f"{[[round(float(diff[i, j]), 4) for j in range(2)] for i in range(2)]}; the derived "
                "covariance uses exp(Sigma_gh) entry by entry. The note's text now says "
                "entrywise"
            ),
        )
    )

    # Monte Carlo
    rng = np.random.default_rng(20261030)
    n = 2_000_000
    s_num = 0.8
    mu_num = np.array([0.5, 1.0, -0.2])
    Sig_mc = np.array([[0.40, 0.15, -0.10], [0.15, 0.30, 0.05], [-0.10, 0.05, 0.25]])
    z = rng.multivariate_normal(mu_num, Sig_mc, size=n)
    x = rng.poisson(s_num * np.exp(z))
    m_th = s_num * np.exp(mu_num + np.diag(Sig_mc) / 2)
    C_th = np.diag(m_th) + np.outer(m_th, m_th) * (np.exp(Sig_mc) - 1)
    xbar = x.mean(axis=0)
    xc = x - xbar
    zs = []
    for g in range(3):
        se = x[:, g].std(ddof=1) / np.sqrt(n)
        zs.append(abs(xbar[g] - m_th[g]) / se)
    for g in range(3):
        for h in range(g, 3):
            prod = xc[:, g] * xc[:, h]
            se = prod.std(ddof=1) / np.sqrt(n)
            zs.append(abs(prod.mean() - C_th[g, h]) / se)
    zmax = float(max(zs))
    tol = 4.5
    out.append(
        _rec(
            "V3c",
            r"Monte Carlo confirmation of the Poisson log-normal moments ($G=3$, fixed $s$)",
            "numeric",
            "pass" if zmax < tol else "fail",
            "NumPy seed 20261030, $n=" + f"{n:,}".replace(",", "{,}") + f"$, $s={s_num}$; largest $|z|$ over 3 means and 6 covariances $={zmax:.2f}$ (tolerance {tol})",
            323,
            detail_full=(
                f"numpy default_rng seed 20261030, n = {n:,}, s = {s_num}; 3 means and 6 "
                f"covariances; largest |MC - formula| / MC standard error = {zmax:.2f} "
                f"(tolerance {tol})"
            ),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V35: population-level dynamics and potency (article app:dynamic-potency).

Derivation route:
- V35a: SCENT, Wang et al. 2026 and Teschendorff & Feinberg 2021 Box 2. For the
  walk P(g->h) = A_gh x_h / (A x)_g on a PPI graph weighted by the expression x
  of one cell, the law Pi_g = x_g (A x)_g / (x' A x) is stationary (Pi P = Pi),
  checked on a 4-gene graph with symbolic positive x.
- V35b: the largest entropy rate of a walk on a connected graph is
  log lambda_max(A) (Demetrius & Manke): the walk P*(g->h) = A_gh v_h /
  (lambda v_g), v the Perron vector, has stationary law v_g^2 / sum v^2 and
  entropy rate log lambda; SCENT entropy rates for random expression vectors
  stay below it. mpmath at 40 digits on two fixed graphs.
- V35c: Gaussian differential entropy -E log f = 1/2 log det(2 pi e Sigma),
  for G = 2, by the change of variables x = mu + L z (Sigma = L L') and an
  exact integral of z'z under the standard normal.
- V35d: Teschendorff & Feinberg Box 1 with net birth rate zero: in one
  dimension nu proportional to exp(-F/D) makes the right-hand side
  d/dzeta (nu F') + D nu'' vanish for any potential F.
- V35e: RNA velocity as the mean of the nascent-mature reaction network
  (constant transcription k_tx, splicing k_spl, degradation k_deg): the
  generator applied to f(N, M) = M gives k_spl N - k_deg M, and to f = N gives
  k_tx - k_spl N.
"""

from __future__ import annotations

import json
import random

import mpmath as mp  # ships with sympy
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


TEX = "article app:dynamic-potency"

# two connected undirected graphs without self-loops (SCENT convention A_gg = 0)
GRAPH_4 = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
GRAPH_6 = [(0, 1), (0, 2), (0, 3), (1, 2), (2, 4), (3, 4), (4, 5), (1, 5)]


def _adjacency(n, edges):
    A = sp.zeros(n, n)
    for g, h in edges:
        A[g, h] = 1
        A[h, g] = 1
    return A


def run():
    out = []

    # V35a: stationarity of the SCENT law
    A = _adjacency(4, GRAPH_4)
    x = sp.Matrix(sp.symbols("x0:4", positive=True))
    Ax = A * x
    P = sp.Matrix(4, 4, lambda g, h: A[g, h] * x[h] / Ax[g])
    Pi = sp.Matrix(1, 4, lambda _, g: x[g] * Ax[g] / (x.T * A * x)[0, 0])
    resid = sp.simplify(Pi * P - Pi)
    rows = sp.simplify(P * sp.ones(4, 1) - sp.ones(4, 1))
    tot = sp.simplify(sum(Pi) - 1)
    ok_a = resid == sp.zeros(1, 4) and rows == sp.zeros(4, 1) and tot == 0
    out.append(
        _rec(
            "V35a",
            r"SCENT: $\Pi_{gc}=x_{gc}(\boldsymbol{A}\boldsymbol{x}_{\cdot c})_{g}/"
            r"\boldsymbol{x}_{\cdot c}^{\top}\boldsymbol{A}\boldsymbol{x}_{\cdot c}$ is stationary for "
            r"$P_{c}(g\to h)=A_{gh}x_{hc}/(\boldsymbol{A}\boldsymbol{x}_{\cdot c})_{g}$",
            "exact",
            "pass" if ok_a else "fail",
            r"$\Pi P-\Pi=" + sp.latex(resid) + r"$ on a 4-gene graph, rows of $P$ sum to 1, $\sum\Pi=1$",
            TEX,
            detail_full=f"resid={resid}; rows={rows.T}; total={tot}",
        )
    )

    # V35b: maximum entropy rate log lambda_max, and SCENT below it
    mp.mp.dps = 40
    rng = random.Random(20261109)
    min_gap = mp.inf
    max_dev = mp.mpf(0)
    for n, edges in ((4, GRAPH_4), (6, GRAPH_6)):
        An = mp.matrix(n, n)
        for g, h in edges:
            An[g, h] = 1
            An[h, g] = 1
        E, ER = mp.eigsy(An)
        k = max(range(n), key=lambda i: E[i])
        lam = E[k]
        v = [abs(ER[i, k]) for i in range(n)]
        s2 = sum(vi**2 for vi in v)
        rate_star = mp.mpf(0)
        for g in range(n):
            pig = v[g] ** 2 / s2
            for h in range(n):
                if An[g, h]:
                    p = v[h] / (lam * v[g])
                    rate_star -= pig * p * mp.log(p)
        max_dev = max(max_dev, abs(rate_star - mp.log(lam)))
        for _ in range(25):
            xv = [mp.mpf(rng.uniform(0.05, 5.0)) for _ in range(n)]
            Axv = [sum(An[g, h] * xv[h] for h in range(n)) for g in range(n)]
            norm = sum(xv[g] * Axv[g] for g in range(n))
            rate = mp.mpf(0)
            for g in range(n):
                pig = xv[g] * Axv[g] / norm
                for h in range(n):
                    if An[g, h]:
                        p = xv[h] / Axv[g]
                        rate -= pig * p * mp.log(p)
            min_gap = min(min_gap, mp.log(lam) - rate)
    ok_b = max_dev < mp.mpf(10) ** -30 and min_gap >= -mp.mpf(10) ** -30
    out.append(
        _rec(
            "V35b",
            r"largest entropy rate of a walk on the graph is $\log\lambda_{\max}(\boldsymbol{A})$, "
            r"and the SCENT rate stays below it",
            "numeric",
            "pass" if ok_b else "fail",
            "mpmath 40 digits, graphs of 4 and 6 genes: |rate(P*) - log lambda| = "
            + mp.nstr(max_dev, 3)
            + "; 50 random expression vectors (seed 20261109), none above log lambda (smallest gap "
            + mp.nstr(min_gap, 3)
            + ", tolerance 1e-30)",
            TEX,
        )
    )

    # V35c: Gaussian entropy for G = 2
    l11, l22 = sp.symbols("l11 l22", positive=True)
    l21 = sp.Symbol("l21", real=True)
    L = sp.Matrix([[l11, 0], [l21, l22]])
    Sig = L * L.T
    z1, z2 = sp.symbols("z1 z2", real=True)
    phi = sp.exp(-(z1**2 + z2**2) / 2) / (2 * sp.pi)
    # -log f(mu + L z) = 1/2 log det(2 pi Sigma) + z'z/2
    Ezz = sp.integrate(sp.integrate((z1**2 + z2**2) * phi, (z1, -sp.oo, sp.oo)), (z2, -sp.oo, sp.oo))
    H = sp.log((2 * sp.pi) ** 2 * Sig.det()) / 2 + Ezz / 2
    target = sp.log((2 * sp.pi * sp.E) ** 2 * Sig.det()) / 2
    diff_c = sp.simplify(sp.expand_log(H - target, force=True))
    ok_c = diff_c == 0 and sp.simplify(Ezz - 2) == 0
    out.append(
        _rec(
            "V35c",
            r"Gaussian entropy $-\mathbb{E}\log f=\tfrac12\log\det(2\pi e\boldsymbol{\Sigma})$ ($G=2$)",
            "exact",
            "pass" if ok_c else "fail",
            r"$\mathbb{E}\,\boldsymbol{z}^{\top}\boldsymbol{z}=" + sp.latex(Ezz) + r"$, difference $" + sp.latex(diff_c) + r"$",
            TEX,
        )
    )

    # V35d: Fokker-Planck stationary law at zero net birth rate (1-D)
    zeta = sp.Symbol("zeta", real=True)
    D = sp.Symbol("D", positive=True)
    F = sp.Function("F")(zeta)
    nu = sp.exp(-F / D)
    rhs = sp.diff(nu * sp.diff(F, zeta), zeta) + D * sp.diff(nu, zeta, 2)
    flux = sp.simplify(nu * sp.diff(F, zeta) + D * sp.diff(nu, zeta))
    ok_d = sp.simplify(rhs) == 0 and flux == 0
    out.append(
        _rec(
            "V35d",
            r"with $k^{\mathrm{net}}=0$, $\nu\propto e^{-\mathcal{F}/D_{\mathrm{diff}}}$ is stationary for "
            r"$\partial_{t}\nu=\nabla\cdot(\nu\nabla\mathcal{F})+D_{\mathrm{diff}}\nabla^{2}\nu$ (1-D)",
            "exact",
            "pass" if ok_d else "fail",
            r"flux $" + sp.latex(flux) + r"$, right-hand side $" + sp.latex(sp.simplify(rhs)) + r"$",
            TEX,
        )
    )

    # V35e: generator of the nascent-mature network
    N, M = sp.symbols("N M", nonnegative=True, integer=True)
    ktx, kspl, kdeg = sp.symbols("k_tx k_spl k_deg", positive=True)

    def gen(f):
        return sp.expand(
            ktx * (f(N + 1, M) - f(N, M))
            + kspl * N * (f(N - 1, M + 1) - f(N, M))
            + kdeg * M * (f(N, M - 1) - f(N, M))
        )

    gM = gen(lambda n, m: m)
    gN = gen(lambda n, m: n)
    ok_e = sp.simplify(gM - (kspl * N - kdeg * M)) == 0 and sp.simplify(gN - (ktx - kspl * N)) == 0
    out.append(
        _rec(
            "V35e",
            r"velocity: $\mathrm{d}\langle x^{\mathrm{M}}\rangle/\mathrm{d}t=k^{\mathrm{spl}}\langle x^{\mathrm{N}}\rangle"
            r"-k^{\mathrm{deg}}\langle x^{\mathrm{M}}\rangle$ from the nascent--mature generator",
            "exact",
            "pass" if ok_e else "fail",
            r"$\mathcal{L}M=" + sp.latex(gM) + r"$, $\mathcal{L}N=" + sp.latex(gN) + r"$",
            TEX,
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

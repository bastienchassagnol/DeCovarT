# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V12: the quantile universal threshold (QUT) of an l1 penalty read as a
multiple-testing threshold (CardamomOT penalty subsection).

Setting (our reading, stated in the note): for one target gene, the l1 KKT
condition at theta = 0 keeps every edge out if and only if lambda >=
max_j |s_j|, with s_j the score (gradient of the loss) of candidate regulator j
at the null model. The QUT at level 1 - alpha is the (1 - alpha) quantile of
max_j |s_j| under the null (Giacobino et al. 2017). For the analysis the m
null scores are taken as independent N(0, sigma^2) (orthogonal design).

Derivation routes:
- V12a (exact): union bound, P(max |s_j| > q) <= m P(|s_1| > q), so the QUT is
  never above the Bonferroni quantile q_B = sigma Phi^-1(1 - alpha/(2m)),
  whatever the dependence. Under independence the QUT is the Sidak quantile,
  per-test level a_S = 1 - (1 - alpha)^(1/m) >= alpha/m (Bernoulli: the
  derivative of (1 - alpha/m)^m - (1 - alpha) is >= 0 and it is 0 at
  alpha = 0, so 1 - alpha/m >= (1 - alpha)^(1/m)), and m a_S -> -log(1 - alpha) as m -> oo, so a_S / (alpha/m) ->
  -log(1 - alpha)/alpha, not 1.
- V12b (numeric): mpmath quantiles for m = 10 ... 10^6 and alpha = 0.05, 0.5:
  q_S / q_B -> 1 and q_B - q_S -> 0, so the two thresholds agree to first
  order, both ~ sqrt(2 log m).
- V12c (exact): k true regulators with scores far above the threshold and
  m - k independent null scores, threshold = QUT computed with all m
  candidates null: per-target FWER = 1 - (1 - alpha)^((m - k)/m) <= alpha,
  expected false edges E V = (m - k) a_S -> -log(1 - alpha) when k = sqrt(m)
  and m -> oo; the false-edge share E V / (E V + k) at m = 100, k = 10; the
  per-test level t = q k / {(1 - q)(m - k)} that gives a share q.
- V12d (exact): the paper's lambda = 0.05 N_c / sqrt(N_g) against the null-score
  scale sigma sqrt(N_c) of a sum of N_c independent per-cell terms: the ratio
  grows as sqrt(N_c).
"""

from __future__ import annotations

import json

import mpmath as mp
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


def run() -> list[dict]:
    out: list[dict] = []
    al = sp.Symbol("alpha", positive=True)
    m = sp.Symbol("m", positive=True, integer=True)

    # ---------------- V12a: QUT <= Bonferroni; Sidak per-test level ----------------
    f = (1 - al / m) ** m - (1 - al)
    df = sp.diff(f, al)
    deriv_ok = sp.simplify(sp.powsimp(sp.expand(df - (1 - (1 - al / m) ** (m - 1))), force=True)) == 0  # >= 0 for 0 <= alpha <= 1, m >= 1
    at0 = sp.simplify(f.subs(al, 0)) == 0
    a_S = 1 - (1 - al) ** (1 / m)
    lim_ma = sp.limit(m * a_S, m, sp.oo)
    lim_ok = sp.simplify(lim_ma + sp.log(1 - al)) == 0
    r05 = sp.N(-sp.log(1 - sp.Rational(1, 20)) / sp.Rational(1, 20), 4)
    r50 = sp.N(-sp.log(sp.Rational(1, 2)) / sp.Rational(1, 2), 4)
    # spot check of a_S > alpha/m at integer m >= 2
    spot = all(
        sp.N(1 - (1 - a) ** sp.Rational(1, mm) - a / mm, 30) > 0
        for a in (sp.Rational(1, 20), sp.Rational(1, 2))
        for mm in (2, 10, 1000)
    )  # m = 1: equality
    out.append(
        _rec(
            "V12a",
            r"QUT $\le$ Bonferroni quantile for any dependence; under independence the QUT is the Sidak "
            r"quantile, with per-test level $1-(1-\alpha)^{1/m}\ge\alpha/m$",
            "exact",
            "pass" if deriv_ok and at0 and lim_ok and spot else "fail",
            r"union bound; $\tfrac{\mathrm d}{\mathrm d\alpha}\{(1-\alpha/m)^m-(1-\alpha)\}="
            + sp.latex(1 - (1 - al / m) ** (m - 1))
            + r"\ge0$, zero at $\alpha=0$; $\lim_{m\to\infty}m\{1-(1-\alpha)^{1/m}\}="
            + sp.latex(lim_ma)
            + r"$, so the per-test levels differ by the factor $-\log(1-\alpha)/\alpha$: $"
            + str(r05)
            + r"$ at $\alpha=0.05$, $"
            + str(r50)
            + r"$ at $\alpha=0.5$",
            0,
            detail_full=(
                f"f'(alpha) identity: {deriv_ok}; f(0) = 0: {at0}; lim m a_S = {lim_ma}; spot checks: {spot}; "
                f"ratio limits {r05} (alpha=0.05), {r50} (alpha=0.5)"
            ),
        )
    )

    # ---------------- V12b: numeric quantiles ----------------
    mp.mp.dps = 50

    def qnorm_upper(p):  # Phi^-1(1 - p)
        return -mp.sqrt(2) * mp.erfinv(2 * p - 1)

    rows = []
    ok = True
    for a in (mp.mpf("0.05"), mp.mpf("0.5")):
        prev_ratio, prev_diff = None, None
        for mm in (10, 100, 1000, 10**4, 10**6):
            aS = 1 - (1 - a) ** (mp.mpf(1) / mm)
            qS = qnorm_upper(aS / 2)
            qB = qnorm_upper(a / (2 * mm))
            # QUT definition: P(max |z| <= qS) = 1 - a under independence
            ok &= abs((1 - mp.erfc(qS / mp.sqrt(2))) ** mm - (1 - a)) < mp.mpf(10) ** -40
            ratio, diff = qS / qB, qB - qS
            ok &= ratio <= 1 and diff >= 0
            if prev_ratio is not None:
                ok &= ratio > prev_ratio and diff < prev_diff
            prev_ratio, prev_diff = ratio, diff
            rows.append((float(a), mm, mp.nstr(qS, 6), mp.nstr(qB, 6), mp.nstr(ratio, 6),
                         mp.nstr(mp.sqrt(2 * mp.log(mm)), 6)))
    last = [r for r in rows if r[1] == 10**6]
    # first-order gap: q_B - q_S ~ log{-log(1 - alpha)/alpha} / q_B (ratio of the two per-test levels)
    for a in (mp.mpf("0.05"), mp.mpf("0.5")):
        aS = 1 - (1 - a) ** (mp.mpf(1) / 10**6)
        qS, qB = qnorm_upper(aS / 2), qnorm_upper(a / (2 * 10**6))
        pred = mp.log(-mp.log(1 - a) / a) / qB
        ok &= abs((qB - qS) / pred - 1) < mp.mpf("0.1")
    small = [r for r in rows if r[1] == 100]
    out.append(
        _rec(
            "V12b",
            r"For independent Gaussian scores the QUT and the Bonferroni quantile agree asymptotically: "
            r"ratio $\to1$, difference $\to0$",
            "numeric",
            "pass" if ok else "fail",
            r"mpmath, 50 digits; $q_{\mathrm{QUT}}/q_{\mathrm B}$ at $m=100$: $"
            + small[0][4] + r"$ ($\alpha=0.05$), $" + small[1][4] + r"$ ($\alpha=0.5$); at $m=10^6$: $"
            + last[0][4] + r"$, $" + last[1][4] + r"$; the ratio increases in $m$ and the gap "
            r"$q_{\mathrm B}-q_{\mathrm{QUT}}$ shrinks as $\log\{-\log(1-\alpha)/\alpha\}/q_{\mathrm B}$ "
            r"(within 10\,\% at $m=10^6$)",
            0,
            detail_full=f"(alpha, m, q_QUT, q_Bonferroni, ratio, sqrt(2 log m)): {rows}",
        )
    )

    # ---------------- V12c: sparsity, per-target FWER and false-edge share ----------------
    k = sp.Symbol("k", positive=True)
    b = sp.Symbol("b", positive=True)  # b = 1 - alpha in (0, 1)
    fwer_b = 1 - (b ** (1 / m)) ** (m - k)
    fwer_ok = sp.simplify(sp.powsimp(fwer_b - (1 - b ** ((m - k) / m)), force=True)) == 0
    aS = 1 - (1 - al) ** (1 / m)
    fwer = 1 - (1 - aS) ** (m - k)
    EV = (m - k) * aS
    lim_EV = sp.limit(EV.subs(k, sp.sqrt(m)), m, sp.oo)
    lim_ok2 = sp.simplify(lim_EV + sp.log(1 - al)) == 0
    # per-test level t giving a false-edge share q: (m - k) t / {(m - k) t + k} = q
    t, qq = sp.symbols("t q", positive=True)
    t_sol = sp.solve(sp.Eq((m - k) * t / ((m - k) * t + k), qq), t)
    t_ok = len(t_sol) == 1 and sp.simplify(t_sol[0] - qq * k / ((1 - qq) * (m - k))) == 0
    vals = {}
    for a in (sp.Rational(1, 20), sp.Rational(1, 2)):
        ev = sp.N(EV.subs({m: 100, k: 10, al: a}), 4)
        vals[str(a)] = (ev, sp.N(ev / (ev + 10), 3), sp.N(fwer.subs({m: 100, k: 10, al: a}), 3))
    out.append(
        _rec(
            "V12c",
            r"With $k$ true and $m-k$ null regulators, the QUT controls the per-target FWER at "
            r"$\le\alpha$; expected false edges $\to-\log(1-\alpha)$ per target",
            "exact",
            "pass" if fwer_ok and lim_ok2 and t_ok else "fail",
            r"FWER $=1-(1-\alpha)^{(m-k)/m}$; $\lim_{m\to\infty}\mathbb E V\,(k=\sqrt m)="
            + sp.latex(lim_EV)
            + r"$; $m=100$, $k=10$: $\mathbb E V=" + str(vals["1/2"][0]) + r"$, share $\mathbb E V/(\mathbb E V+k)="
            + str(vals["1/2"][1]) + r"$ at $\alpha=0.5$; $" + str(vals["1/20"][0]) + r"$ and $"
            + str(vals["1/20"][1]) + r"$ at $\alpha=0.05$; share $q$ needs per-test level $t="
            + sp.latex(qq * k / ((1 - qq) * (m - k))) + r"$ (solved exactly)",
            0,
            detail_full=f"FWER identity: {fwer_ok}; t for share q: {t_sol} ({t_ok}); lim E V = {lim_EV}; (E V, share, FWER) at m=100, k=10: {vals}",
        )
    )

    # ---------------- V12d: the paper's lambda against the null-score scale ----------------
    Nc, Ng, sig, q = sp.symbols("N_c N_g sigma q", positive=True)
    lam_paper = sp.Rational(1, 20) * Nc / sp.sqrt(Ng)
    ratio = sp.simplify(lam_paper / (sig * sp.sqrt(Nc) * q))
    grows = sp.limit(ratio, Nc, sp.oo) == sp.oo
    budget = sp.simplify(lam_paper * sp.sqrt(Ng) - sp.Rational(1, 20) * Nc) == 0
    out.append(
        _rec(
            "V12d",
            r"The paper's $\lambda=0.05N_c/\sqrt{N_g}$ grows as $N_c$, the null score as $\sqrt{N_c}$",
            "exact",
            "pass" if grows and budget else "fail",
            r"$\lambda/(\sigma\sqrt{N_c}\,q)=" + sp.latex(ratio) + r"\to\infty$ as $N_c\to\infty$; "
            r"$\lambda\sqrt{N_g}=0.05N_c$",
            0,
            detail_full=f"ratio = {ratio}; diverges: {grows}; lambda sqrt(N_g) = 0.05 N_c: {budget}",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

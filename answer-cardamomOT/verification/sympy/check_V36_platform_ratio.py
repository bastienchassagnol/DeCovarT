# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V36: bulk discrepancy with a plug-in gene-wise platform ratio
(vignette 03-1, #sec-fs-platform-ratio, @eq-fs-platform-mean,
@eq-fs-discrepancy-b, @eq-fs-var-b).

Setting: (A1)-(A4) of #sec-fs-five-bounds, with A3 relaxed to
E(y_g / a) = b_g sum_j p_j mu_gj: the bulk count of gene g is b_g times a
sum of a cells of the reference law. D^b = y / (a bhat) - sum_j p_j muhat_j.

Derivation route:
- V36a (exact): with bhat = b known, the CGF of y = b S (S the sum of the
  m_j = p_j a bulk cells) is K_y(s) = sum_j m_j K_j(b s); the CGF of
  y / (a b) is then K_y(u / (a b)) = sum_j m_j K_j(u / a), and adding the
  reference term gives the printed K_D of @eq-fs-cgf-d with n = a: the
  thresholds of T1 to T5 are unchanged. Checked with the NB per-cell CGF
  and symbolic J = 3.
- V36b (series): bhat = b exp(eps), eps independent of the bulk sample under
  test (cross-fitting), E eps = 0, Var eps = v. With X = y / (a b) and
  eps = s eta (bookkeeping s), E[eta] = 0, E[eta^2] = w, E[eta^3] = k3,
  E[eta^4] = k4, the moments of X e^{-eps} are expanded to order s^2 (so
  v = s^2 w): Var D^b = Var X + v (E X)^2 + 2 v Var X + O(s^3) and
  E D^b - E D = v E X / 2 + O(s^3). The vignette keeps the first-order term
  (E X)^2 v; the neglected 2 v Var X is relative order 2v. Gaussian eps:
  the exact log-normal moments agree with the expansion.
- V36c (exact): identifiability. (i) a common factor c of all b_g is
  absorbed by a: (c b_g)(a / c) sum p mu = b_g a sum p mu. (ii) one sample,
  b free: for any composition p', b'_g = b_g (sum_j p_j mu_gj) /
  (sum_j p'_j mu_gj) gives the same mean for every gene, so p is not
  identified. (iii) two samples i, k: the ratio of the means of gene g,
  a_i sum_j p_ji mu_gj / (a_k sum_j p_jk mu_gj), does not contain b_g, so
  compositions that differ between samples carry the information on p.
"""

from __future__ import annotations

import json

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
    if sp.cancel(sp.together(sp.expand(e))) == 0:
        return True
    return sp.simplify(e) == 0


psi = sp.Symbol("psi", positive=True)


def K_nb(arg, mu):
    """Per-cell NB CGF of @eq-fs-cgf-d."""
    return -psi * sp.log(1 - mu / psi * (sp.exp(arg) - 1))


def run() -> list[dict]:
    out = []
    J = 3
    u, s = sp.symbols("u s", real=True)
    mu = sp.symbols("mu1:4", positive=True)
    p = sp.symbols("p1:4", positive=True)
    Cs = sp.symbols("C1:4", positive=True)
    a, b, c = sp.symbols("a b c", positive=True)

    # ---- V36a ---------------------------------------------------------------
    # CGF of y = b * S, S = sum of m_j = p_j a cells: K_y(t) = sum_j m_j K_j(b t)
    K_y = lambda t: sum(p[j] * a * K_nb(b * t, mu[j]) for j in range(J))  # noqa: E731
    K_ref = sum(Cs[j] * K_nb(-u * p[j] / Cs[j], mu[j]) for j in range(J))
    K_Db = K_y(u / (a * b)) + K_ref
    K_D = sum(p[j] * a * K_nb(u / a, mu[j]) for j in range(J)) + K_ref
    ok_a = zero(K_Db - K_D)
    out.append(
        _rec(
            "V36a",
            r"With the ratio known, $D^{b}_{gi}$ has the cumulant generating function of $D_{gi}$ under A3, "
            r"so the thresholds of T1 to T5 are unchanged",
            "exact",
            "pass" if ok_a else "fail",
            r"$K_{y}(t)=\sum_j p_ja\,K_{gj}(\tilde b_g t)$ evaluated at $t=u/(a\tilde b_g)$, plus the reference "
            r"term, equals @eq-fs-cgf-d with $n=a$ (NB per-cell CGF, symbolic $J=3$)",
            "",
            detail_full=f"K_Db - K_D == 0: {ok_a}",
        )
    )

    # ---- V36b ---------------------------------------------------------------
    # X: positive random variable with E X = m, E X^2 = q (Var X = q - m^2),
    # independent of eta. E[(X e^{-s eta})^k] = E X^k E e^{-k s eta}.
    m, q, w, k3, k4 = sp.symbols("m q w k3 k4", real=True)
    eta_mom = {1: 0, 2: w, 3: k3, 4: k4}

    def mgf_eta(k):
        # E exp(-k s eta) to order s^3 (moments of eta up to 3 suffice)
        return 1 + sum(((-k * s) ** r) / sp.factorial(r) * eta_mom[r] for r in range(1, 4))

    mean_b = m * mgf_eta(1)
    second_b = q * mgf_eta(2)
    var_b = sp.series(sp.expand(second_b - mean_b**2), s, 0, 3).removeO()
    v = s**2 * w
    var_target = (q - m**2) + v * m**2 + 2 * v * (q - m**2)
    bias = sp.series(sp.expand(mean_b - m), s, 0, 3).removeO()
    ok_var = zero(var_b - var_target)
    ok_bias = zero(bias - v * m / 2)
    # Gaussian eps: exact log-normal moments e^{v/2}, e^{2v}
    vv = sp.Symbol("v", positive=True)
    var_ln = q * sp.exp(2 * vv) - m**2 * sp.exp(vv)
    ser_ln = sp.series(var_ln, vv, 0, 2).removeO()
    ok_ln = zero(ser_ln - ((q - m**2) + vv * m**2 + 2 * vv * (q - m**2)))
    ok_b = ok_var and ok_bias and ok_ln
    out.append(
        _rec(
            "V36b",
            r"With a cross-fitted plug-in $\hat b_g=\tilde b_g e^{\varepsilon_g}$, "
            r"$\operatorname{Var}D^{b}_{gi}\approx\operatorname{Var}D_{gi}+(\sum_jp_{ji}\tilde\mu_{gj})^{2}v_g$",
            "series",
            "pass" if ok_b else "fail",
            r"Expansion in $\varepsilon=s\eta$ to order $s^{2}$: "
            r"$\operatorname{Var}=\operatorname{Var}X+v\,(\operatorname{E}X)^{2}+2v\operatorname{Var}X+O(s^{3})$ "
            r"and bias $v\operatorname{E}X/2$; the term $2v\operatorname{Var}X$ is dropped (relative order $2v$); "
            r"Gaussian $\varepsilon$: log-normal moments agree to order $v$",
            "",
            detail_full=f"variance {ok_var}; bias {ok_bias}; log-normal {ok_ln}",
        )
    )

    # ---- V36c ---------------------------------------------------------------
    bg = sp.symbols("b1:4", positive=True)  # three genes
    mug = [sp.symbols(f"mu{g}_1:4", positive=True) for g in range(1, 4)]
    pp = sp.symbols("q1:4", positive=True)  # an arbitrary other composition
    pk = sp.symbols("r1:4", positive=True)  # composition of a second sample
    ak = sp.Symbol("a_k", positive=True)
    fit = lambda g, comp: sum(comp[j] * mug[g][j] for j in range(J))  # noqa: E731
    ok_scale = all(zero((c * bg[g]) * (a / c) * fit(g, p) - bg[g] * a * fit(g, p)) for g in range(3))
    ok_single = all(
        zero(bg[g] * fit(g, p) / fit(g, pp) * a * fit(g, pp) - bg[g] * a * fit(g, p)) for g in range(3)
    )
    ratio = [bg[g] * a * fit(g, p) / (bg[g] * ak * fit(g, pk)) for g in range(3)]
    ok_pair = all(not ratio[g].has(bg[g]) for g in range(3)) and all(
        zero(ratio[g] - a * fit(g, p) / (ak * fit(g, pk))) for g in range(3)
    )
    ok_c = ok_scale and ok_single and ok_pair
    out.append(
        _rec(
            "V36c",
            r"A diagonal ratio $\tilde{\mathbf B}$ is identified only up to a common factor (absorbed by $a_i$), is "
            r"not identified together with $\boldsymbol p$ from one sample, and cancels from the ratio of two samples",
            "exact",
            "pass" if ok_c else "fail",
            r"Three genes, $J=3$: $(c\tilde b_g)(a/c)=\tilde b_ga$; "
            r"$\tilde b'_g=\tilde b_g\sum_jp_j\tilde\mu_{gj}/\sum_jp'_j\tilde\mu_{gj}$ reproduces every mean; "
            r"$\tilde b_g$ is absent from the ratio of the means of two samples",
            "",
            detail_full=f"scale {ok_scale}; single sample {ok_single}; pair {ok_pair}",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

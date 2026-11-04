# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V27: Gini coefficient of the type shares (vignette 03-1, #sec-gini-entropy,
@eq-gini and @fig-gini-derivation).

Derivation route (nothing typed from the vignette except the two gene
profiles of the figure):
- V27a: the Gini coefficient is defined from the type means as the mean
  absolute difference over the J^2 ordered pairs divided by twice the mean,
  Delta / (2 mu_bar). Sorted means are built as cumulative sums of positive
  increments, so every |mu_j - mu_l| is resolved by SymPy without typing an
  ordering. From that definition we compute (i) the pairwise share form
  (1/(2J)) sum_j sum_l |u_j - u_l| with u_j = mu_j / sum mu, (ii) the sorted
  share form (1/J) sum_k (2k - J - 1) u_(k); all three agree exactly for
  J = 2..10 (symbolic increments). The net weight 2k - J - 1 is recomputed as
  (number of smaller shares) - (number of larger shares) with symbolic J, k.
  The maximum: on the simplex the sorted form is linear in u, so its maximum
  is at a vertex; at u = e_J it is (J - 1)/J (symbolic J), and the correction
  J/(J - 1) maps it to 1. Range of the corrected value: a vertex enumeration
  on the ordered simplex for J = 2..10 (linear objective, the ordered simplex
  has vertices v_k = (0,...,0,1/k,...,1/k)) gives min 0 and max 1.
- V27b: exact rational arithmetic on the two genes of the figure.
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
    e = sp.simplify(e)
    if e == 0:
        return True
    return sp.cancel(sp.together(sp.expand(e))) == 0


def gini_from_means(mu):
    """Definition: mean absolute difference over J^2 pairs / (2 * mean)."""
    J = len(mu)
    delta = sum(sp.Abs(a - b) for a in mu for b in mu) / J**2
    mbar = sum(mu) / J
    return delta / (2 * mbar), delta, mbar


def pairwise_share(u):
    J = len(u)
    return sum(sp.Abs(a - b) for a in u for b in u) / (2 * J)


def sorted_share(u_sorted):
    J = len(u_sorted)
    return sum((2 * k - J - 1) * u_sorted[k - 1] for k in range(1, J + 1)) / J


def run() -> list[dict]:
    out = []

    # ---- V27a ----------------------------------------------------------------
    ok_forms = True
    for J in range(2, 11):
        d = sp.symbols(f"d1:{J + 1}", positive=True)  # mu_(1) = d1 > 0, increments > 0
        mu = [sum(d[: k + 1]) for k in range(J)]  # sorted, strictly increasing
        g_def, _, _ = gini_from_means(mu)
        tot = sum(mu)
        u = [m / tot for m in mu]
        g_pair = pairwise_share(u)
        g_sort = sorted_share(u)
        if not (zero(g_def - g_pair) and zero(g_pair - g_sort)):
            ok_forms = False
    # net weight of the share of rank k: (# smaller) - (# larger)
    Js, ks, i = sp.symbols("J k i", positive=True, integer=True)
    smaller = sp.summation(1, (i, 1, ks - 1))
    larger = sp.summation(1, (i, ks + 1, Js))
    ok_weight = zero(smaller - larger - (2 * ks - Js - 1))
    # maximum at u = e_J (symbolic J): only the top share is 1
    gmax = (2 * Js - Js - 1) * 1 / Js
    ok_max = zero(gmax - (Js - 1) / Js)
    ok_corr = zero(Js / (Js - 1) * gmax - 1)
    # range of the corrected value on the ordered simplex (vertex enumeration)
    ok_range = True
    for J in range(2, 11):
        vals = []
        for k in range(1, J + 1):
            v = [sp.Integer(0)] * (J - k) + [sp.Rational(1, k)] * k
            vals.append(sp.Rational(J, J - 1) * sorted_share(v))
        if min(vals) != 0 or max(vals) != 1:
            ok_range = False
    ok_a = ok_forms and ok_weight and ok_max and ok_corr and ok_range
    out.append(
        _rec(
            "V27a",
            r"Gini of the type shares: pairwise form $=$ sorted form $=$ mean absolute difference of the means over "
            r"$2\bar\mu_g$; net weight $2k-J-1$; maximum $(J-1)/J$ at one expressing type; "
            r"$J\,\mathrm{Gini}_g/(J-1)\in[0,1]$",
            "exact",
            "pass" if ok_a else "fail",
            r"Sorted means as cumulative sums of positive symbols, $J=2,\ldots,10$: the three forms agree; "
            r"$\#\{\text{smaller}\}-\#\{\text{larger}\}=2k-J-1$ for symbolic $J,k$; value $(J-1)/J$ at "
            r"$u=e_J$; vertices of the ordered simplex give corrected range $[0,1]$",
            "",
            detail_full=(
                f"forms {ok_forms}; weight {ok_weight}; max {ok_max}; correction {ok_corr}; range {ok_range}"
            ),
        )
    )

    # ---- V27b ----------------------------------------------------------------
    genes = {
        "low": [12, 15, 9, 14, 11, 13, 10, 16],
        "high": [1, 2, 1, 3, 1, 2, 4, 86],
    }
    res = {}
    for name, mu in genes.items():
        mu = [sp.Integer(x) for x in mu]
        g, delta, mbar = gini_from_means(mu)
        u = [m / sum(mu) for m in mu]
        J = len(mu)
        wsum = sum((2 * k - J - 1) * x for k, x in enumerate(sorted(mu), start=1))
        res[name] = {
            "gini": sp.nsimplify(g),
            "pair": pairwise_share(u),
            "sorted": sorted_share(sorted(u)),
            "corr": sp.Rational(J, J - 1) * g,
            "delta": delta,
            "mbar": mbar,
            "wsum_over_J2mbar": wsum / (J**2 * mbar),
            "min": min(mu),
        }
    lo, hi = res["low"], res["high"]
    ok_vals = (
        lo["gini"] == sp.Rational(21, 200)
        and hi["gini"] == sp.Rational(77, 100)
        and lo["corr"] == sp.Rational(3, 25)
        and hi["corr"] == sp.Rational(22, 25)
    )
    ok_same = all(r["gini"] == r["pair"] == r["sorted"] == r["wsum_over_J2mbar"] for r in res.values())
    ok_fig = (
        lo["mbar"] == hi["mbar"] == sp.Rational(25, 2)
        and lo["delta"] == sp.Rational(21, 8)  # 2.625
        and hi["delta"] == sp.Rational(77, 4)  # 19.25
        and 8**2 * lo["mbar"] == 800
        and lo["min"] >= 1
        and hi["min"] >= 1
    )
    J = 8
    one_type = sorted_share([sp.Integer(0)] * (J - 1) + [sp.Integer(1)])
    ok_one = one_type == sp.Rational(7, 8) and sp.Rational(J, J - 1) * one_type == 1
    ok_b = ok_vals and ok_same and ok_fig and ok_one
    out.append(
        _rec(
            "V27b",
            r"Genes of the Gini figure: $\mathrm{Gini}_g=21/200$ and $77/100$, corrected $3/25$ and $22/25$; "
            r"a gene expressed in one type has $(J-1)/J$",
            "exact",
            "pass" if ok_b else "fail",
            rf"Exact rationals: Gini ${sp.latex(lo['gini'])}$ and ${sp.latex(hi['gini'])}$, corrected "
            rf"${sp.latex(lo['corr'])}$ and ${sp.latex(hi['corr'])}$, $\Delta_g={sp.latex(lo['delta'])}$ and "
            rf"${sp.latex(hi['delta'])}$, $\bar\mu_g={sp.latex(lo['mbar'])}$, $J^2\bar\mu_g=800$; one type: "
            rf"${sp.latex(one_type)}$; all four forms agree",
            "",
            detail_full=json.dumps({k: {kk: str(vv) for kk, vv in v.items()} for k, v in res.items()}),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

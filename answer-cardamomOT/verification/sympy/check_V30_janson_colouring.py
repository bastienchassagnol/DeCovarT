# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V30: the colouring (Janson) bound of T5 (vignette 03-1, #thm-fs-nb-dep,
@eq-fs-colouring and @eq-fs-janson).

Derivation route:
- V30a: Jensen step for a generic number of classes k. With IndexedBase
  y_1..y_k and ybar = k^{-1} sum_c y_c: (i) sum_c (y_c - ybar) = 0 (symbolic
  k); (ii) the tangent-line function f(x) = e^x - 1 - x has f(0) = 0, a single
  critical point x = 0 and f'' = e^x > 0, so f >= 0; (iii) the termwise
  identity e^{y_c} = e^{ybar} {f(y_c - ybar) + 1 + (y_c - ybar)}, summed, gives
  sum_c e^{y_c} - k e^{ybar} = e^{ybar} sum_c f(y_c - ybar) >= 0. With
  y_c = k u S_c and sum_c S_c = S, ybar = u S (symbolic k). Independence
  within a class turns E e^{k u S_c} into the product. Derivatives of
  K_bar_S at 0 for standardised Z_g (K(0) = K'(0) = 0, K''(0) = 1, generic
  functions): K_bar_S(0) = K_bar_S'(0) = 0, K_bar_S''(0) = k |M| (classes 8, 7,
  5: 60). Numeric toy (exact Gaussian MGFs): random graph on 14 genes (seed
  20261101), Welsh-Powell colouring, Z ~ N(0, I + rho A) with rho = 0.15;
  exp(u^2 1'R1 / 2) <= k^{-1} sum_c exp(|I_c| k^2 u^2 / 2) on a grid of u, and
  the pointwise inequality on 2000 pseudo-random draws.
- V30b: (i) algebra: the positive root of chi^2 - chi - 2m = 0 is
  (1 + sqrt(8m + 1))/2 and the other root is <= 0, so chi(chi - 1)/2 <= m
  iff chi <= (1 + sqrt(8m + 1))/2 for chi > 0; (ii) the lemma that a graph of
  chromatic number chi has at least chi(chi - 1)/2 edges, and the greedy
  bounds (Welsh-Powell colouring uses at most min(d_max + 1,
  (1 + sqrt(8m + 1))/2) colours), by exhaustive enumeration of every graph on
  up to 6 labelled vertices; (iii) the simulation numbers m = 30, d_max = 6:
  (1 + sqrt 241)/2 = 8.26 and min = 7.
"""

from __future__ import annotations

import json
import math
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
    if sp.cancel(sp.together(sp.expand(e))) == 0:
        return True
    return sp.simplify(e) == 0


def welsh_powell(nv, adj):
    order = sorted(range(nv), key=lambda i: (-len(adj[i]), i))
    col = {}
    for vtx in order:
        used = {col[w] for w in adj[vtx] if w in col}
        c = 0
        while c in used:
            c += 1
        col[vtx] = c
    return col


def colourable(nv, adj, k):
    col = [-1] * nv

    def bt(i):
        if i == nv:
            return True
        for c in range(k):
            if all(col[w] != c for w in adj[i] if w < i):
                col[i] = c
                if bt(i + 1):
                    return True
        col[i] = -1
        return False

    return bt(0)


def chromatic(nv, adj):
    if nv == 0:
        return 0
    k = 1
    while not colourable(nv, adj, k):
        k += 1
    return k


def run() -> list[dict]:
    out = []

    # ---- V30a ---------------------------------------------------------------
    k = sp.Symbol("k", positive=True, integer=True)
    c = sp.Symbol("c", positive=True, integer=True)
    y = sp.IndexedBase("y")
    S_ = sp.IndexedBase("S")
    uu = sp.Symbol("u", real=True)
    d_ = sp.Symbol("d", positive=True, integer=True)
    ybar = sp.Sum(y[d_], (d_, 1, k)) / k
    centred = sp.expand(sp.Sum(y[c] - ybar, (c, 1, k)).doit())  # ybar is constant in c
    ok_centre = sp.simplify(centred.xreplace({d_: c})) == 0  # rename the bound index
    x = sp.Symbol("x", real=True)
    f = sp.exp(x) - 1 - x
    crit = sp.solve(sp.diff(f, x), x)
    ok_f = f.subs(x, 0) == 0 and crit == [0] and sp.diff(f, x, 2) == sp.exp(x)
    yc, yb = sp.symbols("y_c ybar", real=True)
    ok_term = zero(sp.exp(yb) * (f.subs(x, yc - yb) + 1 + (yc - yb)) - sp.exp(yc))
    # y_c = k u S_c: mean is u S
    mean_y = sp.Sum(k * uu * S_[c], (c, 1, k)) / k
    ok_mean = zero(sp.simplify(mean_y - uu * sp.Sum(S_[c], (c, 1, k))))
    # derivatives of K_bar_S at 0 for standardised Z_g
    ok_L = True
    Ls = {}
    for sizes in ((3, 2, 2), (8, 7, 5)):
        kk = len(sizes)
        # each K_g represented by its 3-jet a_g x + b_g x^2/2 + c_g x^3 (generic coefficients);
        # derivatives of K_bar_S at 0 up to order 2 depend only on the 2-jet
        jets = []
        for cc, sz in enumerate(sizes):
            cls = []
            for g in range(sz):
                a_, b_, c_ = sp.symbols(f"a_{cc}_{g} b_{cc}_{g} c_{cc}_{g}", real=True)
                cls.append((a_, b_, c_))
            jets.append(cls)
        L = sp.log(
            sum(sp.exp(sum(a_ * kk * uu + b_ * (kk * uu) ** 2 / 2 + c_ * (kk * uu) ** 3 for a_, b_, c_ in cls))
                for cls in jets) / kk
        )
        std = {}
        for cls in jets:
            for a_, b_, _ in cls:
                std[a_] = 0  # K'(0) = 0 (centred)
                std[b_] = 1  # K''(0) = 1 (standardised)
        vals = [sp.simplify(sp.diff(L, uu, r).subs(uu, 0).subs(std)) for r in range(3)]
        Ls[sizes] = vals
        if not (vals[0] == 0 and vals[1] == 0 and vals[2] == kk * sum(sizes)):
            ok_L = False
    # numeric toy: exact Gaussian MGFs
    rng = random.Random(20261101)
    N, rho = 14, 0.15
    edges = set()
    while len(edges) < 18:
        a, b = sorted(rng.sample(range(N), 2))
        edges.add((a, b))
    adj = [set() for _ in range(N)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    col = welsh_powell(N, adj)
    kk = max(col.values()) + 1
    proper = all(col[a] != col[b] for a, b in edges)
    sizes = [sum(1 for g in range(N) if col[g] == cc) for cc in range(kk)]
    R = sp.eye(N)
    for a, b in edges:
        R[a, b] = R[b, a] = sp.nsimplify(rho)
    pd = all(ev > 0 for ev in R.eigenvals(multiple=True)) if N <= 14 else False
    one_R_one = sum(R)
    ok_gauss = True
    worst_gap = math.inf
    for i in range(1, 101):
        uv = i / 100
        lhs = math.exp(uv**2 * float(one_R_one) / 2)
        rhs = sum(math.exp(sz * kk**2 * uv**2 / 2) for sz in sizes) / kk
        worst_gap = min(worst_gap, rhs / lhs)
        if lhs > rhs * (1 + 1e-12):
            ok_gauss = False
    ok_point = True
    for _ in range(2000):
        Sc = [rng.gauss(0, 1) * 3 for _ in range(kk)]
        for uv in (-0.2, -0.05, 0.05, 0.2):
            if math.exp(uv * sum(Sc)) > sum(math.exp(kk * uv * s) for s in Sc) / kk * (1 + 1e-12):
                ok_point = False
    ok_a = ok_centre and ok_f and ok_term and ok_mean and ok_L and proper and pd and ok_gauss and ok_point
    out.append(
        _rec(
            "V30a",
            r"Jensen step of T5: $\operatorname{E}e^{uS_M}\le k_G^{-1}\sum_c\prod_{g\in\mathcal M_c}\operatorname{E}"
            r"e^{k_GuZ_g}=\exp\bar K_S(u)$; $\bar K_S(0)=\bar K_S'(0)=0$ and $\bar K_S''(0)=k_G\lvert M\rvert$",
            "exact",
            "pass" if ok_a else "fail",
            r"Symbolic $k$: $\sum_c(y_c-\bar y)=0$, $e^x-1-x\ge0$ (single critical point 0, $f''=e^x$), termwise "
            r"identity gives $\sum_ce^{y_c}-ke^{\bar y}=e^{\bar y}\sum_cf(y_c-\bar y)\ge0$ with $\bar y=uS$; generic "
            rf"CGFs: $\bar K_S''(0)={Ls[(8, 7, 5)][2]}$ for classes 8, 7, 5. Gaussian toy, 14 genes, 18 edges, "
            rf"$k_G={kk}$: inequality holds for $u=0.01,\ldots,1$ (smallest ratio {worst_gap:.3f}); pointwise on 2000 draws",
            "",
            detail_full=(
                f"centre {ok_centre}; tangent {ok_f}; termwise {ok_term}; mean {ok_mean}; L derivs "
                f"{ {str(kx): [str(v) for v in vv] for kx, vv in Ls.items()} }; proper {proper}; pd {pd}; "
                f"gauss {ok_gauss}; pointwise {ok_point}; class sizes {sizes}"
            ),
        )
    )

    # ---- V30b ---------------------------------------------------------------
    chi = sp.Symbol("chi", positive=True)
    m = sp.Symbol("m", nonnegative=True)
    roots = sp.solve(sp.Eq(chi**2 - chi - 2 * m, 0), chi)
    rplus = (1 + sp.sqrt(8 * m + 1)) / 2
    rminus = (1 - sp.sqrt(8 * m + 1)) / 2
    ok_roots = zero((chi**2 - chi - 2 * m).subs(chi, rplus)) and any(zero(r_ - rplus) for r_ in roots)
    ok_fact = zero(sp.expand((chi - rplus) * (chi - rminus)) - (chi**2 - chi - 2 * m))
    ok_neg = sp.simplify(rminus.subs(m, 0)) == 0 and sp.diff(rminus, m).is_negative  # r- <= 0 for m >= 0
    # exhaustive: all graphs on nv <= 6 labelled vertices
    n_graphs = 0
    ok_enum = True
    for nv in range(1, 7):
        pairs = [(a, b) for a in range(nv) for b in range(a + 1, nv)]
        for mask in range(1 << len(pairs)):
            adj = [set() for _ in range(nv)]
            me = 0
            for i, (a, b) in enumerate(pairs):
                if mask >> i & 1:
                    adj[a].add(b)
                    adj[b].add(a)
                    me += 1
            n_graphs += 1
            ch = chromatic(nv, adj)
            gr = max(welsh_powell(nv, adj).values()) + 1
            dmax = max(len(s) for s in adj)
            bound = (1 + math.sqrt(8 * me + 1)) / 2
            if not (me >= ch * (ch - 1) // 2 and ch <= gr and gr <= dmax + 1 and gr <= bound + 1e-12
                    and me >= gr * (gr - 1) // 2):
                ok_enum = False
    val30 = rplus.subs(m, 30)
    ok_num = round(float(val30), 2) == 8.26 and min(sp.Integer(6 + 1), val30) == 7
    ok_b = ok_roots and ok_fact and ok_neg and ok_enum and ok_num
    out.append(
        _rec(
            "V30b",
            r"$k_G\le\min\{d_{\max}+1,(1+\sqrt{8\lvert\mathcal E\rvert+1})/2\}$: a graph of chromatic number $\chi$ has "
            r"at least $\chi(\chi-1)/2$ edges; the Welsh--Powell colouring attains both bounds; $8.26$ and $7$ for "
            r"30 edges and $d_{\max}=6$",
            "exact",
            "pass" if ok_b else "fail",
            rf"$\chi^2-\chi-2m=(\chi-r_+)(\chi-r_-)$ with $r_+={sp.latex(rplus)}$, $r_-\le0$; all {n_graphs} labelled "
            rf"graphs on at most 6 vertices: $m\ge\chi(\chi-1)/2$, and the greedy count $k$ satisfies $k\le d_{{\max}}+1$ "
            rf"and $m\ge k(k-1)/2$; $r_+(30)={float(val30):.4f}$",
            "",
            detail_full=f"roots {ok_roots}; factor {ok_fact}; r- {ok_neg}; enumeration {ok_enum}; numbers {ok_num}",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))

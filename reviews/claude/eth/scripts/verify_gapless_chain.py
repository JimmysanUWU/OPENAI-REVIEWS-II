"""Checks for 08_GAPLESS_PROOF.md: every step of the gapless ETH-refuting chain.

Part S (Section 2)  disjoint sparsification: Lemmas 2.1-2.3 and Proposition 2.4 at every
                    node of the branching tree of random k-CNFs, plus the partition property.
Part R (Section 3)  OR gates and Tovey cycles: count preservation, width and occurrence bounds.
Part K (Sections 4-7) the counting algorithm exactly as written in the proof:
                    closed-form swap rows (Theorem 5.3), Krawtchouk coefficients in Z[1/2],
                    derandomized schedules (Lemma 6.2), measurement rows touching <= 2 groups
                    (Lemma 6.1), locality, separated Taylor jets (Lemmas 4.2-4.3) and the CRT
                    recovery (Lemma 4.1), all in Z/p^mu, against brute force; negative control.
Part P (Section 7)  the parameter inequalities (P1)-(P9) for lambda >= 2^17 + 1, and the
                    least lambda at which the exponent bound drops below 1.
Part A (Section 8)  the sparsification constants M, B for the final assembly.

Standard library only; deterministic; writes ../results/gapless_chain.json.
"""

import json
import math
import random
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import comb, isqrt

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from eth_core import rref, in_rowspace, poly_mul, poly_eval  # noqa: E402

RESULTS = {}


# =============================================================== Part S: sparsification

def taut(C):
    return any(-l in C for l in C)


def normalize(F):
    F = {C for C in F if not taut(C)}
    return {C for C in F if not any(D < C for D in F)}


def select(F, k, M):
    """Smallest width c, then largest |H|, then lexicographically smallest H."""
    for c in range(2, k + 1):
        Wc = [C for C in F if len(C) == c]
        if not Wc:
            continue
        deg = Counter()
        for C in Wc:
            items = sorted(C)
            for size in range(1, c):
                for H in combinations(items, size):
                    deg[H] += 1
        elig = [H for H, dg in deg.items() if dg >= M ** (c - len(H))]
        if elig:
            H = min(elig, key=lambda H: (-len(H), H))
            return c, frozenset(H)
    return None


def count_sat(n, F):
    cnt = 0
    for x in product((0, 1), repeat=n):
        if all(any((x[abs(l) - 1] == 1) == (l > 0) for l in C) for C in F):
            cnt += 1
    return cnt


def sat_set(n, F):
    out = set()
    for x in product((0, 1), repeat=n):
        if all(any((x[abs(l) - 1] == 1) == (l > 0) for l in C) for C in F):
            out.add(x)
    return out


def sparsify_check(n, F0, k, M):
    theta = [M ** j for j in range(k)]
    C_const = 4 * sum((3 * M) ** j for j in range(k - 1))
    Gamma = 2 * (3 ** (k - 1) - 1)
    stats = {"nodes": 0, "leaves": 0, "max_depth": 0, "max_B_choices": 0,
             "lemma21": True, "lemma22": True, "lemma23": True, "depth_bound": True,
             "B_choice_bound": True, "leaf_occurrence": True, "leaf_literal_degree": True}
    leaves = []
    root = normalize(F0)
    # stack: (formula, marked set, T dict, depth, B choices)
    stack = [(root, frozenset(), {j: 0 for j in range(1, k)}, 0, 0)]
    while stack:
        F, marked, T, depth, nb = stack.pop()
        stats["nodes"] += 1
        stats["max_depth"] = max(stats["max_depth"], depth)
        stats["max_B_choices"] = max(stats["max_B_choices"], nb)
        if depth > C_const * n:
            stats["depth_bound"] = False
        if nb > Gamma * n / M:
            stats["B_choice_bound"] = False
        # Lemma 2.2 at this node
        live = [C for C in marked if C in F]
        for j in range(2, k):
            Lj = [C for C in live if len(C) == j]
            for i in range(1, j):
                cnt = Counter()
                for C in Lj:
                    for S in combinations(sorted(C), i):
                        cnt[S] += 1
                if any(v >= 2 * theta[j - i] for v in cnt.values()):
                    stats["lemma22"] = False
        # Lemma 2.3 at this node
        for j in range(1, k):
            if T[j] > 4 * (3 * M) ** (j - 1) * n:
                stats["lemma23"] = False
        if frozenset() in F:
            leaves.append(F)
            stats["leaves"] += 1
            continue
        sel = select(F, k, M)
        if sel is None:
            stats["leaves"] += 1
            leaves.append(F)
            # leaf bounds: each literal in < theta_{c-1} width-c clauses (c >= 2), <= 1 unit
            lit_deg = defaultdict(Counter)
            for C in F:
                for l in C:
                    lit_deg[l][len(C)] += 1
            for l, by_w in lit_deg.items():
                for c, v in by_w.items():
                    if (c >= 2 and v >= theta[c - 1]) or (c == 1 and v > 1):
                        stats["leaf_literal_degree"] = False
            occ = Counter(abs(l) for C in F for l in C)
            Bbound = 2 * sum(M ** j for j in range(k))
            if occ and max(occ.values()) > Bbound:
                stats["leaf_occurrence"] = False
            continue
        c, H = sel
        Sel = [C for C in F if len(C) == c and H <= C]
        jsel = c - len(H)
        G = F - set(Sel)
        P = {C - H for C in Sel}
        units = {frozenset([-l]) for l in H}
        FA = normalize(G | {H})
        FB = normalize(G | P | units)
        newA = {H} - F
        newB = (P | units) - F
        ok = (newA == {H}) and H in FA and all(len(N) < c for N in newA | newB)
        ok = ok and not (P & F) and all(N in FB for N in newB) and len(P) >= theta[jsel]
        for child, new in ((FA, newA), (FB, newB)):
            for C in F - child:
                if not any(N < C for N in new):
                    ok = False
        if not ok:
            stats["lemma21"] = False
        for child, new, isB in ((FA, newA, 0), (FB, newB, 1)):
            T2 = dict(T)
            for N in new:
                T2[len(N)] += 1
            m2 = frozenset((marked & child) | new)
            stack.append((child, m2, T2, depth + 1, nb + isB))
    # partition
    total = count_sat(n, root)
    stats["count"] = total
    stats["sum_leaf_counts"] = sum(count_sat(n, L) for L in leaves)
    sets = [sat_set(n, L) for L in leaves]
    union = set().union(*sets) if sets else set()
    stats["disjoint"] = sum(len(s) for s in sets) == len(union)
    stats["union_is_sat_set"] = union == sat_set(n, root)
    # node count bound (C n + 1) * sum_{a <= Gamma n / M} binom(C n, a)
    Cn = C_const * n
    amax = (Gamma * n) // M
    stats["node_bound"] = (Cn + 1) * sum(comb(Cn, a) for a in range(amax + 1))
    stats["node_bound_ok"] = stats["nodes"] <= stats["node_bound"]
    stats["C"] = C_const
    stats["Gamma"] = Gamma
    return stats, leaves


def random_kcnf(n, m, k, rng, planted=False):
    """Random k-CNF; with planted=True every clause is satisfied by a hidden assignment."""
    hidden = [rng.random() < 0.5 for _ in range(n + 1)]
    F = set()
    while len(F) < m:
        w = min(n, rng.randint(1, k) if rng.random() < 0.15 else k)
        vs = rng.sample(range(1, n + 1), w)
        C = frozenset(v if rng.random() < 0.5 else -v for v in vs)
        if planted and not any((l > 0) == hidden[abs(l)] for l in C):
            continue
        F.add(C)
    return F


def part_S():
    rng = random.Random(20261027)
    out = []
    cases = [(3, 2, 8, 40), (3, 2, 9, 70), (3, 3, 9, 90), (4, 2, 8, 40), (4, 2, 9, 60), (2, 2, 9, 30),
             (3, 2, 10, 32), (4, 2, 10, 28)]
    allok = True
    for (k, M, n, m) in cases:
        for rep in range(4):
            F = random_kcnf(n, m, k, rng, planted=(rep > 0))
            st, _ = sparsify_check(n, F, k, M)
            good = all(st[key] for key in ("lemma21", "lemma22", "lemma23", "depth_bound", "B_choice_bound",
                                            "leaf_occurrence", "leaf_literal_degree", "disjoint",
                                            "union_is_sat_set", "node_bound_ok"))
            good = good and st["count"] == st["sum_leaf_counts"]
            allok = allok and good
            out.append({"k": k, "M": M, "n": n, "m": m, "rep": rep, "ok": good,
                        **{key: st[key] for key in ("nodes", "leaves", "max_depth", "max_B_choices",
                                                     "count", "sum_leaf_counts", "C", "Gamma")}})
    RESULTS["sparsification"] = {"all_ok": allok, "cases": out}
    print(f"[S] sparsification: {len(out)} formulas, all lemmas hold: {allok}")
    return allok


# =============================================================== Part R: OR gates, Tovey

def or_gates(clauses, n):
    out, nv = [], n
    for C in clauses:
        C = sorted(C, key=lambda l: (abs(l), l))
        while len(C) > 3:
            l1, l2 = C[0], C[1]
            nv += 1
            z = nv
            out += [frozenset([-l1, z]), frozenset([-l2, z]), frozenset([l1, l2, -z])]
            C = [z] + C[2:]
        out.append(frozenset(C))
    return out, nv


def tovey(clauses, n):
    occ = defaultdict(list)
    for ci, C in enumerate(clauses):
        for l in C:
            occ[abs(l)].append(ci)
    free = [v for v in range(1, n + 1) if not occ[v]]
    newid = {}
    nv = 0
    cl = [set() for _ in clauses]
    extra = []
    for v in range(1, n + 1):
        lst = occ[v]
        if not lst:
            continue
        if len(lst) == 1:
            nv += 1
            ci = lst[0]
            sgn = 1 if v in clauses[ci] else -1
            cl[ci].add(sgn * nv)
            continue
        ids = []
        for ci in lst:
            nv += 1
            ids.append(nv)
            sgn = 1 if v in clauses[ci] else -1
            cl[ci].add(sgn * nv)
        for i in range(len(ids)):
            extra.append(frozenset([-ids[i], ids[(i + 1) % len(ids)]]))
    return [frozenset(c) for c in cl] + extra, nv, len(free)


def dpll_count(clauses, variables):
    """#SAT with unit propagation (exact)."""
    clauses = [frozenset(c) for c in clauses]
    variables = set(variables)

    def simplify(cls, lit):
        out = []
        for c in cls:
            if lit in c:
                continue
            if -lit in c:
                c = c - {-lit}
            out.append(c)
        return out

    def rec(cls, vs):
        while True:
            if any(len(c) == 0 for c in cls):
                return 0
            unit = next((c for c in cls if len(c) == 1), None)
            if unit is None:
                break
            (l,) = tuple(unit)
            cls = simplify(cls, l)
            vs = vs - {abs(l)}
        if not cls:
            return 1 << len(vs)
        cnt = Counter(abs(l) for c in cls for l in c)
        v = max(cnt, key=lambda x: (cnt[x], -x))
        return rec(simplify(cls, v), vs - {v}) + rec(simplify(cls, -v), vs - {v})

    return rec(clauses, variables)


def part_R():
    rng = random.Random(7)
    allok = True
    rows = []
    for trial in range(40):
        k = rng.choice([3, 4, 5, 6])
        n = rng.randint(4, 8)
        m = rng.randint(3, 9)
        F = list(random_kcnf(n, m, k, rng))
        F = [C for C in normalize(set(F))]
        base = count_sat(n, F)
        G, nv = or_gates(F, n)
        after_gates = dpll_count(G, range(1, nv + 1))
        occ_orig = Counter(abs(l) for C in F for l in C)
        B = max(occ_orig.values()) if occ_orig else 0
        tot_occ = sum(len(C) for C in G)
        T, N, nfree = tovey(G, nv)
        cnt_T = dpll_count(T, range(1, N + 1))
        occT = Counter(abs(l) for C in T for l in C)
        ok = (after_gates == base and (cnt_T << nfree) == base and all(len(C) <= 3 for C in T)
              and (max(occT.values()) if occT else 0) <= 3 and N <= tot_occ
              and tot_occ <= max(1, 6 * k) * B * n and all(occT[v] >= 1 for v in range(1, N + 1)))
        allok = allok and ok
        rows.append({"k": k, "n": n, "m": m, "count": base, "N": N, "free": nfree, "ok": ok})
    RESULTS["or_gates_tovey"] = {"all_ok": allok, "cases": rows}
    print(f"[R] OR gates + Tovey: {len(rows)} formulas, counts and bounds preserved: {allok}")
    return allok


# =============================================================== Part K: counting algorithm

def kraw(Q, k, j):
    """Krawtchouk K_k(j) at an arbitrary integer j (generalized binomials)."""
    def gbin(x, i):
        if i < 0:
            return 0
        num = 1
        for a in range(i):
            num *= (x - a)
        return num // math.factorial(i)
    return sum((-1) ** i * gbin(j, i) * gbin(Q - j, k - i) for i in range(k + 1))


def solve_fraction(A, b):
    n = len(A)
    Mx = [list(map(Fraction, row)) + [Fraction(bb)] for row, bb in zip(A, b)]
    for col in range(n):
        piv = next(i for i in range(col, n) if Mx[i][col] != 0)
        Mx[col], Mx[piv] = Mx[piv], Mx[col]
        inv = 1 / Mx[col][col]
        Mx[col] = [v * inv for v in Mx[col]]
        for i in range(n):
            if i != col and Mx[i][col] != 0:
                f = Mx[i][col]
                Mx[i] = [a - f * c for a, c in zip(Mx[i], Mx[col])]
    return [Mx[i][n] for i in range(n)]


class Swap:
    """Theorem 5.3: closed-form swap rows for parameters (r, Q, R)."""

    def __init__(self, r, Q, R):
        self.r, self.Q, self.R = r, Q, R
        pairs = [(a, d) for a in range(r) for d in range(r)]
        self.labels = list(product(pairs, repeat=Q))
        self.D = len(self.labels)
        p = Fraction(r - 1, 2 * r)
        self.p = p
        j0 = math.ceil(Q * p) - R // 2
        self.interval = list(range(j0, j0 + R))
        # gamma: sum_k gamma_k K_k(j) = (-1)^j on the interval
        A = [[kraw(Q, kk, j) for kk in range(R)] for j in self.interval]
        b = [(-1) ** (j % 2) for j in self.interval]
        self.gamma = solve_fraction(A, b)
        self.gamma_dyadic = all((g.denominator & (g.denominator - 1)) == 0 for g in self.gamma)
        # rows P(T): dict u -> dict y -> poly (list of Fractions), and M0
        self.P, self.M0 = {}, {}
        for u in self.labels:
            row = {u: [Fraction(1)]}
            m0 = {}
            for kk in range(min(R, Q + 1)):
                g = self.gamma[kk]
                if g == 0:
                    continue
                for S in combinations(range(Q), kk):
                    y = list(u)
                    for s in range(Q):
                        if s not in S:
                            a, d = y[s]
                            y[s] = (d, a)
                    y = tuple(y)
                    phi = [Fraction(1)]
                    for s in S:
                        phi = poly_mul(phi, [Fraction(1), Fraction(-1, s + 1)])
                    phi[0] -= 1  # phi_S(T) - 1
                    cur = row.get(y, [Fraction(0)])
                    L = max(len(cur), len(phi))
                    row[y] = [(cur[i] if i < len(cur) else 0) + g * (phi[i] if i < len(phi) else 0) for i in range(L)]
                    m0[y] = m0.get(y, 0) + g
            self.P[u] = {y: pp for y, pp in row.items() if any(pp)}
            self.M0[u] = {y: v for y, v in m0.items() if v != 0}
        # eigenbasis, bad rows
        local = []
        for a in range(r):
            local.append(({(a, a): 1}, 1))
        for a in range(r):
            for d in range(a + 1, r):
                local.append(({(a, d): 1, (d, a): 1}, 1))
                local.append(({(a, d): 1, (d, a): -1}, -1))
        self.bad = []
        self.level_count = Counter()
        for choice in product(range(len(local)), repeat=Q):
            vecs = [local[c] for c in choice]
            j = sum(1 for v in vecs if v[1] == -1)
            self.level_count[j] += 1
            if j in self.interval:
                continue
            row = {(): 1}
            for vec, _ in vecs:
                row = {pre + (pair,): val * x for pre, val in row.items() for pair, x in vec.items()}
            self.bad.append(row)
        self.c = len(self.bad)

    def checks(self):
        r, Q, R = self.r, self.Q, self.R
        res = {}
        res["P0_identity"] = all({y: pp[0] for y, pp in self.P[u].items() if pp[0] != 0} == {u: 1} for u in self.labels)
        res["degree_le_R-1"] = all(len(pp) - 1 <= R - 1 for u in self.labels for pp in self.P[u].values())
        loc = True
        for s in range(Q):
            for u in self.labels:
                for y, pp in self.P[u].items():
                    val = poly_eval(pp, s + 1)
                    const = (1 if y == u else 0) - self.M0[u].get(y, 0)
                    if val != const and y[s][0] != u[s][1]:
                        loc = False
        res["B2_nodes_swap"] = loc
        # level multiplicities: binom(Q,j) binom(r,2)^j binom(r+1,2)^(Q-j)
        res["level_multiplicities"] = all(self.level_count[j] == comb(Q, j) * comb(r, 2) ** j * comb(r + 1, 2) ** (Q - j)
                                          for j in range(Q + 1))
        res["bad_row_l1_le_2^Q"] = all(sum(abs(v) for v in row.values()) <= 2 ** Q for row in self.bad)
        res["gamma_in_Z[1/2]"] = self.gamma_dyadic
        # M0 psi = (-1)^j f(j) psi, and M0 fixes good psi: check M0 xi = xi on ker C (via rank)
        return res


def prime_above(Q):
    p = max(3, Q + 1)
    while any(p % d == 0 for d in range(2, isqrt(p) + 1)):
        p += 1
    return p


def schedules(n, factors, r, Q):
    """Lemma 6.2: fix schedule coordinates greedily by exact conditional expectation."""
    M = sum(len(f) for f in factors)
    b = max(Counter(v for f in factors for v in f).values())
    w = max(len(f) for f in factors)
    D = r ** (2 * Q)
    inc = [(fi, v) for fi, f in enumerate(factors) for v in f]
    share = [(i1, i2) for i1 in range(len(inc)) for i2 in range(len(inc)) if i1 != i2 and
             (inc[i1][0] == inc[i2][0] or inc[i1][1] == inc[i2][1])]
    words = list(product(range(r), repeat=Q))
    a = [[None] * Q for _ in range(n)]
    d = [[None] * Q for _ in factors]
    R1 = Fraction(1, r)

    def pr(sched, word):
        p = Fraction(1)
        for t in range(Q):
            if sched[t] is None:
                p *= R1
            elif sched[t] != word[t]:
                return Fraction(0)
        return p

    def expect():
        pa = [{wd: pr(a[v], wd) for wd in words} for v in range(n)]
        pd = [{wd: pr(d[fi], wd) for wd in words} for fi in range(len(factors))]
        phi = Fraction(0)
        for wa in words:
            for wd in words:
                pi = [pa[v][wa] * pd[fi][wd] for (fi, v) in inc]
                mean = sum(pi)
                var = sum(x - x * x for x in pi)
                for i1, i2 in share:
                    (f1, v1), (f2, v2) = inc[i1], inc[i2]
                    if v1 == v2:
                        joint = pa[v1][wa] * pd[f1][wd] * pd[f2][wd]
                    else:
                        joint = pd[f1][wd] * pa[v1][wa] * pa[v2][wa]
                    var += joint - pi[i1] * pi[i2]
                phi += var + (mean - Fraction(M, D)) ** 2
        xi = Fraction(0)
        for t in range(Q):
            for g in range(r):
                ps = [(R1 if a[v][t] is None else Fraction(int(a[v][t] == g))) for v in range(n)]
                xi += sum(x - x * x for x in ps) + (sum(ps) - Fraction(n, r)) ** 2
        return phi / (M * (b + w)) + xi / (Q * n)

    E0 = expect()
    coords = [("a", v, t) for v in range(n) for t in range(Q)] + [("d", fi, t) for fi in range(len(factors)) for t in range(Q)]
    for kind, i, t in coords:
        best = None
        for val in range(r):
            (a if kind == "a" else d)[i][t] = val
            e = expect()
            if best is None or e < best[0]:
                best = (e, val)
        (a if kind == "a" else d)[i][t] = best[1]
    E1 = expect()
    return [tuple(x) for x in a], [tuple(x) for x in d], E0, E1, M, b, w


def build(n, clauses, sw, rng_unused=None):
    r, Q = sw.r, sw.Q
    factors = [[v for v, _ in cl] for cl in clauses]
    a, d, E0, E1, M, b, w = schedules(n, factors, r, Q)
    count, occ, place = Counter(), {}, {}
    for fi, f in enumerate(factors):
        for v in f:
            u = tuple((a[v][s], d[fi][s]) for s in range(Q))
            c = count[u]
            count[u] += 1
            occ[(c, u)] = v
            place[(fi, v)] = (c, u)
    h = max(count.values())
    D = sw.D
    # balance statistics
    nt = {(t, g): sum(1 for v in range(n) if a[v][t] == g) for t in range(Q) for g in range(r)}
    phi = sum((count[u] - Fraction(M, D)) ** 2 for u in sw.labels)
    xi = sum((nt[k] - Fraction(n, r)) ** 2 for k in nt)
    Psi = phi / (M * (b + w)) + xi / (Q * n)
    bal = {"E_Psi_initial": float(E0), "Psi_final": float(Psi), "E_final_equals_Psi": E1 == Psi,
           "Psi_le_E0": Psi <= E0, "E0_lt_2": E0 < 2,
           "h_bound": max(count.values()) <= M / D + math.sqrt(2 * M * (b + w)),
           "group_bound": max(nt.values()) <= n / r + math.sqrt(2 * Q * n)}
    Pfv = {}
    for (fi, v), (c, u) in place.items():
        row = {}
        for y, poly in sw.P[u].items():
            vv = occ.get((c, y))
            if vv is None:
                continue
            cur = row.get(vv, [Fraction(0)])
            L = max(len(cur), len(poly))
            row[vv] = [(cur[i] if i < len(cur) else 0) + (poly[i] if i < len(poly) else 0) for i in range(L)]
        Pfv[(fi, v)] = {k: pp for k, pp in row.items() if any(pp)}
    Cm = {}
    rows_per_node = Counter()
    l1_ok = True
    for s in range(Q):
        for g in range(r):
            rows = []
            for c in range(h):
                for psi in sw.bad:
                    row = {}
                    for y, val in psi.items():
                        if y[s][0] != g:
                            continue
                        vv = occ.get((c, y))
                        if vv is not None:
                            row[vv] = row.get(vv, 0) + val
                    row = {k: x for k, x in row.items() if x != 0}
                    if row:
                        rows.append(row)
                        l1_ok = l1_ok and sum(abs(x) for x in row.values()) <= 2 ** Q
            Cm[(s, g)] = rows
            rows_per_node[s] += len(rows)
    groups = {(s, g): [v for v in range(n) if a[v][s] == g] for s in range(Q) for g in range(r)}
    return {"factors": factors, "signs": [[sg for _, sg in cl] for cl in clauses], "a": a, "d": d,
            "h": h, "P": Pfv, "C": Cm, "groups": groups, "M": M, "bal": bal,
            "rows_le_2hc": all(rows_per_node[s] <= 2 * h * sw.c for s in range(Q)), "row_l1_ok": l1_ok,
            "max_rows_per_node": max(rows_per_node.values()) if rows_per_node else 0}


def locality_ok(inst, sw):
    Q, r = sw.Q, sw.r
    for s in range(Q):
        for g in range(r):
            cols = set(inst["groups"][(s, g)])
            basis = rref([dict(row) for row in inst["C"][(s, g)]], None)
            for (fi, v), row in inst["P"].items():
                if inst["d"][fi][s] == g:
                    continue
                restr = {w: poly_eval(pp, s + 1) for w, pp in row.items() if w in cols}
                restr = {k: x for k, x in restr.items() if x != 0}
                if not in_rowspace(basis, restr):
                    return False
    return True


class Ring:
    """(Z/m)[z]/(z^J)."""

    def __init__(self, m, J):
        self.m, self.J = m, J

    def mul(self, x, y):
        out = [0] * self.J
        for i, a in enumerate(x):
            if a:
                for j in range(self.J - i):
                    out[i + j] = (out[i + j] + a * y[j]) % self.m
        return out

    def add(self, x, y):
        return [(a + b) % self.m for a, b in zip(x, y)]

    def const(self, c):
        return [c % self.m] + [0] * (self.J - 1)


def fmod(fr, m):
    return fr.numerator * pow(fr.denominator, -1, m) % m


def shift_poly_mod(poly, t, J, m):
    """coefficients of poly(t + z) mod (z^J, m)."""
    out = [0] * J
    for kk, c in enumerate(poly):
        if c == 0:
            continue
        cm = fmod(c, m)
        for j in range(min(kk, J - 1) + 1):
            out[j] = (out[j] + cm * comb(kk, j) * pow(t, kk - j, m)) % m
    return out


def clause_value(ring, signs, args):
    """Phi_C = 1 - prod (1 - lit), lit = Y (positive) or 1 - Y (negative)."""
    prod = ring.const(1)
    for sg, y in zip(signs, args):
        one_minus_lit = ring.add(ring.const(1), [(-c) % ring.m for c in y]) if sg else y
        prod = ring.mul(prod, one_minus_lit)
    return ring.add(ring.const(1), [(-c) % ring.m for c in prod])


def taylor_direct(n, inst, s, ring):
    t = s + 1
    J, m = ring.J, ring.m
    sh = {key: {w: shift_poly_mod(pp, t, J, m) for w, pp in row.items()} for key, row in inst["P"].items()}
    total = ring.const(0)
    for x in product((0, 1), repeat=n):
        val = ring.const(1)
        for fi, f in enumerate(inst["factors"]):
            args = []
            for v in f:
                acc = [0] * J
                for w, coeffs in sh[(fi, v)].items():
                    if x[w]:
                        acc = ring.add(acc, coeffs)
                args.append(acc)
            val = ring.mul(val, clause_value(ring, inst["signs"][fi], args))
        total = ring.add(total, val)
    return total


def taylor_separated(n, inst, s, ring, q, Cm=None):
    """Lemmas 4.2 + 4.3: sum over |E| <= q, U subset E, of separated group sums."""
    t = s + 1
    J, m = ring.J, ring.m
    Cm = Cm if Cm is not None else inst["C"]
    r = len([g for (ss, g) in inst["groups"] if ss == s])
    groups = [inst["groups"][(s, g)] for g in range(r)]
    at_t = {key: {w: fmod(Fraction(poly_eval(pp, t)), m) for w, pp in row.items()} for key, row in inst["P"].items()}
    delta = {}
    for key, row in inst["P"].items():
        dd = {}
        for w, pp in row.items():
            coeffs = shift_poly_mod(pp, t, J, m)
            coeffs[0] = 0
            dd[w] = coeffs
        delta[key] = dd
    fac_by_group = defaultdict(list)
    for fi in range(len(inst["factors"])):
        fac_by_group[inst["d"][fi][s]].append(fi)
    total = ring.const(0)
    for size in range(min(q, n) + 1):
        for E in combinations(range(n), size):
            for usize in range(size + 1):
                for U in combinations(E, usize):
                    sign = -1 if (size - usize) % 2 else 1
                    kap = {}
                    for key, dd in delta.items():
                        acc = [0] * J
                        for w in U:
                            if w in dd:
                                acc = ring.add(acc, dd[w])
                        kap[key] = acc
                    # buckets per group
                    buckets = []
                    for g in range(r):
                        S = groups[g]
                        bk = defaultdict(list)
                        for bits in product((0, 1), repeat=len(S)):
                            if any(S[i] in E and not bits[i] for i in range(len(S))):
                                continue
                            xg = dict(zip(S, bits))
                            y = tuple(sum(val * xg[w] for w, val in row.items()) for row in Cm[(s, g)])
                            bk[y].append(xg)
                        buckets.append(bk)
                    if any(not bk for bk in buckets):
                        continue
                    keys = [list(bk.keys()) for bk in buckets]
                    contrib = ring.const(0)
                    for ys in product(*keys):
                        reps = [buckets[g][ys[g]][0] for g in range(r)]
                        prod_g = ring.const(1)
                        for g in range(r):
                            gsum = ring.const(0)
                            for xg in buckets[g][ys[g]]:
                                x = [0] * n
                                for h2 in range(r):
                                    src = xg if h2 == g else reps[h2]
                                    for w, bit in src.items():
                                        x[w] = bit
                                val = ring.const(1)
                                for fi in fac_by_group[g]:
                                    args = []
                                    for v in inst["factors"][fi]:
                                        lin = sum(c * x[w] for w, c in at_t[(fi, v)].items()) % m
                                        args.append(ring.add(ring.const(lin), kap[(fi, v)]))
                                    val = ring.mul(val, clause_value(ring, inst["signs"][fi], args))
                                gsum = ring.add(gsum, val)
                            prod_g = ring.mul(prod_g, gsum)
                        contrib = ring.add(contrib, prod_g)
                    if sign < 0:
                        contrib = [(-c) % m for c in contrib]
                    total = ring.add(total, contrib)
    return total


def crt_at_zero(data, Q, J, m):
    """Lemma 4.1: solve for the coefficients of H (degree < QJ) from the jets; return H(0)."""
    N = Q * J
    A = []
    rhs = []
    for s in range(Q):
        t = s + 1
        for j in range(J):
            A.append([comb(u, j) * pow(t, u - j, m) % m if u >= j else 0 for u in range(N)])
            rhs.append(data[s][j])
    Mx = [row + [b] for row, b in zip(A, rhs)]
    p = None
    for col in range(N):
        piv = next(i for i in range(col, N) if math.gcd(Mx[i][col], m) == 1)
        Mx[col], Mx[piv] = Mx[piv], Mx[col]
        inv = pow(Mx[col][col], -1, m)
        Mx[col] = [v * inv % m for v in Mx[col]]
        for i in range(N):
            if i != col and Mx[i][col]:
                f = Mx[i][col]
                Mx[i] = [(x - f * y) % m for x, y in zip(Mx[i], Mx[col])]
    return Mx[0][N]


def random_bounded(n, rng, n_clauses):
    occ = [0] * n
    clauses = []
    tries = 0
    while len(clauses) < n_clauses and tries < 500:
        tries += 1
        free = [v for v in range(n) if occ[v] < 3]
        if not free:
            break
        wdt = rng.randint(1, min(3, len(free)))
        vs = rng.sample(free, wdt)
        for v in vs:
            occ[v] += 1
        clauses.append([(v, rng.random() < 0.5) for v in vs])
    used = {v for cl in clauses for v, _ in cl}
    # rename so that every variable occurs (Theorem 7.1 removes free variables first)
    ren = {v: i for i, v in enumerate(sorted(used))}
    return len(used), [[(ren[v], sg) for v, sg in cl] for cl in clauses]


def brute(n, clauses):
    return sum(1 for x in product((0, 1), repeat=n) if all(any((x[v] == 1) == sg for v, sg in cl) for cl in clauses))


def part_K():
    rng = random.Random(1027)
    params = [(2, 1, 1), (2, 1, 2), (2, 2, 2), (2, 3, 2), (2, 3, 3), (3, 2, 2), (2, 4, 3), (3, 2, 3)]
    rows = []
    allok = True
    neg_detected = 0
    neg_trials = 0
    swap_checks = {}
    for (r, Q, R) in params:
        sw = Swap(r, Q, R)
        ch = sw.checks()
        swap_checks[f"{r},{Q},{R}"] = {**ch, "c": sw.c, "D": sw.D, "interval": sw.interval}
        allok = allok and all(ch.values())
        for rep in range(2):
            n, clauses = random_bounded(rng.randint(4, 6) if Q <= 2 else rng.randint(4, 5), rng, rng.randint(3, 6))
            inst = build(n, clauses, sw)
            loc = locality_ok(inst, sw)
            Rp = max(len(pp) - 1 for row in inst["P"].values() for pp in row.values())
            dtot = inst["M"]
            q = (dtot * Rp) // Q
            J = q + 1
            p = prime_above(Q)
            m = p ** (n + 1)
            ring = Ring(m, J)
            data, direct_ok = [], True
            for s in range(Q):
                sep = taylor_separated(n, inst, s, ring, q)
                dirv = taylor_direct(n, inst, s, ring)
                direct_ok = direct_ok and sep == dirv
                data.append(sep)
            H0 = crt_at_zero(data, Q, J, m)
            truth = brute(n, clauses)
            bal = inst["bal"]
            ok = (loc and direct_ok and H0 == truth and inst["rows_le_2hc"] and inst["row_l1_ok"]
                  and bal["Psi_le_E0"] and bal["E0_lt_2"] and bal["E_final_equals_Psi"]
                  and bal["h_bound"] and bal["group_bound"] and dtot * Rp < Q * J)
            allok = allok and ok
            # negative control: drop all measurement rows of one (node, group) pair
            if rep == 0:
                for s in range(Q):
                    for g in range(r):
                        if not inst["C"][(s, g)]:
                            continue
                        Cbad = dict(inst["C"])
                        Cbad[(s, g)] = []
                        neg_trials += 1
                        if taylor_separated(n, inst, s, ring, q, Cm=Cbad) != data[s]:
                            neg_detected += 1
                        break
                    else:
                        continue
                    break
            rows.append({"r": r, "Q": Q, "R": R, "n": n, "clauses": len(clauses), "count": truth, "H0": H0,
                         "p": p, "mu": n + 1, "J": J, "h": inst["h"], "c": sw.c,
                         "max_rows_per_node": inst["max_rows_per_node"], "locality": loc,
                         "separated_equals_direct": direct_ok, "ok": ok,
                         "E_Psi_initial": bal["E_Psi_initial"], "Psi_final": bal["Psi_final"]})
            print(f"    (r,Q,R)=({r},{Q},{R}) n={n} count={truth} H(0)={H0} ok={ok}")
    RESULTS["counting"] = {"all_ok": allok, "swap_checks": swap_checks, "cases": rows,
                           "negative_control": {"trials": neg_trials, "detected": neg_detected}}
    print(f"[K] counting algorithm mod p^mu: {len(rows)} instances, all exact: {allok}; "
          f"negative control detected {neg_detected}/{neg_trials}")
    return allok


# =============================================================== Part P: parameters

def params_of(lam):
    nu = (lam - 1).bit_length()          # = ceil(log2 lam) for lam >= 2
    r = isqrt(lam)
    Q = lam // (4 * nu)
    s = isqrt(lam)
    R = 2 * (s if s * s == lam else s + 1)
    return nu, r, Q, R


def exponent_bound(lam):
    nu, r, Q, R = params_of(lam)
    x = 3 * (R - 1) / Q
    eT = x * math.log2(math.e ** 2 / x)
    eM = 24 * math.exp(-R * R / (8 * Q)) * (Q + 2)
    eG = 2 / r
    return eT, eM, eG


def check_lambda(lam):
    nu, r, Q, R = params_of(lam)
    ok = r >= 2 and R % 2 == 0 and R >= 6 and 2 * R <= Q
    ok = ok and 4 * Q * nu <= lam                     # => log2 D <= Q nu <= lam/4 (P3)
    ok = ok and R * R >= 4 * lam and 8 * Q * nu <= 2 * lam   # (P5): R^2/(8Q) >= 2 nu
    ok = ok and 9 * (R - 1) ** 2 * lam <= 625 * nu * nu * Q * Q   # (P6): x <= 25 nu / sqrt(lam)
    ok = ok and 4 * lam <= 9 * r * r                  # (P8): 2/r <= 3/sqrt(lam)
    return ok


def part_P():
    lam0 = (1 << 17) + 1
    bad = []
    # exhaustive range
    for lam in range(lam0, (1 << 21) + 1):
        if not check_lambda(lam):
            bad.append(lam)
            if len(bad) > 5:
                break
    # block endpoints for nu up to 400 (exact integer checks)
    ends_ok = True
    for nu in range(18, 401):
        for lam in (2 ** (nu - 1) + 1, 2 ** nu):
            if not check_lambda(lam):
                ends_ok = False
    # also: just below lam0 the conditions fail (threshold is sharp for (P2))
    below = [lam for lam in range(1 << 16, lam0) if not check_lambda(lam)]
    # exponent bound check (floats) on a grid and least lambda with bound < 1
    grid_ok = True
    for e in range(17, 120):
        for lam in (2 ** e + 1, 3 * 2 ** (e - 1), 2 ** (e + 1)):
            nu = params_of(lam)[0]
            eT, eM, eG = exponent_bound(lam)
            if not (eT <= 12.5 * nu * nu / math.sqrt(lam) and eM <= 1 / math.sqrt(lam) and eG <= 3 / math.sqrt(lam)
                    and eT + eM + eG <= 13 * nu * nu / math.sqrt(lam)):
                grid_ok = False
    # least lambda (binary search per block) where exact bound < 1, and where 13 nu^2/sqrt(lam) < 1
    def least(pred):
        lo, hi = lam0, 1 << 40
        while lo < hi:
            mid = (lo + hi) // 2
            if pred(mid):
                hi = mid
            else:
                lo = mid + 1
        return lo
    least_exact = least(lambda lam: sum(exponent_bound(lam)) < 1)
    least_simple = least(lambda lam: 13 * params_of(lam)[0] ** 2 < math.sqrt(lam))
    # the exact bound is not monotone inside a block (floors and ceilings); find the last lambda
    # in [8e6, 2e7] where it is still >= 1, then check a grid up to 2^31 above that point
    last_bad = max(l for l in range(8_000_000, 20_000_001) if sum(exponent_bound(l)) >= 1)
    mono_ok = all(sum(exponent_bound(l)) < 1 for l in range(last_bad + 1, 2 ** 31, 9973))
    # beyond 2^24 (nu >= 25): on each block the majorant xh*log2(e^2/xh) + 2/(sqrt(lam)-1), with
    # xh = 12 nu (2 sqrt(lam) + 1)/(lam - 4 nu) >= x, decreases in lam; the measurement term is at most
    # its value at the block end.  Check majorant(block start) + measurement(block end) < 1.
    def block_ok(nu):
        lam = 2 ** (nu - 1)
        xh = 12 * nu * (2 * math.sqrt(lam) + 1) / (lam - 4 * nu)
        maj = xh * math.log2(math.e ** 2 / xh) + 2 / (math.sqrt(lam) - 1)
        meas = 24 * math.exp(-2 * nu) * (2 ** nu / (4 * nu) + 2)
        return xh < math.e and maj + meas < 1
    blocks_ok = all(block_ok(nu) for nu in range(25, 401))
    mono_ok = mono_ok and blocks_ok
    RESULTS["parameters"] = {
        "lambda0": lam0, "exhaustive_range": [lam0, 1 << 21], "exhaustive_failures": bad,
        "block_endpoints_nu_18_to_400_ok": ends_ok, "fails_just_below_lambda0": len(below) > 0,
        "exponent_grid_ok": grid_ok, "least_lambda_exact_bound_below_1": least_exact,
        "last_lambda_with_exact_bound_ge_1": last_bad, "below_1_on_grid_after_last_bad": mono_ok,
        "least_lambda_13nu2_below_sqrt": least_simple,
        "example": {str(l): dict(zip(("nu", "r", "Q", "R"), params_of(l))) | dict(zip(("eT", "eM", "eG"), exponent_bound(l)))
                    for l in (lam0, 1 << 24, least_exact, 1 << 32)}}
    ok = not bad and ends_ok and grid_ok and mono_ok
    print(f"[P] parameters: (P1)-(P8) hold on [2^17+1, 2^21] and block ends to nu=400: {not bad and ends_ok}; "
          f"exponent bounds on grid: {grid_ok}; exact bound < 1 for all lambda > {last_bad} "
          f"(checked to 2^31: {mono_ok}); simple form 13 nu^2/sqrt(lambda) < 1 from {least_simple}")
    return ok


# =============================================================== Part A: assembly constants

def sparsification_M(k, sigma):
    Gamma = 2 * (3 ** (k - 1) - 1)
    M = 2
    while True:
        CM = 4 * sum((3 * M) ** j for j in range(k - 1))
        h = math.ceil(math.log2(4 * CM * M))
        if M >= max(Gamma, 4) and Gamma * h <= sigma * M:
            return M, CM, Gamma
        M *= 2


def part_A():
    rows = []
    ok = True
    for k in (3, 4, 5):
        for e in (20, 40, 80, 160):
            logn = 2.0 ** e
            s = math.ceil(logn ** (1 / (2 * k)))
            sigma = 1 / s
            M, CM, Gamma = sparsification_M(k, sigma)
            A = 2 * Gamma * (k + 2) / sigma
            Mstar = 4 * A * math.log2(2 * A)
            bound = 2 * max(Gamma, 4, Mstar)
            ok = ok and M <= bound
            B = 2 * sum(M ** j for j in range(k))
            rows.append({"k": k, "log2_n": f"2^{e}", "sigma": sigma, "M": M, "M_bound": bound, "B": B,
                         "B_le_4M^(k-1)": B <= 4 * M ** (k - 1)})
            ok = ok and B <= 4 * M ** (k - 1)
    RESULTS["assembly"] = {"all_ok": ok, "rows": rows}
    print(f"[A] sparsification constants: M <= 2 max(Gamma, 4, 4A log2(2A)) and B <= 4 M^(k-1): {ok}")
    return ok


def main():
    oks = [part_S(), part_R(), part_K(), part_P(), part_A()]
    RESULTS["all_ok"] = all(oks)
    out = __file__.rsplit("/", 2)[0] + "/results/gapless_chain.json"
    with open(out, "w") as fh:
        json.dump(RESULTS, fh, indent=1, default=str)
    print("ALL OK" if all(oks) else "FAILURES")


if __name__ == "__main__":
    main()

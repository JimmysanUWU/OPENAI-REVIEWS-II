"""Checks for 18_AGGREGATE_TRANSFER.md: GPT's transfer finding, repaired under the algorithm's own schedules.

Part T  [Lemma 18.1, Theorem 18.2] On actual incidence structures (n variables, clause width <= 3, at most
        3 occurrences, r = 2, Q = 2 and 3), with schedules whose label occupancies satisfy 08 Lemma 6.2's
        variance bound V1 < 2M(b+w) and 2(b+w)L^2 <= M:
          - no label is empty and at most 3 labels have c_u < M/(2L);
          - support injection: inside one copy no variable meets two labels of one eigenvector's support,
            so every filled group row has l1-norm = number of occupied support labels in that group;
          - every eigenvector avoiding the under-filled labels gives, in each of the first ceil(M/(2L))
            copies and at every node, a nonzero group row of norm >= 2^max(j-1,0);
          - the aggregate count c_*(S) >= ceil(M/(2L)) * (|S| - #{psi in S meeting an under-filled label})
            for the measured sets S = {j >= j0}.
        Run on balanced random schedules and on adversarial schedules with three labels pushed below M/(2L)
        while the variance bound still holds.  GPT's stronger hypothesis M >= 8(b+w)L^2 fails on all of them.
        Controls: (a) four under-filled labels violate the variance bound (so "3" is the right count);
        (b) on the adversarial schedules some eigenvector does lose rows in late copies (the per-copy premise
        of 07 E4 really fails inside the hypotheses, as GPT's assessment says).
Part F  [Corollary 18.3, Theorem 18.4] the arithmetic that lets the additive loss 3 * 2^Q / L be absorbed in
        14 Theorem 3.2 and in GPT 16 Theorem 3.1.
Part S  [independent re-check of GPT 15 and 17, fresh code, arithmetic mod a large prime and over Q]
          - full-cube spline module and trace dimensions (GPT 15 Lemma 1.2), Q <= 4;
          - parity completion on greedy maximal good-edge graphs: #deleted >= 2^Q tau_D(Q), and the degree
            d = floor((Q+D)/2) trace space is everything (GPT 15 Theorem 2.1); control at degree d - 1;
          - the single-bad-edge design at degree Q-1 (GPT 15 Proposition 2.4), exactly, Q = 2..6;
          - graded designs at degree Q-1 need >= Q bad edges and Q is attained (so ungraded designs beat
            graded ones by the factor Q in total count, while the largest direction is the same);
          - the uniform defect energy rho/(1-rho) (GPT 17 Theorem 4.1) by exact normal equations, Q <= 12.

Standard library only; deterministic; writes ../results/aggregate_transfer.json (about 20 s).
"""

import json
import math
import os
import random
import time
from fractions import Fraction
from itertools import product

RESULTS = {}
P = 2 ** 61 - 1          # prime for modular ranks


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)


# ======================================================================= Part T: transfer

def local_vectors(r):
    """Local vectors on [r]x[r]: (dict pair->coeff, kind, indices)."""
    out = []
    for a in range(r):
        out.append(({(a, a): 1}, "diag", (a,)))
    for a in range(r):
        for b in range(a + 1, r):
            out.append(({(a, b): 1, (b, a): 1}, "sym", (a, b)))
            out.append(({(a, b): 1, (b, a): -1}, "anti", (a, b)))
    return out


def eigenvectors(r, Q):
    """All r^(2Q) tensor eigenvectors: support dict label->coeff, j (antisymmetric count), local data."""
    locs = local_vectors(r)
    eig = []
    for combo in product(range(len(locs)), repeat=Q):
        support = {(): 1}
        for k in combo:
            vec = locs[k][0]
            support = {lab + (pair,): c * cv for lab, c in support.items() for pair, cv in vec.items()}
        j = sum(locs[k][1] == "anti" for k in combo)
        eig.append(dict(support=support, j=j, idx=[locs[k][2] for k in combo]))
    return eig


def make_instance(n, seed):
    """Incidences of a CNF with clause width <= 3 and <= 3 occurrences per variable."""
    rng = random.Random(seed)
    slots = [v for v in range(n) for _ in range(3)]
    rng.shuffle(slots)
    clauses = []
    for i in range(0, len(slots), 3):
        clauses.append(sorted(set(slots[i:i + 3])))
    inc = [(f, v) for f, cl in enumerate(clauses) for v in cl]
    occ = [0] * n
    for _, v in inc:
        occ[v] += 1
    width = max(len(c) for c in clauses)
    return clauses, inc, max(occ), width


def labels_of(inc, a, d, Q):
    return [tuple((a[v][s], d[f][s]) for s in range(Q)) for f, v in inc]


def occupancy(labs, all_labels):
    c = {u: 0 for u in all_labels}
    for u in labs:
        c[u] += 1
    return c


def variance(c, M, L):
    return sum((Fraction(x) - Fraction(M, L)) ** 2 for x in c.values())


def designed_instance(r, Q, n_targets, seed):
    """A CNF (width <= 3, <= 3 occurrences) whose label occupancy is designed: n_targets labels sit just
    below M/(2L), the others are balanced, and M = 2(b+w)L^2 with b = w = 3 (the extreme allowed by 08
    Lemma 6.2).  Returns incidences with explicit schedules."""
    L = r ** (2 * Q)
    all_labels = list(product(product(range(r), repeat=2), repeat=Q))
    M = 2 * 6 * L * L
    rng = random.Random(seed)
    targets = sorted(all_labels, key=lambda u: rng.random())[:n_targets]
    low = math.ceil(Fraction(M, 2 * L)) - 1
    c = {u: low for u in targets}
    rest = M - low * n_targets
    others = [u for u in all_labels if u not in c]
    for i, u in enumerate(others):
        c[u] = rest // len(others) + (1 if i < rest % len(others) else 0)
    check(sum(c.values()) == M, "occupancy sums to M")
    # incidences (label, copy); variables group by first schedule, clauses by second schedule
    A = lambda u: tuple(p[0] for p in u)
    Dd = lambda u: tuple(p[1] for p in u)
    incid = [(u, k) for u in all_labels for k in range(c[u])]
    by_a, by_d = {}, {}
    for idx, (u, k) in enumerate(incid):
        by_a.setdefault(A(u), []).append(idx)
        by_d.setdefault(Dd(u), []).append(idx)
    var = [None] * len(incid)
    nvar = 0
    for key in sorted(by_a):
        lst = sorted(by_a[key], key=lambda i: (incid[i][1], incid[i][0]))   # cycle through labels
        for s in range(0, len(lst), 3):
            for i in lst[s:s + 3]:
                var[i] = nvar
            nvar += 1
    clause = [None] * len(incid)
    ncl = 0
    for key in sorted(by_d):
        lst = sorted(by_d[key], key=lambda i: (incid[i][0], incid[i][1]))   # by label, then copy
        for s in range(0, len(lst), 3):
            for i in lst[s:s + 3]:
                clause[i] = ncl
            ncl += 1
    members = {}
    for i in range(len(incid)):
        members.setdefault(clause[i], []).append(var[i])
    check(all(len(set(v)) == len(v) for v in members.values()), "a clause repeats a variable")
    occ = [0] * nvar
    for v in var:
        occ[v] += 1
    b, w = max(occ), max(len(v) for v in members.values())
    check(b <= 3 and w <= 3, "occurrence and width bounds")
    labs = [u for u, k in incid]
    inc = [(clause[i], var[i]) for i in range(len(incid))]
    return inc, labs, nvar, b, w, c, all_labels


def transfer_check(r, Q, rho_M, seed, n_targets):
    """Balanced run (n_targets = 0): a random CNF with M ~ rho_M * 2(b+w)L^2 and random schedules satisfying
    the variance bound.  Adversarial run: a designed CNF at M = 2(b+w)L^2 with n_targets under-filled labels."""
    L = r ** (2 * Q)
    all_labels = list(product(product(range(r), repeat=2), repeat=Q))
    if n_targets == 0:
        n = math.ceil(rho_M * 2 * 6 * L * L / 3)
        rng = random.Random(seed + 1)
        while True:
            clauses, inc, b, w = make_instance(n, seed)
            if len(inc) >= 2 * (b + w) * L * L:
                break
            n += 1
        M = len(inc)
        for attempt in range(200):
            a = [[rng.randrange(r) for _ in range(Q)] for _ in range(n)]
            d = [[rng.randrange(r) for _ in range(Q)] for _ in range(len(clauses))]
            labs = labels_of(inc, a, d, Q)
            c = occupancy(labs, all_labels)
            if variance(c, M, L) < 2 * M * (b + w):
                break
        else:
            raise AssertionError("no balanced schedule found")
    else:
        inc, labs, n, b, w, c, _ = designed_instance(r, Q, n_targets, seed)
        M = len(inc)
    check(2 * (b + w) * L * L <= M and M >= n, "hypotheses 2(b+w)L^2 <= M and M >= n")
    V = variance(c, M, L)
    check(V < 2 * M * (b + w), "variance bound of 08 Lemma 6.2")
    # copies: number incidences with a given label consecutively
    seen = {u: 0 for u in all_labels}
    pos = {}                      # (copy, label) -> variable
    for (f, v), u in zip(inc, labs):
        pos[(seen[u], u)] = v
        seen[u] += 1
    # support injection precondition: a variable's labels share one first-coordinate schedule
    first = {}
    for (f, v), u in zip(inc, labs):
        sched = tuple(p[0] for p in u)
        check(first.setdefault(v, sched) == sched, "a variable with two first schedules")
    h = max(c.values())
    half = Fraction(M, 2 * L)
    deficient = [u for u in all_labels if c[u] < half]
    check(min(c.values()) >= 1, "an empty label under the variance bound")
    check(len(deficient) <= 3, "more than three under-filled labels")
    check(8 * (b + w) * L * L > M, "this instance was meant to violate GPT's stronger hypothesis")
    k_good = math.ceil(M / (2 * L))
    eig = eigenvectors(r, Q)
    check(len(eig) == L, "eigenbasis size")
    deficient_set = set(deficient)
    partial_late_rows = 0
    lost_row_from_underfilled = False
    per_j = {}
    for t in range(Q):
        counts = {j0: 0 for j0 in range(Q + 1)}
        affected = {j0: 0 for j0 in range(Q + 1)}
        sizes = {j0: 0 for j0 in range(Q + 1)}
        for psi in eig:
            meets = any(y in deficient_set for y in psi["support"])
            j = psi["j"]
            o = sum(len(ix) == 2 for ix in psi["idx"])
            groups = {g: [(y, cf) for y, cf in psi["support"].items() if y[t][0] == g] for g in range(r)}
            nonzero_rows = 0
            full_rows_in_good_copies = 0
            for kappa in range(h):
                for g in range(r):
                    row = {}
                    for y, coef in groups[g]:
                        v = pos.get((kappa, y))
                        if v is None:
                            continue
                        check(v not in row, "support injection fails: one variable meets two support labels")
                        row[v] = coef
                    norm = sum(abs(x) for x in row.values())
                    check(norm == len(row), "cancellation inside a filled row")
                    if norm:
                        nonzero_rows += 1
                    if groups[g] and norm < len(groups[g]):
                        partial_late_rows += 1
                        if meets and kappa < k_good:
                            lost_row_from_underfilled = True
                    if kappa < k_good and g == psi["idx"][t][0] and not meets:
                        check(norm == len(groups[g]), "a good copy misses a support label")
                        check(norm >= 2 ** max(j - 1, 0), "full-copy row norm below 2^(j-1)")
                        check(norm in (2 ** o, 2 ** max(o - 1, 0)), "full-copy row norm is 2^o or 2^(o-1)")
                        full_rows_in_good_copies += 1
            if not meets:
                check(full_rows_in_good_copies == k_good, "a clean eigenvector misses a good copy")
            for j0 in range(Q + 1):
                if j >= j0:
                    counts[j0] += nonzero_rows
                    sizes[j0] += 1
                    affected[j0] += meets
        for j0 in range(Q + 1):
            lower = k_good * (sizes[j0] - affected[j0])
            check(counts[j0] >= lower, "aggregate count below the lemma's bound")
            per_j.setdefault(str(t), {})[str(j0)] = dict(c_star=counts[j0], lemma_lower=lower, size=sizes[j0],
                                                          meeting_underfilled=affected[j0])
    check(partial_late_rows > 0, "control: some late-copy row should be partial (per-copy premise fails)")
    if n_targets:
        check(len(deficient) == n_targets, "adversarial run should have the requested under-filled labels")
        check(lost_row_from_underfilled, "an eigenvector meeting an under-filled label should lose a row early")
    return dict(r=r, Q=Q, n=n, M=M, labels=L, b=b, w=w, h=h, copies_guaranteed=k_good,
                rho_M=M / (2 * (b + w) * L * L), variance=float(V), variance_bound=2 * M * (b + w),
                gpt_stronger_hypothesis=8 * (b + w) * L * L <= M,
                min_occupancy=min(c.values()), underfilled=len(deficient),
                partial_rows_in_late_copies=partial_late_rows, aggregate=per_j)


def four_underfilled_control():
    """Four labels below M/(2L) force V1 >= M^2/L^2 >= 2M(b+w) whenever 2(b+w)L^2 <= M."""
    out = []
    for L, b, w in [(16, 3, 3), (64, 3, 3), (81, 2, 3)]:
        M = 2 * (b + w) * L * L
        half = Fraction(M, 2 * L)
        below = math.ceil(half) - 1
        lower = 4 * (Fraction(M, L) - below) ** 2
        check(lower >= 2 * M * (b + w), "four under-filled labels compatible with the variance bound")
        three = 3 * (Fraction(M, L) - below) ** 2
        check(three < 2 * M * (b + w), "three under-filled labels should be compatible")
        out.append(dict(L=L, M=M, four_labels_min_variance=float(lower), three_labels_variance=float(three),
                        bound=2 * M * (b + w)))
    return out


def part_T():
    t0 = time.time()
    runs = [transfer_check(2, 2, 2.0, 11, 0), transfer_check(2, 2, 1.0, 12, 3),
            transfer_check(2, 3, 1.0, 14, 3)]
    RESULTS["T"] = dict(runs=runs, four_underfilled_control=four_underfilled_control(), seconds=time.time() - t0)
    print("Part T: transfer lemma on 3 instances (2 designed with 3 under-filled labels): OK")


# ======================================================================= Part F: arithmetic

def part_F():
    rows = []
    # 14 Theorem 3.2 regime: Q >= 64, 2 <= R <= Q/8 + 1, eps' = e^{-2(R+1)^2/Q}/(4(Q+1)) >= 6 * 2^-Q
    for Q in range(64, 4001):
        R = Q // 8 + 1
        lhs = -2 * (R + 1) ** 2 / Q - math.log(4 * (Q + 1))
        rhs = math.log(6) - Q * math.log(2)
        check(lhs - rhs > 0, "eps' >= 6 * 2^-Q fails at Q=%d" % Q)
        check(Q / (512 * (Q + 1)) >= 1 / 520, "1/520 constant")
    rows.append("14 regime: eps'(Q, Q/8+1) >= 6*2^-Q for 64 <= Q <= 4000; margin grows linearly in Q")
    # GPT 16 Theorem 3.1 regime: Q >= 512, D <= Q/16: (3/(8Q)) e^{-25-4D^2/Q} >= 6 * 2^-Q
    for Q in range(512, 6001):
        D = Q // 16
        lhs = math.log(3 / (8 * Q)) - 25 - 4 * D * D / Q
        check(lhs - (math.log(6) - Q * math.log(2)) > 0, "heavy-defect mass >= 6*2^-Q fails at Q=%d" % Q)
    rows.append("GPT 16 regime: (3/(8Q)) e^{-25-4D^2/Q} >= 6*2^-Q for 512 <= Q <= 6000, D = Q/16")
    # step 4 of 14 Theorem 3.2 with the new constant: ln(1/(520 delta)) >= (1/2) ln(1/delta) iff delta <= 520^-2
    check(520 ** 2 <= 2 ** 19, "520^2 <= 2^19")
    # the R = 1 case: (1/2)(1/4 - 6/r^2) >= 1/16 once r >= 7
    check(Fraction(1, 2) * (Fraction(1, 4) - Fraction(6, 49)) >= Fraction(1, 16), "R = 1 case")
    # (2/r^2)^Q <= min(2^-Q, 2/r^2)
    for r in range(2, 12):
        for Q in range(1, 30):
            x = Fraction(2, r * r) ** Q
            check(x <= Fraction(1, 2 ** Q) and x <= Fraction(2, r * r), "2^Q/L bounds")
    RESULTS["F"] = rows
    print("Part F: constants of Corollary 18.3 and Theorem 18.4: OK")


# ======================================================================= Part S: independent re-checks

def rank_mod(rows, ncols):
    rows = [[x % P for x in row] for row in rows]
    rk = 0
    col = 0
    while rk < len(rows) and col < ncols:
        piv = next((i for i in range(rk, len(rows)) if rows[i][col]), None)
        if piv is None:
            col += 1
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][col], P - 2, P)
        rows[rk] = [x * inv % P for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][col]:
                f = rows[i][col]
                rows[i] = [(x - f * y) % P for x, y in zip(rows[i], rows[rk])]
        rk += 1
        col += 1
    return rk


def cube_edges(Q):
    return [(x, x | (1 << t), t) for t in range(Q) for x in range(1 << Q) if not (x >> t) & 1]


def edge_row(Q, d, nodes, e):
    x, y, t = e
    row = [0] * ((1 << Q) * (d + 1))
    for k in range(d + 1):
        row[x * (d + 1) + k] = pow(nodes[t], k, P)
        row[y * (d + 1) + k] = -pow(nodes[t], k, P)
    return row


def trace_dim(Q, d, nodes, kept):
    n = (1 << Q) * (d + 1)
    rows = [edge_row(Q, d, nodes, e) for e in kept]
    r0 = rank_mod(rows, n)
    tr = []
    for x in range(1 << Q):
        row = [0] * n
        row[x * (d + 1)] = 1
        tr.append(row)
    return rank_mod(rows + tr, n) - r0, n - r0


def tau_count(Q, D):
    return sum(math.comb(Q, j) for j in range(Q + 1) if 2 * j > Q + D)


def parity_feasible(Q, D, nodes, kept):
    """Is there a degree-<=D spline on the kept edges with f_x(0) = parity(x)?  (consistency mod P)"""
    n = (1 << Q) * (D + 1)
    rows = [edge_row(Q, D, nodes, e) + [0] for e in kept]
    for x in range(1 << Q):
        row = [0] * (n + 1)
        row[x * (D + 1)] = 1
        row[n] = (-1) ** bin(x).count("1") % P
        rows.append(row)
    return rank_mod([row[:n] for row in rows], n) == rank_mod(rows, n + 1)


def part_S():
    out = {}
    nodes_all = [3, 5, 7, 11, 13, 17]
    # (1) full-cube module and trace dimensions
    dims = []
    for Q in range(1, 5):
        nodes = nodes_all[:Q]
        for d in range(Q + 2):
            tdim, mdim = trace_dim(Q, d, nodes, cube_edges(Q))
            check(tdim == sum(math.comb(Q, j) for j in range(min(d, Q) + 1)), "trace dim")
            check(mdim == sum(math.comb(Q, s) * max(0, d - s + 1) for s in range(Q + 1)), "module dim")
            dims.append((Q, d, mdim, tdim))
    out["full_cube_dims"] = len(dims)
    # (2) parity completion on greedy maximal good-edge graphs
    rng = random.Random(7)
    comp = []
    for Q in range(2, 5):
        nodes = nodes_all[:Q]
        for D in range(0, Q):
            best = None
            for trial in range(6 if Q < 4 else 3):
                E = cube_edges(Q)
                rng.shuffle(E)
                kept = []
                for e in E:
                    if parity_feasible(Q, D, nodes, kept + [e]):
                        kept.append(e)
                deleted = len(E) - len(kept)
                check(deleted >= tau_count(Q, D), "parity completion bound violated")
                d = (Q + D) // 2
                tdim, _ = trace_dim(Q, d, nodes, kept)
                check(tdim == 1 << Q, "trace space at degree floor((Q+D)/2) is not everything")
                if d - 1 >= 0 and tau_count(Q, D) > 0:
                    tdim2, _ = trace_dim(Q, d - 1, nodes, kept)
                    # control: at one degree less, completion need not hold (record, do not require)
                    comp.append(dict(Q=Q, D=D, deleted=deleted, bound=tau_count(Q, D), trace_d=tdim,
                                     trace_d_minus_1=tdim2))
                best = deleted if best is None else min(best, deleted)
    check(any(c["trace_d_minus_1"] < (1 << c["Q"]) for c in comp), "control: degree d-1 never fails")
    out["parity_completion_cases"] = comp
    # (3) single bad edge at degree Q - 1 (exact over Q)
    single = []
    for Q in range(2, 7):
        nodes = [Fraction(i + 1) for i in range(Q)]
        k = Q - 1

        def poly_mul(p, q):
            res = [Fraction(0)] * (len(p) + len(q) - 1)
            for i, x in enumerate(p):
                for j, y in enumerate(q):
                    res[i + j] += x * y
            return res

        def ev(p, t):
            return sum(c * t ** i for i, c in enumerate(p))

        def g(S, x):
            poly = [Fraction(1)]
            for i in range(Q):
                if (S >> i) & 1:
                    poly = poly_mul(poly, [Fraction(1), -1 / nodes[i]])
            sign = (-1) ** bin(S & x).count("1")
            return [sign * c for c in poly]

        hpoly = [Fraction(1)]
        for i in range(Q):
            if i != k:
                hpoly = poly_mul(hpoly, [Fraction(1), -1 / nodes[i]])
        full = (1 << Q) - 1
        n_x = {}
        for x in range(1 << Q):
            f = [Fraction(0)] * Q
            if x == 0:
                f = [Fraction(2 ** Q) * c for c in hpoly] + [Fraction(0)] * (Q - len(hpoly))
            for S in range(full):            # proper subsets
                gs = g(S, x)
                for i, c in enumerate(gs):
                    f[i] -= c
            par = (-1) ** bin(x).count("1")
            n_x[x] = [par * c for c in f]
            check(n_x[x][0] == 1, "n_x(0) = 1")
            check(len(n_x[x]) <= Q, "degree <= Q-1")
        bad = [(x, y, t) for x, y, t in cube_edges(Q) if ev(n_x[x], nodes[t]) + ev(n_x[y], nodes[t]) != 0]
        check(bad == [(0, 1 << k, k)], "single bad edge at {0, e_k}")
        c_empty = [sum(n_x[x][i] for x in range(1 << Q)) / 2 ** Q for i in range(Q)]
        check(any(c_empty[1:]), "identity coefficient nonconstant (ungraded)")
        single.append(Q)
    out["single_bad_edge_Q"] = single
    # (4) graded designs at degree Q-1: each direction needs a bad edge; Q is attained
    graded = []
    for Q in range(2, 6):
        total = 0
        for t in range(Q):
            # G_t = parity except on the edge {0, e_t}, where G_t * parity = -(2^{Q-1} - 1)
            G = {}
            for x in range(1 << Q):
                par = (-1) ** bin(x).count("1")
                G[x] = Fraction(par)
                if x in (0, 1 << t):
                    G[x] = Fraction(-(2 ** (Q - 1) - 1) * par)
            top = sum(G[x] * (-1) ** bin(x).count("1") for x in range(1 << Q))
            check(top == 0, "G_t has degree <= Q-1")
            check(all(G[x ^ (1 << t)] == -G[x] for x in range(1 << Q)), "G_t odd in x_t")
            E_t = [x for x in range(1 << Q) if G[x] != (-1) ** bin(x).count("1")]
            check(len(E_t) == 2, "one bad edge in direction t")
            total += len(E_t) // 2
        # lower bound: Smolensky with D = Q-1 gives |E_t| >= 1, flip invariance gives >= 2 points = 1 edge
        check(tau_count(Q, Q - 1) == 1, "tau count at D = Q-1")
        graded.append(dict(Q=Q, graded_min_total_edges=total, ungraded_min_total_edges=1))
    out["graded_vs_ungraded_top_degree"] = graded
    # (5) uniform defect energy (GPT 17 Theorem 4.1), exact normal equations
    energy = []
    for Q in range(1, 13):
        for D in range(0, Q):
            # minimize sum_{i=1}^Q q(i)^2 with q(0) = 1, q = 1 + sum_{k=1}^D c_k T^k
            A = [[Fraction(sum(i ** (a + b) for i in range(1, Q + 1))) for b in range(1, D + 1)]
                 for a in range(1, D + 1)]
            rhs = [-Fraction(sum(i ** a for i in range(1, Q + 1))) for a in range(1, D + 1)]
            # solve A c = rhs
            m = len(A)
            M_ = [row[:] + [rhs[i]] for i, row in enumerate(A)]
            for col in range(m):
                piv = next(i for i in range(col, m) if M_[i][col] != 0)
                M_[col], M_[piv] = M_[piv], M_[col]
                pv = M_[col][col]
                M_[col] = [x / pv for x in M_[col]]
                for i in range(m):
                    if i != col and M_[i][col] != 0:
                        f = M_[i][col]
                        M_[i] = [x - f * y for x, y in zip(M_[i], M_[col])]
            c = [M_[i][m] for i in range(m)]
            val = sum((1 + sum(ck * Fraction(i) ** (k + 1) for k, ck in enumerate(c))) ** 2 for i in range(1, Q + 1))
            rho = Fraction(1)
            for j in range(D + 1):
                rho *= Fraction(Q - j, Q + j + 1)
            check(val == rho / (1 - rho), "energy formula rho/(1-rho)")
            energy.append((Q, D))
    out["energy_cases"] = len(energy)
    RESULTS["S"] = out
    print("Part S: independent re-check of GPT 15 (dims, completion, single edge, graded vs ungraded) "
          "and GPT 17 (energy): OK")


def main():
    t0 = time.time()
    part_T()
    part_F()
    part_S()
    RESULTS["seconds"] = round(time.time() - t0, 1)
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "results", "aggregate_transfer.json"), "w") as fh:
        json.dump(RESULTS, fh, indent=1, default=str)
    print("ALL OK")


if __name__ == "__main__":
    main()

"""Checks for the distilled 'swap interpolation' formula (DISTILLATION.md, Theorem B).

With Krawtchouk coefficients gamma_k of f (f = sum_k gamma_k K_k, deg f < R) put

    M0'  = V0 * sum_k gamma_k e_k(Pi_1, ..., Pi_Q),
    P'(T) = I - M0' + V0 * sum_k gamma_k e_k(B_1(T), ..., B_Q(T)),   B_s(T) = (1 - T/s) Pi_s,

where e_k is the k-th elementary symmetric polynomial of commuting matrices.  We verify
    (a) M0' = M0 = V0 f(W0) of the paper (entrywise, exact);
    (b) P'(0) = I and deg_T P' <= R - 1;
    (c) for every t in [Q] and every row u, the entries of P'(t) - (I - M0') are t-local:
        nonzero only at y with a_y(t) = d_u(t);
    (d) consequently the kernel condition of Lemma 4.3 holds for P' (exact rank test);
    (e) the symmetric lower bound: for every f of degree < R, every window of R+1
        consecutive j in {0..Q} contains j with (-1)^j f(j) != 1 (finite-difference argument),
        tested on random rational f; and the interval construction's exact rank fraction is
        compared with the best window bound min_{j in window} Pr[W = j].
"""

import json
import random
import sys
from fractions import Fraction
from itertools import combinations, product
from math import comb

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from eth_core import Aux, rref, in_rowspace, poly_eval, poly_mul, binom_tail_outside  # noqa: E402


def krawtchouk(Q, k, j):
    return sum((-1) ** i * comb(j, i) * comb(Q - j, k - i) for i in range(0, k + 1))


def gammas(aux):
    """Solve f(j) = sum_{k<R} gamma_k K_k(j) for j = 0..R-1."""
    Q, R = aux.Q, aux.R
    fc = aux.f_coeffs()
    A = [[Fraction(krawtchouk(Q, k, j)) for k in range(R)] + [poly_eval(fc, j)] for j in range(R)]
    n = R
    for c in range(n):
        p = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[p] = A[p], A[c]
        A[c] = [x / A[c][c] for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[c])]
    return [A[k][n] for k in range(R)]


def distilled(aux):
    Q, R = aux.Q, aux.R
    g = gammas(aux)
    # rows: dict y -> polynomial in T (list of Fractions)
    M0p = {u: {} for u in aux.labels}
    Pp = {u: {} for u in aux.labels}
    for u in aux.labels:
        for k in range(R):
            if g[k] == 0:
                continue
            for S in combinations(range(Q), k):
                # V0 * prod_{s in S} Pi_s maps u to theta_{S^c}(u): (V0 Pi_S zeta)_u = zeta_{theta_all theta_S u}
                y = list(aux.swap_all(u))
                for s in S:
                    a, d = y[s]
                    y[s] = (d, a)
                y = tuple(y)
                M0p[u][y] = M0p[u].get(y, 0) + g[k]
                phi = [Fraction(1)]
                for s in S:
                    phi = poly_mul(phi, [Fraction(1), Fraction(-1, s + 1)])
                cur = Pp[u].get(y, [Fraction(0)])
                m = max(len(cur), len(phi))
                Pp[u][y] = [(cur[i] if i < len(cur) else 0) + g[k] * (phi[i] if i < len(phi) else 0) for i in range(m)]
        # I - M0'
        cur = Pp[u].get(u, [Fraction(0)])
        Pp[u][u] = [cur[0] + 1] + cur[1:]
        for y, m0 in M0p[u].items():
            cur = Pp[u].get(y, [Fraction(0)])
            Pp[u][y] = [cur[0] - m0] + cur[1:]
        Pp[u] = {y: p for y, p in Pp[u].items() if any(p)}
        M0p[u] = {y: v for y, v in M0p[u].items() if v != 0}
    return M0p, Pp


def check(r, Q, R):
    aux = Aux(r, Q, R)
    M0p, Pp = distilled(aux)
    res = {"r": r, "Q": Q, "R": R}
    res["M0_equal"] = all(M0p[u] == {y: v for y, v in aux.M0[u].items() if v != 0} for u in aux.labels)
    res["P0_identity"] = all({y: p[0] for y, p in Pp[u].items() if p[0] != 0} == {u: 1} for u in aux.labels)
    res["degree_le_R-1"] = all(len(p) - 1 <= R - 1 for u in aux.labels for p in Pp[u].values())
    loc = True
    for s in range(Q):
        t = s + 1
        for u in aux.labels:
            for y, p in Pp[u].items():
                val = poly_eval(p, t)
                const = (1 if y == u else 0) - M0p[u].get(y, 0)
                if val - const != 0 and y[s][0] != u[s][1]:
                    loc = False
    res["t_local_beyond_I-M0"] = loc
    ok_k = True
    for s in range(Q):
        t = s + 1
        for g in range(r):
            cols = set(aux.S_aux(s, g))
            basis = rref([{c: v for c, v in row.items() if c in cols} for row in aux.CB], aux.D)
            for u in aux.labels:
                if u[s][1] == g:
                    continue
                restr = {y: poly_eval(p, t) for y, p in Pp[u].items() if y in cols}
                restr = {k: v for k, v in restr.items() if v != 0}
                if not in_rowspace(basis, restr):
                    ok_k = False
    res["kernel_condition"] = ok_k
    res["differs_from_paper_rows"] = any(Pp[u] != aux.P[u] for u in aux.labels)
    return res


def window_lower_bound(Q, r, R, trials, rng):
    """(e): no degree<R polynomial matches (-1)^j on R+1 consecutive integers."""
    ok = True
    for _ in range(trials):
        A = rng.randint(0, Q - R)
        # interpolate (-1)^j on R points chosen inside the window, test the (R+1)-th
        pts = list(range(A, A + R + 1))
        drop = rng.randrange(R + 1)
        use = [p for i, p in enumerate(pts) if i != drop]
        # Lagrange value at the dropped point
        x0 = pts[drop]
        val = Fraction(0)
        for i, xi in enumerate(use):
            li = Fraction(1)
            for k, xk in enumerate(use):
                if k != i:
                    li *= Fraction(x0 - xk, xi - xk)
            val += (-1) ** xi * li
        ok &= val != (-1) ** x0
    p = Fraction(r - 1, 2 * r)
    pmf = [comb(Q, j) * p ** j * (1 - p) ** (Q - j) for j in range(Q + 1)]
    best_window = max(min(pmf[A:A + R + 1]) for A in range(0, Q - R + 1))
    aux_interval_tail = None
    mu = Fraction(Q * (r - 1), 2 * r)
    A0 = max(0, min(int(mu - Fraction(R - 1, 2)), Q - R + 1))
    tail = binom_tail_outside(Q, p, set(range(A0, A0 + R)))
    return {"Q": Q, "r": r, "R": R, "no_R+1_run": ok, "lower_bound_rank_fraction": float(best_window),
            "interval_construction_rank_fraction": float(tail),
            "ratio_upper_over_lower": float(tail / best_window) if best_window else None}


def main():
    out = {"cases": [check(r, Q, R) for r, Q, R in [(2, 2, 2), (2, 3, 2), (2, 3, 3), (2, 4, 2), (2, 4, 3), (3, 2, 2)]]}
    rng = random.Random(3)
    out["symmetric_optimality"] = [window_lower_bound(Q, r, R, 200, rng)
                                   for (Q, r, R) in [(20, 2, 6), (40, 2, 10), (80, 4, 14), (200, 8, 30), (400, 16, 44)]]
    ok = all(all(v for k, v in c.items() if isinstance(v, bool) and k != "differs_from_paper_rows") for c in out["cases"])
    ok &= all(x["no_R+1_run"] for x in out["symmetric_optimality"])
    out["all_ok"] = ok
    print(json.dumps(out, indent=1))
    print("ALL PASS" if ok else "SOME CHECK FAILED")
    path = __file__.rsplit("/", 2)[0] + "/results/distilled_swap_formula.json"
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()

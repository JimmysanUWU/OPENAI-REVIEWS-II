"""Independent re-check of GPT's ETH 21 (product-bias energy) and ETH 22 (minimal contradictions), written
without GPT's code: exact rational arithmetic, direct quadratic programs, brute-force circuit enumeration.

  T  [ETH 21 Theorem 21.1] the exact biased defect energy, computed as a direct quadratic program in all
     D * 2^Q vertex coefficients, equals E_S C_D(S; w/theta), with C_D evaluated by the Cauchy-Binet ratio
     Delta_{D+1}/Gamma_D (ETH 21 Lemma 21.2); theta^D U <= E <= U (Theorem 21.3); Example 21.6 = 3/20.
     Control: the wrong thinning law (theta_i = p_i) gives a different value.
  C  [ETH 22] random affine systems (M <= 6 equations, d <= 3 parameters, exact rationals):
     - Lemma 22.1: a subsystem is inconsistent iff it contains an inconsistent circuit (all subsets);
     - (22.4): the least number of failed equations equals the least hitting set of inconsistent circuits;
     - Lemma 22.2: a sign-compatible circuit decomposition of the optimal dual vector exists, its weights sum
       to at most w, and the circuit losses add up exactly to the full loss (squared loss);
     - Theorems 22.3 and 22.6: under independent retention, theta_(1)...theta_(c-1) L <= E L(V) <= L and
       kappa L <= E L(V);
     - Proposition 22.4: the Vandermonde system gives exactly theta^d L, and the parity law gives 0;
     - Example 22.8 (5/16 < 17/32) and the two rank-one values 503/1820 and 22/65.

Standard library only; deterministic; writes ../results/gpt_round2.json (about 5 s).
"""

import json
import os
import random
from fractions import Fraction as F
from itertools import combinations, product
from math import prod

RESULTS = {}


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)


# ------------------------------------------------------------------ exact linear algebra

def solve_psd(G, h):
    """Minimize c^T G c + 2 h^T c for PSD G: solve G c = -h (consistent), free variables set to 0."""
    n = len(G)
    M = [row[:] + [-h[i]] for i, row in enumerate(G)]
    piv_cols, r = [], 0
    for col in range(n):
        p = next((i for i in range(r, n) if M[i][col] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][col]
        M[r] = [x / pv for x in M[r]]
        for i in range(n):
            if i != r and M[i][col] != 0:
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv_cols.append(col)
        r += 1
    check(all(M[i][n] == 0 for i in range(r, n)), "normal equations consistent")
    c = [F(0)] * n
    for i, col in enumerate(piv_cols):
        c[col] = M[i][n]
    return c


def rank(rows):
    rows = [list(r) for r in rows]
    if not rows:
        return 0
    n = len(rows[0])
    rk = 0
    for col in range(n):
        p = next((i for i in range(rk, len(rows)) if rows[i][col] != 0), None)
        if p is None:
            continue
        rows[rk], rows[p] = rows[p], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][col] != 0:
                f = rows[i][col] / rows[rk][col]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk


def det(M):
    M = [list(map(F, r)) for r in M]
    n = len(M)
    d = F(1)
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c] != 0), None)
        if p is None:
            return F(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return d


# ------------------------------------------------------------------ Part T: ETH 21

def capacity_CB(S, D, t, v):
    """C_D(S; v) by Cauchy-Binet: Delta_{D+1}/Gamma_D; zero if |S| <= D."""
    S = list(S)
    if len(S) <= D:
        return F(0)
    vand = lambda J: prod((t[j] - t[i]) ** 2 for a, i in enumerate(J) for j in J[a + 1:])
    Delta = sum(prod(v[i] for i in J) * vand(J) for J in combinations(S, D + 1))
    Gamma = sum(prod(v[i] * t[i] ** 2 for i in J) * vand(J) for J in combinations(S, D))
    return Delta / Gamma


def biased_energy_qp(Q, D, t, p, w):
    """Direct minimization over all vertex coefficients c_{x,k} (k = 1..D)."""
    n = (2 ** Q) * D
    idx = lambda x, k: x * D + (k - 1)
    mu = lambda x: prod(p[i] if (x >> i) & 1 else 1 - p[i] for i in range(Q))
    G = [[F(0)] * n for _ in range(n)]
    h = [F(0)] * n
    const = F(0)
    for i in range(Q):
        for x in range(2 ** Q):
            y = x ^ (1 << i)
            g = {}
            for k in range(1, D + 1):
                g[idx(x, k)] = g.get(idx(x, k), 0) + F(t[i] ** k, 2)
                g[idx(y, k)] = g.get(idx(y, k), 0) + F(t[i] ** k, 2)
            wt = w[i] * mu(x)
            for a, ga in g.items():
                h[a] += wt * ga
                for b, gb in g.items():
                    G[a][b] += wt * ga * gb
            const += wt
    c = solve_psd(G, h)
    val = const + 2 * sum(hi * ci for hi, ci in zip(h, c)) + sum(c[a] * G[a][b] * c[b] for a in range(n) for b in range(n))
    return val


def part_T():
    rng = random.Random(21)
    cases = []
    for Q in (2, 3):
        for D in range(1, Q):
            for trial in range(3):
                t = [F(v) for v in rng.sample([1, 2, 3, -1, -2, 5], Q)]
                p = [F(rng.choice([1, 2, 3]), rng.choice([4, 5, 7])) for _ in range(Q)]
                w = [F(rng.choice([1, 2, 3])) for _ in range(Q)]
                theta = [4 * pi * (1 - pi) for pi in p]
                E_qp = biased_energy_qp(Q, D, t, p, w)
                E_thin = sum(prod(theta[i] if (S >> i) & 1 else 1 - theta[i] for i in range(Q))
                             * capacity_CB([i for i in range(Q) if (S >> i) & 1], D, t, [w[i] / theta[i] for i in range(Q)])
                             for S in range(2 ** Q))
                check(E_qp == E_thin, "Theorem 21.1: QP optimum = E_S C_D(S; w/theta)")
                U = capacity_CB(range(Q), D, t, w)
                tmin = min(theta)
                check(tmin ** D * U <= E_qp <= U, "Theorem 21.3: theta_min^D U <= E <= U")
                # control: the wrong thinning law (inclusion probability p_i instead of theta_i)
                E_wrong = sum(prod(p[i] if (S >> i) & 1 else 1 - p[i] for i in range(Q))
                              * capacity_CB([i for i in range(Q) if (S >> i) & 1], D, t, [w[i] / p[i] for i in range(Q)])
                              for S in range(2 ** Q))
                cases.append(dict(Q=Q, D=D, t=[str(x) for x in t], p=[str(x) for x in p], w=[str(x) for x in w],
                                  energy=str(E_qp), uniform=str(U), wrong_law_differs=E_wrong != E_qp))
    check(any(c["wrong_law_differs"] for c in cases), "control: the wrong thinning law must differ")
    ex = biased_energy_qp(2, 1, [F(1), F(2)], [F(1, 4), F(1, 4)], [F(1), F(1)])
    check(ex == F(3, 20), "Example 21.6: 3/20")
    check(capacity_CB([0, 1], 1, [F(1), F(2)], [F(1), F(1)]) == F(1, 5), "scalar optimum 1/5")
    RESULTS["T"] = dict(cases=cases, example_21_6=str(ex))
    print("Part T: ETH 21 Theorem 21.1 (direct QP = spectral thinning), Theorem 21.3, Example 21.6 = 3/20: OK")


# ------------------------------------------------------------------ Part C: ETH 22

def fit_L2(A, b, v):
    """min_beta sum v_i (b_i - a_i beta)^2 and the residual at a minimizer (exact)."""
    d = len(A[0])
    G = [[sum(v[i] * A[i][j] * A[i][k] for i in range(len(A))) for k in range(d)] for j in range(d)]
    h = [-sum(v[i] * A[i][j] * b[i] for i in range(len(A))) for j in range(d)]
    beta = solve_psd(G, h)
    r = [b[i] - sum(A[i][j] * beta[j] for j in range(d)) for i in range(len(A))]
    return sum(v[i] * r[i] ** 2 for i in range(len(A))), r


def consistent(A, b, T):
    if not T:
        return True
    return rank([A[i] for i in T]) == rank([A[i] + [b[i]] for i in T])


def circuits(A, b):
    """All row circuits with their kernel vectors and consistency flag."""
    M = len(A)
    out = []
    for size in range(1, len(A[0]) + 2):
        for J in combinations(range(M), size):
            if rank([A[i] for i in J]) != size - 1:
                continue
            if any(rank([A[i] for i in J if i != j]) != size - 1 for j in J):
                continue
            # kernel vector of A_J^T: solve sum g_i a_i = 0 with g_{J[0]} = 1
            d = len(A[0])
            rows = [[A[i][k] for i in J[1:]] + [-A[J[0]][k]] for k in range(d)]
            # least-squares-free exact solve: pick independent equations
            sol = solve_kernel(rows, size - 1)
            g = {J[0]: F(1)}
            for idx_, i in enumerate(J[1:]):
                g[i] = sol[idx_]
            check(all(sum(g[i] * A[i][k] for i in J) == 0 for k in range(d)), "kernel vector")
            out.append(dict(J=J, g=g, inconsistent=sum(g[i] * b[i] for i in J) != 0))
    return out


def solve_kernel(rows, n):
    M = [list(map(F, r)) for r in rows]
    piv, rk = [], 0
    for col in range(n):
        p = next((i for i in range(rk, len(M)) if M[i][col] != 0), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        pv = M[rk][col]
        M[rk] = [x / pv for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][col] != 0:
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
        piv.append(col)
        rk += 1
    x = [F(0)] * n
    for i, col in enumerate(piv):
        x[col] = M[i][n]
    return x


def part_C():
    rng = random.Random(22)
    stats = dict(systems=0, subsets=0, splittings=0)
    for trial in range(40):
        d = rng.choice([1, 2, 3])
        M = rng.randint(d + 1, min(6, d + 3))
        A = [[F(rng.randint(-3, 3)) for _ in range(d)] for _ in range(M)]
        if rank(A) < d:
            continue
        b = [F(rng.randint(-3, 3)) for _ in range(M)]
        w = [F(rng.randint(1, 3)) for _ in range(M)]
        circ = circuits(A, b)
        bad = [set(c["J"]) for c in circ if c["inconsistent"]]
        stats["systems"] += 1
        # Lemma 22.1 on every subset
        for mask in range(1, 2 ** M):
            T = [i for i in range(M) if (mask >> i) & 1]
            stats["subsets"] += 1
            check((not consistent(A, b, T)) == any(J <= set(T) for J in bad), "Lemma 22.1")
        # (22.4) with F = cardinality
        best_fail = min(M - len(T) for T in (tuple(i for i in range(M) if (mask >> i) & 1) for mask in range(2 ** M))
                        if consistent(A, b, list(T)))
        best_hit = min(len(D_) for k in range(M + 1) for D_ in combinations(range(M), k)
                       if all(set(D_) & J for J in bad))
        check(best_fail == best_hit, "(22.4): failures = hitting sets")
        L, r = fit_L2(A, b, w)
        if L == 0:
            check(not bad, "consistent system has no inconsistent circuit")
            continue
        # Lemma 22.2: sign-compatible decomposition of lam = w * r (squared loss: lam_i = w_i r_i up to factor)
        lam = [w[i] * r[i] for i in range(M)]
        rest = lam[:]
        comps = []
        for _ in range(4 * M):
            if all(x == 0 for x in rest):
                break
            cand = None
            for c in circ:
                g = c["g"]
                for sgn in (1, -1):
                    if all(rest[i] != 0 and (sgn * g[i]) * rest[i] > 0 for i in g):
                        cand = {i: sgn * g[i] for i in g}
                        break
                if cand:
                    break
            check(cand is not None, "a sign-compatible circuit exists (conformal decomposition)")
            s = min(rest[i] / cand[i] for i in cand)
            comp = {i: s * cand[i] for i in cand}
            comps.append(comp)
            for i in comp:
                rest[i] -= comp[i]
        check(all(x == 0 for x in rest), "decomposition exhausts lambda")
        wk = [[(w[i] * comp.get(i, 0) / lam[i]) if lam[i] != 0 else F(0) for i in range(M)] for comp in comps]
        check(all(sum(wk[k][i] for k in range(len(wk))) <= w[i] for i in range(M)), "split weights <= w")
        check(all(all(x >= 0 for x in row) for row in wk), "split weights nonnegative")
        total = sum(fit_L2(A, b, row)[0] for row in wk)
        check(total == L, "Lemma 22.2: circuit losses sum exactly to the full loss")
        for comp in comps:
            J = set(comp)
            check(any(set(c["J"]) == J and c["inconsistent"] for c in circ), "components are inconsistent circuits")
        stats["splittings"] += 1
        # Theorems 22.3 and 22.6 under independent retention
        theta = [F(rng.choice([1, 2, 3]), 4) for _ in range(M)]
        EL = F(0)
        for mask in range(2 ** M):
            pr = prod(theta[i] if (mask >> i) & 1 else 1 - theta[i] for i in range(M))
            V = [w[i] / theta[i] if (mask >> i) & 1 else F(0) for i in range(M)]
            EL += pr * fit_L2(A, b, V)[0]
        c_max = max(len(J) for J in bad)
        ths = sorted(theta)
        check(prod(ths[:c_max - 1]) * L <= EL <= L, "Theorem 22.3")
        Gam = lambda v: det([[sum(v[i] * A[i][j] * A[i][k] for i in range(M)) for k in range(d)] for j in range(d)])
        kappa = Gam(w) / Gam([w[i] / theta[i] for i in range(M)])
        check(kappa * L <= EL, "Theorem 22.6 (determinant bound)")
    # Proposition 22.4: Vandermonde sharpness and the parity law
    for d in (1, 2, 3):
        ts = [F(i + 1) for i in range(d + 1)]
        A = [[t ** k for k in range(d)] for t in ts]
        b = [t ** d for t in ts]
        w = [F(1)] * (d + 1)
        L, _ = fit_L2(A, b, w)
        th = F(1, 2)
        EL = F(0)
        for mask in range(2 ** (d + 1)):
            k = bin(mask).count("1")
            pr = th ** k * (1 - th) ** (d + 1 - k)
            V = [w[i] / th if (mask >> i) & 1 else F(0) for i in range(d + 1)]
            EL += pr * fit_L2(A, b, V)[0]
        check(L > 0 and EL == th ** d * L, "Proposition 22.4(1): exactly theta^d L")
        # parity law: masks with parity opposite to the all-ones mask, uniform
        masks = [m for m in range(2 ** (d + 1)) if bin(m).count("1") % 2 != (d + 1) % 2]
        EP = sum(fit_L2(A, b, [w[i] * 2 if (m >> i) & 1 else F(0) for i in range(d + 1)])[0] for m in masks) / len(masks)
        check(EP == 0, "Proposition 22.4(2): the d-wise independent parity law gives zero")
    # rank-one example values
    A = [[F(1)], [F(2)], [F(3)]]
    b = [F(1)] * 3
    ind = sum(fit_L2(A, b, [F(2) if (m >> i) & 1 else F(0) for i in range(3)])[0] for m in range(8)) / 8
    par = sum(fit_L2(A, b, [F(2) if (m >> i) & 1 else F(0) for i in range(3)])[0] for m in (0, 3, 5, 6)) / 4
    check(ind == F(503, 1820) and par == F(22, 65), "rank-one example values")
    # Example 22.8: nonlinear six-state model
    from itertools import permutations
    states = list(set(permutations([F(0), F(1, 4), F(1)])))
    lossN = lambda v: min(sum(v[i] * s[i] ** 2 for i in range(3)) for s in states)
    full = lossN([F(1)] * 3)
    EN = sum(lossN([F(2) if (m >> i) & 1 else F(0) for i in range(3)]) for m in range(8)) / 8
    check(full == F(17, 16) and EN == F(5, 16) and EN < F(1, 2) * full, "Example 22.8: affine structure is needed")
    RESULTS["C"] = dict(stats=stats, rank_one=[str(ind), str(par)], nonlinear=[str(full), str(EN)])
    print("Part C: ETH 22 Lemma 22.1 (%d subsets), (22.4), Lemma 22.2 (%d exact splittings), Theorems 22.3/22.6, "
          "Proposition 22.4 (d = 1..3), Example 22.8: OK" % (stats["subsets"], stats["splittings"]))


def main():
    part_T()
    part_C()
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "results", "gpt_round2.json"), "w") as fh:
        json.dump(RESULTS, fh, indent=1)
    print("ALL OK")


if __name__ == "__main__":
    main()

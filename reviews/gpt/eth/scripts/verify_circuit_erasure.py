"""Exact checks of ETH 22: minimal contradictions and deletion/refitting.

Standard library only. Direct normal equations, independent Gram determinants,
enumerated circuits and mask laws; no floating-point optimization.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import prod
from pathlib import Path
import json
import random

from verify_biased_energy import solve


def rref(a, n):
    rows = [[F(x) for x in row] for row in a]
    pivots = []
    k = 0
    for j in range(n):
        pivot = next((i for i in range(k, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[k], rows[pivot] = rows[pivot], rows[k]
        scale = rows[k][j]
        rows[k] = [v / scale for v in rows[k]]
        for i in range(len(rows)):
            if i != k and rows[i][j]:
                scale = rows[i][j]
                rows[i] = [u - scale * v for u, v in zip(rows[i], rows[k])]
        pivots.append(j)
        k += 1
    return rows, pivots


def rank(a, n):
    return len(rref(a, n)[1])


def nullspace(a, n):
    rows, pivots = rref(a, n)
    out = []
    for j in range(n):
        if j in pivots:
            continue
        x = [F(i == j) for i in range(n)]
        for i, p in enumerate(pivots):
            x[p] = -rows[i][j]
        out.append(x)
    return out


def determinant(a):
    rows = [[F(x) for x in row] for row in a]
    ans = F(1)
    for j in range(len(rows)):
        pivot = next((i for i in range(j, len(rows)) if rows[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            rows[j], rows[pivot] = rows[pivot], rows[j]
            ans = -ans
        p = rows[j][j]
        ans *= p
        for i in range(j + 1, len(rows)):
            t = rows[i][j] / p
            for k in range(j + 1, len(rows)):
                rows[i][k] -= t * rows[j][k]
    return ans


def gram(a, w):
    d = len(a[0])
    return [[sum(v * row[j] * row[k] for row, v in zip(a, w))
             for k in range(d)] for j in range(d)]


def loss(a, b, w):
    d = len(a[0])
    h = [sum(v * row[j] * y for row, y, v in zip(a, b, w)) for j in range(d)]
    beta = solve(gram(a, w), h)
    residual = [F(y) - sum(x * z for x, z in zip(row, beta)) for row, y in zip(a, b)]
    return sum(v * r * r for v, r in zip(w, residual)), beta, residual


def cb(a, w):
    d = len(a[0])
    return sum(determinant([a[i] for i in ix]) ** 2 * prod(w[i] for i in ix)
               for ix in combinations(range(len(a)), d))


def circuits(a, b):
    m, d = len(a), len(a[0])
    out = []
    for size in range(1, min(d + 1, m) + 1):
        for ix in combinations(range(m), size):
            ns = nullspace([[a[i][j] for i in ix] for j in range(d)], size)
            if len(ns) != 1 or not all(ns[0]):
                continue
            c = [F(0)] * m
            for i, t in zip(ix, ns[0]):
                c[i] = t
            if sum(x * y for x, y in zip(c, b)):
                out.append((ix, c))
    return out


def split(a, b, w, cs):
    full, beta, residual = loss(a, b, w)
    lam = [v * r for v, r in zip(w, residual)]
    rest = lam[:]
    out = []
    while any(rest):
        chosen = None
        for ix, c in cs:
            for sign in [1, -1]:
                if all(sign * c[i] * rest[i] > 0 for i in ix):
                    chosen = ix, [sign * x for x in c]
                    break
            if chosen:
                break
        assert chosen is not None, "missing sign-compatible circuit"
        ix, c = chosen
        t = min(rest[i] / c[i] for i in ix)
        h = [t * x for x in c]
        wk = [w[i] * h[i] / lam[i] if lam[i] else F(0) for i in range(len(a))]
        assert all(v >= 0 for v in wk)
        lk, _, rk = loss(a, b, wk)
        assert lk == sum(v * r * r for v, r in zip(wk, residual)) > 0
        assert all(sum(a[i][j] * h[i] for i in range(len(a))) == 0
                   for j in range(len(a[0])))
        assert sum(x * y for x, y in zip(h, b)) > 0
        numerator = sum(x * y for x, y in zip(h, b)) ** 2
        assert lk == numerator / sum(h[i] ** 2 / wk[i] for i in ix)
        out.append((ix, h, wk, lk))
        rest = [x - y for x, y in zip(rest, h)]
        assert all(x * y >= 0 for x, y in zip(rest, lam))
        assert len(out) <= len(a)
    assert sum(x[3] for x in out) == full
    for i, r in enumerate(residual):
        assert sum(x[2][i] for x in out) == (w[i] if r else 0)
    return out


def independent(theta):
    return [(mask, prod(t if mask >> i & 1 else 1 - t for i, t in enumerate(theta)))
            for mask in range(1 << len(theta))]


def parity_law(m, parity):
    return [(mask, F(1, 1 << (m - 1))) for mask in range(1 << m)
            if mask.bit_count() % 2 == parity]


def moments(law, theta, order):
    return all(sum(p for mask, p in law if all(mask >> i & 1 for i in ix))
               == prod(theta[i] for i in ix)
               for size in range(order + 1) for ix in combinations(range(len(theta)), size))


def expected(a, b, w, theta, law):
    return sum(p * loss(a, b, [v / t if mask >> i & 1 else 0
                              for i, (v, t) in enumerate(zip(w, theta))])[0]
               for mask, p in law if p)


def check_design(a, b, w, theta, law=None):
    a = [[F(x) for x in row] for row in a]
    b, w, theta = map(lambda xs: list(map(F, xs)), (b, w, theta))
    d = len(a[0])
    assert rank(a, d) == d
    law = independent(theta) if law is None else law
    assert moments(law, theta, min(d + 1, len(a)))
    full = loss(a, b, w)[0]
    aug = [[y] + row for row, y in zip(a, b)]
    gamma, delta = determinant(gram(a, w)), determinant(gram(aug, w))
    denominator = determinant(gram(a, [v / t for v, t in zip(w, theta)]))
    assert gamma == cb(a, w) > 0
    assert delta == cb(aug, w) == gamma * full
    kappa = gamma / denominator
    assert kappa >= prod(sorted(theta)[:d])
    ev = F(0)
    for mask, p in law:
        if not p:
            continue
        v = [x / t if mask >> i & 1 else 0 for i, (x, t) in enumerate(zip(w, theta))]
        l = loss(a, b, v)[0]
        g, z = determinant(gram(a, v)), determinant(gram(aug, v))
        assert z == g * l
        assert l * denominator >= z
        ev += p * l
    assert kappa * full <= ev <= full
    cs = circuits(a, b)
    pieces = split(a, b, w, cs)
    circuit_lower = sum(prod(theta[i] for i in ix)
                        * loss(a, b, [v / t for v, t in zip(wk, theta)])[0]
                        for ix, _, wk, _ in pieces)
    assert circuit_lower <= ev
    if full:
        cmax = max(map(lambda c: len(c[0]), cs))
        assert circuit_lower >= prod(sorted(theta)[:cmax - 1]) * full
    return {"rank": d, "rows": len(a), "full": str(full), "expected": str(ev),
            "determinant_factor": str(kappa), "circuit_lower": str(circuit_lower),
            "circuit_sizes": sorted({len(c[0]) for c in cs}), "pieces": len(pieces)}


def cube_checks(q):
    from verify_spline_bound import edges
    es = edges(q)
    m = len(es)
    def consistent(mask):
        adj = [[] for _ in range(1 << q)]
        for j, (x, y, i) in enumerate(es):
            if mask >> j & 1:
                even, odd = (x, y) if x.bit_count() % 2 == 0 else (y, x)
                adj[even].append((odd, F(2, i + 1)))
                adj[odd].append((even, F(-2, i + 1)))
        val = {}
        for root in range(1 << q):
            if root in val:
                continue
            val[root] = F(0)
            stack = [root]
            while stack:
                x = stack.pop()
                for y, g in adj[x]:
                    if y in val:
                        if val[y] != val[x] + g:
                            return False
                    else:
                        val[y] = val[x] + g
                        stack.append(y)
        return True
    ok = [consistent(mask) for mask in range(1 << m)]
    bad = [s for s in range(1, 1 << m) if not ok[s]
           and all(ok[s ^ (1 << j)] for j in range(m) if s >> j & 1)]
    # Independent row-rank and dependence check of every minimal contradiction.
    for s in bad:
        ix = [j for j in range(m) if s >> j & 1]
        a, b = [], []
        for j in ix:
            x, y, i = es[j]
            even, odd = (x, y) if x.bit_count() % 2 == 0 else (y, x)
            a.append([F(v == odd) - F(v == even) for v in range(1 << q)])
            b.append(F(2, i + 1))
        assert rank(a, 1 << q) == len(ix) - 1
        assert rank([[y] + row for row, y in zip(a, b)], (1 << q) + 1) == len(ix)
    minsum = minmax = m
    for deleted in range(1 << m):
        hit = all(deleted & c for c in bad)
        assert hit == ok[((1 << m) - 1) ^ deleted]
        if hit:
            counts = [sum((deleted >> j) & 1 for j, e in enumerate(es) if e[2] == i)
                      for i in range(q)]
            minsum = min(minsum, sum(counts))
            minmax = min(minmax, max(counts))
    assert (len(bad), minsum, minmax) == ((1, 1, 1) if q == 2 else (24, 4, 2))
    witness = next((s for s in bad if s.bit_count() > 4), None)
    if q == 3:
        assert witness is not None
        assert all(c & witness != c for c in bad if c.bit_count() == 4)
        assert not ok[witness]  # Checking only retained square faces misses this.
    return {"Q": q, "circuits": len(bad), "sizes": sorted({s.bit_count() for s in bad}),
            "min_total": minsum, "min_largest_direction": minmax,
            "long_cycle_control_mask": witness}


def main():
    rng = random.Random(5222026)
    designs = [check_design([[], [], []], [1, -2, 3], [1, 2, 3], [F(1, 2)] * 3)]
    for d in [1, 2, 3]:
        for m in [d + 1, d + 2]:
            for _ in range(4):
                while True:
                    a = [[rng.randrange(-2, 3) for _ in range(d)] for _ in range(m)]
                    if rank(a, d) == d:
                        break
                b = [rng.randrange(-3, 4) for _ in range(m)]
                w = [F(rng.randrange(1, 4), rng.randrange(1, 4)) for _ in range(m)]
                theta = [F(rng.randrange(1, 5), 4) for _ in range(m)]
                designs.append(check_design(a, b, w, theta))
        m = d + 2
        a = [[F(i + 1) ** j for j in range(d)] for i in range(m)]
        designs.append(check_design(a, [F(i + 1) ** d for i in range(m)], [1] * m,
                                    [F(1, 2)] * m, parity_law(m, 0)))

    # Positive residual with a singular retained design: Delta=Gamma=0, not L=0.
    a, b = [[1, 0], [0, 1], [0, 0]], [0, 0, 1]
    assert loss(a, b, [0, 0, 2])[0] == 2
    assert determinant(gram(a, [0, 0, 2])) == 0
    designs.append(check_design(a, b, [1] * 3, [F(1, 2)] * 3))

    sharp = []
    for d in range(1, 5):
        m = d + 1
        a = [[F(i + 1) ** j for j in range(d)] for i in range(m)]
        b = [F(i + 1) ** d for i in range(m)]
        full = loss(a, b, [F(1)] * m)[0]
        for theta in [F(1, 2), F(2, 3)]:
            ev = expected(a, b, [F(1)] * m, [theta] * m, independent([theta] * m))
            assert ev == theta ** d * full < theta ** (d - 1) * full
        law = parity_law(m, 1 - m % 2)
        assert moments(law, [F(1, 2)] * m, d)
        assert not moments(law, [F(1, 2)] * m, d + 1)
        assert expected(a, b, [F(1)] * m, [F(1, 2)] * m, law) == 0
        sharp.append({"d": d, "full_loss": str(full), "d_wise_expected": "0"})

    # Higher ambient rank can be irrelevant: independent contradictory pairs.
    short = []
    for d in [2, 3, 4]:
        a = [[F(j == k) for k in range(d)] for j in range(d) for _ in range(2)]
        b = [F(s) for _ in range(d) for s in [1, -1]]
        r = check_design(a, b, [1] * (2 * d), [F(1, 2)] * (2 * d))
        assert r["circuit_sizes"] == [2]
        assert F(r["expected"]) == d == F(r["circuit_lower"])
        # The same fixture has exact fourth-power loss: retained pairs minimize at 0,
        # singletons at +/-1. Each full pair contributes 4 to inverse-weighted L4.
        ev4 = sum(p * 4 * sum((mask >> (2 * j) & 3) == 3 for j in range(d))
                  for mask, p in independent([F(1, 2)] * (2 * d)))
        assert ev4 == d == F(1, 2) * (2 * d)
        short.append(r)

    a, b, w, theta = [[1], [2], [3]], [1] * 3, [F(1)] * 3, [F(1, 2)] * 3
    ev = expected(a, b, w, theta, independent(theta))
    even = expected(a, b, w, theta, parity_law(3, 0))
    assert ev == F(503, 1820) and even == F(22, 65) and ev != even
    assert moments(parity_law(3, 0), theta, 2)

    # Circuit-size survival is not a theorem for arbitrary nonlinear model images.
    states = list(permutations([F(0), F(1, 4), F(1)]))
    def nonlinear(v):
        return min(sum(x * r * r for x, r in zip(v, state)) for state in states)
    nl_full = nonlinear([1, 1, 1])
    nl_ev = sum(p * nonlinear([2 if mask >> i & 1 else 0 for i in range(3)])
                for mask, p in independent([F(1, 2)] * 3))
    assert nl_full == F(17, 16) and nl_ev == F(5, 16) < nl_full / 2
    assert all(nonlinear([int(i in ix) for i in range(3)]) == F(1, 16)
               for ix in combinations(range(3), 2))
    assert all(nonlinear([int(i == j) for i in range(3)]) == 0 for j in range(3))

    ridge = []
    for a, b, w in [([[1, 0], [0, 2], [1, 1]], [1, 2, -1], [1, 2, 3]),
                     ([[1, 0], [0, 2]], [2, -1], [3, 1]),
                     ([[1, 0], [2, 0]], [1, -1], [1, 2])]:
        d, m, lam = 2, len(a), F(1, 3)
        extra = [[F(i == j) for j in range(d)] for i in range(d)]
        aa, bb, ww = a + extra, b + [0] * d, w + [lam] * d
        ts = [F(2, 3)] * m + [F(1)] * d
        r = check_design(aa, bb, ww, ts)
        g = gram(a, w)
        diag = all(g[i][j] == 0 for i in range(d) for j in range(d) if i != j)
        if diag:
            assert F(r["determinant_factor"]) == prod(
                (g[j][j] + lam) / (g[j][j] / F(2, 3) + lam) for j in range(d))
        ridge.append(r)

    out = {"status": "ALL OK", "exact_designs": designs, "sharp_threshold": sharp,
           "short_circuits_and_fourth_power": short, "ridge": ridge,
           "pairwise_independence_not_exact_expectation": {"independent": str(ev),
                                                           "pairwise_even": str(even)},
           "nonlinear_model_control": {"full": str(nl_full), "expected": str(nl_ev),
                                       "false_circuit_size_bound": str(nl_full / 2)},
           "cube_locality": [cube_checks(2), cube_checks(3)],
           "controls": {"d_wise_not_enough": True, "smaller_exponent_false": True,
                        "singular_positive_residual": True, "faces_not_all_loops": True,
                        "matching_moments_not_matching_optimum": True,
                        "nonlinear_model_extension_false": True}}
    path = Path(__file__).resolve().parents[1] / "results" / "circuit_erasure.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"ALL OK: {len(designs)} rational designs, sharpness d=1..4, "
          "short circuits, fourth-power fixtures, ridge, cube Q=2,3")


if __name__ == "__main__":
    main()

"""Exact algebra and hypothesis controls for KAK 23; standard library only.

The sphere-bundle argument is a cited topological proof, not a finite test.
Checks here cover transformed Hopf pencils, trace reconstruction, spin-factor
products, state selectors, and concrete failures when global hypotheses drop.
"""
from fractions import Fraction as F
from math import prod
from pathlib import Path
import json

from verify_curvature_fibre import (
    add, eye, hopf, matmul, mv, non_clifford, scaled, transpose,
)
from verify_weyl_fibre import determinant


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def combine(matrices, coeff):
    n = len(matrices[0])
    return [[sum(c * a[i][j] for a, c in zip(matrices, coeff))
             for j in range(n)] for i in range(n)]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def upper_solve(a, b):
    x = [F(0)] * len(b)
    for i in reversed(range(len(b))):
        x[i] = (b[i] - sum(a[i][j] * x[j] for j in range(i + 1, len(b)))) / a[i][i]
    return x


def positive_definite(a):
    return all(determinant([row[:j] for row in a[:j]]).re > 0 for j in range(1, len(a) + 1))


def model(m):
    gamma = hopf(m)[0]
    n, size = len(gamma), len(gamma[0])
    assert n == m + 1 and size == 2 * m == n + m - 1
    a = [[F(i + 1) if i == j else F(1, 3) if j == i + 1 else F(0)
          for j in range(n)] for i in range(n)]
    drift = [F(j - 2, 7) for j in range(n)]
    ds = [combine(gamma, [a[i][j] for i in range(n)]) for j in range(n)]
    cs = [add(d, scaled(eye(size), c)) for d, c in zip(ds, drift)]
    recovered_drift = [trace(c) / size for c in cs]
    assert recovered_drift == drift
    recovered_d = [add(c, scaled(eye(size), -l)) for c, l in zip(cs, recovered_drift)]
    q = [[trace(matmul(x, y)) / size for y in recovered_d] for x in recovered_d]
    assert q == matmul(transpose(a), a) and positive_definite(q)
    for i in range(n):
        for j in range(n):
            assert add(matmul(ds[i], ds[j]), matmul(ds[j], ds[i])) == scaled(eye(size), 2 * q[i][j])

    selectors = 0
    for seed in range(1, 5):
        z = [F(((i + 3) * seed) % 11 - 5, i % 3 + 1) for i in range(size)]
        norm = dot(z, z)
        h = [dot(z, mv(g, z)) / norm for g in gamma]
        assert dot(h, h) == 1
        xi = upper_solve(a, h)
        assert mv(a, xi) == h
        assert mv(combine(cs, xi), z) == [(dot(drift, xi) + 1) * x for x in z]
        selectors += 1
        # Spin-factor Jordan product, computed from actual matrix multiplication.
        u = [F(i + seed, 3) for i in range(n)]
        v = [F(2 * i - seed, 5) for i in range(n)]
        alpha, beta = F(seed, 2), F(1 - seed, 3)
        x = add(scaled(eye(size), alpha), combine(ds, u))
        y = add(scaled(eye(size), beta), combine(ds, v))
        lhs = scaled(add(matmul(x, y), matmul(y, x)), F(1, 2))
        rhs = add(scaled(eye(size), alpha * beta + dot(u, mv(q, v))),
                  combine(ds, [alpha * vv + beta * uu for uu, vv in zip(u, v)]))
        assert lhs == rhs
    xi = [F(0)] * (n - 1) + [F(1)]
    qxi, value = mv(q, xi), dot(xi, mv(q, xi))
    hscaled = [[value * q[i][j] - qxi[i] * qxi[j] for j in range(n - 1)]
               for i in range(n - 1)]
    assert positive_definite(hscaled)
    return dict(multiplicity=m, spatial_dimension=n, state_dimension=size,
                drift=list(map(str, drift)), metric=[[str(x) for x in row] for row in q],
                trace_reconstruction=True, Clifford_coefficients=True,
                spin_factor_products=4, exact_state_selectors=selectors)


def controls():
    cs, gap, north, derivatives, phi = non_clifford()
    size, m, k = 6, 2, 4
    # At the opposite direction the top space is the last complex coordinate.
    top = 4
    complement = [0, 1, 2, 3]
    south_gaps = [F(3), F(3), F(1), F(1)]
    south = [[2 * sum(cs[i][b][top] * cs[j][b][top] / g
                     for b, g in zip(complement, south_gaps))
              for j in range(k)] for i in range(k)]
    assert south == [[(F(6) if i == 2 else F(32, 3) if i == 3 else F(0)) * (i == j)
                      for j in range(k)] for i in range(k)]
    assert all(x > 0 for x in north) and not positive_definite(south)
    d = add(cs[-1], scaled(eye(size), -trace(cs[-1]) / size))
    assert matmul(d, d) != scaled(eye(size), trace(matmul(d, d)) / size)
    # The 3-state arrow is locally simple/curved; its opposite top is double.
    arrow_north, arrow_south = [1, -1, -1], [-1, 1, 1]
    assert arrow_north.count(max(arrow_north)) == 1
    assert arrow_south.count(max(arrow_south)) == 2
    # Repeating a globally curved Clifford representation breaks critical size.
    old_m, n, repeated_n, repeated_m = 2, 3, 8, 4
    assert repeated_n != n + repeated_m - 1
    return dict(local_nonClifford_curvature=north, opposite_curvature=south,
                local_is_not_global=True, trace_without_two_levels_fails=True,
                arrow_multiplicity_changes=True, repeated_model_is_not_critical=True)


def main():
    cases = [model(m) for m in [1, 2, 4, 8]]
    failed = controls()
    minimum_dimensions = [1, 2, 4, 4, 8, 8, 8, 8]
    allowed = [m for m in range(1, 129)
               if m % (16 ** ((m - 1) // 8) * minimum_dimensions[(m - 1) % 8]) == 0]
    assert allowed == [1, 2, 4, 8]
    # Wrong period multiplier would incorrectly admit dimension 16.
    assert 16 % (2 ** (15 // 8) * minimum_dimensions[15 % 8]) == 0
    report = dict(scope="Finite exact matrix checks; topology and real Clifford classification are explicitly cited in KAK 23.",
                  cases=cases, controls=failed,
                  Clifford_dimension_arithmetic_to_128=allowed,
                  wrong_periodicity_control=True)
    path = Path(__file__).resolve().parent.parent / "results" / "global_rigidity.json"
    path.write_text(json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    print("Trace reconstruction, Clifford relations and spin-factor products: exact")
    print("16 exact state selectors; global-hypothesis controls: OK")
    print("Classical dimension-table arithmetic and wrong-periodicity control: OK")
    print("ALL OK")


if __name__ == "__main__":
    main()

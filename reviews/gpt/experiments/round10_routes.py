"""Exact probes retained from round 10's abstraction search.

Not part of verify_all.py. No asymptotic or matrix-classification conclusion
is inferred from these finite examples.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "reviews/gpt/kakeya/scripts"))
sys.path.insert(0, str(ROOT / "reviews/gpt/eth/scripts"))
from verify_curvature_fibre import hopf
from verify_circuit_erasure import determinant, expected, independent, loss


def commutant_dimension(cs):
    """Solve X C_j=C_j X exactly for signed-permutation involutions.

    Entry equations identify two unknowns up to sign. A negative cycle kills
    that whole component. Unforced components are a basis of the commutant.
    """
    n = len(cs[0])
    adj = [[] for _ in range(n * n)]
    for c in cs:
        perm, signs = [], []
        for row in c:
            ix = [j for j, x in enumerate(row) if x]
            assert len(ix) == 1 and abs(row[ix[0]]) == 1
            perm.append(ix[0])
            signs.append(row[ix[0]])
        assert all(perm[perm[i]] == i and signs[i] * signs[perm[i]] == 1
                   for i in range(n))
        for i in range(n):
            for j in range(n):
                x, y = i * n + j, perm[i] * n + perm[j]
                t = signs[i] * signs[j]
                adj[x].append((y, t))
                adj[y].append((x, t))
    labels = {}
    dimension = 0
    for root in range(n * n):
        if root in labels:
            continue
        labels[root] = F(1)
        stack = [root]
        forced_zero = False
        while stack:
            x = stack.pop()
            for y, s in adj[x]:
                if y in labels:
                    if labels[y] != s * labels[x]:
                        forced_zero = True
                else:
                    labels[y] = s * labels[x]
                    stack.append(y)
        dimension += not forced_zero
    return dimension


def arrow(xi):
    x, y, z = map(F, xi)
    return [[z, x, y], [x, -z, 0], [y, 0, -z]]


def main():
    commute = []
    for m, dim in [(1, 1), (2, 2), (4, 4), (8, 1)]:
        cs = hopf(m)[0]
        actual = commutant_dimension(cs)
        assert actual == dim
        commute.append({"m": m, "N": 2 * m, "real_commutant_dimension": actual})
    assert commute[-1]["real_commutant_dimension"] < commute[-1]["m"]

    # Quadratic phase factor does not imply a fixed two-state polarization model.
    probes = 0
    for xi in [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, -1),
               (F(2, 3), F(-1, 2), F(3, 4))]:
        c = arrow(xi)
        q = sum(F(x) ** 2 for x in xi)
        for t in map(F, range(-3, 4)):
            matrix = [[t * (i == j) - c[i][j] for j in range(3)] for i in range(3)]
            assert determinant(matrix) == (t + xi[2]) * (t * t - q)
            probes += 1
    vectors = [[1, 1, 0], [1, 0, 1], [1, 0, 0]]
    directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    for a, xi in zip(vectors, directions):
        c = arrow(xi)
        assert [sum(x * y for x, y in zip(row, a)) for row in c] == a
    assert determinant(vectors) != 0

    # Same inconsistent-circuit hypergraph, different energy ratios.
    ratios = []
    for a, b in [([[1], [2], [3]], [1, 1, 1]),
                 ([[1], [1], [1]], [-1, 0, 1])]:
        full = loss(a, b, [F(1)] * 3)[0]
        ev = expected(a, b, [F(1)] * 3, [F(1, 2)] * 3, independent([F(1, 2)] * 3))
        ratios.append(ev / full)
    assert ratios == [F(503, 780), F(5, 8)] and ratios[0] != ratios[1]

    out = {"status": "exact probes passed", "commutants": commute,
           "arrow": {"characteristic_probes": probes, "polarization_span_rank": 3,
                     "curved_phase_factor_degree": 2,
                     "global_constant_multiplicity": False},
           "same_circuit_hypergraph_ratios": list(map(str, ratios))}
    path = Path(__file__).resolve().parent / "results" / "round10_routes.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print("Exploration complete: commutants, arrow polarization, circuit geometry")


if __name__ == "__main__":
    main()

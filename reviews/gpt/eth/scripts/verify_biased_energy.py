"""Exact checks of ETH 21: product-bias defect energy; standard library only.

Independent checks: full vertex-coefficient normal equations versus tensor
spectral sectors; explicit minimizing polynomials; Cauchy--Binet capacities;
bias comparison, sharp high-degree case, and falsified shortcut controls.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb, prod
from pathlib import Path
import json


def solve(a, b):
    """One solution of a rational, possibly singular consistent system."""
    n = len(b)
    if not n:
        return []
    rows = [[F(x) for x in row] + [F(y)] for row, y in zip(a, b)]
    pivots = []
    k = 0
    for j in range(n):
        pivot = next((i for i in range(k, n) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[k], rows[pivot] = rows[pivot], rows[k]
        scale = rows[k][j]
        rows[k] = [v / scale for v in rows[k]]
        for i in range(k + 1, n):
            if rows[i][j]:
                scale = rows[i][j]
                rows[i] = [u - scale * v for u, v in zip(rows[i], rows[k])]
        pivots.append(j)
        k += 1
    assert all(any(row[:-1]) or row[-1] == 0 for row in rows)
    x = [F(0)] * n
    for i, j in reversed(list(enumerate(pivots))):
        x[j] = rows[i][-1] - sum(rows[i][r] * x[r] for r in range(j + 1, n))
    assert all(sum(u * v for u, v in zip(row, x)) == y for row, y in zip(a, b))
    return x


def value(poly, t):
    ans = F(0)
    for c in reversed(poly):
        ans = ans * t + c
    return ans


def capacity(nodes, weights, degree):
    """Scalar constrained least squares, and one explicit minimizer."""
    if len(nodes) <= degree:
        coeff = [F(1)]
        for t in nodes:
            nxt = [F(0)] * (len(coeff) + 1)
            for k, c in enumerate(coeff):
                nxt[k] += c
                nxt[k + 1] -= c / t
            coeff = nxt
        return F(0), coeff + [F(0)] * (degree + 1 - len(coeff))
    a = [[sum(w * t ** (r + s) for t, w in zip(nodes, weights))
          for s in range(1, degree + 1)] for r in range(1, degree + 1)]
    b = [-sum(w * t ** r for t, w in zip(nodes, weights))
         for r in range(1, degree + 1)]
    coeff = [F(1)] + solve(a, b)
    energy = sum(w * value(coeff, t) ** 2 for t, w in zip(nodes, weights))
    return energy, coeff


def vandermonde(nodes):
    return prod((b - a) ** 2 for a, b in combinations(nodes, 2))


def determinants(nodes, weights, degree):
    numerator = sum(prod(weights[i] for i in ix) *
                    vandermonde([nodes[i] for i in ix])
                    for ix in combinations(range(len(nodes)), degree + 1))
    denominator = sum(prod(weights[i] * nodes[i] ** 2 for i in ix) *
                      vandermonde([nodes[i] for i in ix])
                      for ix in combinations(range(len(nodes)), degree))
    return numerator, denominator


def spectral(nodes, biases, degree, direction_weights=None):
    q = len(nodes)
    direction_weights = direction_weights or [F(1)] * q
    theta = [4 * p * (1 - p) for p in biases]
    coeffs = [[F(0)] * (degree + 1) for _ in range(1 << q)]
    optimum = F(0)
    for mask in range(1 << q):
        ix = [i for i in range(q) if mask >> i & 1]
        probability = prod(theta[i] if mask >> i & 1 else 1 - theta[i] for i in range(q))
        cost, poly = capacity([nodes[i] for i in ix],
                              [direction_weights[i] / theta[i] for i in ix], degree)
        num, den = determinants([nodes[i] for i in ix],
                                [direction_weights[i] / theta[i] for i in ix], degree)
        if den:
            assert cost == num / den
        else:
            assert cost == num == 0
        optimum += probability * cost
        for x in range(1 << q):
            factor = F(1)
            for i, p in enumerate(biases):
                active = 2 * (1 - p) if x >> i & 1 else 2 * p
                factor *= active if mask >> i & 1 else 1 - active
            for k, c in enumerate(poly):
                coeffs[x][k] += factor * c
    assert all(poly[0] == 1 for poly in coeffs)
    return optimum, coeffs


def direct_energy(coeffs, nodes, biases, direction_weights=None):
    q = len(nodes)
    direction_weights = direction_weights or [F(1)] * q
    out = F(0)
    for i, t in enumerate(nodes):
        for x in range(1 << q):
            if x >> i & 1:
                continue
            y = x | (1 << i)
            mass = prod(p if x >> j & 1 else 1 - p
                        for j, p in enumerate(biases) if j != i)
            defect = (value(coeffs[x], t) + value(coeffs[y], t)) / 2
            out += direction_weights[i] * mass * defect ** 2
    return out


def vertex_optimum(nodes, biases, degree, direction_weights=None):
    """Full original QP: no spectral decomposition or scalar capacities."""
    q = len(nodes)
    direction_weights = direction_weights or [F(1)] * q
    size = degree * (1 << q)
    a = [[F(0)] * size for _ in range(size)]
    b = [F(0)] * size
    constant = sum(direction_weights)
    for i, t in enumerate(nodes):
        for x in range(1 << q):
            if x >> i & 1:
                continue
            y = x | (1 << i)
            mass = direction_weights[i] * prod(p if x >> j & 1 else 1 - p
                                               for j, p in enumerate(biases) if j != i)
            row = [(z * degree + k - 1, t ** k / 2)
                   for z in (x, y) for k in range(1, degree + 1)]
            for r, v in row:
                b[r] -= mass * v
                for s, w in row:
                    a[r][s] += mass * v * w
    coeff = solve(a, b)
    return constant - sum(v * w for v, w in zip(b, coeff))


def local_operator_check(p):
    theta = 4 * p * (1 - p)
    mass = [1 - p, p]
    # K=P^*P, computed in value coordinates with the biased inner product.
    k = [[F(1, 4) / mass[r] for _ in range(2)] for r in range(2)]
    active = [2 * p, 2 * (1 - p)]
    inactive = [1 - c for c in active]
    assert sum(m * a * b for m, a, b in zip(mass, active, inactive)) == 0
    assert sum(m * a * a for m, a in zip(mass, active)) == theta
    assert sum(m * b * b for m, b in zip(mass, inactive)) == 1 - theta
    assert [sum(k[r][s] * active[s] for s in range(2)) for r in range(2)] == [v / theta for v in active]
    assert all(sum(k[r][s] * inactive[s] for s in range(2)) == 0 for r in range(2))


def main():
    for p in [F(1, 7), F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(6, 7)]:
        local_operator_check(p)
    cases = []
    tests = [([F(i) for i in range(1, q + 1)], [p] * q, d, None)
             for q in range(1, 5) for d in range(q)
             for p in [F(1, 4), F(1, 2)]]
    tests += [([F(-2), F(1), F(3), F(7)],
               [F(1, 5), F(1, 2), F(2, 3), F(3, 4)], d,
               [F(1), F(2), F(3), F(5)]) for d in range(4)]
    for nodes, biases, degree, weights in tests:
        exact, polynomials = spectral(nodes, biases, degree, weights)
        direct = direct_energy(polynomials, nodes, biases, weights)
        full = vertex_optimum(nodes, biases, degree, weights)
        assert exact == direct == full
        actual_weights = weights or [F(1)] * len(nodes)
        uniform, _ = capacity(nodes, actual_weights, degree)
        theta = [4 * p * (1 - p) for p in biases]
        num, _ = determinants(nodes, actual_weights, degree)
        _, den = determinants(nodes, [w / t for w, t in zip(actual_weights, theta)], degree)
        lower = num / den
        assert lower <= exact <= uniform
        if degree == len(nodes) - 1:
            assert lower == exact
        if all(p == F(1, 2) for p in biases):
            rho = prod(F(len(nodes) - j, len(nodes) + j + 1) for j in range(degree + 1))
            assert exact == rho / (1 - rho)
        cases.append(dict(nodes=list(map(str, nodes)), biases=list(map(str, biases)), degree=degree,
                          weights=list(map(str, actual_weights)), minimum=str(exact),
                          determinant_lower_bound=str(lower), uniform_minimum=str(uniform)))

    # Larger scalar-sector calculations; no full-QP shortcut is used above.
    for q in range(2, 8):
        nodes = list(map(F, range(1, q + 1)))
        for p in [F(1, 7), F(1, 4), F(1, 3), F(1, 2)]:
            theta = 4 * p * (1 - p)
            exact, poly = spectral(nodes, [p] * q, q - 1)
            assert exact == theta ** (q - 1) / (comb(2 * q, q) - 1)
            assert direct_energy(poly, nodes, [p] * q) == exact

    # Coordinatewise monotonicity in theta, including a bias crossing 1/2.
    monotone_cases = 0
    for q in range(2, 5):
        nodes = list(map(F, range(1, q + 1)))
        for d in range(1, q):
            for coordinate in range(q):
                energies = []
                for p in [F(1, 10), F(1, 4), F(1, 3), F(1, 2)]:
                    biases = [F(1, 4)] * q
                    biases[coordinate] = p
                    energies.append(spectral(nodes, biases, d)[0])
                    reflected = biases.copy()
                    reflected[coordinate] = 1 - p
                    assert spectral(nodes, reflected, d)[0] == energies[-1]
                assert energies == sorted(energies)
                monotone_cases += 1

    example, family = spectral([F(1), F(2)], [F(1, 4)] * 2, 1)
    assert example == F(3, 20)
    # Controls: ordinary Parseval/scalar optimality, wrong thinning law,
    # a stronger exponent D-1 comparison, and omission of 1/theta must fail.
    uniform = F(1, 5)
    theta = F(3, 4)
    controls = dict(scalar_optimality_is_false=example != uniform,
                    thinning_at_2p_is_false=example != F(1, 2) * uniform,
                    stronger_exponent_D_minus_1_bound_is_false=example < uniform,
                    omitted_inverse_theta_is_false=example != theta * example)
    assert all(controls.values())
    output = dict(scope="Finite exact checks support the separate prose proofs; no SAT algorithm or Lean build.",
                  local_operator_biases=6, independent_quadratic_programs=len(cases),
                  monotonicity_cases=monotone_cases, sharp_high_degree_cases=24,
                  cases=cases, negative_controls=controls,
                  example=dict(Q=2, D=1, p="1/4", minimum=str(example),
                               coefficients=[[str(c) for c in row] for row in family]))
    path = Path(__file__).resolve().parent.parent / "results" / "biased_energy.json"
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"Independent vertex QPs: {len(cases)}; explicit families and determinants: exact")
    print(f"Coordinate monotonicity: {monotone_cases}; sharp high-degree cases: 24")
    print("Negative controls: all fail as required")
    print("ALL OK")


if __name__ == "__main__":
    main()

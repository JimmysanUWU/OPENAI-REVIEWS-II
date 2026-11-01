"""Exact finite checks for ETH 17. Standard library, Python 3.10+.

The formulas are checked against independently constructed linear systems,
direct sums, and finite-field designs. This is not a formal proof.
"""
from copy import deepcopy
from fractions import Fraction as F
from math import comb, factorial, log
from pathlib import Path
import json
import random

from verify_spline_bound import (LinearSystem, rank, parity, edges, constraint,
                                trace_rows, evaluate, biased_mass, expectation_tail)


class ModSystem:
    def __init__(self, n, prime):
        self.n, self.p, self.rows = n, prime, {}

    def add(self, row, rhs=0):
        row = [int(x) % self.p for x in row] + [rhs % self.p]
        for pivot in sorted(self.rows):
            a = row[pivot]
            row = [(x-a*y) % self.p for x, y in zip(row, self.rows[pivot])]
        pivot = next((i for i in range(self.n) if row[i]), None)
        if pivot is None:
            return row[-1] == 0
        inv = pow(row[pivot], -1, self.p)
        self.rows[pivot] = [x*inv % self.p for x in row]
        return True

    def solution(self):
        answer = [0]*self.n
        for pivot in sorted(self.rows, reverse=True):
            row = self.rows[pivot]
            answer[pivot] = (row[-1]-sum(row[j]*answer[j]
                                          for j in range(pivot+1, self.n))) % self.p
        return answer


def mod_rank(rows, n, prime):
    s = ModSystem(n, prime)
    for row in rows:
        assert s.add(row)
    return len(s.rows)


def field_checks():
    dimensions, designs = [], []
    for prime in [3, 5, 7]:
        for q in range(1, min(prime, 5)):
            nodes = list(range(1, q+1))
            for d in range(q+1):
                n = (1 << q)*(d+1)
                rows = [constraint(q, d, nodes, e) for e in edges(q)]
                r = mod_rank(rows, n, prime)
                trace_dim = mod_rank(rows+trace_rows(q, d), n, prime)-r
                assert n-r == sum(max(0, d-s.bit_count()+1) for s in range(1 << q))
                assert trace_dim == sum(comb(q, k) for k in range(d+1))
                dimensions.append(dict(prime=prime, q=q, degree=d, trace_dimension=trace_dim))
                if d == q:
                    continue
                for seed in range(2):
                    system = ModSystem(n, prime)
                    for x, row in enumerate(trace_rows(q, d)):
                        assert system.add(row, parity(x))
                    order = edges(q)
                    random.Random(7300+seed).shuffle(order)
                    for e in order:
                        system.add(constraint(q, d, nodes, e))
                    solution = system.solution()
                    coeff = [solution[x*(d+1):(x+1)*(d+1)] for x in range(1 << q)]
                    assert all(coeff[x][0] == parity(x) % prime for x in range(1 << q))
                    bad = [e for e in edges(q)
                           if (evaluate(coeff[e[0]], nodes[e[2]])-
                               evaluate(coeff[e[1]], nodes[e[2]])) % prime]
                    lower = sum(comb(q, k) for k in range(q+1) if 2*k > q+d)
                    assert len(bad) >= lower
                    probabilities = [F(1, 4), F(1, 3), F(2, 5), F(1, 2)][:q]
                    masses = biased_mass(q, bad, probabilities)
                    lhs = sum(p*m for p, m in zip(probabilities, masses))
                    rhs = expectation_tail(q, d, probabilities)
                    assert lhs >= rhs
                    designs.append(dict(prime=prime, q=q, degree=d, bad_edges=len(bad),
                                        lower_bound=lower, biased_lhs=str(lhs), biased_rhs=str(rhs)))
    # Distinct integer nodes can collide in the residue field.
    p, q = 5, 4
    nodes = [1+i*p for i in range(1, q+1)]
    assert len(set(nodes)) == q and len({x % p for x in nodes}) == 1
    assert all((1-t) % p == 0 for t in nodes)
    assert sum(comb(q, k) for k in range(q+1) if 2*k > q+1) > 0
    # Characteristic two: scalar sum vanishes, but the nilpotent matrix acts nontrivially.
    assert (1+1) % 2 == 0
    N = [[1, 1], [1, 1]]
    assert [[sum(N[i][k]*N[k][j] for k in range(2)) % 2
             for j in range(2)] for i in range(2)] == [[0, 0], [0, 0]]
    assert [N[i][0] for i in range(2)] == [1, 1]
    return dict(dimensions=dimensions, feasible_designs=designs,
                controls=dict(residue_separation_needed=True, odd_characteristic_needed=True,
                              nilpotence_does_not_annihilate_readouts=True))


def character_union(q, generators, degrees, d):
    answer = set()
    for j in range(1 << len(generators)):
        centre, cost = 0, 0
        for k, a in enumerate(generators):
            if j & (1 << k):
                centre ^= a
                cost += degrees[k]
        answer.update(centre ^ s for s in range(1 << q) if s.bit_count() <= d-cost)
    return answer


def trace_checks():
    # The concrete two-generator example is realized by a jointly consistent graph.
    q, nodes, generators, degrees = 4, [1, 2, 3, 4], [3, 12], [1, 1]
    systems = []
    for a, degree in zip(generators, degrees):
        system = LinearSystem((1 << q)*(degree+1))
        for x, row in enumerate(trace_rows(q, degree)):
            assert system.add(row, parity(a & x))
        systems.append(system)
    order = edges(q)
    random.Random(419).shuffle(order)
    retained = []
    for e in order:
        trial = deepcopy(systems)
        if all(s.add(constraint(q, degree, nodes, e)) for s, degree in zip(trial, degrees)):
            systems = trial
            retained.append(e)
    families = []
    for s, degree, a in zip(systems, degrees, generators):
        flat = s.solution()
        coeff = [flat[x*(degree+1):(x+1)*(degree+1)] for x in range(1 << q)]
        assert all(coeff[x][0] == parity(a & x) for x in range(1 << q))
        families.append(coeff)
    bad = [e for e in edges(q) if any(evaluate(c[e[0]], nodes[e[2]]) !=
                                    evaluate(c[e[1]], nodes[e[2]]) for c in families)]
    assert len(bad) == len(edges(q))-len(retained)
    union = character_union(q, generators, degrees, 2)
    assert union == set(range(16))
    assert len(bad) >= len(union)-11 == 5
    direct = {s for s in range(16) if s.bit_count() <= 2}
    for a in generators:
        direct.update(a ^ s for s in range(16) if s.bit_count() <= 1)
    assert len(direct) == 15  # The product contributes the remaining character.
    joint = dict(q=q, generators=generators, degrees=degrees, bad_edges=len(bad),
                 joint_bound=5, without_product=4, single_bounds=[2, 2])

    # General noncharacter traces: compare the multiplication rank with actual erasures.
    cases = []
    q, degree, nodes = 3, 1, [1, 3, 7]
    for seed in range(6):
        rng = random.Random(1780+seed)
        h = [F(rng.randrange(-3, 4)) for _ in range(1 << q)]
        system = LinearSystem((1 << q)*(degree+1))
        for x, row in enumerate(trace_rows(q, degree)):
            assert system.add(row, h[x])
        order = edges(q)
        rng.shuffle(order)
        for e in order:
            system.add(constraint(q, degree, nodes, e))
        flat = system.solution()
        coeff = [flat[x*(degree+1):(x+1)*(degree+1)] for x in range(1 << q)]
        bad = [e for e in edges(q) if evaluate(coeff[e[0]], nodes[e[2]]) !=
               evaluate(coeff[e[1]], nodes[e[2]])]
        fourier = [sum(h[x]*parity(s & x) for x in range(1 << q))/(1 << q)
                   for s in range(1 << q)]
        for d in range(degree, q+1):
            low = [s for s in range(1 << q) if s.bit_count() <= d-degree]
            high = [s for s in range(1 << q) if s.bit_count() > d]
            r = rank([[fourier[t ^ s] for s in low] for t in high], len(low))
            assert r <= len(bad)
            cases.append(dict(seed=seed, degree=d, rank_bound=r, bad_edges=len(bad)))
    return dict(joint_example=joint, arbitrary_traces=cases)


def falling(x, j):
    value = 1
    for k in range(j):
        value *= x-k
    return value


def H(q, j, x):
    def R(y):
        return falling(y, j)*falling(y-q-1, j)
    numerator = sum((-1)**(j-a)*comb(j, a)*R(x+a) for a in range(j+1))
    denominator = (-1)**j*factorial(j)*falling(q, j)
    return F(numerator, denominator)


def rho(q, d):
    answer = F(1)
    for j in range(d+1):
        answer *= F(q-j, q+j+1)
    return answer


def gram_minimum(nodes, weights, degree):
    """Direct normal equations after fixing the constant coefficient to one."""
    system = LinearSystem(degree)
    for r in range(1, degree+1):
        row = [sum(w*t**(r+s) for t, w in zip(nodes, weights))
               for s in range(1, degree+1)]
        assert system.add(row, -sum(w*t**r for t, w in zip(nodes, weights)))
    assert len(system.rows) == degree
    coefficients = [F(1)] + system.solution()
    energy = sum(w*evaluate(coefficients, t)**2 for t, w in zip(nodes, weights))
    return coefficients, energy


def orthogonal_and_gram_checks():
    cases, orthogonality_pairs = [], 0
    for q in range(1, 13):
        values = [[H(q, j, x) for x in range(q+1)] for j in range(q+1)]
        norms = []
        for j in range(q+1):
            expected = F(factorial(q+j+1)*factorial(q-j),
                         (2*j+1)*factorial(q)**2)
            assert values[j][0] == 1
            assert sum(v*v for v in values[j]) == expected
            for k in range(j):
                assert sum(v*w for v, w in zip(values[j], values[k])) == 0
                orthogonality_pairs += 1
            norms.append(expected)
        for d in range(q):
            r = rho(q, d)
            assert r == F(factorial(q)**2, factorial(q-d-1)*factorial(q+d+1))
            kernel = [sum(values[j][x]/norms[j] for j in range(d+1)) for x in range(q+1)]
            assert kernel[0] == 1-r
            optimal = [v/kernel[0] for v in kernel]
            assert optimal[0] == 1
            energy = sum(v*v for v in optimal[1:])
            assert energy == r/(1-r)
            coefficients, independent = gram_minimum(list(range(1, q+1)), [1]*q, d)
            assert energy == independent
            assert optimal == [evaluate(coefficients, x) for x in range(q+1)]
            if d == 0:
                assert energy == q
            if d == 1:
                assert energy == F(q*(q-1), 2*(2*q+1))
            if d == q-1:
                assert energy == F(1, comb(2*q, q)-1)
            cases.append(dict(q=q, degree=d, minimum_energy=str(energy)))
    weighted = []
    for nodes, weights in [([-3, 1, 5], [F(1, 2), F(2), F(3, 2)]),
                            ([1, 2, 4, 9], [1, 3, 2, 1])]:
        for d in range(len(nodes)):
            coeff, optimum = gram_minimum(nodes, weights, d)
            system = LinearSystem(d+1)
            for r in range(d+1):
                assert system.add([sum(w*t**(r+s) for t, w in zip(nodes, weights))
                                   for s in range(d+1)], int(r == 0))
            inverse_column = system.solution()
            assert optimum == 1/inverse_column[0]
            weighted.append(dict(nodes=nodes, degree=d, minimum_energy=str(optimum)))
    return dict(orthogonality_pairs=orthogonality_pairs, exact_optima=cases, weighted_nodes=weighted)


def energy_identity_checks():
    cases = []
    for q in range(1, 6):
        for d in range(q):
            rng = random.Random(9200+17*q+d)
            coeff = [[F(1)] + [F(rng.randrange(-3, 4)) for _ in range(d)]
                     for x in range(1 << q)]
            direct = sum(((evaluate(coeff[x], i+1)+evaluate(coeff[x ^ (1 << i)], i+1))/2)**2
                         for i in range(q) for x in range(1 << q))/(1 << q)
            fourier = [[sum(parity(s & x)*coeff[x][k] for x in range(1 << q))/(1 << q)
                        for k in range(d+1)] for s in range(1 << q)]
            parseval = sum(evaluate(fourier[s], i+1)**2 for s in range(1 << q)
                           for i in range(q) if not s & (1 << i))
            assert direct == parseval
            r = rho(q, d)
            assert direct >= r/(1-r)
            cases.append(dict(q=q, degree=d, energy=str(direct)))
    # Degree Q removes the obstruction: all nodes can be roots.
    for q in range(1, 6):
        vals = []
        for t in range(q+1):
            value = F(1)
            for i in range(1, q+1):
                value *= 1-F(t, i)
            vals.append(value)
        assert vals == [1]+[0]*q
    # Energy alone does not count support: one failed edge can be made arbitrarily large.
    # On the two-cube, perturb n only at x=0 by c*T*(T-2). Only direction 1 is affected.
    q, c = 2, 100
    base = [[F(1), -F(3, 2), F(1, 2)] for _ in range(4)]  # roots 1,2
    base[0][1] -= 2*c
    base[0][2] += c
    defects = [(evaluate(base[x], i+1)+evaluate(base[y], i+1))/2
               for x, y, i in edges(q)]
    assert sum(v != 0 for v in defects) == 1 and max(v*v for v in defects) == F(c*c, 4)
    return dict(parseval_cases=cases, degree_threshold_control=True,
                energy_does_not_count_support_without_amplitude_control=True)


def scale_checks():
    count = 0
    for q in [2, 8, 31, 100, 500]:
        for d in range(q):
            logarithm = sum(log(q-j)-log(q+j+1) for j in range(d+1))
            assert -(d+1)**2/(q-d)-1e-10 <= logarithm
            assert logarithm <= -(d+1)**2/(q+d+1)+1e-10
            count += 1
    return dict(cases=count, note="Numerical finite check of proved logarithmic bounds; all other checks exact.")


def main():
    results = {}
    for name, check in [("trace_multiplication", trace_checks), ("residue_fields", field_checks),
                        ("orthogonal_polynomials", orthogonal_and_gram_checks),
                        ("energy_identity", energy_identity_checks), ("scales", scale_checks)]:
        results[name] = check()
        print(name+": OK", flush=True)
    results["all_ok"] = True
    output = Path(__file__).resolve().parent.parent / "results" / "defect_capacity.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(results, indent=2)+"\n", encoding="utf-8")
    print("ALL OK")


if __name__ == "__main__":
    main()

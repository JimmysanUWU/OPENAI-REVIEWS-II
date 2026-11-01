"""Exact finite checks for 15_UNGRADED_SPLINE_BOUND.md.

Standard library only. These checks test the mechanisms independently of the
asymptotic proof, including genuinely ungraded feasible designs. They do not
constitute a formal proof. Run from any directory with Python 3.10+.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, exp
from pathlib import Path
import json
import random


class LinearSystem:
    """Incremental exact echelon form. Rejected equations do not mutate it."""
    def __init__(self, n):
        self.n = n
        self.rows = {}

    def add(self, row, rhs=0):
        row = [F(x) for x in row] + [F(rhs)]
        for pivot in sorted(self.rows):
            if row[pivot]:
                a = row[pivot]
                row = [x-a*y for x, y in zip(row, self.rows[pivot])]
        pivot = next((i for i in range(self.n) if row[i]), None)
        if pivot is None:
            return row[-1] == 0
        a = row[pivot]
        self.rows[pivot] = [x/a for x in row]
        return True

    def solution(self):
        answer = [F(0)] * self.n
        for pivot in sorted(self.rows, reverse=True):
            row = self.rows[pivot]
            answer[pivot] = row[-1]-sum(row[j]*answer[j]
                                        for j in range(pivot+1, self.n))
        return answer


def rank(rows, n):
    system = LinearSystem(n)
    for row in rows:
        assert system.add(row)
    return len(system.rows)


def parity(x):
    return (-1)**x.bit_count()


def edges(q):
    return [(x, x ^ (1 << t), t) for t in range(q)
            for x in range(1 << q) if not x & (1 << t)]


def constraint(q, d, nodes, edge, jet=0):
    x, y, t = edge
    row = [0]*((1 << q)*(d+1))
    for k in range(jet, d+1):
        v = comb(k, jet)*nodes[t]**(k-jet)
        row[x*(d+1)+k] = v
        row[y*(d+1)+k] = -v
    return row


def trace_rows(q, d):
    rows = []
    for x in range(1 << q):
        row = [0]*((1 << q)*(d+1))
        row[x*(d+1)] = 1
        rows.append(row)
    return rows


def evaluate(coeff, t, jet=0):
    return sum(comb(k, jet)*a*t**(k-jet)
               for k, a in enumerate(coeff) if k >= jet)


def multiply(a, b):
    result = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def spline_basis(q, nodes, weights, subset):
    poly = [F(1)]
    for t in range(q):
        if subset & (1 << t):
            for _ in range(weights[t]):
                poly = multiply(poly, [1, -F(1, nodes[t])])
    return [[parity(subset & x)*a for a in poly] for x in range(1 << q)]


def tail(q, degree):
    return F(sum(comb(q, j) for j in range(q+1)
                 if 2*j > q+degree), 1 << q)


def full_cube_checks():
    cases = []
    for q, weights in [(1, [1]), (2, [1, 1]), (3, [1]*3),
                       (4, [1]*4), (2, [1, 2]), (3, [1, 2, 1])]:
        nodes = [2, -3, 5, 11][:q]
        subset_weights = [sum(weights[t] for t in range(q) if s & (1 << t))
                          for s in range(1 << q)]
        for s in range(1 << q):
            g = spline_basis(q, nodes, weights, s)
            assert [a[0] for a in g] == [parity(s & x) for x in range(1 << q)]
            for x, y, t in edges(q):
                for j in range(weights[t]):
                    assert evaluate(g[x], nodes[t], j) == evaluate(g[y], nodes[t], j)
        for d in range(sum(weights)+1):
            rows = [constraint(q, d, nodes, e, j) for e in edges(q)
                    for j in range(weights[e[2]])]
            n = (1 << q)*(d+1)
            r = rank(rows, n)
            module_dim = n-r
            trace_dim = rank(rows+trace_rows(q, d), n)-r
            assert module_dim == sum(max(0, d-w+1) for w in subset_weights)
            assert trace_dim == sum(w <= d for w in subset_weights)
            cases.append(dict(q=q, weights=weights, degree=d,
                              module_dimension=module_dim, trace_dimension=trace_dim))
    return cases


def deletion_checks():
    rng = random.Random(20261101)
    cases = []
    for q, weights in [(2, [1, 1]), (3, [1, 1, 1]), (3, [1, 2, 1])]:
        nodes = [1, 3, 7][:q]
        all_constraints = [(e, j) for e in edges(q) for j in range(weights[e[2]])]
        for d in range(sum(weights)+1):
            n = (1 << q)*(d+1)
            full_trace = sum(sum(weights[t] for t in range(q) if s & (1 << t)) <= d
                             for s in range(1 << q))
            for _ in range(4):
                # Retain an initial jet segment on each edge: multiplication remains valid.
                keep = {e: rng.randrange(weights[e[2]]+1) for e in edges(q)}
                kept = [(e, j) for e, j in all_constraints if j < keep[e]]
                rows = [constraint(q, d, nodes, e, j) for e, j in kept]
                r = rank(rows, n)
                trace = rank(rows+trace_rows(q, d), n)-r
                erased = len(all_constraints)-len(kept)
                assert trace <= full_trace+erased
                cases.append(dict(q=q, degree=d, erased=erased, trace_dimension=trace))
    return cases


def greedy_design(q, degree, seed, weights=None):
    """Choose a consistent edge-jet system with prescribed parity at T=0."""
    weights = weights or [1]*q
    nodes = [1, 3, 7, 11, 17][:q]
    n = (1 << q)*(degree+1)
    system = LinearSystem(n)
    for x, row in enumerate(trace_rows(q, degree)):
        assert system.add(row, parity(x))
    order = edges(q)
    random.Random(seed).shuffle(order)
    for e in order:
        for j in range(weights[e[2]]):
            if not system.add(constraint(q, degree, nodes, e, j)):
                break
    flat = system.solution()
    coeff = [flat[x*(degree+1):(x+1)*(degree+1)] for x in range(1 << q)]
    assert all(coeff[x][0] == parity(x) for x in range(1 << q))
    bad = [e for e in edges(q)
           if evaluate(coeff[e[0]], nodes[e[2]]) != evaluate(coeff[e[1]], nodes[e[2]])]
    return nodes, coeff, bad


def biased_mass(q, bad, probabilities):
    events = [{x for a, b, t2 in bad if t2 == t for x in (a, b)} for t in range(q)]
    masses = [sum(product_mass(x, probabilities) for x in event) for event in events]
    return masses


def product_mass(x, probabilities):
    result = F(1)
    for t, p in enumerate(probabilities):
        result *= p if x & (1 << t) else 1-p
    return result


def expectation_tail(q, degree, probabilities):
    return sum(product_mass(z, [2*p for p in probabilities])*tail(z.bit_count(), degree)
               for z in range(1 << q))


def ungraded(coeff, q, degree):
    # n_x = parity(x) f_x. Grading forbids nonconstant c_U for |U^c| >= R.
    return any(sum(parity(u & x)*parity(x)*coeff[x][k] for x in range(1 << q)) != 0
               for u in range(1 << q) if q-u.bit_count() >= degree+1
               for k in range(1, degree+1))


def design_checks():
    cases = []
    genuinely_ungraded = 0
    for q in range(2, 6):
        for degree in range(q):
            for seed in range(3):
                nodes, coeff, bad = greedy_design(q, degree, seed)
                lower = (1 << q)*tail(q, degree)
                assert len(bad) >= lower
                is_ungraded = ungraded(coeff, q, degree)
                genuinely_ungraded += is_ungraded
                d = (q+degree)//2
                # Explicit character completion, including every high Fourier character.
                for s in range(1 << q):
                    if s.bit_count() <= d:
                        completed = spline_basis(q, nodes, [1]*q, s)
                    else:
                        g = spline_basis(q, nodes, [1]*q, ((1 << q)-1) ^ s)
                        completed = [multiply(coeff[x], g[x]) for x in range(1 << q)]
                    assert all(len(p) <= d+1 for p in completed)
                    assert [p[0] for p in completed] == [parity(s & x) for x in range(1 << q)]
                    for x, y, t in edges(q):
                        if (x, y, t) not in bad:
                            assert evaluate(completed[x], nodes[t]) == evaluate(completed[y], nodes[t])
                biases = [[F(1, 4)]*q, [F(1, 3)]*q, [F(1, 2)]*q,
                          [F(1, 4) if t % 2 else F(2, 5) for t in range(q)]]
                for probs in biases:
                    masses = biased_mass(q, bad, probs)
                    rhs = expectation_tail(q, degree, probs)
                    assert sum(p*m for p, m in zip(probs, masses)) >= rhs
                    # Independently enumerate the coupling to catch its factor of 2.
                    coupling = F(0)
                    for z in range(1 << q):
                        active = [t for t in range(q) if z & (1 << t)]
                        vertices = [x for x in range(1 << q) if not x & ~z]
                        restricted = sum(1 for a, b, t in bad
                                         if t in active and a in vertices and b in vertices)
                        coupling += product_mass(z, [2*p for p in probs])*F(2*restricted, len(vertices))
                    assert coupling == 2*sum(p*m for p, m in zip(probs, masses))
                cases.append(dict(q=q, degree=degree, seed=seed, bad_edges=len(bad),
                                  lower_bound=int(lower), genuinely_ungraded=is_ungraded))
    assert genuinely_ungraded > 0
    return dict(cases=cases, genuinely_ungraded=genuinely_ungraded)


def weighted_checks():
    cases = []
    for q, weights in [(2, [1, 2]), (3, [1, 2, 1]), (3, [2, 2, 2])]:
        total = sum(weights)
        for degree in range(total):
            nodes, coeff, _ = greedy_design(q, degree, 47+degree, weights)
            erased = 0
            for x, y, t in edges(q):
                matched = 0
                while matched < weights[t] and evaluate(coeff[x], nodes[t], matched) == evaluate(coeff[y], nodes[t], matched):
                    matched += 1
                erased += weights[t]-matched
            d = (total+degree)//2
            lower = sum(sum(weights[t] for t in range(q) if s & (1 << t)) > d
                        for s in range(1 << q))
            assert erased >= lower
            cases.append(dict(weights=weights, degree=degree, erased=erased, lower_bound=lower))
    return cases


def single_edge_checks():
    cases = []
    for q in range(1, 7):
        nodes = list(range(1, q+1))
        k = q-1
        poly = [F(1)]
        for t in range(q):
            if t != k:
                poly = multiply(poly, [1, -F(1, nodes[t])])
        f = [[F(0)]*q for _ in range(1 << q)]
        f[0] = [(1 << q)*a for a in poly]
        for s in range((1 << q)-1):
            g = spline_basis(q, nodes, [1]*q, s)
            for x in range(1 << q):
                for j, a in enumerate(g[x]):
                    f[x][j] -= a
        assert [a[0] for a in f] == [parity(x) for x in range(1 << q)]
        bad = [e for e in edges(q) if evaluate(f[e[0]], nodes[e[2]]) != evaluate(f[e[1]], nodes[e[2]])]
        assert bad == [(0, 1 << k, k)]
        assert (1 << q)*tail(q, q-1) == 1
        assert q == 1 or ungraded(f, q, q-1)
        cases.append(dict(q=q, degree=q-1, bad_edge=list(bad[0]), sharp_total_bound=1))
    return cases


def negative_controls():
    # Coincident edge labels: f_x(T)=parity(x)(1-T) obeys every edge at T=1.
    assert tail(3, 1) > 0
    bad_collision = [(x, y, t) for x, y, t in edges(3)
                     if evaluate([parity(x), -parity(x)], 1) != evaluate([parity(y), -parity(y)], 1)]
    assert not bad_collision
    # With no parity boundary value, the zero spline violates the claimed obstruction.
    assert 0 < (1 << 3)*tail(3, 1)
    # Ignoring degree: the full character spline has parity trace and no bad edges.
    f = spline_basis(3, [1, 2, 3], [1]*3, 7)
    assert all(evaluate(f[x], t+1) == evaluate(f[y], t+1) for x, y, t in edges(3))
    assert len(f[0])-1 == 3
    return {"distinct_nodes_needed": True, "parity_trace_needed": True,
            "degree_bound_needed": True}


def tail_checks():
    cases = 0
    for m in [256, 257, 320, 512, 1024]:
        for degree in range(m//8+1):
            assert float(tail(m, degree)) >= .75*exp(-25-2*degree*degree/m)
            cases += 1
    return {"cases": cases, "note": "Finite numerical check of an independently proved tail bound."}


def main():
    results = {}
    for name, fn in [("full_cube", full_cube_checks), ("deleted_constraints", deletion_checks),
                     ("ungraded_designs", design_checks), ("weighted_jets", weighted_checks),
                     ("single_edge_extremizers", single_edge_checks),
                     ("negative_controls", negative_controls), ("tail_bound", tail_checks)]:
        results[name] = fn()
        print(name+": OK", flush=True)
    results["all_ok"] = True
    path = Path(__file__).resolve().parent.parent / "results" / "spline_bound.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(results, indent=2)+"\n", encoding="utf-8")
    print("ALL OK")


if __name__ == "__main__":
    main()

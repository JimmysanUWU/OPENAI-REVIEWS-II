"""Exhaustive tiny-instance exploration of T.ETH.7; standard library only.

For each cap, enumerate every permitted deleted-edge set, then decide the
remaining polynomial equations over the rationals. Rational linear systems
are feasible over R iff they are feasible over Q. This does not settle an
asymptotic directional bound. Not part of verify_all's verification suite.
"""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import json

from verify_spline_bound import LinearSystem, edges, evaluate


def search(q, degree):
    all_edges = edges(q)
    directions = [[e for e in all_edges if e[2] == i] for i in range(q)]
    tested = 0
    for cap in range(1+(1 << (q-1))):
        choices = [[c for j in range(cap+1) for c in combinations(es, j)] for es in directions]
        for parts in product(*choices):
            deleted = {e for part in parts for e in part}
            system = LinearSystem(degree*(1 << q))
            consistent = True
            for x, y, i in all_edges:
                if (x, y, i) in deleted:
                    continue
                row = [0]*(degree*(1 << q))
                for k in range(1, degree+1):
                    row[x*degree+k-1] = row[y*degree+k-1] = (i+1)**k
                if not system.add(row, -2):
                    consistent = False
                    break
            tested += 1
            if not consistent:
                continue
            flat = system.solution()
            polynomials = [[F(1)]+flat[x*degree:(x+1)*degree] for x in range(1 << q)]
            bad = [e for e in all_edges if evaluate(polynomials[e[0]], e[2]+1)+
                   evaluate(polynomials[e[1]], e[2]+1) != 0]
            counts = [sum(e[2] == i for e in bad) for i in range(q)]
            assert max(counts) == cap and set(bad) <= deleted
            return dict(q=q, degree=degree, nodes=list(range(1, q+1)), minimum_max_edges=cap,
                        minimum_max_uniform_mass=str(F(cap, 1 << (q-1))),
                        systems_examined=tested, bad_edges=bad,
                        coefficients=[[str(c) for c in p] for p in polynomials])
    raise AssertionError("The empty constraint system must be feasible")


def main():
    cases = [search(q, d) for q, d in [(2, 1), (3, 1), (3, 2)]]
    assert [x["minimum_max_edges"] for x in cases] == [1, 2, 1]
    output = Path(__file__).resolve().parent.parent / "results" / "directional_extremizers.json"
    output.write_text(json.dumps(dict(cases=cases, scope="Exhaustive only for the three stated tiny instances."),
                                 indent=2)+"\n", encoding="utf-8")
    for case in cases:
        print("Q={q}, D={degree}: min maximum directional mass = {minimum_max_uniform_mass}".format(**case))
    print("Exploration complete")


if __name__ == "__main__":
    main()

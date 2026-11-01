"""Exact counterexample to the per-row premise used in Toolkit E4.

This does not disprove the asymptotic floor; it shows a transfer lemma is needed.
"""
from itertools import product
from fractions import Fraction
import json
from pathlib import Path

def main():
    r, Q = 2, 6
    pairs = list(product(range(r), repeat=2))
    labels = list(product(pairs, repeat=Q))
    D = len(labels)
    full_copies = 4 * D
    # One unary clause on a different variable for each occupied position.
    # Every label occurs full_copies times, with one additional position at u0.
    u0 = ((0, 1),) * Q
    n = M = full_copies * D + 1
    h = full_copies + 1
    assert 2 * (1+1) * D**2 <= M  # balanced-schedule lemma, b=w=1
    assert 2 * Q * r**2 <= n
    assert h <= Fraction(2*M, D)
    assert full_copies*D//r + 1 <= Fraction(2*n, r)

    def anti(pair):
        return 1 if pair == (0, 1) else -1 if pair == (1, 0) else 0

    def sym(pair):
        return int(pair in {(0, 1), (1, 0)})

    psi = {u: anti(u[0])*anti(u[1])*anti(u[2])*sym(u[3])*sym(u[4])*sym(u[5]) for u in labels}
    label_norm = sum(abs(a) for a in psi.values())
    assert label_norm == 64
    # At node 1 either nonzero abstract group restriction has norm 32.
    abstract_group_norm = sum(abs(psi[u]) for u in labels if u[0][0] == 0)
    assert abstract_group_norm == 32
    # The final copy contains just u0, hence the actual row is x_v, of norm 1.
    transferred_group_norm = abs(psi[u0])
    assert transferred_group_norm == 1
    assert transferred_group_norm < 2**(3-1)  # j=3 antisymmetric factors

    # A different measured eigenvector (j=3) vanishes on this occupied position.
    phi = {u: int(u[0] == (0, 0))*anti(u[1])*anti(u[2])*anti(u[3])*sym(u[4])*sym(u[5]) for u in labels}
    assert sum(abs(a) for a in phi.values()) == 32
    assert phi[u0] == 0
    # With R=2 and the central agreement window {1,2}, j=3 is measured.
    # R=2 <= Q/e, as required in the nonconstant-cost part of E4's argument.
    # The actual nonzero row has 2 possible values, versus 33 for the abstract
    # norm-32 restriction. Thus per-row measurement bit counts do not transfer.
    result = dict(r=r, Q=Q, D=D, n=n, h=h, balanced_hypotheses=True,
                  abstract_row_l1=label_norm, abstract_group_l1=abstract_group_norm,
                  actual_group_l1=transferred_group_norm,
                  other_measured_row_vanishes=True,
                  conclusion='Per-row premise fails; aggregate floor remains unresolved by this example.')
    print(json.dumps(result, indent=2))
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print('ALL OK')

if __name__ == '__main__':
    main()

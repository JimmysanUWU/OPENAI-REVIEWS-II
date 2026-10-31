"""Exact checks for the complete-copy interface; no external dependencies."""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import json


def local_basis(r):
    for a in range(r):
        yield 0, {(a, a): 1}
    for a in range(r):
        for b in range(a+1, r):
            yield 0, {(a, b): 1, (b, a): 1}
            yield 1, {(a, b): 1, (b, a): -1}


def tensor_row(factors):
    row = {}
    for entries in product(*(list(f[1].items()) for f in factors)):
        label = tuple(e[0] for e in entries)
        coefficient = 1
        for _, v in entries:
            coefficient *= v
        row[label] = coefficient
    return row


def transfer(row, occupied, t):
    # Identify ALL positions with the same first-coordinate schedule as one variable.
    # This deliberately introduces repeated incidences, which must not cause cancellation.
    groups = {}
    for label, value in row.items():
        if label not in occupied:
            continue
        variable = tuple(pair[0] for pair in label)
        group = variable[t]
        groups.setdefault(group, {})
        groups[group][variable] = groups[group].get(variable, 0)+value
    return {g: {v: c for v, c in a.items() if c} for g, a in groups.items()
            if any(a.values())}


def main():
    cases = []
    for r, q in [(2, 1), (2, 3), (2, 4), (3, 2), (3, 3)]:
        labels = set(product(list(product(range(r), repeat=2)), repeat=q))
        spectral = [0]*(q+1)
        row_checks = 0
        heavy_rows = [0]*q
        max_norm = 1
        for factors in product(list(local_basis(r)), repeat=q):
            j = sum(f[0] for f in factors)
            spectral[j] += 1
            row = tensor_row(factors)
            schedules = [tuple(pair[0] for pair in label) for label in row]
            assert len(schedules) == len(set(schedules))
            for t in range(q):
                groups = transfer(row, labels, t)
                assert 1 <= len(groups) <= 2
                for group, values in groups.items():
                    abstract_norm = sum(abs(a) for u, a in row.items() if u[t][0] == group)
                    norm = sum(abs(a) for a in values.values())
                    assert norm == abstract_norm
                    assert norm >= 2**max(j-1, 0)
                    max_norm = max(max_norm, norm)
                if j >= max(1, q//2):
                    heavy_rows[t] += len(groups)
                row_checks += 1
        p = F(r-1, 2*r)
        label_count = r**(2*q)
        assert spectral == [label_count*comb(q, j)*p**j*(1-p)**(q-j) for j in range(q+1)]
        measured = sum(spectral[max(1, q//2):])
        copies = 3
        assert all(copies*c >= copies*measured for c in heavy_rows)
        cases.append(dict(r=r, q=q, eigenvectors=label_count, row_node_checks=row_checks,
                          spectral_multiplicities=spectral, heavy_eigenvectors=measured,
                          group_rows_per_complete_copy=heavy_rows, max_norm=max_norm))

    # Same sparse-copy failure isolated in the previous audit: no full-copy premise.
    anti = (1, {(0, 1): 1, (1, 0): -1})
    sym = (0, {(0, 1): 1, (1, 0): 1})
    diag = (0, {(0, 0): 1})
    u0 = ((0, 1),)*6
    psi = tensor_row([anti]*3+[sym]*3)
    sparse = transfer(psi, {u0}, 0)
    assert sum(abs(a) for values in sparse.values() for a in values.values()) == 1 < 2**(3-1)
    phi = tensor_row([diag]+[anti]*3+[sym]*2)
    assert not transfer(phi, {u0}, 0)

    # Variance-slack implication checked with integer arithmetic (no square roots).
    variance_cases = 0
    for label_count in [4, 8, 16, 32]:
        bw = 2
        mean = 8*bw*label_count+1
        counts = [mean+(i % 3)-1 for i in range(label_count)]
        m = sum(counts)
        variance = sum((F(c)-F(m, label_count))**2 for c in counts)
        if variance < 2*m*bw and m >= 8*bw*label_count**2:
            assert min(counts)*2*label_count >= m
            variance_cases += 1
    assert variance_cases > 0
    result = dict(cases=cases, variance_cases=variance_cases,
                  negative_controls={"sparse_copy_norm_loss": True, "sparse_copy_zero_row": True},
                  all_ok=True)
    path = Path(__file__).resolve().parent.parent / "results" / "complete_copies.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print("support injection, norms, multiplicities, aggregate rows, variance: OK")
    print("sparse-copy negative controls: OK")
    print("ALL OK")


if __name__ == "__main__":
    main()

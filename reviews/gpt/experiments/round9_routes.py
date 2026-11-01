"""Deterministic probes retained from abstraction round 9; standard library.

Exploration, not an asymptotic proof or part of the verification runner.
"""
from fractions import Fraction as F
from math import prod
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "reviews/gpt/eth/scripts"))
from verify_biased_energy import capacity


def energy(theta):
    q, result = len(theta), F(0)
    for mask in range(1 << q):
        ix = [i for i in range(q) if mask >> i & 1]
        probability = prod(theta[i] if mask >> i & 1 else 1 - theta[i] for i in range(q))
        result += probability * capacity([F(i + 1) for i in ix],
                                         [1 / theta[i] for i in ix], 1)[0]
    return result


def main():
    a = [F(3, 4), F(1, 4), F(3, 4)]
    b = [F(3, 4), F(3, 4), F(1, 4)]
    midpoint = [(x + y) / 2 for x, y in zip(a, b)]
    ea, eb, em = map(energy, [a, b, midpoint])
    gap = em - (ea + eb) / 2
    assert gap == F(-9374301, 11806674880) < 0
    # Local coordinate concavity survives: changing just one coordinate does
    # not give this counterexample. The false step was joint Jensen.
    c = [F(1, 4), a[1], a[2]]
    assert energy([(x + y) / 2 for x, y in zip(a, c)]) >= (energy(a) + energy(c)) / 2

    def even(t):
        return (1 + t / 3) / (1 - t / 3)

    def odd(t):
        return 1 - 3 * t

    assert even(F(0)) == odd(F(0)) == 1
    assert all(even(F(t)) + odd(F(t)) == 0 for t in [1, 2])
    for t in [F(j, 100) for j in range(201)]:
        assert 1 - t / 3 > 0
    # Common denominator exposes the degree hidden by the rational form.
    for t in map(F, [-1, 0, 1, 2, 4]):
        assert (even(t) + odd(t)) * (1 - t / 3) == 2 - 3 * t + t ** 2
    assert not all(2 - 3 * F(t) == 0 for t in [1, 2])
    assert all(2 - 3 * F(t) + F(t) ** 2 == 0 for t in [1, 2])

    height_cases = []
    for q in range(2, 9):
        slopes = [(-1) ** x.bit_count() * 2 * sum(F(1, i + 1) for i in range(q) if x >> i & 1)
                  for x in range(1 << q)]
        bad = []
        for i in range(q):
            count = sum(2 + (slopes[x] + slopes[x | (1 << i)]) * (i + 1) != 0
                        for x in range(1 << q) if not (x >> i & 1))
            assert count == 1 << (q - 2)
            bad.append(count)
        height_cases.append(dict(Q=q, bad_edges_by_direction=bad))

    quenched = energy([F(1, 2)] * 3)
    annealed = capacity(list(map(F, [1, 2, 3])), [F(1)] * 3, 1)[0]
    assert quenched == F(503, 1820) and annealed == F(3, 7) and quenched < annealed
    # Weyl central amplitude (lambda + lambda^-1)/2: uncharged period is
    # its Laurent coefficient at zero, which vanishes. A charged period
    # survives but has not been turned into conserved quotient data.
    amplitude = {1: F(1, 2), -1: F(1, 2)}
    assert amplitude.get(0, 0) == 0
    assert amplitude.get(1, 0) == F(1, 2)

    result = dict(
        joint_concavity=dict(theta_a=list(map(str, a)), theta_b=list(map(str, b)),
                             E_a=str(ea), E_b=str(eb), E_midpoint=str(em), Jensen_gap=str(gap)),
        rational_pair=dict(even="(1+T/3)/(1-T/3)", odd="1-3T", nodes=[1, 2],
                           common_numerator="2-3T+T^2", no_poles_on="[0,2]"),
        degree_one_heights=height_cases,
        annealed_shortcut=dict(exact=str(quenched), ratio_of_expected_determinants=str(annealed)),
        central_period=dict(uncharged="0", charged_coefficient="1/2",
                            conserved_charged_quotient="not constructed"),
        scope="Exact route probes only; no sharp directional theorem or new counting algorithm.")
    path = Path(__file__).resolve().parent / "results/round9_routes.json"
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("Exact joint-Jensen counterexample, rational degree accounting, height designs, and period probes: OK")
    print("Exploration complete")


if __name__ == "__main__":
    main()

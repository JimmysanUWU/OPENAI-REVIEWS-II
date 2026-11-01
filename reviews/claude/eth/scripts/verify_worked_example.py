"""A worked example for 19_HOW_TO_READ_THE_PROOF.md: Part I's counting algorithm on one tiny formula, with
every intermediate quantity printed.

The formula (4 variables, 4 clauses, every variable in at most 3 clauses):
    (x1 or not x2) and (x2 or x3 or not x4) and (not x1 or x4) and (x3 or x4)
Parameters (r, Q, R) = (2, 2, 2): 2 groups, 2 nodes, deformation degree R - 1 = 1.

Steps printed (the numbering is that of the guide's Level 1):
  1. #F by brute force.
  2. The swap design (08 Theorem 5.4): labels, the agreement window, the Krawtchouk coefficients, the
     measured eigenvectors.
  3. Schedules (08 Lemma 6.2), labels of incidences, copies.
  4. The global polynomial H(T) = sum_x prod_f Phi_f(<P_{f,v}(T), x>): degree bound d*R' < Q*J.
  5. At each node t = 1, 2: groups, measurement rows, the Taylor jet H(t + z) mod z^J computed two ways:
     directly (2^n terms) and by separation into group sums (08 Lemmas 4.4-4.5); they agree.
  6. Recovery of H(0) = #F from the jets by the Chinese remainder theorem (08 Lemma 4.3), mod p^mu.
Controls: dropping one group's measurement rows changes a jet; recovering from too few jet coefficients
(J - 1) gives a wrong count.

Reuses the implementation of verify_gapless_chain.py (Part K).  Standard library only; deterministic;
writes ../results/worked_example.json (about 2 s).
"""

import json
import os
import sys
from fractions import Fraction
from itertools import product
from math import comb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_gapless_chain import (Swap, build, locality_ok, Ring, taylor_separated, taylor_direct,  # noqa: E402
                                  crt_at_zero, prime_above, brute)
from eth_core import poly_eval  # noqa: E402

OUT = {}


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)


def lit(v, sign):
    return ("x%d" % (v + 1)) if sign else ("not x%d" % (v + 1))


def main():
    n = 4
    clauses = [[(0, True), (1, False)], [(1, True), (2, True), (3, False)], [(0, False), (3, True)],
               [(2, True), (3, True)]]
    text = " and ".join("(" + " or ".join(lit(v, s) for v, s in cl) + ")" for cl in clauses)
    count = brute(n, clauses)
    print("1. The formula:", text)
    print("   #F by brute force over all 2^%d assignments: %d" % (n, count))
    OUT["formula"] = text
    OUT["count"] = count

    r, Q, R = 2, 2, 2
    sw = Swap(r, Q, R)
    ch = sw.checks()
    check(all(ch.values()), "swap design checks")
    levels = {j: sw.level_count[j] for j in sorted(sw.level_count)}
    print("\n2. Swap design (r, Q, R) = (%d, %d, %d)" % (r, Q, R))
    print("   labels: ([r] x [r])^Q, D = %d;  p_r = (r-1)/(2r) = %s" % (sw.D, sw.p))
    print("   eigenvectors per antisymmetric level j:", levels, " (= D * Pr[Bin(Q, p_r) = j])")
    print("   agreement window I =", sw.interval, "; Krawtchouk coefficients gamma =", [str(g) for g in sw.gamma])
    print("   measured eigenvectors (levels outside I): c = %d of %d" % (sw.c, sw.D))
    OUT["design"] = dict(r=r, Q=Q, R=R, D=sw.D, p_r=str(sw.p), levels=levels, window=sw.interval,
                         gamma=[str(g) for g in sw.gamma], measured=sw.c, checks=ch)

    inst = build(n, clauses, sw)
    check(locality_ok(inst, sw), "locality at every node")
    print("\n3. Schedules (derandomized, 08 Lemma 6.2)")
    print("   variable schedules a_v:", {("x%d" % (v + 1)): inst["a"][v] for v in range(n)})
    print("   clause schedules d_f:", {("C%d" % (f + 1)): inst["d"][f] for f in range(len(clauses))})
    print("   incidences M = %d, copies h = %d, potential Psi = %.3f (< 2 as the lemma promises)"
          % (inst["M"], inst["h"], inst["bal"]["Psi_final"]))
    print("   locality: at every node t, P_{f,v}(t) is local for the partition into groups, modulo measurements")
    OUT["schedules"] = dict(a=inst["a"], d=inst["d"], M=inst["M"], h=inst["h"], Psi=inst["bal"]["Psi_final"])

    Rp = max(len(pp) - 1 for row in inst["P"].values() for pp in row.values())
    dtot = inst["M"]
    q = (dtot * Rp) // Q
    J = q + 1
    p = prime_above(Q)
    mu = n + 1
    m = p ** mu
    print("\n4. The global polynomial H(T) = sum_x prod_f Phi_f(<P_{f,v}(T), x>), with H(0) = #F")
    print("   deformation degree R' = %d, total clause degree d = %d: deg H <= d R' = %d" % (Rp, dtot, dtot * Rp))
    print("   Taylor depth q = floor(d R'/Q) = %d, J = q + 1 = %d jets per node; Q J = %d > deg H" % (q, J, Q * J))
    print("   arithmetic modulo p^mu = %d^%d = %d > 2^n (p = least prime > Q, so the nodes' differences are units)"
          % (p, mu, m))
    check(dtot * Rp < Q * J, "degree bound")
    check(m > 2 ** n, "modulus exceeds the count")
    OUT["polynomial"] = dict(R_prime=Rp, d=dtot, q=q, J=J, p=p, mu=mu, modulus=m)

    ring = Ring(m, J)
    jets = []
    node_info = []
    for s in range(Q):
        t = s + 1
        sep = taylor_separated(n, inst, s, ring, q)
        dirv = taylor_direct(n, inst, s, ring)
        check(sep == dirv, "separated jet equals direct jet at node %d" % t)
        groups = {g: ["x%d" % (v + 1) for v in inst["groups"][(s, g)]] for g in range(r)}
        rows = {g: len(inst["C"][(s, g)]) for g in range(r)}
        # number of distinct measurement-value tuples per group (the enumeration over y in 08 Lemma 4.5)
        buckets = {}
        for g in range(r):
            S = inst["groups"][(s, g)]
            vals = set()
            for bits in product((0, 1), repeat=len(S)):
                xg = dict(zip(S, bits))
                vals.add(tuple(sum(c * xg[w] for w, c in row.items()) for row in inst["C"][(s, g)]))
            buckets[g] = len(vals)
        corrections = sum(comb(n, j) * 2 ** j for j in range(min(q, n) + 1))
        print("\n5. Node t = %d" % t)
        print("   groups S_{t,g}:", groups, "; measurement rows per group:", rows)
        print("   distinct measurement values per group:", buckets, "(each group sum is computed once per value)")
        print("   correction terms: sum_{j<=q} C(n,j) 2^j = %d (Boolean substitution, 08 Lemma 4.4)" % corrections)
        print("   jet H(t+z) mod (z^J, p^mu), separated == direct:", sep)
        jets.append(sep)
        node_info.append(dict(t=t, groups=groups, rows=rows, values=buckets, corrections=corrections, jet=sep))
    OUT["nodes"] = node_info

    H0 = crt_at_zero(jets, Q, J, m)
    print("\n6. Recovery (Chinese remainder at the nodes, 08 Lemma 4.3): H(0) mod p^mu = %d" % H0)
    check(H0 == count, "recovered count equals brute force")
    print("   equals #F = %d" % count)
    OUT["recovered"] = H0

    # controls
    controls = {}
    s = next(s for s in range(Q) for g in range(r) if inst["C"][(s, g)])
    g = next(g for g in range(r) if inst["C"][(s, g)])
    Cbad = dict(inst["C"])
    Cbad[(s, g)] = []
    bad = taylor_separated(n, inst, s, ring, q, Cm=Cbad)
    controls["drop_measurements_changes_jet"] = bad != jets[s]
    check(controls["drop_measurements_changes_jet"], "control: dropping measurements must change the jet")
    short = [jt[:J - 1] for jt in jets]
    H0_short = crt_at_zero(short, Q, J - 1, m)
    controls["too_few_jets_wrong"] = H0_short != count
    check(controls["too_few_jets_wrong"], "control: J - 1 jets must give a wrong count")
    print("\nControls: dropping the measurement rows of node %d, group %d changes the jet: %s; "
          "recovering from J - 1 jets gives %d != %d: %s"
          % (s + 1, g, controls["drop_measurements_changes_jet"], H0_short, count, controls["too_few_jets_wrong"]))
    OUT["controls"] = controls
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "results", "worked_example.json"), "w") as fh:
        json.dump(OUT, fh, indent=1, default=str)
    print("ALL OK")


if __name__ == "__main__":
    main()

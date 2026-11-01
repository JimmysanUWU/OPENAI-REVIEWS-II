"""Checks for 20_DIRECTIONAL_DEFECTS.md (target T.ETH.7): subcube averaging, and small-cube data.

Setting (GPT 15): vertex polynomials n_x(T) of degree <= D on {0,1}^Q with n_x(0) = 1, distinct nonzero nodes
t_i; the edge {x, x+e_i} is bad if n_x(t_i) + n_{x+e_i}(t_i) != 0; b_i = number of bad edges in direction i.
Equivalently (GPT 15 Corollary 2.2): a spline with parity trace exists on the good-edge graph.

  A  [Proposition 20.1] on designs produced by randomized greedy retention (exact feasibility mod 2^61-1),
     for Q = 3, 4, 5 and every D < Q:
         total bad edges >= (Q/m) 2^{Q-m} N(m, D)   for every m in (D, Q],   N(m,D) = sum_{j>(m+D)/2} C(m,j).
     The case m = Q is GPT 15 Theorem 2.1, m = D+1 is the face lemma (14 Lemma 5.1), m = 1 at D = 0 is tight.
     Control: the same inequality with the right side multiplied by 3/2 fails (tightness at D = 0).
  B  [Corollary 20.2] max_i mu(B_i) >= max_m (2/m) tau_D(m): the gain over GPT's 2 tau_D(Q)/Q, computed exactly
     with fractions for Q up to 400 and D = floor(sqrt(Q))/2, floor(sqrt(Q)), 2 floor(sqrt(Q)); and the constant
     tau_D(max(D+1, D^2)) >= 1/16 for D <= 30 (minimum at D = 2).
  C  [data for Conjecture 20.3] the least max_i b_i found, against N(Q, D) (the total bound) and N(Q, D)/2
     (what graded designs pay per direction, by Smolensky): every design found has max_i b_i >= N(Q, D)/2.

Standard library only; deterministic; writes ../results/directional_defects.json (about 5 s).
"""

import json
import math
import os
import random
from fractions import Fraction

P = (1 << 61) - 1
RESULTS = {}


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)


def N_high(m, D):
    return sum(math.comb(m, j) for j in range(m + 1) if 2 * j > m + D)


def cube_edges(Q):
    return [(x, x | (1 << t), t) for t in range(Q) for x in range(1 << Q) if not (x >> t) & 1]


class Sys:
    def __init__(self, n):
        self.n = n
        self.rows = {}

    def add(self, row):
        row = row[:]
        for piv in sorted(self.rows):
            c = row[piv]
            if c:
                pr = self.rows[piv]
                row = [(a - c * b) % P for a, b in zip(row, pr)]
        piv = next((i for i in range(self.n) if row[i]), None)
        if piv is None:
            return row[self.n] == 0
        inv = pow(row[piv], P - 2, P)
        self.rows[piv] = [a * inv % P for a in row]
        return True

    def copy(self):
        s = Sys(self.n)
        s.rows = dict(self.rows)
        return s


def greedy_design(Q, D, rng, roundrobin):
    """Retain edges greedily while a degree-D spline with parity trace exists; return bad edges per direction."""
    nodes = [3 + 2 * i for i in range(Q)]
    n = (1 << Q) * (D + 1)
    S = Sys(n)
    for x in range(1 << Q):
        row = [0] * (n + 1)
        row[x * (D + 1)] = 1
        row[n] = (-1) ** bin(x).count("1") % P
        check(S.add(row), "trace rows are independent")
    order = cube_edges(Q)
    rng.shuffle(order)
    if roundrobin:
        byd = {t: [e for e in order if e[2] == t] for t in range(Q)}
        order = []
        while any(byd.values()):
            for t in rng.sample(range(Q), Q):
                if byd[t]:
                    order.append(byd[t].pop())
    bad = [0] * Q
    for (x, y, t) in order:
        row = [0] * (n + 1)
        for k in range(D + 1):
            v = pow(nodes[t], k, P)
            row[x * (D + 1) + k] = v
            row[y * (D + 1) + k] = (-v) % P
        T = S.copy()
        if T.add(row):
            S = T
        else:
            bad[t] += 1
    return bad


def part_A_C():
    rng = random.Random(20261101)
    cases = []
    control_violations = 0
    for Q in (3, 4, 5):
        for D in range(Q):
            trials = 40 if Q < 5 else 20
            best = None
            totals = []
            for tr in range(trials):
                bad = greedy_design(Q, D, rng, roundrobin=(tr % 2 == 1))
                tot = sum(bad)
                totals.append(tot)
                for m in range(D + 1, Q + 1):
                    rhs = Fraction(Q, m) * 2 ** (Q - m) * N_high(m, D)
                    check(tot >= rhs, "Proposition 20.1 violated at Q=%d D=%d m=%d" % (Q, D, m))
                    if tot < Fraction(3, 2) * rhs:
                        control_violations += 1
                # data for the conjecture
                check(2 * max(bad) >= N_high(Q, D), "a found design beats graded designs per direction")
                if best is None or max(bad) < max(best):
                    best = bad
            bound = max(Fraction(Q, m) * 2 ** (Q - m) * N_high(m, D) for m in range(D + 1, Q + 1))
            cases.append(dict(Q=Q, D=D, N_high=N_high(Q, D), best_total_bound=str(bound), min_total_found=min(totals),
                              least_max_direction_found=max(best), bad_per_direction=best))
    check(control_violations > 0, "control: the inequality with factor 3/2 should fail somewhere")
    RESULTS["A_C"] = dict(cases=cases, control_violations_of_3_over_2=control_violations)
    print("Part A: subcube averaging holds on every greedy design (Q = 3..5, all D); control (x3/2) fails %d times: OK"
          % control_violations)
    print("Part C: least largest-direction defect found vs N(Q,D):")
    for c in cases:
        print("   Q=%d D=%d  N(Q,D)=%d  least max_i b_i found=%d  per direction %s"
              % (c["Q"], c["D"], c["N_high"], c["least_max_direction_found"], c["bad_per_direction"]))


def tau(m, D):
    return Fraction(N_high(m, D), 2 ** m)


def part_B():
    rows = []
    for Q in (16, 64, 144, 256, 400):
        s = math.isqrt(Q)
        for D in sorted({max(1, s // 2), s, 2 * s}):
            if D >= Q:
                continue
            gpt = 2 * tau(Q, D) / Q
            best = max((Fraction(2, m) * tau(m, D), m) for m in range(D + 1, Q + 1))
            check(best[0] >= gpt, "the family includes m = Q")
            rows.append(dict(Q=Q, D=D, gpt_bound=float(gpt), averaged_bound=float(best[0]), best_m=best[1],
                             gain=float(best[0] / gpt)))
    # the constant in Corollary 20.2: tau_D(max(D+1, D^2)) >= 1/16 (exact for D <= 30; minimum at D = 2)
    consts = {}
    for D in range(1, 31):
        m = max(D + 1, D * D)
        consts[D] = tau(m, D)
        check(consts[D] >= Fraction(1, 16), "tau_D(max(D+1,D^2)) >= 1/16")
    check(min(consts.values()) == Fraction(1, 16) and consts[2] == Fraction(1, 16), "minimum 1/16 at D = 2")
    RESULTS["B"] = dict(rows=rows, tau_at_m_D2={str(D): float(v) for D, v in consts.items()})
    print("Part B: max_m (2/m) tau_D(m) versus 2 tau_D(Q)/Q:")
    for r in rows:
        print("   Q=%3d D=%2d  gain x%.2f at m=%d" % (r["Q"], r["D"], r["gain"], r["best_m"]))


def main():
    part_A_C()
    part_B()
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "results", "directional_defects.json"), "w") as fh:
        json.dump(RESULTS, fh, indent=1)
    print("ALL OK")


if __name__ == "__main__":
    main()

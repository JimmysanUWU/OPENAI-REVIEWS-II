"""Checks for 14_DESIGN_FLOOR.md: every single-polynomial swap design pays the cost floor.

Part A  [Lemma 2.1, Smolensky over Q] a real polynomial of degree <= D on {0,1}^Q agrees with parity on at
        most N_D(Q) = sum_{i <= (Q+D)/2} binom(Q,i) points.  For Q <= 6 and every D: random sets of N_D + 1
        points admit no degree-<=D interpolant of parity (exact rational elimination), while the central
        windows of D+1 levels do.  Control: the same random sets are often interpolable in degree D + 2.
Part B  [Lemma 2.2, restriction] for p <= 1/2 and every F of degree <= D,
            mu_p(F != parity) >= E_{m ~ Bin(Q,2p)} tau_D(m),   tau_D(m) = P[Bin(m,1/2) > (m+D)/2].
        Checked exactly for window designs and for greedy designs that interpolate parity on maximal random
        sets and on maximal mu_p-heavy sets (Q <= 6, p in {1/4, 1/3, 3/8, 1/2}).  The coupling identity
        mu_p(E) = E_Z[bad fraction on the subcube Z] is checked exactly.  Control: window designs of degree
        D + 2 violate the degree-D bound in some cases, so the bound depends on the degree.
Part C  [Lemma 2.3] the explicit forms: tau_D(m) >= e^{-(D+2)^2/m}/(m+1) for D + 2 <= m <= 300, and the
        low-weight tail P[Bin(m,1/2) < Q/32] <= e^{-0.229 Q} for Q/2 <= m <= Q <= 2000.
Part D  [Theorem 3.2] the floor model for all single-polynomial designs: minimize
        max(x log2(1/x), 1/r, Q eps'(Q,R)/64) over integers (r, Q, R) with r^(2Q) <= n, where eps' is the
        all-design lower bound of Lemma 2.4; compare with (log2 log2 n)^2/sqrt(log2 n), log2 n = 2^10..2^300.
Part E  [§5, evidence on the open ungraded question] for (Q, R) = (6, 3), (7, 3) (and (8, 4) with --full,
        about 3 minutes): general (ungraded) group-algebra designs meet the node conditions on exactly the
        same central windows as single-polynomial designs (exact elimination).

Standard library only; deterministic; writes ../results/design_floor.json (about 1 minute).
"""

import json
import math
import random
import time
from fractions import Fraction
from itertools import combinations, product

RESULTS = {}


# =============================================================== exact incremental elimination

class Eliminator:
    """A linear system over Q kept in reduced row echelon form; the right-hand side is the last entry of a row.

    add(row) returns False, and leaves the state unchanged, when the row makes the system inconsistent."""

    def __init__(self, ncols):
        self.n = ncols
        self.rows = []      # (pivot column, row as a list of n + 1 Fractions)

    def add(self, row):
        r = [Fraction(v) for v in row]
        for c, b in self.rows:
            f = r[c]
            if f:
                r = [x - f * y for x, y in zip(r, b)]
        c = next((i for i in range(self.n) if r[i]), None)
        if c is None:
            return r[self.n] == 0
        pv = r[c]
        r = [x / pv for x in r]
        for k, (c2, b) in enumerate(self.rows):
            f = b[c]
            if f:
                self.rows[k] = (c2, [x - f * y for x, y in zip(b, r)])
        self.rows.append((c, r))
        return True

    def solution(self):
        """The solution with all free variables set to zero."""
        coef = [Fraction(0)] * self.n
        for c, b in self.rows:
            coef[c] = b[self.n]
        return coef


def monomials(Q, D):
    return [A for k in range(D + 1) for A in combinations(range(Q), k)]


def chi(A, x):
    return -1 if sum(x[a] for a in A) % 2 else 1


def parity(x):
    return -1 if sum(x) % 2 else 1


def interp_row(ms, x):
    return [chi(A, x) for A in ms] + [parity(x)]


def interpolable(S, Q, D):
    """Is there a real polynomial of degree <= D equal to parity on S?"""
    ms = monomials(Q, D)
    el = Eliminator(len(ms))
    return all(el.add(interp_row(ms, x)) for x in S)


def N_D(Q, D):
    return sum(math.comb(Q, i) for i in range((Q + D) // 2 + 1))


# =============================================================== Part A: Smolensky over Q

def part_A():
    rng = random.Random(1401)
    tests = violations = control = 0
    windows_ok = True
    for Q in range(3, 7):
        cube = list(product((0, 1), repeat=Q))
        for D in range(0, Q):
            nd = N_D(Q, D)
            lo = math.ceil((Q - D) / 2)
            win = [x for x in cube if lo <= sum(x) <= lo + D]
            windows_ok &= interpolable(win, Q, D) and len(win) <= nd
            if nd + 1 > 2 ** Q:
                continue
            for _ in range(10):
                S = rng.sample(cube, nd + 1)
                tests += 1
                violations += int(interpolable(S, Q, D))
                if D + 2 <= Q:
                    control += int(interpolable(S, Q, D + 2))
    ok = violations == 0 and windows_ok and control > 0
    RESULTS["A_smolensky"] = {"random_sets_of_N_D_plus_1": tests, "interpolable_(must_be_0)": violations,
                              "central_windows_interpolable": windows_ok,
                              "control_interpolable_in_degree_D_plus_2": control, "ok": ok}
    print(f"[A] Smolensky over Q: {tests} random sets of N_D+1 points, {violations} interpolable (must be 0); "
          f"central windows interpolable: {windows_ok}; control (degree D+2) interpolable: {control}")
    return ok


# =============================================================== Part B: the p-biased bound by restriction

def tau(D, m):
    """P[Bin(m,1/2) > (m+D)/2] exactly."""
    k = (m + D) // 2
    return Fraction(sum(math.comb(m, i) for i in range(k + 1, m + 1)), 2 ** m)


def mu(x, p):
    w = sum(x)
    return p ** w * (1 - p) ** (len(x) - w)


def restriction_bound(Q, D, p):
    q = 2 * p
    return sum(math.comb(Q, m) * q ** m * (1 - q) ** (Q - m) * tau(D, m) for m in range(Q + 1))


def design_from(order, Q, D):
    """Interpolate parity in degree <= D on a maximal subset of `order`, taken greedily in that order."""
    ms = monomials(Q, D)
    el = Eliminator(len(ms))
    for x in order:
        el.add(interp_row(ms, x))
    return ms, el.solution()


def evaluate(poly, x):
    ms, coef = poly
    return sum(c * chi(A, x) for A, c in zip(ms, coef) if c)


def bad_measure(poly, cube, p):
    return sum((mu(x, p) for x in cube if evaluate(poly, x) != parity(x)), Fraction(0))


def coupling_side(poly, Q, p):
    """E_Z[fraction of y in {0,1}^Z with F(Z and y) != parity(y)], Z ~ Bern(2p)^Q, computed exactly."""
    tot = Fraction(0)
    for z in product((0, 1), repeat=Q):
        pz = mu(z, 2 * p)
        if not pz:
            continue
        A = [i for i in range(Q) if z[i]]
        bad = 0
        for y in product((0, 1), repeat=len(A)):
            xx = [0] * Q
            for i, a in enumerate(A):
                xx[a] = y[i]
            bad += evaluate(poly, xx) != parity(xx)
        tot += pz * Fraction(bad, 2 ** len(A))
    return tot


def part_B():
    rng = random.Random(1402)
    rows, coupling_ok = [], True
    control_violations = tight = 0
    best_over_window = []
    for Q in (4, 5, 6):
        cube = list(product((0, 1), repeat=Q))
        for D in range(0, Q - 1):
            greedy_random = [design_from(rng.sample(cube, len(cube)), Q, D) for _ in range(2)]
            for p in (Fraction(1, 4), Fraction(1, 3), Fraction(3, 8), Fraction(1, 2)):
                bound = restriction_bound(Q, D, p)
                levels = sorted(range(Q + 1), key=lambda j: (abs(j - Q * p), j))
                lo, hi = min(levels[:D + 1]), max(levels[:D + 1])
                window = design_from([x for x in cube if lo <= sum(x) <= hi], Q, D)
                heavy = design_from(sorted(cube, key=lambda x: (-mu(x, p), rng.random())), Q, D)
                designs = [("window", window), ("greedy_random_0", greedy_random[0]),
                           ("greedy_random_1", greedy_random[1]), ("greedy_mu_heavy", heavy)]
                measures = {}
                for name, poly in designs:
                    mE = bad_measure(poly, cube, p)
                    measures[name] = mE
                    tight += mE == bound
                    rows.append({"Q": Q, "D": D, "p": str(p), "design": name, "mu_p(E)": float(mE),
                                 "bound": float(bound), "holds": mE >= bound})
                    if Q <= 5 and name == "greedy_random_0":
                        coupling_ok &= coupling_side(poly, Q, p) == mE
                best_over_window.append(float(min(measures.values()) / measures["window"]))
                # control: a window design of degree D + 2 can beat the degree-D bound
                lo2, hi2 = min(levels[:D + 3]), max(levels[:D + 3])
                window2 = design_from([x for x in cube if lo2 <= sum(x) <= hi2], Q, D + 2)
                control_violations += int(bad_measure(window2, cube, p) < bound)
    holds = all(r["holds"] for r in rows)
    ratios = [r["mu_p(E)"] / r["bound"] for r in rows if r["bound"] > 0]
    ok = holds and coupling_ok and control_violations > 0
    RESULTS["B_p_biased"] = {"cases": len(rows), "all_hold": holds, "coupling_identity_exact": coupling_ok,
                             "min_ratio_mu_over_bound": min(ratios), "cases_with_equality": tight,
                             "best_found_over_window_min": min(best_over_window),
                             "control_degree_D_plus_2_violates_degree_D_bound": control_violations,
                             "sample": rows[:16], "ok": ok}
    print(f"[B] restriction bound: {len(rows)} (Q, D, p, design) cases hold: {holds}; coupling identity exact: "
          f"{coupling_ok}; min mu_p(E)/bound = {min(ratios):.3f}; best design found / window >= "
          f"{min(best_over_window):.3f}; control violations: {control_violations}")
    return ok


# =============================================================== Part C: explicit forms

def part_C():
    worst = math.inf
    for m in range(2, 301):
        suffix = [0] * (m + 2)
        for i in range(m, -1, -1):
            suffix[i] = suffix[i + 1] + math.comb(m, i)
        for D in range(0, m - 1):
            log_tau = math.log(suffix[(m + D) // 2 + 1]) - m * math.log(2)
            worst = min(worst, log_tau - (-(D + 2) ** 2 / m - math.log(m + 1)))
    tail_margin = math.inf
    for Q in range(32, 2001, 7):
        k = -(-Q // 32)                      # weights below Q/32 are 0, ..., ceil(Q/32) - 1
        for m in range(-(-Q // 2), Q + 1, max(1, Q // 16)):
            low = sum(math.comb(m, i) for i in range(k))
            tail_margin = min(tail_margin, -0.229 * Q - (math.log(low) - m * math.log(2)))
    ok = worst > 0 and tail_margin > 0
    RESULTS["C_explicit"] = {"min_log_tau_minus_log_lower_bound": worst, "min_tail_log_margin": tail_margin,
                             "ok": ok}
    print(f"[C] tau_D(m) >= e^(-(D+2)^2/m)/(m+1) for m <= 300: min log-margin {worst:.3f}; "
          f"low-weight tail <= e^(-0.229Q) for Q <= 2000: min log-margin {tail_margin:.3f}")
    return ok


# =============================================================== Part D: the floor for all designs

def eps_all(Q, R):
    """Lemma 2.4: mu_p(E and weight >= Q/32) >= this, for every design of degree < R and p in [1/4, 1/2]."""
    if Q < 4 or R > Q / 4:
        return None
    val = 0.5 * (math.exp(-2 * (R + 1) ** 2 / Q) / (Q + 1) - math.exp(-0.229 * Q))
    return val if val > 0 else None


def taylor(R, Q):
    x = R / Q
    if x >= 1 / math.e:
        return 1 / math.e * math.log2(math.e)
    return x * math.log2(1 / x)


def best(lam):
    bestv = (float("inf"), None)
    r = 2
    while r <= 2 ** 200:
        Qmax = int(lam / (2 * math.log2(r)))
        if Qmax < 2:
            break
        Qs = sorted({max(2, int(Qmax * f)) for f in (1, 0.8, 0.6, 0.4, 0.25, 0.1, 0.05)})
        for Q in Qs:
            R = 1
            while R <= Q:
                e = eps_all(Q, R)
                meas = Q * e / 64 if e is not None else 0.0
                d = max(taylor(R, Q), 1 / r, meas)
                if d < bestv[0]:
                    bestv = (d, (r, Q, R))
                R = max(R + 1, int(R * 1.15))
        r = max(r + 1, int(r * 1.3))
    return bestv


def part_D():
    rows = []
    for e in (10, 20, 30, 40, 60, 80, 100, 150, 200, 300):
        lam = 2.0 ** e
        d, arg = best(lam)
        ref = math.log2(lam) ** 2 / math.sqrt(lam)
        rows.append({"log2_n": f"2^{e}", "min_delta": d, "argmin_(r,Q,R)": arg, "reference": ref,
                     "ratio": d / ref})
    ratios = [r["ratio"] for r in rows if r["min_delta"] < 1]
    ok = bool(ratios) and min(ratios) > 0.01
    RESULTS["D_floor_all_designs"] = {"rows": rows, "ratio_range": [min(ratios), max(ratios)] if ratios else None,
                                      "ok": ok}
    print(f"[D] floor for all single-polynomial designs: ratio to (log log n)^2/sqrt(log n) in "
          f"[{min(ratios):.3f}, {max(ratios):.3f}] for log2 n = 2^10..2^300")
    return ok


# =============================================================== Part E: node-dependent designs on windows

def node_window_feasible(Q, R, W, node_dependent):
    """Can the node conditions  sum_{S containing t} psi_S(t) chi_S(X) = parity(X)  hold for all nodes t = 1..Q
    and all X with |X| or |X + e_t| in W?

    node_dependent: psi_S(t) is free for |S| < R and t in S; psi_S(T) = sum_{1 <= j < R} beta_{S,j} T^j for |S| >= R.
    single polynomial: psi_S(t) = g_S for t in S when |S| < R, and psi_S = 0 on S when |S| >= R."""
    subsets = [S for k in range(1, Q + 1) for S in combinations(range(Q), k)]
    idx = {}
    for S in subsets:
        if len(S) < R:
            for t in (S if node_dependent else (None,)):
                idx[(S, t)] = len(idx)
        elif node_dependent:
            for j in range(1, R):
                idx[(S, -j)] = len(idx)
    el = Eliminator(len(idx))
    for t in range(Q):
        for x in product((0, 1), repeat=Q):
            # both sides are odd in x_t, so X and X + e_t give the same condition: use x_t = 0
            if x[t] or not (sum(x) in W or sum(x) + 1 in W):
                continue
            row = [0] * (len(idx) + 1)
            for S in subsets:
                if t not in S:
                    continue
                c = chi(S, x)
                if len(S) < R:
                    row[idx[(S, t if node_dependent else None)]] += c
                elif node_dependent:
                    for j in range(1, R):
                        row[idx[(S, -j)]] += c * (t + 1) ** j
            row[-1] = parity(x)
            if not el.add(row):
                return False
    return True


def part_E(full=False):
    cases = [(6, 3, [{3}, {2, 3}, {2, 3, 4}]), (7, 3, [{3}, {4}, {3, 4}])]
    if full:
        cases.append((8, 4, [{3, 4}, {3, 4, 5}]))          # about 3 minutes
    out, ok = [], True
    for Q, R, windows in cases:
        rows = [{"Q": Q, "R": R, "window": sorted(W),
                 "node_dependent_feasible": node_window_feasible(Q, R, W, True),
                 "single_polynomial_feasible": node_window_feasible(Q, R, W, False)} for W in windows]
        ok &= all(r["node_dependent_feasible"] == r["single_polynomial_feasible"] for r in rows)
        ok &= any(r["node_dependent_feasible"] for r in rows) and not all(r["node_dependent_feasible"] for r in rows)
        out += rows
    RESULTS["E_node_dependent_windows"] = {"rows": out, "same_feasible_windows": ok, "ok": ok}
    print("[E] central windows, feasible for node-dependent / single-polynomial designs: " +
          "; ".join(f"Q={r['Q']} R={r['R']} {r['window']}: {r['node_dependent_feasible']}/"
                    f"{r['single_polynomial_feasible']}" for r in out))
    return ok


def main():
    import sys
    full = "--full" in sys.argv
    oks = []
    for name, fn in (("A", part_A), ("B", part_B), ("C", part_C), ("D", part_D),
                     ("E", lambda: part_E(full))):
        t0 = time.time()
        oks.append(fn())
        RESULTS.setdefault("seconds", {})[name] = round(time.time() - t0, 1)
    RESULTS["all_ok"] = all(oks)
    with open(__file__.rsplit("/", 1)[0] + "/../results/design_floor.json", "w") as fh:
        json.dump(RESULTS, fh, indent=1, default=str)
    print("ALL OK" if all(oks) else "FAILURES")


if __name__ == "__main__":
    main()

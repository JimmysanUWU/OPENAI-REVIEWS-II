"""End-to-end run of the bounded-occurrence counting algorithm (Sections 2-4) on small CNFs.

Pipeline, following the paper step by step:

  1. auxiliary rows P_u(T) and C^B (Lemmas 4.1-4.3), built by eth_core.Aux;
  2. schedules a_v, d_f (random here; Lemma 4.5 only affects running time, not correctness);
  3. copies, iota, rows P_{f,v}(T) and matrices C_{t,g} (Section 4.5, (4.14)-(4.15));
  4. exact check of Definition 3.1(4) for every t, f, v and g != sigma_f(t);
  5. for every t in [Q] and every pair U subset E, |E| <= q, the restricted sum
     Z_{t,E,U}(z) evaluated by the measurement-separation procedure of Lemma 3.2;
  6. the Taylor polynomial H(t+z) mod z^J by the signed sum (3.13) (Lemma 3.3);
  7. Hermite interpolation of H from the QJ Taylor coefficients, returning H(0).

The answer H(0) is compared with brute-force #F.  Independently, H(T) is expanded
symbolically by brute force to confirm deg H, H(0) = #F and the Taylor data of step 6.
Arithmetic is exact over Q, or over F_p (p = 2^61-1) for the larger runs.
"""

import json
import random
import sys
import time
from fractions import Fraction
from itertools import combinations, product
from math import comb

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from eth_core import Aux, rref, in_rowspace, poly_eval, poly_mul, random_bounded_cnf, brute_count  # noqa: E402

P61 = (1 << 61) - 1


class F:
    """Either exact rationals (p=None) or integers modulo a prime p."""

    def __init__(self, p=None):
        self.p = p

    def c(self, x):
        x = Fraction(x)
        if self.p is None:
            return x
        return x.numerator % self.p * pow(x.denominator, -1, self.p) % self.p

    def red(self, x):
        return x if self.p is None else x % self.p

    def pmul(self, a, b, J):
        out = [0] * J
        for i, x in enumerate(a):
            if x:
                for j in range(J - i):
                    y = b[j]
                    if y:
                        out[i + j] += x * y
        return [self.red(v) for v in out]

    def solve(self, A, b):
        """Solve the square system A h = b (Gaussian elimination)."""
        n = len(A)
        M = [list(row) + [bb] for row, bb in zip(A, b)]
        for col in range(n):
            piv = next(i for i in range(col, n) if M[i][col] != 0)
            M[col], M[piv] = M[piv], M[col]
            inv = (1 / M[col][col]) if self.p is None else pow(M[col][col], -1, self.p)
            M[col] = [self.red(v * inv) for v in M[col]]
            for i in range(n):
                if i != col and M[i][col] != 0:
                    f = M[i][col]
                    M[i] = [self.red(a - f * b2) for a, b2 in zip(M[i], M[col])]
        return [M[i][n] for i in range(n)]


def build_instance(n, clauses, r, Q, R, rng, aux=None):
    aux = aux or Aux(r, Q, R)
    used = {v for cl in clauses for v, _ in cl}
    factors = [{"scope": [v for v, _ in cl], "signs": [s for _, s in cl], "dummy": False} for cl in clauses]
    for v in range(n):
        if v not in used:
            factors.append({"scope": [v], "signs": [True], "dummy": True})
    a = [tuple(rng.randrange(r) for _ in range(Q)) for _ in range(n)]
    dsch = [tuple(rng.randrange(r) for _ in range(Q)) for _ in factors]
    # labels and copies
    count = {}
    occ = {}
    place = {}
    for fi, f in enumerate(factors):
        for v in f["scope"]:
            u = tuple((a[v][s], dsch[fi][s]) for s in range(Q))
            c = count.get(u, 0)
            count[u] = c + 1
            occ[(c, u)] = v
            place[(fi, v)] = (c, u)
    h = max(count.values())
    # rows P_{f,v}(T) in Q[T]^{1 x n}: P_u(T) pi_c iota
    Pfv = {}
    for (fi, v), (c, u) in place.items():
        row = {}
        for y, poly in aux.P[u].items():
            vv = occ.get((c, y))
            if vv is None:
                continue
            cur = row.get(vv, [Fraction(0)])
            m = max(len(cur), len(poly))
            cur = cur + [Fraction(0)] * (m - len(cur))
            for i, x in enumerate(poly):
                cur[i] += x
            row[vv] = cur
        Pfv[(fi, v)] = {k: p for k, p in row.items() if any(p)}
    # measurement matrices C_{t,g}: rows indexed by (copy c, C^B row i)
    C = {}
    for s in range(Q):
        for g in range(r):
            rows = []
            for c in range(h):
                for crow in aux.CB:
                    row = {}
                    for y, val in crow.items():
                        if y[s][0] != g:
                            continue
                        vv = occ.get((c, y))
                        if vv is not None:
                            row[vv] = row.get(vv, 0) + val
                    rows.append({k: x for k, x in row.items() if x != 0})
            C[(s, g)] = rows
    return {"aux": aux, "factors": factors, "a": a, "d": dsch, "h": h, "P": Pfv, "C": C,
            "n": n, "r": r, "Q": Q, "R": R, "label_multiplicity": max(count.values())}


def check_design(inst):
    """Definition 3.1 (1)-(4) in exact arithmetic (group sizes reported, not bounded)."""
    aux, r, Q, n = inst["aux"], inst["r"], inst["Q"], inst["n"]
    ok_p0 = all({k: p[0] for k, p in row.items() if p[0] != 0} == {v: 1} for (fi, v), row in inst["P"].items())
    ok_deg = all(len(p) - 1 <= inst["R"] for row in inst["P"].values() for p in row.values())
    ok_kernel = True
    max_abs = 0
    total_rows = []
    for s in range(Q):
        t = s + 1
        groups = {g: {v for v in range(n) if inst["a"][v][s] == g} for g in range(r)}
        bases = {g: rref(inst["C"][(s, g)], n) for g in range(r)}
        total_rows.append(sum(len(inst["C"][(s, g)]) for g in range(r)))
        for g in range(r):
            for row in inst["C"][(s, g)]:
                max_abs = max(max_abs, sum(abs(x) for x in row.values()))
                assert set(row) <= groups[g]
        for (fi, v), row in inst["P"].items():
            sig = inst["d"][fi][s]
            val = {k: poly_eval(p, t) for k, p in row.items()}
            for g in range(r):
                if g == sig:
                    continue
                restr = {k: x for k, x in val.items() if k in groups[g] and x != 0}
                if not in_rowspace(bases[g], restr):
                    ok_kernel = False
    return {"P0_delta": ok_p0, "deg_le_R": ok_deg, "kernel_condition": ok_kernel,
            "max_meas_row_abs_sum": max_abs, "meas_abs_sum_le_D": max_abs <= aux.D,
            "rows_per_point": total_rows, "rows_equal_r_h_cB": all(x == r * inst["h"] * aux.cB for x in total_rows)}


def phi_value(fld, f, inputs, J):
    """Clause polynomial 1 - prod(1 - L) on truncated power series inputs; dummy -> 1."""
    one = [fld.c(1)] + [0] * (J - 1)
    if f["dummy"]:
        return one
    prod_ = one
    for inp, pos in zip(inputs, f["signs"]):
        L = inp if pos else [fld.red((1 if i == 0 else 0) - x) for i, x in enumerate(inp)]
        prod_ = fld.pmul(prod_, [fld.red((1 if i == 0 else 0) - x) for i, x in enumerate(L)], J)
    return [fld.red((1 if i == 0 else 0) - x) for i, x in enumerate(prod_)]


def run_algorithm(inst, fld, tight=False):
    """Proposition 3.4 with Lemma 3.2 separation and Lemma 3.3 Taylor recovery."""
    n, r, Q, R = inst["n"], inst["r"], inst["Q"], inst["R"]
    factors = inst["factors"]
    d = sum(len(f["scope"]) for f in factors if not f["dummy"])
    dH = d * ((R - 1) if tight else R)
    q = dH // Q
    J = q + 1
    stats = {"d": d, "d_H": dH, "q": q, "J": J, "pairs_per_point": sum(2 ** s * comb(n, s) for s in range(min(q, n) + 1))}
    taylor = {}
    max_bucket = 0
    tuples_seen = 0
    for s in range(Q):
        t = s + 1
        # P(t) rows and shift polynomials Delta P(z) = P(t+z) - P(t) (z^1..z^{J-1})
        rowt, shift = {}, {}
        for key, row in inst["P"].items():
            rowt[key] = {k: fld.c(poly_eval(p, t)) for k, p in row.items()}
            sh = {}
            for k, p in row.items():
                co = [Fraction(0)] * J
                for deg, cf in enumerate(p):
                    if cf:
                        for j in range(1, min(deg, J - 1) + 1):
                            co[j] += cf * comb(deg, j) * Fraction(t) ** (deg - j)
                sh[k] = [fld.c(x) for x in co]
            shift[key] = sh
        groups = [[v for v in range(n) if inst["a"][v][s] == g] for g in range(r)]
        sigma = [inst["d"][fi][s] for fi in range(len(factors))]
        fac_by_group = [[fi for fi in range(len(factors)) if sigma[fi] == g] for g in range(r)]
        # buckets of group assignments by integer measurement vector (Lemma 3.2)
        buckets = []
        for g in range(r):
            bk = {}
            for bits in product((0, 1), repeat=len(groups[g])):
                xg = dict(zip(groups[g], bits))
                key = tuple(sum(cf * xg[v] for v, cf in row.items()) for row in inst["C"][(s, g)])
                bk.setdefault(key, []).append(xg)
            buckets.append(bk)
            max_bucket = max(max_bucket, max(len(b) for b in bk.values()))
        keys = [list(b.keys()) for b in buckets]
        total = [0] * J
        for size in range(min(q, n) + 1):
            for E in combinations(range(n), size):
                for usize in range(size + 1):
                    for U in combinations(E, usize):
                        sign = -1 if (size - usize) % 2 else 1
                        Uset = set(U)
                        # fixed corrections <Delta P_{f,v}(z), 1_U>
                        corr = {}
                        for key, sh in shift.items():
                            co = [0] * J
                            for k, poly in sh.items():
                                if k in Uset:
                                    co = [a_ + b_ for a_, b_ in zip(co, poly)]
                            corr[key] = [fld.red(x) for x in co]
                        Eset = set(E)
                        Z = [0] * J
                        for ytuple in product(*keys):
                            tuples_seen += 1
                            reps = [buckets[g][ytuple[g]][0] for g in range(r)]
                            prod_total = [fld.c(1)] + [0] * (J - 1)
                            for g in range(r):
                                sg = [0] * J
                                for xg in buckets[g][ytuple[g]]:
                                    if any(xg[v] != 1 for v in Eset if v in xg):
                                        continue
                                    x = {}
                                    for hh in range(r):
                                        x.update(xg if hh == g else reps[hh])
                                    val = [fld.c(1)] + [0] * (J - 1)
                                    for fi in fac_by_group[g]:
                                        f = factors[fi]
                                        inputs = []
                                        for v in f["scope"]:
                                            base = sum(cf * x[k] for k, cf in rowt[(fi, v)].items())
                                            inp = list(corr[(fi, v)])
                                            inp[0] = fld.red(inp[0] + base)
                                            inputs.append(inp)
                                        val = fld.pmul(val, phi_value(fld, f, inputs, J), J)
                                    sg = [fld.red(a_ + b_) for a_, b_ in zip(sg, val)]
                                prod_total = fld.pmul(prod_total, sg, J)
                                if not any(prod_total):
                                    break
                            Z = [fld.red(a_ + b_) for a_, b_ in zip(Z, prod_total)]
                        total = [fld.red(a_ + sign * b_) for a_, b_ in zip(total, Z)]
        taylor[t] = total
    # Hermite interpolation: unknown coefficients h_0..h_{QJ-1}
    N = Q * J
    A, b = [], []
    for t in range(1, Q + 1):
        for j in range(J):
            A.append([fld.c(comb(u, j) * t ** (u - j)) if u >= j else 0 for u in range(N)])
            b.append(taylor[t][j])
    hcoef = fld.solve(A, b)
    stats["max_bucket_size"] = max_bucket
    stats["measurement_tuples_visited"] = tuples_seen
    return hcoef[0], taylor, stats


def symbolic_H(inst):
    """Brute-force H(T) = sum_x prod_f Phi_f(<P_{f,v}(T), x>) as an exact polynomial in T."""
    n, factors = inst["n"], inst["factors"]
    H = [Fraction(0)]
    for x in product((0, 1), repeat=n):
        val = [Fraction(1)]
        for fi, f in enumerate(factors):
            if f["dummy"]:
                continue
            prod_ = [Fraction(1)]
            for v, pos in zip(f["scope"], f["signs"]):
                inp = [Fraction(0)]
                for k, p in inst["P"][(fi, v)].items():
                    if x[k]:
                        m = max(len(inp), len(p))
                        inp = [(inp[i] if i < len(inp) else 0) + (p[i] if i < len(p) else 0) for i in range(m)]
                L = inp if pos else [(1 if i == 0 else 0) - c for i, c in enumerate(inp)]
                prod_ = poly_mul(prod_, [(1 if i == 0 else 0) - c for i, c in enumerate(L)])
            phi = [(1 if i == 0 else 0) - c for i, c in enumerate(prod_)]
            val = poly_mul(val, phi)
        m = max(len(H), len(val))
        H = [(H[i] if i < len(H) else 0) + (val[i] if i < len(val) else 0) for i in range(m)]
    while len(H) > 1 and H[-1] == 0:
        H.pop()
    return H


def one_run(n, r, Q, R, seed, fld, tight, n_clauses=None):
    rng = random.Random(seed)
    clauses = random_bounded_cnf(n, rng, n_clauses=n_clauses)
    inst = build_instance(n, clauses, r, Q, R, rng)
    design = check_design(inst)
    t0 = time.time()
    h0, taylor, stats = run_algorithm(inst, fld, tight=tight)
    elapsed = time.time() - t0
    truth = brute_count(n, clauses)
    H = symbolic_H(inst)
    d = stats["d"]
    # Taylor data agreement with the symbolic polynomial
    agree = True
    for t, data in taylor.items():
        exact = []
        for j in range(stats["J"]):
            exact.append(fld.c(sum(c * comb(k, j) * Fraction(t) ** (k - j) for k, c in enumerate(H) if k >= j)))
        agree &= [fld.red(x) for x in data] == exact
    res = {"n": n, "r": r, "Q": Q, "R": R, "seed": seed, "field": "Q" if fld.p is None else "F_(2^61-1)",
           "tight_degree": tight, "clauses": len(clauses), "brute_count": truth,
           "algorithm_H0": str(h0), "match": fld.c(truth) == h0,
           "h_copies": inst["h"], "cB": inst["aux"].cB, "D": inst["aux"].D,
           "symbolic_deg_H": len(H) - 1, "deg_bound_dR": d * R, "deg_le_d(R-1)": len(H) - 1 <= d * (R - 1),
           "symbolic_H0_equals_count": H[0] == truth, "taylor_data_matches_symbolic": agree,
           "design": design, "stats": stats, "seconds": round(elapsed, 2)}
    return res


def main():
    out = []
    regimes = [
        # (n, r, Q, R, tight, field, seeds): see module docstring
        (5, 2, 3, 3, False, None, [1, 2]),        # c_B = 1: nontrivial measurement buckets, exact Q
        (6, 2, 3, 3, False, P61, [3, 4, 5]),      # c_B = 1, F_p
        (6, 2, 2, 2, False, P61, [6, 7]),         # D = 16: many copies (h >= 2)
        (6, 2, 4, 2, True, P61, [8, 9, 11, 12]),  # q < n: truncated inclusion-exclusion
        (6, 2, 4, 2, False, P61, [13]),           # paper's d_H = dR
        (7, 2, 4, 2, True, P61, [10, 14]),        # q < n, larger
        (8, 2, 4, 2, True, P61, [15]),            # q < n, n = 8
        (6, 3, 2, 2, False, P61, [16, 17]),       # r = 3 groups
    ]
    allok = True
    for n, r, Q, R, tight, p, seeds in regimes:
        for seed in seeds:
            res = one_run(n, r, Q, R, seed, F(p), tight)
            ok = (res["match"] and res["symbolic_H0_equals_count"] and res["taylor_data_matches_symbolic"]
                  and res["deg_le_d(R-1)"] and all(v for k, v in res["design"].items() if isinstance(v, bool)))
            res["all_ok"] = ok
            allok &= ok
            out.append(res)
            print(json.dumps({k: res[k] for k in ("n", "r", "Q", "R", "seed", "field", "brute_count", "algorithm_H0",
                                                   "match", "h_copies", "cB", "all_ok", "seconds")}),
                  json.dumps(res["stats"]), flush=True)
    # negative control: delete one measurement row from every C_{t,g}; the kernel condition of
    # Definition 3.1(4) must then fail, and the separated sum is no longer the true count.
    neg = []
    for seed in (2, 9, 16):
        n, r, Q, R = (5, 2, 3, 3) if seed == 2 else ((6, 2, 4, 2) if seed == 9 else (6, 3, 2, 2))
        rng = random.Random(seed)
        clauses = random_bounded_cnf(n, rng)
        inst = build_instance(n, clauses, r, Q, R, rng)
        for key in inst["C"]:
            nonzero = [i for i, row in enumerate(inst["C"][key]) if row]
            if nonzero:
                inst["C"][key] = [row for i, row in enumerate(inst["C"][key]) if i != nonzero[0]]
        design = check_design(inst)
        h0, _, _ = run_algorithm(inst, F(P61), tight=(Q == 4))
        truth = brute_count(n, clauses)
        neg.append({"seed": seed, "kernel_condition_after_deletion": design["kernel_condition"],
                    "algorithm_H0": str(h0), "brute_count": truth, "still_matches": F(P61).c(truth) == h0})
        print("negative control", neg[-1], flush=True)
    neg_ok = all(not x["kernel_condition_after_deletion"] for x in neg) and any(not x["still_matches"] for x in neg)
    print("ALL PASS" if allok and neg_ok else "SOME CHECK FAILED")
    path = __file__.rsplit("/", 2)[0] + "/results/counting_pipeline.json"
    with open(path, "w") as fh:
        json.dump({"runs": out, "negative_controls": neg, "all_ok": allok, "negative_controls_detect": neg_ok}, fh, indent=1)


if __name__ == "__main__":
    main()

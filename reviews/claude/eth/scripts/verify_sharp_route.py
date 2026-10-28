"""Checks for 09_GAPLESS_SHARP_EXPONENT.md: the paper's Sections 5-6 route (exponent rho(n)).

Labels in brackets are the doc's numbering; the paper's follows.

Part L  [Lemma 2.1] The h-adic digits of the Beigel-Tarui polynomial in closed form:
        e_0 = 1 - Y, e_j = (-1)^(j+1) C(2j-1, j-1) (2Y - 1).  For every ell <= 128 they
        equal the digits from repeated division by h = Y^2 - Y; the congruence
        p_ell(y) = [y even] mod 2^ell holds on |y| <= 40; ||e_j||_1 <= 2 * 4^j.
Part Z  [Theorem 3.1] Counting common zeros of a sparse F_2 system (paper Theorem 2.2):
        one factor B_i = sum_j e_j(L_i) h(L_i)^j V^j per equation, product truncated in
        A_ell = (Z/p^mu)[V]/(V^ell), fed to the deformation-counting algorithm of
        08_GAPLESS_PROOF.md (closed-form swap rows, derandomized schedules, separated
        Taylor jets, CRT at zero); the V-coefficients of H(0) are recovered exactly and
        S_ell(H(0)) mod 2^ell is compared with brute force.  Negative control: the
        untruncated degree would need far more Taylor coefficients.
Part T  [Lemma 4.1] Paper Lemma 6.1's mechanism on toy parameters (epsilon-biased shifts,
        the mixing bound, a local code, exhaustive minimum weight of the Tanner code), plus
        the exact arithmetic of the lemma's constants on a grid of (B, m).
Part G  [Lemma 4.2] Paper Lemma 6.2: greedy rows chosen by exact counts on small CNFs; each
        round keeps at most 3/4 of the surviving nonsatisfying assignments; after 3n rounds
        N(A) = #F.

Writes ../results/sharp_route.json.
"""

import json
import math
import random
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import comb

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from eth_core import poly_eval  # noqa: E402
from verify_gapless_chain import Swap, schedules, prime_above, shift_poly_mod, fmod  # noqa: E402
from verify_parity_lift import p_ell, divide_by_h  # noqa: E402

RESULTS = {}


# =============================================================== Part L

def closed_digits(ell):
    """(a_j, b_j) with e_j = a_j + b_j Y, from the closed form of Lemma 2.1."""
    out = [(1, -1)]
    for j in range(1, ell):
        c = (-1) ** (j + 1) * comb(2 * j - 1, j - 1)
        out.append((-c, 2 * c))
    return out


def part_L():
    ok, rows = True, []
    for ell in range(1, 129):
        div = [tuple(e) for e in divide_by_h(p_ell(ell))]
        div += [(0, 0)] * (ell - len(div))
        cl = closed_digits(ell)
        same = div == cl
        norm_ok = all(abs(a) + abs(b) <= 2 * 4 ** j for j, (a, b) in enumerate(cl))
        cong_ok = True
        if ell <= 40:
            for y in range(-40, 41):
                hy = y * y - y
                val = sum((a + b * y) * hy ** j for j, (a, b) in enumerate(cl))
                cong_ok &= (val - (1 if y % 2 == 0 else 0)) % (2 ** ell) == 0
        ok &= same and norm_ok and cong_ok
        if ell in (1, 2, 3, 8, 32, 128):
            rows.append({"ell": ell, "closed_equals_division": same, "norm_bound": norm_ok,
                         "congruence_on_|y|<=40": cong_ok if ell <= 40 else "not run",
                         "log2_max_digit_norm": round(math.log2(max(abs(a) + abs(b) for a, b in cl)), 3)})
    # negative control: the unhalved central binomial C(2j, j) in place of C(2j-1, j-1)
    bad = [(1, -1)] + [(-(-1) ** (j + 1) * comb(2 * j, j), 2 * (-1) ** (j + 1) * comb(2 * j, j))
                       for j in range(1, 6)]
    neg = any((sum((a + b * y) * (y * y - y) ** j for j, (a, b) in enumerate(bad))
               - (1 if y % 2 == 0 else 0)) % 2 ** 6 for y in range(-10, 11))
    ok &= neg
    RESULTS["parity_lift_closed_form"] = {"all_ok": ok, "ell_range": [1, 128], "samples": rows,
                                          "negative_control_detects": neg}
    print(f"[L] closed-form h-adic digits = division digits for ell <= 128; congruence; norms <= 2*4^j: {ok}")
    return ok


# =============================================================== Part Z

class Ring2:
    """(Z/m)[z, V]/(z^J, V^L): elements are J x L lists."""

    def __init__(self, m, J, L):
        self.m, self.J, self.L = m, J, L

    def zero(self):
        return [[0] * self.L for _ in range(self.J)]

    def const(self, c):
        x = self.zero()
        x[0][0] = c % self.m
        return x

    def from_z(self, coeffs):
        x = self.zero()
        for i, c in enumerate(coeffs[:self.J]):
            x[i][0] = c % self.m
        return x

    def add(self, a, b):
        return [[(a[i][j] + b[i][j]) % self.m for j in range(self.L)] for i in range(self.J)]

    def scale(self, a, c):
        return [[(a[i][j] * c) % self.m for j in range(self.L)] for i in range(self.J)]

    def mul(self, a, b):
        J, L, m = self.J, self.L, self.m
        out = [[0] * L for _ in range(J)]
        for i1 in range(J):
            for j1 in range(L):
                x = a[i1][j1]
                if not x:
                    continue
                for i2 in range(J - i1):
                    row = b[i2]
                    o = out[i1 + i2]
                    for j2 in range(L - j1):
                        y = row[j2]
                        if y:
                            o[j1 + j2] = (o[j1 + j2] + x * y) % m
        return out

    def shiftV(self, a, k):
        out = self.zero()
        for i in range(self.J):
            for j in range(self.L - k):
                out[i][j + k] = a[i][j]
        return out


def zmul(a, b, J, m):
    out = [0] * J
    for i, x in enumerate(a):
        if x:
            for j in range(J - i):
                if b[j]:
                    out[i + j] = (out[i + j] + x * b[j]) % m
    return out


def eq_factor(ring, monos, args, es):
    """B = sum_j e_j(L) h(L)^j V^j with L = sum_{M in monos} prod_{v in M} Y_v.

    The arguments are z-polynomials (no V), so every coefficient of V^j is the z-polynomial
    e_j(L) h(L)^j; it is placed in column j of the (z, V) array."""
    J, m = ring.J, ring.m
    L = [0] * J
    for M in monos:
        t = [1] + [0] * (J - 1)
        for v in M:
            t = zmul(t, args[v], J, m)
        L = [(p + q) % m for p, q in zip(L, t)]
    Lm1 = list(L)
    Lm1[0] = (Lm1[0] - 1) % m
    h = zmul(L, Lm1, J, m)
    out = ring.zero()
    hp = [1] + [0] * (J - 1)
    for j, (a0, a1) in enumerate(es):
        if j >= ring.L:
            break
        ej = [(a1 * c) % m for c in L]
        ej[0] = (ej[0] + a0) % m
        col = zmul(ej, hp, J, m)
        for i in range(J):
            out[i][j] = col[i]
        hp = zmul(hp, h, J, m)
    return out


def build_general(n, scopes, sw):
    r, Q = sw.r, sw.Q
    a, d, E0, E1, M, b, w = schedules(n, scopes, r, Q)
    count, occ, place = Counter(), {}, {}
    for fi, f in enumerate(scopes):
        for v in f:
            u = tuple((a[v][s], d[fi][s]) for s in range(Q))
            c = count[u]
            count[u] += 1
            occ[(c, u)] = v
            place[(fi, v)] = (c, u)
    h = max(count.values())
    Pfv = {}
    for (fi, v), (c, u) in place.items():
        row = {}
        for y, poly in sw.P[u].items():
            vv = occ.get((c, y))
            if vv is None:
                continue
            cur = row.get(vv, [Fraction(0)])
            Ln = max(len(cur), len(poly))
            row[vv] = [(cur[i] if i < len(cur) else 0) + (poly[i] if i < len(poly) else 0) for i in range(Ln)]
        Pfv[(fi, v)] = {k: pp for k, pp in row.items() if any(pp)}
    Cm = {}
    for s in range(Q):
        for g in range(r):
            rows = []
            for c in range(h):
                for psi in sw.bad:
                    row = {}
                    for y, val in psi.items():
                        if y[s][0] != g:
                            continue
                        vv = occ.get((c, y))
                        if vv is not None:
                            row[vv] = row.get(vv, 0) + val
                    row = {k: x for k, x in row.items() if x != 0}
                    if row:
                        rows.append(row)
            Cm[(s, g)] = rows
    groups = {(s, g): [v for v in range(n) if a[v][s] == g] for s in range(Q) for g in range(r)}
    return {"scopes": scopes, "d": d, "P": Pfv, "C": Cm, "groups": groups, "r": r, "Q": Q}


def factor_value(ring, fi, inst, args, eqs, es):
    if fi < len(eqs):
        return eq_factor(ring, eqs[fi], args, es)
    return ring.const(1)                      # neutral factor


def taylor_direct(n, inst, s, ring, eqs, es):
    t = s + 1
    sh = {key: {w: shift_poly_mod(pp, t, ring.J, ring.m) for w, pp in row.items()} for key, row in inst["P"].items()}
    total = ring.zero()
    for x in product((0, 1), repeat=n):
        val = ring.const(1)
        for fi, f in enumerate(inst["scopes"]):
            args = {}
            for v in f:
                acc = [0] * ring.J
                for w, coeffs in sh[(fi, v)].items():
                    if x[w]:
                        acc = [(p + q) % ring.m for p, q in zip(acc, coeffs)]
                args[v] = acc
            if fi < len(eqs):
                val = ring.mul(val, eq_factor(ring, eqs[fi], args, es))
        total = ring.add(total, val)
    return total


def taylor_separated(n, inst, s, ring, q, eqs, es):
    t = s + 1
    J, m, r = ring.J, ring.m, inst["r"]
    groups = [inst["groups"][(s, g)] for g in range(r)]
    at_t = {key: {w: fmod(Fraction(poly_eval(pp, t)), m) for w, pp in row.items()} for key, row in inst["P"].items()}
    delta = {}
    for key, row in inst["P"].items():
        dd = {}
        for w, pp in row.items():
            co = shift_poly_mod(pp, t, J, m)
            co[0] = 0
            dd[w] = co
        delta[key] = dd
    fac_by_group = defaultdict(list)
    for fi in range(len(inst["scopes"])):
        fac_by_group[inst["d"][fi][s]].append(fi)
    total = ring.zero()
    for size in range(min(q, n) + 1):
        for E in combinations(range(n), size):
            for usize in range(size + 1):
                for U in combinations(E, usize):
                    sign = -1 if (size - usize) % 2 else 1
                    kap = {}
                    for key, dd in delta.items():
                        acc = [0] * J
                        for w in U:
                            if w in dd:
                                acc = [(p + qq) % m for p, qq in zip(acc, dd[w])]
                        kap[key] = acc
                    buckets = []
                    for g in range(r):
                        S = groups[g]
                        bk = defaultdict(list)
                        for bits in product((0, 1), repeat=len(S)):
                            if any(S[i] in E and not bits[i] for i in range(len(S))):
                                continue
                            xg = dict(zip(S, bits))
                            y = tuple(sum(val * xg[w] for w, val in row.items()) for row in inst["C"][(s, g)])
                            bk[y].append(xg)
                        buckets.append(bk)
                    if any(not bk for bk in buckets):
                        continue
                    contrib = ring.zero()
                    for ys in product(*[list(bk.keys()) for bk in buckets]):
                        reps = [buckets[g][ys[g]][0] for g in range(r)]
                        prod_g = ring.const(1)
                        for g in range(r):
                            gsum = ring.zero()
                            for xg in buckets[g][ys[g]]:
                                x = [0] * n
                                for h2 in range(r):
                                    for w, bit in (xg if h2 == g else reps[h2]).items():
                                        x[w] = bit
                                val = ring.const(1)
                                for fi in fac_by_group[g]:
                                    if fi >= len(eqs):
                                        continue
                                    args = {}
                                    for v in inst["scopes"][fi]:
                                        lin = sum(c * x[w] for w, c in at_t[(fi, v)].items()) % m
                                        co = list(kap[(fi, v)])
                                        co[0] = (co[0] + lin) % m
                                        args[v] = co
                                    val = ring.mul(val, eq_factor(ring, eqs[fi], args, es))
                                gsum = ring.add(gsum, val)
                            prod_g = ring.mul(prod_g, gsum)
                        contrib = ring.add(contrib, prod_g)
                    if sign < 0:
                        contrib = ring.scale(contrib, -1)
                    total = ring.add(total, contrib)
    return total


def crt_at_zero_vec(data, Q, J, L, m):
    """Solve V h = data coefficientwise in V (Lemma 4.3); return the V-coefficients of H(0)."""
    N = Q * J
    A = []
    for s in range(Q):
        t = s + 1
        for j in range(J):
            A.append([comb(u, j) * pow(t, u - j, m) % m if u >= j else 0 for u in range(N)])
    rhs = [[data[s][j][k] for k in range(L)] for s in range(Q) for j in range(J)]
    Mx = [A[i] + rhs[i] for i in range(N)]
    for col in range(N):
        piv = next(i for i in range(col, N) if math.gcd(Mx[i][col], m) == 1)
        Mx[col], Mx[piv] = Mx[piv], Mx[col]
        inv = pow(Mx[col][col], -1, m)
        Mx[col] = [v * inv % m for v in Mx[col]]
        for i in range(N):
            if i != col and Mx[i][col]:
                f = Mx[i][col]
                Mx[i] = [(x - f * y) % m for x, y in zip(Mx[i], Mx[col])]
    return Mx[0][N:N + L]


def brute_zeros(n, eqs):
    cnt = 0
    for x in product((0, 1), repeat=n):
        if all(sum(all(x[v] for v in M) for M in monos) % 2 == 0 for monos in eqs):
            cnt += 1
    return cnt


def random_system(n, m_test, k, rng, max_monos=3):
    eqs = []
    for _ in range(m_test):
        monos = set()
        for _ in range(rng.randint(1, max_monos)):
            deg = rng.randint(1, k)
            monos.add(frozenset(rng.sample(range(n), deg)))
        eqs.append([sorted(M) for M in monos])
    return eqs


def part_Z():
    rng = random.Random(55)
    rows, ok = [], True
    for (r, Q, R) in [(2, 1, 1), (2, 2, 2), (2, 3, 2), (3, 2, 2)]:
        sw = Swap(r, Q, R)
        for rep in range(2):
            n = 3 if rep == 0 else 4
            k = 2
            eqs = random_system(n, rng.randint(2, 3), k, rng)
            used = {v for monos in eqs for M in monos for v in M}
            scopes = [sorted({v for M in monos for v in M}) for monos in eqs] + [[v] for v in range(n)]
            ell = n + 1
            es = divide_by_h(p_ell(ell))
            es = [(e[0], e[1]) for e in es][:ell]
            inst = build_general(n, scopes, sw)
            Rp = max(len(pp) - 1 for row in inst["P"].values() for pp in row.values())
            d = k * (len(eqs) + 2 * (ell - 1))                     # Lemma 5.2 degree bound
            q = (d * Rp) // Q
            J = q + 1
            # bound for the integer V-coefficients of H(0)
            ymax = max(len(monos) for monos in eqs)
            normB = max(sum(abs(a0) + abs(a1) * y for a0, a1 in es[:ell]) * (y * (y + 1)) ** ell + 1
                        for y in range(ymax + 1))
            bound = 2 ** n * normB ** len(eqs)
            p = prime_above(Q)
            mu = 1
            while p ** mu <= 4 * bound:
                mu += 1
            mod = p ** mu
            ring = Ring2(mod, J, ell)
            data, sep_ok = [], True
            for s in range(Q):
                sep = taylor_separated(n, inst, s, ring, q, eqs, es)
                dirv = taylor_direct(n, inst, s, ring, eqs, es)
                sep_ok = sep_ok and sep == dirv
                data.append(sep)
            H0 = crt_at_zero_vec(data, Q, J, ell, mod)
            bs = [x if x <= mod // 2 else x - mod for x in H0]
            count = sum(bs) % (2 ** ell)
            truth = brute_zeros(n, eqs)
            good = sep_ok and count == truth
            # negative control: interpolate from too few Taylor coefficients (J/3 per node)
            neg = None
            if Q * J > 3:
                Js = max(1, J // 3)
                short = [[row[:] for row in dat[:Js]] for dat in data]
                H0s = crt_at_zero_vec(short, Q, Js, ell, mod)
                bss = [x if x <= mod // 2 else x - mod for x in H0s]
                neg = (sum(bss) % (2 ** ell)) != truth or bss != bs
            ok = ok and good
            rows.append({"r": r, "Q": Q, "R": R, "n": n, "equations": len(eqs), "ell": ell, "J": J,
                         "untruncated_degree": k * len(eqs) * (2 * ell - 1), "truncated_degree": d,
                         "V_coeffs_of_H0": bs, "S_ell_mod_2^ell": count, "brute": truth,
                         "separated_equals_direct": sep_ok, "ok": good, "negative_control_detects": neg})
            print(f"    Z (r,Q,R)=({r},{Q},{R}) n={n} eqs={len(eqs)} N={truth} got={count} ok={good}")
    RESULTS["zeros_counting"] = {"all_ok": ok, "cases": rows}
    print(f"[Z] common-zero counting with V-truncation, mod p^mu: all exact: {ok}")
    return ok


# =============================================================== Part T

def tanner_toy(seed, loc_rows=5):
    """Lemma 6.1's mechanism on a toy: q_G = 16 vertices per side, d_G = 8 ports."""
    rng = random.Random(seed)
    lG, dG = 4, 8
    qG = 2 ** lG

    def bias(sh):
        return max(abs(sum((-1) ** bin(b & s).count("1") for s in sh)) for b in range(1, qG)) / dG

    # epsilon-biased shifts: best of several random tuples (the lemma enumerates exhaustively)
    shifts = min((tuple(rng.randrange(qG) for _ in range(dG)) for _ in range(400)), key=bias)
    eps = bias(shifts)
    # mixing bound on all pairs of small sets and on random larger ones
    mix_ok = True
    pairs = [([u], {v}) for u in range(qG) for v in range(qG)]
    for _ in range(4000):
        U = [u for u in range(qG) if rng.random() < 0.5] or [0]
        V = {v for v in range(qG) if rng.random() < 0.5} or {0}
        pairs.append((U, V))
    for U, V in pairs:
        e = sum(1 for u in U for s in shifts if (u ^ s) in V)
        if e > dG * len(U) * len(V) / qG + eps * dG * math.sqrt(len(U) * len(V)) + 1e-9:
            mix_ok = False
    # local code with large minimum distance
    best = None
    for _ in range(300):
        Aloc = [[rng.randrange(2) for _ in range(dG)] for _ in range(loc_rows)]
        ker = [w for w in product((0, 1), repeat=dG) if any(w)
               and all(sum(a * b for a, b in zip(row, w)) % 2 == 0 for row in Aloc)]
        if not ker:
            continue
        dist = min(sum(w) for w in ker)
        if best is None or dist > best[0]:
            best = (dist, Aloc)
    dloc, Aloc = best
    delta = dloc / dG
    E = qG * dG
    eidx = {(u, j): u * dG + j for u in range(qG) for j in range(dG)}
    rows = []
    for u in range(qG):
        for row in Aloc:
            rows.append({eidx[(u, j)] for j in range(dG) if row[j]})
    for v in range(qG):
        for row in Aloc:
            rows.append({eidx[(v ^ shifts[j], j)] for j in range(dG) if row[j]})
    M = [sum(1 << e for e in r) for r in rows]
    piv_rows, pivots = [], []
    for r in M:
        for pr, pc in zip(piv_rows, pivots):
            if r >> pc & 1:
                r ^= pr
        if r:
            pc = r.bit_length() - 1
            for i in range(len(piv_rows)):
                if piv_rows[i] >> pc & 1:
                    piv_rows[i] ^= r
            piv_rows.append(r)
            pivots.append(pc)
    free = [c for c in range(E) if c not in pivots]
    basis = []
    for f in free:
        vec = 1 << f
        for pr, pc in zip(piv_rows, pivots):
            if pr >> f & 1:
                vec |= 1 << pc
        basis.append(vec)
    minw = None
    if 0 < len(basis) <= 22:
        for mask in range(1, 2 ** len(basis)):
            v = 0
            for i, bvec in enumerate(basis):
                if mask >> i & 1:
                    v ^= bvec
            w = bin(v).count("1")
            minw = w if minw is None else min(minw, w)
    bound = delta * (delta - eps) * qG * dG if delta > eps else 0
    return {"seed": seed, "eps": eps, "delta_local": delta, "mixing_ok": mix_ok, "code_dim": len(basis),
            "min_weight": minw, "bound_delta(delta-eps)qd": bound, "bound_nontrivial": bound > 0,
            "ok": mix_ok and (minw is None or minw >= bound - 1e-9)}


def lemma61_constants():
    ok = True
    for B in (1, 2, 5, 17, 100):
        for m in (10, 1000, 10 ** 6, 10 ** 12):
            mu = int(math.floor(math.log2(m))) + 1
            P = 20000 * B * B * mu * mu
            D0 = 320000 * B * B * mu ** 3
            eps = 1 / (200 * B * mu)
            delta = 1 / (100 * B * mu)
            dG = min(m, D0)
            ok &= abs(eps * eps * D0 / 2 - 4 * mu) < 1e-6 * mu
            ok &= 2 * 2 ** (mu + 1) * math.exp(-4 * mu) < 1
            ok &= abs(delta * (delta - eps) - 1 / (20000 * B * B * mu * mu)) < 1e-18
            # word-count bound: (2 dG)^(2h) < 2^(dG/(25B)) with h < delta dG, 2dG <= 2^(mu+1)
            h = math.ceil(delta * dG) - 1
            if h > 0:
                ok &= 2 * h * math.log2(2 * dG) < dG / (25 * B)
            cG = math.ceil(dG / (8 * B))
            ok &= cG >= dG / (8 * B)
    return ok


def part_T():
    toys = [tanner_toy(s, 5) for s in range(3)] + [tanner_toy(s, 6) for s in range(3, 7)] + [tanner_toy(7, 7)]
    const_ok = lemma61_constants()
    ok = all(t["ok"] for t in toys) and const_ok and any(t["bound_nontrivial"] for t in toys)
    RESULTS["tanner"] = {"all_ok": ok, "toys": toys, "constants_ok": const_ok}
    print(f"[T] Tanner mechanism on toys (mixing, min weight >= delta(delta-eps)qd): "
          f"{all(t['ok'] for t in toys)}; Lemma 6.1 constant arithmetic: {const_ok}")
    return ok


# =============================================================== Part G

def violations(x, clauses):
    return tuple(0 if any((x[abs(l) - 1] == 1) == (l > 0) for l in C) else 1 for C in clauses)


def part_G():
    rng = random.Random(77)
    ok = True
    rows_out = []
    for trial in range(6):
        n = rng.randint(4, 6)
        m = 4
        clauses = []
        while len(clauses) < m:
            w = rng.randint(1, 3)
            vs = rng.sample(range(1, n + 1), w)
            clauses.append([v if rng.random() < 0.5 else -v for v in vs])
        P = m                                           # every nonzero violation vector has weight >= m/P
        cands = list(product([(e, c) for e in range(m) for c in (0, 1)], repeat=P))
        A = []
        nsat = sum(1 for x in product((0, 1), repeat=n) if not any(violations(x, clauses)))

        def N_of(Arows):
            cnt = 0
            for x in product((0, 1), repeat=n):
                vv = violations(x, clauses)
                if all(sum(vv[e] for e in row) % 2 == 0 for row in Arows):
                    cnt += 1
            return cnt

        def as_row(gamma):
            row = Counter()
            for e, c in gamma:
                if c:
                    row[e] ^= 1
            return [e for e, b in row.items() if b]

        Ds = [N_of(A) - nsat]
        exp_ok = True
        for rnd in range(3 * n):
            Ns = [N_of(A + [as_row(g)]) for g in cands]
            Dmean = sum(Nv - nsat for Nv in Ns) / len(Ns)
            exp_ok &= Dmean <= 0.75 * Ds[-1] + 1e-12
            best = min(range(len(cands)), key=lambda i: Ns[i])
            A.append(as_row(cands[best]))
            Ds.append(Ns[best] - nsat)
            if Ds[-1] == 0:
                break
        good = exp_ok and Ds[-1] == 0 and all(Ds[i + 1] <= 0.75 * Ds[i] for i in range(len(Ds) - 1))
        ok = ok and good
        rows_out.append({"n": n, "m": m, "#F": nsat, "D_sequence": Ds, "rounds": len(Ds) - 1, "ok": good})
    RESULTS["greedy_rows"] = {"all_ok": ok, "cases": rows_out}
    print(f"[G] Lemma 6.2 greedy rows on small CNFs: D shrinks by >= 1/4 per round to 0: {ok}")
    return ok


def main():
    oks = [part_L(), part_Z(), part_T(), part_G()]
    RESULTS["all_ok"] = all(oks)
    with open(__file__.rsplit("/", 2)[0] + "/results/sharp_route.json", "w") as fh:
        json.dump(RESULTS, fh, indent=1, default=str)
    print("ALL OK" if all(oks) else "FAILURES")


if __name__ == "__main__":
    main()

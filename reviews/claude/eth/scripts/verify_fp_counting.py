"""Checks for 12_FP_COUNTING.md: counting common zeros of sparse polynomial systems over F_p.

Part L  [Lemma 2.1, Remark 2.3] The p-adic lift.  (i) The lift p_ell = F_ell(Y^(p-1)) of toolkit E2
        and its h-adic digits (h = Y^p - Y, deg e_j < p) satisfy S_ell(prod B_i mod V^ell) =
        prod [p | y_i] mod p^ell on random integer tuples.  (ii) The closed-form idempotent
        E~_ell = (1 - z) + sum_{1<=j<ell} (-1)^(j+1) C(2j-1, j-1) (2z - 1) Y^((p-2) j) h^j,
        z = Y^(p-1), satisfies E~_ell(y) = [p | y] mod p^ell, agrees with p_ell modulo h^ell
        (exact polynomial division), and for p = 2 has the digits of Part II's Lemma 2.1.
        p in {2, 3, 5, 7}, ell <= 8.  Negative control: the unhalved binomial fails.
Part Z  [Theorem 1.1] The whole pipeline over F_3 and F_5.  Each F_p-variable is encoded by
        t = ceil(log2 p) bits, a validity factor [value < p] is added per variable, each
        equation becomes the p-adic factor B_i of its integer lift, and the product is
        truncated in (Z/P^mu)[V]/(V^ell) with ell = n + 1.  The deformation of Part I
        (closed-form swap rows, derandomized schedules, measured partitions) gives the
        Taylor jets at the Q nodes by separation, which agree with direct evaluation; CRT at
        zero gives the V-coefficients exactly and S_ell mod p^ell is compared with brute
        force over F_p^n.  Negative controls: without the validity factors, and with too few
        Taylor coefficients, the count is wrong.

Writes ../results/fp_counting.json.
"""

import json
import math
import random
import sys
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import comb

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from eth_core import poly_eval  # noqa: E402
from verify_gapless_chain import Swap, prime_above, shift_poly_mod, fmod  # noqa: E402
from verify_sharp_route import Ring2, zmul, build_general, crt_at_zero_vec  # noqa: E402
from verify_padic_truncation import p_l, divide_by_h, pmul, padd, trim, ev  # noqa: E402

RESULTS = {}


# =============================================================== Part L

def closed_idempotent(p, ell):
    """E~_ell as an integer polynomial (list, low -> high)."""
    z = [0] * (p - 1) + [1]                           # Y^(p-1)
    h = [0, -1] + [0] * (p - 2) + [1]                 # Y^p - Y
    out = padd([1], [-c for c in z])                  # 1 - z
    hp = [1]
    for j in range(1, ell):
        hp = pmul(hp, h)
        c = (-1) ** (j + 1) * comb(2 * j - 1, j - 1)
        term = pmul([c * x for x in padd([-1], [2 * x for x in z])], [0] * ((p - 2) * j) + [1])
        out = padd(out, pmul(term, hp))
    return trim(out)


def poly_rem_monic(num, den):
    """Remainder of num modulo the monic integer polynomial den."""
    num = list(num)
    dn = len(den) - 1
    for k in range(len(num) - 1, dn - 1, -1):
        c = num[k]
        if c:
            for i, dv in enumerate(den):
                num[k - dn + i] -= c * dv
    return trim(num[:dn]) if len(num) > dn else trim(num)


def part_L():
    rng = random.Random(3)
    ok = True
    rows = []
    for p in (2, 3, 5, 7):
        for ell in range(1, 9):
            pl = p_l(p, ell)
            es = divide_by_h(pl, p)[:ell]
            h = [0, -1] + [0] * (p - 2) + [1]
            # (i) truncation identity on random tuples
            trunc_ok = True
            for _ in range(40):
                ys = [rng.randint(-30, 30) for _ in range(rng.randint(1, 4))]
                prod_v = [1] + [0] * (ell - 1)
                for y in ys:
                    hy = ev(h, y)
                    Bi = [ev(es[j], y) * hy ** j for j in range(ell)]
                    newv = [0] * ell
                    for a_, x in enumerate(prod_v):
                        for b_, yv in enumerate(Bi):
                            if a_ + b_ < ell:
                                newv[a_ + b_] += x * yv
                    prod_v = newv
                want = 1 if all(y % p == 0 for y in ys) else 0
                trunc_ok &= (sum(prod_v) - want) % p ** ell == 0
            # (ii) closed-form idempotent
            Et = closed_idempotent(p, ell)
            cong_ok = all((ev(Et, y) - (1 if y % p == 0 else 0)) % p ** ell == 0 for y in range(-60, 61))
            hl = [1]
            for _ in range(ell):
                hl = pmul(hl, h)
            diff = padd(Et, [-c for c in pl])
            rem = poly_rem_monic(diff, hl)
            same_mod = not any(rem)
            ok &= trunc_ok and cong_ok and same_mod
            if ell in (1, 4, 8):
                rows.append({"p": p, "ell": ell, "E2_truncation": trunc_ok, "closed_form_congruence": cong_ok,
                             "closed_form = p_ell mod h^ell": same_mod})
    # p = 2: the closed form is Part II's Lemma 2.1 (digits (2Y-1) * signed central binomials)
    p2 = closed_idempotent(2, 6)
    lemma21 = [1, -1]
    h2 = [0, -1, 1]
    hp = [1]
    for j in range(1, 6):
        hp = pmul(hp, h2)
        c = (-1) ** (j + 1) * comb(2 * j - 1, j - 1)
        lemma21 = padd(lemma21, pmul([-c, 2 * c], hp))
    ok &= trim(p2) == trim(lemma21)
    # negative control: binom(2j, j) instead of binom(2j-1, j-1) for p = 3
    bad = [1, 0, -1]
    h3 = [0, -1, 0, 1]
    hp = [1]
    for j in range(1, 5):
        hp = pmul(hp, h3)
        c = (-1) ** (j + 1) * comb(2 * j, j)
        bad = padd(bad, pmul(pmul([-c, 0, 2 * c], [0] * j + [1]), hp))
    neg = any((ev(bad, y) - (1 if y % 3 == 0 else 0)) % 3 ** 5 for y in range(-20, 21))
    ok &= neg
    RESULTS["lift"] = {"all_ok": ok, "rows": rows, "p2_equals_part_II_lemma_2_1": trim(p2) == trim(lemma21),
                       "negative_control_detects": neg}
    print(f"[L] p-adic lift (E2 digits and the closed-form idempotent), p = 2,3,5,7, ell <= 8: {ok}")
    return ok


# =============================================================== Part Z

def int_bits(p):
    return max(1, math.ceil(math.log2(p)))


def validity_coeffs(p, t):
    """Moebius coefficients of [sum 2^b y_b < p] on t bits: dict frozenset(bits) -> integer."""
    vals = {}
    for bits in product((0, 1), repeat=t):
        vals[bits] = 1 if sum(b << i for i, b in enumerate(bits)) < p else 0
    co = {}
    for S in range(1 << t):
        Sset = [i for i in range(t) if S >> i & 1]
        c = 0
        for Tm in range(1 << t):
            if Tm & ~S:
                continue
            Tbits = tuple(1 if (Tm >> i & 1) else 0 for i in range(t))
            c += (-1) ** (len(Sset) - bin(Tm).count("1")) * vals[Tbits]
        if c:
            co[frozenset(Sset)] = c
    return co


def zpow(a, e, J, m):
    out = [1] + [0] * (J - 1)
    for _ in range(e):
        out = zmul(out, a, J, m)
    return out


def make_factor_fn(ring, p, t, eqs, nF, es):
    """Factors: equations first, then one validity factor per F_p-variable.
    Bit variable index of (v, b) is v * t + b."""
    J, m = ring.J, ring.m
    vco = validity_coeffs(p, t)

    def eq_fac(eq, args):
        # x_v = sum_b 2^b Y_{v,b}  as z-polynomials
        xs = {}
        for (c, mono) in eq:
            for v in mono:
                if v not in xs:
                    acc = [0] * J
                    for b in range(t):
                        acc = [(q + (1 << b) * w) % m for q, w in zip(acc, args[v * t + b])]
                    xs[v] = acc
        L = [0] * J
        for (c, mono) in eq:
            term = [c % m] + [0] * (J - 1)
            for v, e in mono.items():
                term = zmul(term, zpow(xs[v], e, J, m), J, m)
            L = [(q + w) % m for q, w in zip(L, term)]
        Lp = zpow(L, p, J, m)
        h = [(q - w) % m for q, w in zip(Lp, L)]
        out = ring.zero()
        hp = [1] + [0] * (J - 1)
        for j in range(ring.L):
            ej = [0] * J
            Lpow = [1] + [0] * (J - 1)
            for coef in es[j]:
                ej = [(q + coef * w) % m for q, w in zip(ej, Lpow)]
                Lpow = zmul(Lpow, L, J, m)
            col = zmul(ej, hp, J, m)
            for i in range(J):
                out[i][j] = col[i]
            hp = zmul(hp, h, J, m)
        return out

    def val_fac(v, args):
        acc = [0] * J
        for S, c in vco.items():
            term = [c % m] + [0] * (J - 1)
            for b in S:
                term = zmul(term, args[v * t + b], J, m)
            acc = [(q + w) % m for q, w in zip(acc, term)]
        return ring.from_z(acc)

    cache = {}

    def factor_fn(fi, args, with_validity=True):
        if fi >= len(eqs) and not with_validity:
            return ring.const(1)
        key = (fi, tuple(sorted((v, tuple(a)) for v, a in args.items())))
        hit = cache.get(key)
        if hit is None:
            hit = eq_fac(eqs[fi], args) if fi < len(eqs) else val_fac(fi - len(eqs), args)
            cache[key] = hit
        return hit
    return factor_fn


def taylor_direct_g(nb, inst, s, ring, factor_fn):
    t = s + 1
    sh = {key: {w: shift_poly_mod(pp, t, ring.J, ring.m) for w, pp in row.items()} for key, row in inst["P"].items()}
    total = ring.zero()
    for x in product((0, 1), repeat=nb):
        val = ring.const(1)
        for fi, f in enumerate(inst["scopes"]):
            args = {}
            for v in f:
                acc = [0] * ring.J
                for w, coeffs in sh[(fi, v)].items():
                    if x[w]:
                        acc = [(a + b) % ring.m for a, b in zip(acc, coeffs)]
                args[v] = acc
            val = ring.mul(val, factor_fn(fi, args))
        total = ring.add(total, val)
    return total


def taylor_separated_g(nb, inst, s, ring, q, factor_fn):
    """Lemma 4.4 (Moebius over (E, U)) and Lemma 4.5 (separation by the measured partition)."""
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
    for size in range(min(q, nb) + 1):
        for E in combinations(range(nb), size):
            for usize in range(size + 1):
                for U in combinations(E, usize):
                    sign = -1 if (size - usize) % 2 else 1
                    kap = {}
                    for key, dd in delta.items():
                        acc = [0] * J
                        for w in U:
                            if w in dd:
                                acc = [(a + b) % m for a, b in zip(acc, dd[w])]
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
                                x = [0] * nb
                                for h2 in range(r):
                                    for w, bit in (xg if h2 == g else reps[h2]).items():
                                        x[w] = bit
                                val = ring.const(1)
                                for fi in fac_by_group[g]:
                                    args = {}
                                    for v in inst["scopes"][fi]:
                                        lin = sum(c * x[w] for w, c in at_t[(fi, v)].items()) % m
                                        co = list(kap[(fi, v)])
                                        co[0] = (co[0] + lin) % m
                                        args[v] = co
                                    val = ring.mul(val, factor_fn(fi, args))
                                gsum = ring.add(gsum, val)
                            prod_g = ring.mul(prod_g, gsum)
                        contrib = ring.add(contrib, prod_g)
                    if sign < 0:
                        contrib = ring.scale(contrib, -1)
                    total = ring.add(total, contrib)
    return total


def random_fp_system(p, nF, m_eq, k, rng):
    eqs = []
    for _ in range(m_eq):
        eq = []
        for _ in range(rng.randint(1, 2)):
            deg = rng.randint(1, k)
            mono = {}
            for _ in range(deg):
                v = rng.randrange(nF)
                mono[v] = mono.get(v, 0) + 1
            eq.append((rng.randint(1, p - 1), mono))
        eqs.append(eq)
    return eqs


def brute_fp(p, nF, eqs):
    cnt = 0
    for x in product(range(p), repeat=nF):
        if all(sum(c * math.prod(x[v] ** e for v, e in mono.items()) for c, mono in eq) % p == 0 for eq in eqs):
            cnt += 1
    return cnt


def part_Z():
    rng = random.Random(1212)
    ok, rows = True, []
    neg_valid = 0
    neg_jets = 0
    jets_cases = 0
    cases = [(3, 2, 2, (2, 1, 1)), (3, 2, 2, (2, 2, 2)), (3, 2, 2, (3, 2, 2)), (3, 3, 2, (2, 1, 1)),
             (3, 3, 2, (2, 2, 2)), (5, 2, 2, (2, 1, 1)), (5, 2, 1, (2, 2, 2)), (7, 2, 1, (2, 1, 1))]
    for (p, nF, m_eq, (r, Q, R)) in cases:
        k = 2
        t = int_bits(p)
        nb = nF * t
        while True:
            eqs = random_fp_system(p, nF, m_eq, k, rng)
            truth = brute_fp(p, nF, eqs)
            if 0 < truth < p ** nF:
                break
        ell = nF + 1
        es = divide_by_h(p_l(p, ell), p)[:ell]
        es += [[0] * p] * (ell - len(es))
        scopes = [sorted({v * t + b for _, mono in eq for v in mono for b in range(t)}) for eq in eqs] \
            + [[v * t + b for b in range(t)] for v in range(nF)]
        sw = Swap(r, Q, R)
        inst = build_general(nb, scopes, sw)
        Rp = max(len(pp) - 1 for row in inst["P"].values() for pp in row.values())
        d = nF * t + k * (len(eqs) * (p - 1) + p * (ell - 1))          # degree in A[Y]
        q = (d * Rp) // Q
        J = q + 1
        # bound for the integer V-coefficients of H(0)
        Lmax = sum(c * math.prod(((1 << t) - 1) ** e for e in mono.values()) for eq in eqs for c, mono in eq)
        normB = max(sum(abs(ev(es[j], y)) * abs(y ** p - y) ** j for j in range(ell)) for y in range(Lmax + 1))
        bound = 2 ** nb * normB ** len(eqs) + 1
        P = prime_above(Q)
        mu = 1
        while P ** mu <= 4 * bound:
            mu += 1
        mod = P ** mu
        ring = Ring2(mod, J, ell)
        fac = make_factor_fn(ring, p, t, eqs, nF, es)
        data, sep_ok = [], True
        for s in range(Q):
            sep = taylor_separated_g(nb, inst, s, ring, q, fac)
            if s == 0:
                dirv = taylor_direct_g(nb, inst, s, ring, fac)
                sep_ok = sep_ok and sep == dirv
            data.append(sep)
        H0 = crt_at_zero_vec(data, Q, J, ell, mod)
        bs = [x if x <= mod // 2 else x - mod for x in H0]
        count = sum(bs) % (p ** ell)
        good = sep_ok and count == truth
        ok &= good
        # negative control 1: no validity factors (direct evaluation suffices for the identity)
        facnv = lambda fi, args: fac(fi, args, with_validity=False)
        data_nv = [taylor_direct_g(nb, inst, s, ring, facnv) for s in range(Q)]
        H0nv = crt_at_zero_vec(data_nv, Q, J, ell, mod)
        cnt_nv = sum(x if x <= mod // 2 else x - mod for x in H0nv) % (p ** ell)
        neg_valid += cnt_nv != truth
        # negative control 2: too few Taylor coefficients
        jets_applicable = Q * J > 3
        if jets_applicable:
            Js = max(1, J // 3)
            short = [[row[:] for row in dat[:Js]] for dat in data]
            H0s = crt_at_zero_vec(short, Q, Js, ell, mod)
            cnt_s = sum(x if x <= mod // 2 else x - mod for x in H0s) % (p ** ell)
            neg_jets += cnt_s != truth
            jets_cases += 1
        rows.append({"p": p, "F_p variables": nF, "bits": nb, "equations": [[(c, {str(v): e for v, e in mono.items()})
                     for c, mono in eq] for eq in eqs], "(r,Q,R)": (r, Q, R), "ell": ell, "degree_d": d, "J": J,
                     "V_coeffs_of_H0": bs, "count": count, "brute": truth, "separated_equals_direct(node 1)": sep_ok,
                     "without_validity_factors": cnt_nv, "ok": good})
        print(f"    F_{p}: n = {nF} ({nb} bits), {len(eqs)} equations, (r,Q,R) = ({r},{Q},{R}), N = {truth}, "
              f"pipeline = {count}, ok = {good}")
    ok &= neg_valid >= 1 and jets_cases > 0 and neg_jets == jets_cases
    RESULTS["pipeline"] = {"all_ok": ok, "cases": rows, "negative_validity_detects": f"{neg_valid}/{len(rows)}",
                           "negative_jets_detects": f"{neg_jets}/{jets_cases} (cases with Q*J > 3)"}
    print(f"[Z] F_p pipeline exact: {all(r['ok'] for r in rows)}; without validity factors wrong in "
          f"{neg_valid}/{len(rows)}; too few jets wrong in {neg_jets}/{jets_cases} applicable cases")
    return ok


def main():
    oks = [part_L(), part_Z()]
    RESULTS["all_ok"] = all(oks)
    with open(__file__.rsplit("/", 2)[0] + "/results/fp_counting.json", "w") as fh:
        json.dump(RESULTS, fh, indent=1, default=str)
    print("ALL OK" if all(oks) else "FAILURES")


if __name__ == "__main__":
    main()

"""Checks for 10_GAPLESS_TENSOR_CONTRACTION.md (the paper's Appendix D, Theorem 1.1).

Labels in brackets are the doc's numbering; the paper's follows.

Part C  [Lemma 6.3] The child identity (D.16), end to end and exactly.  Random bipartite tensor
        networks are packed as in Lemma D.3 (labels in {0,1}^a, packets of L = 2^a slots with
        orientation bits).  At nodes t with v2(t) = i and values zeta in O = Z_(2)[U]/(f), the
        deformed contraction F(t + alpha zeta) is computed (i) directly over all 2^(2M) port
        assignments and (ii) as the sum over mode tuples beta of the product over groups g of
        the child contractions Z_{g,beta}, each child contracted by brute force.  They agree
        modulo 2^p' for every case.  Negative controls: dropping the crossing terms (all modes),
        or crediting each mode to its source group instead of its target, breaks the identity.
        The tile expansion (D.12) behind the modes is checked separately modulo 2^40, where the
        modes of weight >= 2 matter: keeping only b = 0 fails.
Part F  [Lemmas 3.2, 4.1] One level of the recursion at the honest precision: the jets of
        F(t + alpha z) mod 2^(p_t), p_t = p + 8 h D_t + 8 M (a+1), are obtained by evaluating at
        M(2L-1)+1 binary-coordinate points of O and Lagrange interpolation over
        (Z/2^(p_t))[U]/(f); tapered Hermite recovery then returns the contraction mod 2^p.
        Negative control: recovery with p_t = p fails.
Part M  [Lemma 6.4] The mode-tuple count: exact counts by dynamic programming against the
        bound exp(2 pi sqrt(c p / 6)) for c variables of every weight below p.
Part L  [Lemma 9.1] Paper Lemma D.2's encoding of a sparse F_2 system: tests and variables in
        groups, incidence subsets indexed by ceil(log2 binom(N_gh, d_gh)) wires per group
        pair, one network per degree array; the sum of the contractions, after S_ell, equals
        the number of common zeros modulo 2^(n+1).
Part D  [Section 8] Paper D.2's clause network (satisfaction tables, equality trees,
        subdivision): its contraction equals #F.
Part P  [Lemma 7.1] The recursion's parameter inequalities, in logarithmic form at every active
        depth: (A1)-(A2) (precision absorption), (D.6), the number of interpolation points and
        the packet-entry cost.  Reports the least lambda = ceil(log2 s) from which all hold
        (checked interval by interval between depth transitions, then on a grid to
        log2(lambda) ~ 260), and where the paper's floor u >= lambda^(4/5) would start to hold.

Writes ../results/tensor_contraction.json.
"""

import json
import math
import random
import sys
from fractions import Fraction
from itertools import product
from math import comb

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from verify_parity_lift import p_ell, divide_by_h  # noqa: E402

RESULTS = {}


# =============================================================== the rings K, O/2^P

class K:
    """K = Q[U]/(f) with f monic, irreducible mod 2.  Elements: tuples of Fractions."""

    def __init__(self, f):
        self.f = [int(c) for c in f]          # low -> high, monic
        self.e = len(f) - 1

    def red(self, c):
        c = [Fraction(x) for x in c]
        e = self.e
        for k in range(len(c) - 1, e - 1, -1):
            ck = c[k]
            if ck:
                for j in range(e + 1):
                    c[k - e + j] -= ck * self.f[j]
        c = c[:e] + [Fraction(0)] * max(0, e - len(c))
        return tuple(c[:e])

    def const(self, n):
        return tuple([Fraction(n)] + [Fraction(0)] * (self.e - 1))

    def add(self, a, b):
        return tuple(x + y for x, y in zip(a, b))

    def sub(self, a, b):
        return tuple(x - y for x, y in zip(a, b))

    def scal(self, c, a):
        return tuple(Fraction(c) * x for x in a)

    def mul(self, a, b):
        out = [Fraction(0)] * (2 * self.e - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y:
                        out[i + j] += x * y
        return self.red(out)

    def inv(self, a):
        e = self.e
        # multiplication matrix: column j = a * U^j
        cols = []
        for j in range(e):
            uj = [Fraction(0)] * e
            uj[j] = Fraction(1)
            cols.append(self.mul(a, tuple(uj)))
        A = [[cols[j][i] for j in range(e)] + [Fraction(1 if i == 0 else 0)] for i in range(e)]
        for c in range(e):
            piv = next(r for r in range(c, e) if A[r][c] != 0)
            A[c], A[piv] = A[piv], A[c]
            pv = A[c][c]
            A[c] = [x / pv for x in A[c]]
            for r in range(e):
                if r != c and A[r][c] != 0:
                    fac = A[r][c]
                    A[r] = [x - fac * y for x, y in zip(A[r], A[c])]
        return tuple(A[i][e] for i in range(e))

    def v2(self, a):
        best = None
        for x in a:
            if x != 0:
                n, d = x.numerator, x.denominator
                v = 0
                while n % 2 == 0:
                    n //= 2
                    v += 1
                while d % 2 == 0:
                    d //= 2
                    v -= 1
                best = v if best is None else min(best, v)
        return 10 ** 9 if best is None else best


class R:
    """(Z/2^P)[U]/(f): elements are tuples of ints in [0, 2^P)."""

    def __init__(self, f, P):
        self.f = [int(c) for c in f]
        self.e = len(f) - 1
        self.P = P
        self.m = 1 << P

    def zero(self):
        return (0,) * self.e

    def one(self):
        return (1,) + (0,) * (self.e - 1)

    def const(self, n):
        return (n % self.m,) + (0,) * (self.e - 1)

    def add(self, a, b):
        return tuple((x + y) % self.m for x, y in zip(a, b))

    def sub(self, a, b):
        return tuple((x - y) % self.m for x, y in zip(a, b))

    def scal(self, c, a):
        return tuple((c * x) % self.m for x in a)

    def mul(self, a, b):
        e = self.e
        out = [0] * (2 * e - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y:
                        out[i + j] += x * y
        for k in range(len(out) - 1, e - 1, -1):
            ck = out[k]
            if ck:
                for j in range(e + 1):
                    out[k - e + j] -= ck * self.f[j]
        return tuple(x % self.m for x in out[:e])

    def from_K(self, a):
        out = []
        for x in a:
            n, d = x.numerator, x.denominator
            assert d % 2 == 1, "not 2-adically integral"
            out.append((n * pow(d, -1, self.m)) % self.m)
        return tuple(out)

    def reduce(self, a, P2):
        m2 = 1 << P2
        return tuple(x % m2 for x in a)


def irreducible_mod2(e):
    """Least monic f of degree e with 0/1 coefficients, irreducible over F_2 (by trial division)."""
    def pmod(a, b):
        a = a[:]
        while len(a) >= len(b) and any(a):
            if a[-1] == 0:
                a.pop()
                continue
            sh = len(a) - len(b)
            for i, c in enumerate(b):
                a[sh + i] ^= c
            a.pop()
        while a and a[-1] == 0:
            a.pop()
        return a
    for mask in range(1 << e):
        f = [(mask >> k) & 1 for k in range(e)] + [1]
        if f[0] == 0:
            continue
        ok = True
        for dg in range(1, e // 2 + 1):
            for gm in range(1 << dg):
                g = [(gm >> k) & 1 for k in range(dg)] + [1]
                if not pmod(f, g):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return f
    raise ValueError


def binary_points(Kf, count):
    """Elements of O with coordinates in {0,1}; pairwise differences are units."""
    pts = []
    for mask in range(1 << Kf.e):
        pts.append(tuple(Fraction((mask >> k) & 1) for k in range(Kf.e)))
        if len(pts) == count:
            return pts
    raise ValueError("need a larger residue field")


# =============================================================== networks and packing

def random_network(rng, n_rows, n_cols, M, max_deg, modulus):
    """Bipartite multigraph: edges (row vertex, col vertex); tensors on each vertex.  Resampled
    until the contraction is nonzero modulo `modulus`, so that the identities are not vacuous."""
    while True:
        net = _random_network(rng, n_rows, n_cols, M, max_deg, modulus)
        if direct_contraction(net, modulus) % 2:
            return net


def _random_network(rng, n_rows, n_cols, M, max_deg, modulus):
    while True:
        edges = [(rng.randrange(n_rows), rng.randrange(n_cols)) for _ in range(M)]
        rdeg = [sum(1 for e in edges if e[0] == v) for v in range(n_rows)]
        cdeg = [sum(1 for e in edges if e[1] == v) for v in range(n_cols)]
        if max(rdeg) <= max_deg and max(cdeg) <= max_deg:
            break
    rports = [[k for k, e in enumerate(edges) if e[0] == v] for v in range(n_rows)]
    cports = [[k for k, e in enumerate(edges) if e[1] == v] for v in range(n_cols)]
    # 0/1 tables (constraint-like), except one row vertex with arbitrary entries: the contraction
    # is then a generic residue rather than a sum of many products of even numbers
    rten = [{bits: (rng.randrange(modulus) if v == 0 else int(rng.random() < 0.75))
             for bits in product((0, 1), repeat=len(p))} for v, p in enumerate(rports)]
    cten = [{bits: int(rng.random() < 0.75) for bits in product((0, 1), repeat=len(p))} for p in cports]
    return {"edges": edges, "rports": rports, "cports": cports, "rten": rten, "cten": cten}


def pack(rng, net, a):
    """Lemma D.3's structure with r = 2: labels in {0,1}^a; a packet has L = 2^a slots; in a
    coordinate with unequal symbols bit i of the slot is the row symbol; equal coordinates index
    copies.  Returns slot[edge] = (packet, slot) and the labels."""
    nr, nc = len(net["rports"]), len(net["cports"])
    rlab = [tuple(rng.randrange(2) for _ in range(a)) for _ in range(nr)]
    clab = [tuple(rng.randrange(2) for _ in range(a)) for _ in range(nc)]
    classes = {}
    for k, (u, v) in enumerate(net["edges"]):
        rho, kap = rlab[u], clab[v]
        pattern = tuple((rho[i], kap[i]) if rho[i] == kap[i] else ("x",) for i in range(a))
        orient = tuple(rho[i] for i in range(a) if rho[i] != kap[i])
        classes.setdefault(pattern, {}).setdefault(orient, []).append(k)
    slot_of = {}
    packets = []
    for pattern, byorient in classes.items():
        unequal = [i for i in range(a) if pattern[i] == ("x",)]
        equal = [i for i in range(a) if pattern[i] != ("x",)]
        ncopy = 1 << len(equal)
        nbatch = max(math.ceil(len(v) / ncopy) for v in byorient.values())
        for bidx in range(nbatch):
            pid = len(packets)
            packets.append(pattern)
            for orient, lst in byorient.items():
                chunk = lst[bidx * ncopy:(bidx + 1) * ncopy]
                for c, k in enumerate(chunk):
                    s = 0
                    for idx, i in enumerate(unequal):
                        s |= orient[idx] << i
                    for idx, i in enumerate(equal):
                        s |= ((c >> idx) & 1) << i
                    slot_of[k] = (pid, s)
    return slot_of, packets, rlab, clab


# =============================================================== the Cauchy kernel

def cauchy(Kf, a):
    L = 1 << a
    S = list(range(-L, L))
    Cp0 = 1
    for q in S:
        if q != 0:
            Cp0 *= -q
    return L, S, Cp0


def W_entry(Kf, S, Cp0, T, d):
    """W_xy(T) with y - x = d: prod_{q in S, q != d} (T - q) / C'(0)."""
    out = Kf.const(1)
    for q in S:
        if q != d:
            out = Kf.mul(out, Kf.sub(T, Kf.const(q)))
    return Kf.scal(Fraction(1, Cp0), out)


def deformed_contraction(Kf, Rp, net, slot_of, a, T):
    """F(T) = sum_{X,Y} A(Y) B(X) prod_rows l_{Y_x}((W(T) X)_x), brute force, in Rp."""
    L, S, Cp0 = cauchy(Kf, a)
    M = len(net["edges"])
    at = {v: k for k, v in slot_of.items()}           # (packet, slot) -> edge
    Wc = {}
    rows = []                                          # for each edge as a row port: list of (col edge, W)
    for k in range(M):
        pid, x = slot_of[k]
        terms = []
        for y in range(L):
            k2 = at.get((pid, y))
            if k2 is None:
                continue
            d = y - x
            if d not in Wc:
                Wc[d] = Rp.from_K(W_entry(Kf, S, Cp0, T, d))
            terms.append((k2, Wc[d]))
        rows.append(terms)
    total = Rp.zero()
    for Y in product((0, 1), repeat=M):
        Aval = Rp.one()
        for v, ports in enumerate(net["rports"]):
            Aval = Rp.scal(net["rten"][v][tuple(Y[k] for k in ports)], Aval)
        if not any(Aval):
            continue
        for X in product((0, 1), repeat=M):
            Bval = Aval
            for v, ports in enumerate(net["cports"]):
                Bval = Rp.scal(net["cten"][v][tuple(X[k] for k in ports)], Bval)
            if not any(Bval):
                continue
            val = Bval
            for k in range(M):
                w = Rp.zero()
                for k2, Wv in rows[k]:
                    if X[k2]:
                        w = Rp.add(w, Wv)
                f = w if Y[k] else Rp.sub(Rp.one(), w)
                val = Rp.mul(val, f)
                if not any(val):
                    break
            total = Rp.add(total, val)
    return total


def direct_contraction(net, modulus):
    M = len(net["edges"])
    tot = 0
    for x in product((0, 1), repeat=M):
        val = 1
        for v, ports in enumerate(net["rports"]):
            val = val * net["rten"][v][tuple(x[k] for k in ports)]
        for v, ports in enumerate(net["cports"]):
            val = val * net["cten"][v][tuple(x[k] for k in ports)]
        tot += val
    return tot % modulus


# =============================================================== Part C: the child identity

def modes_and_children(Kf, Rp, net, slot_of, rlab, clab, a, t, T, pprime, drop_modes=False,
                       swap_roles=False):
    """Sum over mode tuples beta of prod_g Z_{g,beta}, with the children contracted by brute force."""
    L, S, Cp0 = cauchy(Kf, a)
    i = 0
    while t % (2 ** (i + 1)) == 0:
        i += 1
    M = len(net["edges"])
    at = {v: k for k, v in slot_of.items()}
    rgroup = [rlab[u][i] for (u, v) in net["edges"]]   # group of the row port of edge k
    cgroup = [clab[v][i] for (u, v) in net["edges"]]   # group of the column port of edge k
    # local terms and modes
    local = {}                                          # row edge -> list of (col edge, W)
    modes = []                                          # (target, source, weight, packet, terms)
    row_modes = {k: [] for k in range(M)}               # row edge -> list of (mode index, kappa)
    tiles_done = {}
    for k in range(M):
        pid, x = slot_of[k]
        local[k] = []
        for y in range(L):
            k2 = at.get((pid, y))
            if k2 is None:
                continue
            if cgroup[k2] == rgroup[k]:
                local[k].append((k2, Rp.from_K(W_entry(Kf, S, Cp0, T, y - x))))
        for y0 in range(2 ** (i + 1)):
            tile = [y for y in range(y0, L, 2 ** (i + 1))]
            real = [(at[(pid, y)], j) for j, y in enumerate(tile) if (pid, y) in at]
            crossing = [cgroup[k2] != rgroup[k] for k2, _ in real]
            if not real:
                # a tile with no real column: modes would have L_Z = 0; skip (it cannot contribute)
                continue
            assert all(crossing) or not any(crossing), "tile not homogeneous"
            if not any(crossing):
                continue
            # bit i of y0 equals bit i of x (paper's characterization of crossing tiles)
            assert ((y0 >> i) & 1) == ((x >> i) & 1)
            dxy = y0 - x
            assert dxy in S
            Tm = Kf.sub(T, Kf.const(dxy))
            assert Kf.v2(Kf.sub(Kf.const(t), Kf.const(dxy))) <= i
            base = W_entry(Kf, S, Cp0, T, dxy)
            ratio = Kf.scal(2 ** (i + 1), Kf.inv(Tm))
            key = (pid, y0)
            src = cgroup[real[0][0]]
            tgt = rgroup[k]
            assert all(cgroup[k2] == src for k2, _ in real)
            powr = Kf.const(1)
            for b in range(pprime - 1):
                if drop_modes:
                    break
                kap = Kf.mul(base, powr)
                assert Kf.v2(kap) >= b + 1
                if (key, b) in tiles_done:
                    assert modes[tiles_done[(key, b)]]["target" if not swap_roles else "source"] == tgt
                else:
                    tiles_done[(key, b)] = len(modes)
                    terms = [(k2, j ** b) for k2, j in real]
                    s_, t_ = (tgt, src) if swap_roles else (src, tgt)
                    modes.append({"target": t_, "source": s_, "w": b + 1, "pid": pid, "terms": terms})
                row_modes[k].append((tiles_done[(key, b)], Rp.from_K(kap)))
                powr = Kf.mul(powr, ratio)
    # enumerate beta with sum w beta < p'
    nm = len(modes)
    betas = []

    def gen(idx, budget, cur):
        if idx == nm:
            betas.append(tuple(cur))
            return
        w = modes[idx]["w"]
        c = 0
        while c * w <= budget:
            cur.append(c)
            gen(idx + 1, budget - c * w, cur)
            cur.pop()
            c += 1
    gen(0, pprime - 1, [])
    rvg = [rlab[v][i] for v in range(len(net["rports"]))]      # vertex groups (isolated ones too)
    cvg = [clab[v][i] for v in range(len(net["cports"]))]
    groups = sorted(set(rvg) | set(cvg))
    tables = {g: child_tables(Rp, net, g, modes, local, row_modes, rgroup, cgroup, slot_of, pprime,
                              rvg, cvg)
              for g in groups}
    total = Rp.zero()
    for beta in betas:
        prod_g = Rp.one()
        for g in groups:
            prod_g = Rp.mul(prod_g, child_value(Rp, tables[g], beta, modes, g))
            if not any(prod_g):
                break
        total = Rp.add(total, prod_g)
    return total, len(betas), nm


def poly_mul(Rp, P1, P2, maxw, wt):
    out = {}
    for m1, c1 in P1.items():
        for m2, c2 in P2.items():
            m = tuple(x + y for x, y in zip(m1, m2))
            if sum(w * e for w, e in zip(wt, m)) > maxw:
                continue
            out[m] = Rp.add(out.get(m, Rp.zero()), Rp.mul(c1, c2))
    return out


def child_tables(Rp, net, g, modes, local, row_modes, rgroup, cgroup, slot_of, pprime, rvg, cvg):
    """For the child of group g: for every assignment of its port bits, the product of its old
    vertex tensors, and per packet the row-factor product as a polynomial in the modes and the
    values L_Z(X) of the modes sourced at g."""
    M = len(net["edges"])
    rows_g = [k for k in range(M) if rgroup[k] == g]
    cols_g = [k for k in range(M) if cgroup[k] == g]
    rverts = [v for v in range(len(net["rports"])) if rvg[v] == g]
    cverts = [v for v in range(len(net["cports"])) if cvg[v] == g]
    for v in rverts:
        assert all(rgroup[k] == g for k in net["rports"][v])
    for v in cverts:
        assert all(cgroup[k] == g for k in net["cports"][v])
    nm = len(modes)
    wt = [modes[z]["w"] for z in range(nm)]
    packets = sorted({slot_of[k][0] for k in rows_g + cols_g} | {m["pid"] for m in modes})
    out = []
    for Yb in product((0, 1), repeat=len(rows_g)):
        Y = dict(zip(rows_g, Yb))
        Aval = Rp.one()
        for v in rverts:
            Aval = Rp.scal(net["rten"][v][tuple(Y[k] for k in net["rports"][v])], Aval)
        if not any(Aval):
            continue
        for Xb in product((0, 1), repeat=len(cols_g)):
            X = dict(zip(cols_g, Xb))
            base = Aval
            for v in cverts:
                base = Rp.scal(net["cten"][v][tuple(X[k] for k in net["cports"][v])], base)
            if not any(base):
                continue
            per_packet = []
            for pid in packets:
                poly = {(0,) * nm: Rp.one()}
                for k in rows_g:
                    if slot_of[k][0] != pid:
                        continue
                    loc = Rp.zero()
                    for k2, Wv in local[k]:
                        if X[k2]:
                            loc = Rp.add(loc, Wv)
                    fac = {(0,) * nm: loc if Y[k] else Rp.sub(Rp.one(), loc)}
                    for z, kap in row_modes[k]:
                        mono = tuple(1 if zz == z else 0 for zz in range(nm))
                        fac[mono] = Rp.add(fac.get(mono, Rp.zero()), kap if Y[k] else Rp.sub(Rp.zero(), kap))
                    poly = poly_mul(Rp, poly, fac, pprime - 1, wt)
                Lz = {z: sum(c for k2, c in modes[z]["terms"] if X.get(k2, 0))
                      for z in range(nm) if modes[z]["source"] == g and modes[z]["pid"] == pid}
                per_packet.append((pid, poly, Lz))
            out.append((base, per_packet))
    return out


def child_value(Rp, table, beta, modes, g):
    """Z_{g,beta}: the child's contraction, from its precomputed table."""
    nm = len(modes)
    total = Rp.zero()
    for base, per_packet in table:
        val = base
        for pid, poly, Lz in per_packet:
            want = tuple(beta[z] if (modes[z]["target"] == g and modes[z]["pid"] == pid) else 0
                         for z in range(nm))
            coef = poly.get(want)
            if coef is None:
                val = None
                break
            for z, lz in Lz.items():
                if beta[z]:
                    coef = Rp.scal(pow(lz, beta[z], Rp.m), coef)
            val = Rp.mul(val, coef)
            if not any(val):
                val = None
                break
        if val is not None:
            total = Rp.add(total, val)
    return total


def row_form_check(Kf, Rp, net, slot_of, rlab, clab, a, t, T, bmax, rng, samples=40):
    """(D.12): for every row port, (W(T) X)_x = local part + sum over crossing tiles and b < bmax of
    kappa_{x,y0,b} L_{y0,b}(X), modulo 2^P, on sampled X."""
    L, S, Cp0 = cauchy(Kf, a)
    i = 0
    while t % (2 ** (i + 1)) == 0:
        i += 1
    M = len(net["edges"])
    at = {v: k for k, v in slot_of.items()}
    rgroup = [rlab[u][i] for (u, v) in net["edges"]]
    cgroup = [clab[v][i] for (u, v) in net["edges"]]
    for _ in range(samples):
        X = [rng.randrange(2) for _ in range(M)]
        for k in range(M):
            pid, x = slot_of[k]
            exact = Rp.zero()
            expd = Rp.zero()
            for y in range(L):
                k2 = at.get((pid, y))
                if k2 is None or not X[k2]:
                    continue
                Wv = Rp.from_K(W_entry(Kf, S, Cp0, T, y - x))
                exact = Rp.add(exact, Wv)
                if cgroup[k2] == rgroup[k]:
                    expd = Rp.add(expd, Wv)
            for y0 in range(2 ** (i + 1)):
                real = [(at[(pid, y)], j) for j, y in enumerate(range(y0, L, 2 ** (i + 1))) if (pid, y) in at]
                if not real or cgroup[real[0][0]] == rgroup[k]:
                    continue
                Tm = Kf.sub(T, Kf.const(y0 - x))
                base = W_entry(Kf, S, Cp0, T, y0 - x)
                ratio = Kf.scal(2 ** (i + 1), Kf.inv(Tm))
                powr = Kf.const(1)
                for b in range(bmax):
                    kap = Rp.from_K(Kf.mul(base, powr))
                    Lz = sum(j ** b for k2, j in real if X[k2])
                    expd = Rp.add(expd, Rp.scal(Lz, kap))
                    powr = Kf.mul(powr, ratio)
            if exact != expd:
                return False
    return True


def part_C():
    rng = random.Random(1010)
    f = irreducible_mod2(3)
    Kf = K(f)
    ok, rows = True, []
    cases = [(1, 5, 3), (1, 6, 4), (2, 5, 3), (2, 6, 4), (2, 6, 5)]
    for case_id, (a, M, pprime) in enumerate(cases):
        P = pprime
        Rp = R(f, P)
        net = random_network(rng, 3, 3, M, 3, 1 << P)
        slot_of, packets, rlab, clab = pack(rng, net, a)
        L = 1 << a
        ts = [t for t in range(-L, L) if t != 0 and (t & -t).bit_length() - 1 < a]
        for t in ts:
            for zeta in [Kf.const(0), (Fraction(1), Fraction(1), Fraction(0)), (Fraction(0), Fraction(1), Fraction(1))]:
                T = Kf.add(Kf.const(t), Kf.scal(2 ** (a + 1), zeta))
                direct = deformed_contraction(Kf, Rp, net, slot_of, a, T)
                via, nb, nm = modes_and_children(Kf, Rp, net, slot_of, rlab, clab, a, t, T, pprime)
                good = direct == via
                ok &= good
                rows.append({"case": case_id, "a": a, "M": M, "p'": pprime, "t": t,
                             "zeta": [str(c) for c in zeta], "modes": nm, "beta_tuples": nb, "agree": good})
    # negative controls on one nontrivial case
    a, M, pprime = 2, 6, 4
    Rp = R(f, pprime)
    net = random_network(random.Random(77), 3, 3, M, 3, 1 << pprime)
    slot_of, packets, rlab, clab = pack(random.Random(78), net, a)
    neg = {}
    for label, kw in (("drop the crossing terms (all modes)", {"drop_modes": True}),
                      ("credit modes to the source group", {"swap_roles": True})):
        detects = False
        for t in (1, 3, -1, 2, -2):
            for zeta in [(Fraction(1), Fraction(0), Fraction(1)), (Fraction(1), Fraction(1), Fraction(1))]:
                T = Kf.add(Kf.const(t), Kf.scal(2 ** (a + 1), zeta))
                direct = deformed_contraction(Kf, Rp, net, slot_of, a, T)
                via, _, _ = modes_and_children(Kf, Rp, net, slot_of, rlab, clab, a, t, T, pprime, **kw)
                if direct != via:
                    detects = True
        neg[label] = detects
        ok &= detects
    # (D.12) at high precision, where the modes of weight >= 2 matter
    f5 = irreducible_mod2(3)
    K5 = K(f5)
    Rbig = R(f5, 40)
    rng2 = random.Random(31)
    d12 = []
    for a_ in (2, 3):
        net = random_network(rng2, 3, 3, 8, 4, 1 << 8)
        slot_of, packets, rlab, clab = pack(rng2, net, a_)
        L_ = 1 << a_
        for t in [t for t in range(-L_, L_) if t != 0 and (t & -t).bit_length() - 1 < a_]:
            zeta = (Fraction(1), Fraction(0), Fraction(1))
            T = K5.add(K5.const(t), K5.scal(2 ** (a_ + 1), zeta))
            full = row_form_check(K5, Rbig, net, slot_of, rlab, clab, a_, t, T, 39, rng2)
            trunc = row_form_check(K5, Rbig, net, slot_of, rlab, clab, a_, t, T, 1, rng2)
            d12.append({"a": a_, "t": t, "expansion_exact_mod_2^40": full, "b=0_only_exact": trunc})
    d12_ok = all(r["expansion_exact_mod_2^40"] for r in d12) and any(not r["b=0_only_exact"] for r in d12)
    ok &= d12_ok
    RESULTS["child_identity"] = {"all_ok": ok, "cases": rows, "negative_controls_detect": neg,
                                 "tile_expansion_D12": d12}
    print(f"[C] child identity (D.16), {len(rows)} (network, t, zeta) cases: "
          f"{all(r['agree'] for r in rows)}; negative controls detect: {neg}; tile expansion (D.12) "
          f"mod 2^40 at {len(d12)} nodes: {all(r['expansion_exact_mod_2^40'] for r in d12)}, "
          f"b = 0 alone fails at {sum(1 for r in d12 if not r['b=0_only_exact'])}")
    return ok


# =============================================================== Part F: one honest level

def newton_interpolate(Kf, Rp, pts, vals):
    """Coefficients (in z) of the polynomial through (pts[k], vals[k]) over Rp; the pairwise
    differences of the points are units of O, so every division is exact."""
    n = len(pts)
    Pz = [Rp.from_K(x) for x in pts]
    coef = list(vals)
    for j in range(1, n):
        for k in range(n - 1, j - 1, -1):
            dinv = Rp.from_K(Kf.inv(Kf.sub(pts[k], pts[k - j])))
            coef[k] = Rp.mul(Rp.sub(coef[k], coef[k - 1]), dinv)
    poly = [coef[n - 1]]
    for k in range(n - 2, -1, -1):
        new = [Rp.zero() for _ in range(len(poly) + 1)]
        for i2, c in enumerate(poly):
            new[i2 + 1] = Rp.add(new[i2 + 1], c)
            new[i2] = Rp.sub(new[i2], Rp.mul(c, Pz[k]))
        new[0] = Rp.add(new[0], coef[k])
        poly = new
    return poly


def v2int(t):
    v = 0
    while t % 2 == 0:
        t //= 2
        v += 1
    return v


def recover_one_level(Kf, f, net, slot_of, a, u, p, honest=True):
    """Tapered Hermite recovery of F(0) mod 2^p from jets of F(t + alpha z) mod 2^(p_t)."""
    L, S, Cp0 = cauchy(Kf, a)
    M = len(net["edges"])
    alpha = 2 ** (a + 1)
    h = -(-3 * M // u)
    mult = {t: (h * max(0, u - v2int(t)) if t != 0 else 0) for t in S}
    assert sum(mult.values()) > M * (2 * L - 1)
    deg = M * (2 * L - 1)
    pts = binary_points(Kf, deg + 1)
    total = [Fraction(0)] * Kf.e
    for t in S:
        qt = mult[t]
        if qt == 0:
            continue
        Dt = L // 2 ** v2int(t)
        pt = p + 8 * h * Dt + 8 * M * (a + 1) if honest else p
        Rp = R(f, pt)
        vals = [deformed_contraction(Kf, Rp, net, slot_of, a, Kf.add(Kf.const(t), Kf.scal(alpha, z)))
                for z in pts]
        Fz = newton_interpolate(Kf, Rp, pts, vals)[:qt]
        # E_t(t + alpha z) = prod_{q != 0, t} (1 + alpha z / (t - q))^{q_q}, truncated at z^qt
        Ez = [Fraction(1)] + [Fraction(0)] * (qt - 1)
        for q in S:
            if q in (0, t) or mult[q] == 0:
                continue
            c = Fraction(alpha, t - q)
            fac = [Fraction(comb(mult[q], k)) * c ** k for k in range(min(qt, mult[q] + 1))]
            fac += [Fraction(0)] * (qt - len(fac))
            Ez = [sum(Ez[i2] * fac[k - i2] for i2 in range(k + 1)) for k in range(qt)]
        Ez_mod = [Rp.from_K(Kf.const(x)) for x in Ez]
        # inverse power series of Ez (constant term 1)
        inv = [Rp.one()] + [Rp.zero() for _ in range(qt - 1)]
        for k in range(1, qt):
            acc = Rp.zero()
            for i2 in range(1, k + 1):
                acc = Rp.add(acc, Rp.mul(Ez_mod[i2], inv[k - i2]))
            inv[k] = Rp.sub(Rp.zero(), acc)
        Fz += [Rp.zero()] * (qt - len(Fz))
        Qt = [Rp.zero() for _ in range(qt)]
        for k in range(qt):
            acc = Rp.zero()
            for i2 in range(k + 1):
                acc = Rp.add(acc, Rp.mul(Fz[i2], inv[k - i2]))
            Qt[k] = acc
        Et0 = Fraction(1)
        for q in S:
            if q in (0, t) or mult[q] == 0:
                continue
            Et0 *= Fraction(-q, t - q) ** mult[q]
        x0 = Fraction(-t, alpha)
        for j in range(qt):
            w = Et0 * x0 ** j
            for c in range(Kf.e):
                total[c] += w * Qt[j][c]
    out = []
    for x in total:
        if x.denominator % 2 == 0:
            return None                       # not even 2-adically integral
        out.append((x.numerator * pow(x.denominator, -1, 1 << p)) % (1 << p))
    return tuple(out)


def part_F():
    rng = random.Random(2027)
    ok, rows = True, []
    for (a, u, M, p, e) in [(1, 1, 3, 6, 4), (1, 1, 3, 10, 4), (2, 2, 3, 6, 5), (2, 1, 3, 8, 5)]:
        f = irreducible_mod2(e)
        Kf = K(f)
        for rep in range(2):
            net = random_network(rng, 2, 2, M, 3, 1 << p)
            slot_of, packets, rlab, clab = pack(rng, net, a)
            truth = direct_contraction(net, 1 << p)
            got = recover_one_level(Kf, f, net, slot_of, a, u, p)
            bad = recover_one_level(Kf, f, net, slot_of, a, u, p, honest=False)
            want = (truth,) + (0,) * (e - 1)
            good = got == want
            ok &= good
            rows.append({"a": a, "u": u, "M": M, "p": p, "e": e, "truth": truth,
                         "recovered": list(got) if got else None, "exact": good,
                         "without_extra_precision": "fails" if bad != want else "succeeds"})
    neg_fail = sum(1 for r in rows if r["without_extra_precision"] == "fails")
    ok &= neg_fail >= len(rows) // 2
    RESULTS["one_level_recovery"] = {"all_ok": ok, "cases": rows,
                                     "negative_control_failures": f"{neg_fail}/{len(rows)}"}
    print(f"[F] one level at honest precision (Lagrange over O, tapered Hermite): "
          f"{all(r['exact'] for r in rows)}; without extra precision fails in {neg_fail}/{len(rows)}")
    return ok


# =============================================================== Part M: mode-tuple count

def count_tuples(c, p):
    """Number of tuples beta over c variables of each weight w = 1, ..., p-1 with
    sum w * beta < p (exact, by the generating function prod_w (1 - x^w)^(-c))."""
    coeff = [0] * p
    coeff[0] = 1
    for w in range(1, p):
        for _ in range(c):
            for sidx in range(w, p):
                coeff[sidx] += coeff[sidx - w]
    return sum(coeff)


def part_M():
    ok, rows = True, []
    for c in (1, 2, 5, 12, 40):
        for p in (5, 10, 20, 40, 80):
            exact = count_tuples(c, p)
            bound = 2 * math.pi * math.sqrt(c * p / 6)
            good = math.log(exact) <= bound
            ok &= good
            rows.append({"c": c, "p": p, "ln_count": round(math.log(exact), 3),
                         "bound_2pi_sqrt(cp/6)": round(bound, 3), "ok": good})
    RESULTS["mode_count"] = {"all_ok": ok, "cases": rows}
    print(f"[M] mode-tuple count <= exp(2 pi sqrt(c p / 6)) on {len(rows)} (c, p): {ok}")
    return ok


# =============================================================== Part L: Lemma D.2's encoding

def part_L():
    rng = random.Random(4242)
    ok, rows = True, []
    neg_detects = 0
    for trial in range(6):
        n, ntests, k = (4, 3, 2) if trial < 2 else (5, 4, 2)
        while True:
            tests = []
            for _ in range(ntests):
                monos = set()
                for _ in range(rng.randint(1, 3)):
                    monos.add(frozenset(rng.sample(range(n), rng.randint(1, k))))
                tests.append(sorted(monos, key=sorted))
            if all(tests):
                break
        ell = pmod = n + 1
        mod = 1 << pmod
        es = [tuple(e) for e in divide_by_h(p_ell(ell))][:ell]
        scopes = [sorted({v for m in T for v in m}) for T in tests]
        inc = [(i, v) for i in range(ntests) for v in scopes[i]]          # incidences
        K_deg = k * (ntests + 2 * (ell - 1))

        def Bvec(y):                                                       # B_i at an integer y, in B
            hy = y * y - y
            return [((a0 + a1 * y) * hy ** j) % mod for j, (a0, a1) in enumerate(es)]

        def vmul(x, y):
            out = [0] * ell
            for i1, c1 in enumerate(x):
                if c1:
                    for i2 in range(ell - i1):
                        out[i1 + i2] = (out[i1 + i2] + c1 * y[i2]) % mod
            return out

        # coefficient of the multilinear monomial z_{i,S} in B-bar_i (Moebius inversion)
        coefcache = {}

        def coef(i, Ssel):
            key = (i, Ssel)
            if key not in coefcache:
                tot = [0] * ell
                Sl = sorted(Ssel)
                for r_ in range(len(Sl) + 1):
                    for Tsub in __import__("itertools").combinations(Sl, r_):
                        y = sum(1 for m in tests[i] if m <= set(Tsub))
                        sign = (-1) ** (len(Sl) - r_)
                        tot = [(x + sign * b) % mod for x, b in zip(tot, Bvec(y))]
                coefcache[key] = tot
            return coefcache[key]

        # groups: J = 2 on each side
        J = 2
        tg = [i % J for i in range(ntests)]
        vg = [v % J for v in range(n)]
        pair_inc = {(g, hh): [c for c in inc if tg[c[0]] == g and vg[c[1]] == hh] for g in range(J) for hh in range(J)}
        total = [0] * ell
        total_neg = [0] * ell
        arrays = 0
        for darr in product(*[range(len(pair_inc[(g, hh)]) + 1) for g in range(J) for hh in range(J)]):
            if sum(darr) > K_deg:
                continue
            arrays += 1
            dmap = dict(zip([(g, hh) for g in range(J) for hh in range(J)], darr))
            subsets = {gh: list(__import__("itertools").combinations(pair_inc[gh], dmap[gh])) for gh in dmap}
            wires = {gh: (math.ceil(math.log2(len(subsets[gh]))) if len(subsets[gh]) > 1 else 0) for gh in dmap}
            # contract: sum over all wire bits; decode; invalid index -> entry 0
            names = [gh for gh in dmap for _ in range(wires[gh])]
            for bits in product((0, 1), repeat=len(names)):
                idx = {gh: 0 for gh in dmap}
                pos = {gh: 0 for gh in dmap}
                for gh, b in zip(names, bits):
                    idx[gh] |= b << pos[gh]
                    pos[gh] += 1
                if any(idx[gh] >= len(subsets[gh]) for gh in dmap):
                    continue
                chosen = set()
                for gh in dmap:
                    chosen |= set(subsets[gh][idx[gh]])
                # test-group entries: product over tests of the coefficients (in B)
                val = [1] + [0] * (ell - 1)
                for i in range(ntests):
                    Si = frozenset(v for (i2, v) in chosen if i2 == i)
                    val = vmul(val, coef(i, Si))
                # variable-group entries: 2^{# untouched variables}
                touched = {v for (_, v) in chosen}
                total_neg = [(x + y) % mod for x, y in zip(total_neg, val)]
                val = [(x * 2 ** (n - len(touched))) % mod for x in val]
                total = [(x + y) % mod for x, y in zip(total, val)]
        count = sum(total) % (1 << ell)
        count_neg = sum(total_neg) % (1 << ell)
        brute = sum(1 for x in product((0, 1), repeat=n)
                    if all(sum(1 for m in T if all(x[v] for v in m)) % 2 == 0 for T in tests))
        good = count == brute
        ok &= good
        neg_detects += count_neg != brute
        rows.append({"n": n, "tests": [[sorted(m) for m in T] for T in tests], "degree_arrays": arrays,
                     "count_from_networks": count, "brute": brute, "ok": good,
                     "without_variable_side_weights": count_neg})
    ok &= neg_detects >= len(rows) // 2
    RESULTS["sparse_system_encoding"] = {"all_ok": ok, "cases": rows,
                                         "negative_control_detects": f"{neg_detects}/{len(rows)}"}
    print(f"[L] Lemma D.2 encoding (incidence subsets on group wires) = common-zero count: "
          f"{all(r['ok'] for r in rows)}; dropping the variable-side weights detected in "
          f"{neg_detects}/{len(rows)}")
    return ok


# =============================================================== Part D: the clause network

def contract_generic(vertices, n_edges):
    """vertices: list of (ports, table) where table maps a tuple of port bits to an integer."""
    tot = 0
    for x in product((0, 1), repeat=n_edges):
        val = 1
        for ports, table in vertices:
            val *= table[tuple(x[e] for e in ports)]
            if not val:
                break
        tot += val
    return tot


def clause_network(n, clauses):
    """Paper D.2: clause tables, equality gadgets, subdivision; returns (vertices, edges, unused)."""
    edges = []                      # (endpoint vertex A, endpoint vertex B)
    verts = []                      # (kind, data)
    var_ports = {v: [] for v in range(n)}
    for C in clauses:
        vs = sorted({abs(l) - 1 for l in C})
        cid = len(verts)
        verts.append(("clause", (C, vs)))
        for v in vs:
            var_ports[v].append(cid)
    for v, cl in var_ports.items():
        d = len(cl)
        if d == 0:
            continue
        if d <= 2:
            gid = len(verts)
            verts.append(("eq", d))
            for c in cl:
                edges.append((c, gid, v))
        else:
            # a path-shaped tree of d-2 trivalent equality vertices with d external ports
            inner = []
            for _ in range(d - 2):
                inner.append(len(verts))
                verts.append(("eq", 3))
            edges.append((cl[0], inner[0], v))
            for j in range(d - 2):
                edges.append((cl[j + 1], inner[j], v))
                if j + 1 < d - 2:
                    edges.append((inner[j], inner[j + 1], v))
            edges.append((cl[d - 1], inner[-1], v))
    unused = sum(1 for v in range(n) if not var_ports[v])
    # subdivide every edge by a binary equality vertex
    ports = {vi: [] for vi in range(len(verts))}
    sub_vertices = []
    e_id = 0
    var_of_edge = {}
    for (A, B, v) in edges:
        mid = ("eq", 2)
        e1, e2 = e_id, e_id + 1
        e_id += 2
        ports[A].append((e1, v))
        ports[B].append((e2, v))
        sub_vertices.append(([e1, e2], {(0, 0): 1, (1, 1): 1, (0, 1): 0, (1, 0): 0}))
    out = []
    for vi, (kind, data) in enumerate(verts):
        pl = [e for e, _ in ports[vi]]
        if kind == "clause":
            C, vs = data
            table = {}
            for bits in product((0, 1), repeat=len(pl)):
                val = {}
                for (e, v), b in zip(ports[vi], bits):
                    val[v] = b
                table[bits] = int(any((val[abs(l) - 1] == 1) == (l > 0) for l in C))
            out.append((pl, table))
        else:
            table = {bits: int(len(set(bits)) <= 1) for bits in product((0, 1), repeat=len(pl))}
            out.append((pl, table))
    return out + sub_vertices, e_id, unused


def part_D():
    rng = random.Random(99)
    ok, rows = True, []
    max_ports = 0
    neg = 0
    fixed = [(3, [[1], [1, 2], [-1, 3], [1, -2]]),         # variable 1 has 4 ports (two inner vertices)
             (2, [[1], [-1], [1, 2]])]                      # unsatisfiable
    while len(rows) < 10:
        if fixed:
            n, clauses = fixed.pop(0)
        else:
            n = rng.randint(3, 5)
            m = rng.randint(2, 4)
            clauses = []
            for _ in range(m):
                w = rng.randint(1, 3)
                vs = rng.sample(range(1, n + 1), min(w, n))
                clauses.append([v if rng.random() < 0.5 else -v for v in vs])
        verts, ne, unused = clause_network(n, clauses)
        if ne > 18:
            continue
        got = contract_generic(verts, ne) * 2 ** unused
        brute = sum(1 for x in product((0, 1), repeat=n)
                    if all(any((x[abs(l) - 1] == 1) == (l > 0) for l in C) for C in clauses))
        good = got == brute
        ok &= good
        ports = max(sum(1 for C in clauses if any(abs(l) == v for l in C)) for v in range(1, n + 1))
        max_ports = max(max_ports, ports)
        # negative control: let every equality vertex accept all bit patterns
        broken = [(pl, {k: 1 for k in tb}) if (len(pl) >= 2 and all(tb.get((x,) * len(pl)) == 1 for x in (0, 1))
                                                and sum(tb.values()) == 2) else (pl, tb) for pl, tb in verts]
        got_neg = contract_generic(broken, ne) * 2 ** unused
        neg += got_neg != brute
        rows.append({"n": n, "clauses": clauses, "edges": ne, "max_ports_per_variable": ports,
                     "contraction_times_2^unused": got, "#F": brute, "ok": good,
                     "with_equalities_removed": got_neg})
    ok &= max_ports >= 4 and neg >= len(rows) // 2
    RESULTS["clause_network"] = {"all_ok": ok, "cases": rows,
                                 "negative_control_detects": f"{neg}/{len(rows)}"}
    print(f"[D] clause network contraction = #F on {len(rows)} CNFs (variables with up to "
          f"{max_ports} ports): {all(r['ok'] for r in rows)}; removing the equalities detected in "
          f"{neg}/{len(rows)}")
    return ok


# =============================================================== Part P: parameters

def recursion_conditions(lam, C, eps, r, use_paper_u_floor=False):
    """Check, in logarithmic form, every parameter inequality the recursion uses, at every active
    depth, for s in (2^(lam-1), 2^lam].  Returns (all_ok, details)."""
    lg = math.log2
    rl = lg(r)
    sq = math.isqrt(lam)
    logH = (lam - 1) - lg(sq)                    # H >= s / floor(sqrt(lam)) > 2^(lam-1) / sq
    logDelta = 0.9 * lam + 1                     # Delta = ceil(s^0.9) <= 2^(0.9 lam) + 1 <= 2^(0.9 lam + 1)
    # active depths: C s 8^-d > H, so 8^d < C sq
    dmax = 0
    while 8 ** (dmax + 1) < C * sq:
        dmax += 1
    u = [math.ceil(eps * lam / 2 ** d) for d in range(dmax + 1)]
    a = []
    acc = 0
    for d in range(dmax + 1):
        acc += 2 * u[d]
        a.append(acc)
    fails = []
    for d in range(dmax + 1):
        ud, ad = u[d], a[d]
        logDarr = 2 * ad * rl
        if not (1 <= ud <= ad):
            fails.append((d, "1 <= u <= a"))
        # (A1), (A2): precision absorption
        if ud < lg(20.5 * C * (d + 1) * math.sqrt(lam) * ud):
            fails.append((d, "(A1) 2^u >= 20.5 C (d+1) sqrt(lam) u"))
        if ud < lg(4 * (ad + 1) * ud):
            fails.append((d, "(A2) 2^u >= 4 (a+1) u"))
        if use_paper_u_floor and ud < lam ** 0.8:
            fails.append((d, "paper (D.5): u >= lam^(4/5)"))
        # (D.6)
        if logH < 1 + 2 * logDarr + logDelta:
            fails.append((d, "M >= 2 Darr^2 Delta"))
        if logH < ad + logDarr:
            fails.append((d, "M >= L Darr"))
        if ad + 1 > 0.9 * (lam - 1):
            fails.append((d, "2L <= Delta"))
        # interpolation points: M(2L-1)+1 <= 2^e, e = 4 lam
        if lg(C) + lam + ad + 1 > 4 * lam:
            fails.append((d, "points"))
        # packet-tensor entries: (L P)^L terms with P <= 41 C s (d+2) 2^a; need L log(L P) <= s^(9/10)
        logP = lg(41 * C * (d + 2)) + lam + ad
        if ad + lg(ad + logP) > 0.9 * (lam - 1):
            fails.append((d, "packet entry cost"))
    return (not fails), {"dmax": dmax, "u": u, "a": a, "fails": fails[:3]}


def least_lambda(C, eps, r, dmax_check=60):
    """Least lambda_0 such that all inequalities hold for every lambda >= lambda_0 (checked up to
    the interval in which depth dmax_check becomes active).  A new depth d becomes active when
    floor(sqrt(lambda)) first exceeds 8^d / C; inside an interval between two such transitions
    every condition only improves as lambda grows, so each interval is settled by checking its
    start and, if needed, binary-searching inside it."""
    starts = [64]
    for d in range(1, dmax_check + 1):
        q = 8 ** d // C + 1                      # least floor(sqrt(lambda)) with 8^d < C * floor(sqrt)
        starts.append(max(64, q * q))
    starts = sorted(set(starts))
    lam0 = 64
    for j, lo in enumerate(starts):
        hi = starts[j + 1] - 1 if j + 1 < len(starts) else lo * 64
        if recursion_conditions(lo, C, eps, r)[0]:
            continue
        if not recursion_conditions(hi, C, eps, r)[0]:
            lam0 = hi + 1
            continue
        a_, b_ = lo, hi                           # fails at a_, holds at b_
        while b_ - a_ > 1:
            mid = (a_ + b_) // 2
            if recursion_conditions(mid, C, eps, r)[0]:
                b_ = mid
            else:
                a_ = mid
        lam0 = max(lam0, b_)
    return lam0


def part_P():
    rows = []
    ok = True
    for C in (1, 4, 16):
        for label, eps in (("paper eps = 2^-17", 2.0 ** -17),):
            lam0 = least_lambda(C, eps, 64)
            # independent spot-check beyond the threshold on a geometric grid
            beyond = all(recursion_conditions(int(lam0 * 1.5 ** j), C, eps, 64)[0] for j in range(0, 400))
            ok &= beyond
            rows.append({"C": C, "constants": label, "least_lambda_all_inequalities": lam0,
                         "log2_of_that_lambda": round(math.log2(lam0), 2),
                         "grid_check_up_to_log2_lambda": round(math.log2(lam0 * 1.5 ** 399), 1),
                         "grid_ok": beyond})
    # the paper's own floor u >= lambda^(4/5) is far stronger than needed
    paper_floor = {}
    for C in (1,):
        lam = 2.0 ** 500
        while True:
            dmax = 0
            sq = math.sqrt(lam)
            while 8 ** (dmax + 1) < C * sq:
                dmax += 1
            if 2.0 ** -17 * lam / 2 ** dmax >= lam ** 0.8:
                break
            lam *= 1.05
        paper_floor[C] = round(math.log2(lam), 1)
    RESULTS["parameters"] = {"all_ok": ok, "thresholds": rows,
                             "log2_lambda_where_paper_floor_u>=lambda^(4/5)_starts_to_hold": paper_floor}
    print(f"[P] recursion inequalities (A1),(A2),(D.6), points: thresholds log2(lambda) = "
          f"{[r['log2_of_that_lambda'] for r in rows]} for C = 1, 4, 16; paper's u >= lambda^(4/5) needs "
          f"log2(lambda) ~ {paper_floor[1]}")
    return ok


def main():
    oks = [part_C(), part_F(), part_M(), part_L(), part_D(), part_P()]
    RESULTS["all_ok"] = all(oks)
    with open(__file__.rsplit("/", 2)[0] + "/results/tensor_contraction.json", "w") as fh:
        json.dump(RESULTS, fh, indent=1, default=str)
    print("ALL OK" if all(oks) else "FAILURES")


if __name__ == "__main__":
    main()

"""Exact reference implementation of the replacement-form construction (Sections 3-4).

This module rebuilds, from the definitions in the paper and nothing else, the objects
used by the bounded-occurrence counting algorithm:

  * the auxiliary label set B = ([r] x [r])^Q and the swap involutions Pi_s;
  * W0 = (1/2) sum_s (I - Pi_s), V0 = prod_s Pi_s, f_I and M0 = V0 f_I(W0)   (4.6)-(4.7);
  * the eigenbasis psi and the measurement matrix C^B of Lemma 4.1;
  * the auxiliary rows P_u(T) of (4.10);
  * the transfer to formula incidences of Section 4.5 (copies, iota, C_{t,g}).

All arithmetic is exact (fractions.Fraction) unless a prime modulus is requested.
Positions s in code are 0..Q-1 and stand for the evaluation points t = s+1.
"""

from fractions import Fraction
from itertools import product
from math import comb, floor
import random


# ---------------------------------------------------------------- linear algebra

def rref(rows, ncols):
    """Row-reduce a list of dict rows {col: Fraction}. Returns list of (pivot, row)."""
    basis = []  # list of (pivot_col, row dict normalised so row[pivot] == 1)
    for row in rows:
        r = {c: Fraction(v) for c, v in row.items() if v != 0}
        for piv, b in basis:
            if piv in r:
                coef = r[piv]
                for c, v in b.items():
                    nv = r.get(c, 0) - coef * v
                    if nv == 0:
                        r.pop(c, None)
                    else:
                        r[c] = nv
        if r:
            piv = min(r)
            inv = 1 / r[piv]
            r = {c: v * inv for c, v in r.items()}
            # keep basis fully reduced on the new pivot
            new_basis = []
            for p2, b in basis:
                if piv in b:
                    coef = b[piv]
                    b = dict(b)
                    for c, v in r.items():
                        nv = b.get(c, 0) - coef * v
                        if nv == 0:
                            b.pop(c, None)
                        else:
                            b[c] = nv
                new_basis.append((p2, b))
            basis = new_basis + [(piv, r)]
    return basis


def reduce_against(basis, row):
    r = {c: Fraction(v) for c, v in row.items() if v != 0}
    for piv, b in basis:
        if piv in r:
            coef = r[piv]
            for c, v in b.items():
                nv = r.get(c, 0) - coef * v
                if nv == 0:
                    r.pop(c, None)
                else:
                    r[c] = nv
    return r


def in_rowspace(basis, row):
    return not reduce_against(basis, row)


# ---------------------------------------------------------------- polynomials in T

def poly_mul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def poly_eval(p, t):
    acc = Fraction(0)
    for c in reversed(p):
        acc = acc * t + c
    return acc


def poly_shift_coeffs(p, t, J):
    """Coefficients of z^0..z^{J-1} in p(t+z)."""
    out = [Fraction(0)] * J
    for k, c in enumerate(p):
        if c:
            for j in range(min(k, J - 1) + 1):
                out[j] += c * comb(k, j) * Fraction(t) ** (k - j)
    return out


# ---------------------------------------------------------------- auxiliary array

class Aux:
    """Auxiliary construction of Section 4.1-4.3 for parameters r, Q, R and interval I."""

    def __init__(self, r, Q, R, A0=None):
        assert r >= 2 and Q >= 1 and R >= 1
        self.r, self.Q, self.R = r, Q, R
        pairs = [(a, d) for a in range(r) for d in range(r)]
        self.labels = list(product(pairs, repeat=Q))
        self.index = {u: i for i, u in enumerate(self.labels)}
        self.D = len(self.labels)
        assert self.D == r ** (2 * Q)
        if A0 is None:
            mu = Fraction(Q * (r - 1), 2 * r)
            A0 = floor(mu - Fraction(R - 1, 2))
            A0 = max(0, min(A0, Q - R + 1))
        self.A0 = A0
        self.interval = set(range(A0, A0 + R))
        assert self.interval <= set(range(Q + 1)), "interval must lie in {0..Q}"
        self._build_M0()
        self._build_CB()
        self._build_rows()

    # theta_s exchanges a_u(s) and d_u(s)
    def swap(self, u, s):
        u = list(u)
        a, d = u[s]
        u[s] = (d, a)
        return tuple(u)

    def swap_all(self, u):
        return tuple((d, a) for (a, d) in u)

    def E(self, u, y):
        """E_{u,y} = {s : a_y(s) != d_u(s)} as positions 0..Q-1."""
        return [s for s in range(self.Q) if y[s][0] != u[s][1]]

    def f_coeffs(self):
        """Monomial coefficients of f_I(z) = (-1)^A0 sum_{j<R} (-2)^j binom(z-A0, j)."""
        A0, R = self.A0, self.R
        total = [Fraction(0)]
        for j in range(R):
            # binom(z - A0, j) = prod_{a<j} (z - A0 - a) / j!
            p = [Fraction(1)]
            for a in range(j):
                p = poly_mul(p, [Fraction(-(A0 + a)), Fraction(1)])
            p = [c * Fraction((-2) ** j, _fact(j)) for c in p]
            if len(p) > len(total):
                total += [Fraction(0)] * (len(p) - len(total))
            for i, c in enumerate(p):
                total[i] += c
        sign = -1 if A0 % 2 else 1
        return [sign * c for c in total]

    def _build_M0(self):
        # W0 as sparse rows: W0 = (1/2) sum_s (I - Pi_s); (Pi_s zeta)_u = zeta_{theta_s(u)}
        Q = self.Q
        half = Fraction(1, 2)
        W = {}
        for u in self.labels:
            row = {u: Fraction(Q, 2)}
            for s in range(Q):
                w = self.swap(u, s)
                row[w] = row.get(w, 0) - half
            W[u] = {k: v for k, v in row.items() if v != 0}
        self.W0 = W
        coeffs = self.f_coeffs()
        # Horner: F = c_top I; F = F W + c I ...
        F = {u: {u: coeffs[-1]} for u in self.labels}
        for c in reversed(coeffs[:-1]):
            G = {}
            for u in self.labels:
                acc = {}
                for w, fv in F[u].items():
                    for y, wv in W[w].items():
                        acc[y] = acc.get(y, 0) + fv * wv
                if c:
                    acc[u] = acc.get(u, 0) + c
                G[u] = {k: v for k, v in acc.items() if v != 0}
            F = G
        self.fW = F
        # M0 = V0 f(W0): M0[u, y] = F[theta_all(u), y]
        self.M0 = {u: dict(F[self.swap_all(u)]) for u in self.labels}

    def _build_CB(self):
        r, Q = self.r, self.Q
        local = []  # (vector on [r]^2 as dict, eigenvalue)
        for a in range(r):
            local.append(({(a, a): 1}, 1))
        for a in range(r):
            for b in range(a + 1, r):
                local.append(({(a, b): 1, (b, a): 1}, 1))
                local.append(({(a, b): 1, (b, a): -1}, -1))
        assert len(local) == r * r
        self.local = local
        self.psi = []  # (row dict, j)
        for choice in product(range(len(local)), repeat=Q):
            vecs = [local[c] for c in choice]
            j = sum(1 for v in vecs if v[1] == -1)
            row = {(): 1}
            for vec, _ in vecs:
                new = {}
                for prefix, val in row.items():
                    for pair, x in vec.items():
                        new[prefix + (pair,)] = val * x
                row = new
            self.psi.append((row, j))
        self.CB = [row for row, j in self.psi if j not in self.interval]
        self.cB = len(self.CB)

    def _build_rows(self):
        """P_u(T) = delta_u (I - M0) + sum_{y: M0[u,y] != 0} M0[u,y] phi_{u,y}(T) delta_y."""
        self.P = {}
        for u in self.labels:
            row = {}
            # constant part delta_u (I - M0)
            row[u] = [Fraction(1)]
            for y, m in self.M0[u].items():
                cur = row.get(y, [Fraction(0)])
                row[y] = _padd(cur, [-m])
            for y, m in self.M0[u].items():
                phi = [Fraction(1)]
                for s in self.E(u, y):
                    phi = poly_mul(phi, [Fraction(1), Fraction(-1, s + 1)])
                row[y] = _padd(row.get(y, [Fraction(0)]), [m * c for c in phi])
            self.P[u] = {y: _trim(p) for y, p in row.items() if any(p)}

    def S_aux(self, s, g):
        """S^B_{t,g} = {u : a_u(t) = g} with t = s+1."""
        return [u for u in self.labels if u[s][0] == g]


def _fact(j):
    out = 1
    for i in range(2, j + 1):
        out *= i
    return out


def _padd(a, b):
    n = max(len(a), len(b))
    a = list(a) + [Fraction(0)] * (n - len(a))
    for i, c in enumerate(b):
        a[i] += c
    return a


def _trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def binom_tail_outside(Q, p, interval):
    """Exact Pr[Bin(Q,p) not in interval] as a Fraction."""
    return sum(comb(Q, j) * p ** j * (1 - p) ** (Q - j) for j in range(Q + 1) if j not in interval)


# ---------------------------------------------------------------- formulas

def random_bounded_cnf(n, rng, max_occ=3, max_width=3, n_clauses=None):
    """Random CNF with every variable in at most max_occ clauses; clauses as [(var, positive)]."""
    occ = [0] * n
    clauses = []
    if n_clauses is None:
        n_clauses = n
    attempts = 0
    while len(clauses) < n_clauses and attempts < 1000:
        attempts += 1
        free = [v for v in range(n) if occ[v] < max_occ]
        if not free:
            break
        w = rng.randint(1, min(max_width, len(free)))
        vs = rng.sample(free, w)
        cl = [(v, rng.random() < 0.5) for v in vs]
        for v in vs:
            occ[v] += 1
        clauses.append(cl)
    return clauses


def brute_count(n, clauses):
    cnt = 0
    for x in product((0, 1), repeat=n):
        if all(any((x[v] == 1) == pos for v, pos in cl) for cl in clauses):
            cnt += 1
    return cnt

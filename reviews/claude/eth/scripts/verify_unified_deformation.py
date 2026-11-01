"""Checks for 13_UNIFIED_DEFORMATION.md: one deformation theorem for Parts I-III.

The general lemmas are implemented once, independently of the instance code, and tested on
random readout systems; then the same generic code is run on genuine Part I and Part III data.

Part A  [Lemma 2.3, targeted modes] the node lemma in a p-adic filtration (p = 2, 3), with
        measured partitions: Z(rho) = sum_y sum_{sum w beta < P} prod_g Z_g(y, beta) mod p^P on
        random readout systems.  Negative controls: a coefficient one valuation too low; the
        truncation one weight too early; a mode credited to the wrong target group.
Part B  [Lemma 2.4, Boolean coordinate modes] the node lemma in the z-adic filtration of
        (Z/p^mu)[z]/(z^J), with weights 1 and 2 and measured partitions, on random readout
        systems.  Negative controls: truncation one weight too early; separation without the
        measurements.
Part C  [Lemma 3.1] Hermite recovery with valuation loss.  (C1) p-adic unit nodes (Part I's
        setting): zero loss and exact recovery mod p^mu; a node set with a non-unit difference
        loses precision.  (C2) the complete residue system mod 2^(a+1) with arbitrary depth
        profiles: the loss formula -v(E_t(0)) = D_i Delta(i) holds exactly; recovery of Cauchy-
        family polynomials succeeds at P_t = p + l_t + (m_t - 1)(a + 1 - i) for random lifts, and
        fails when one node gets one bit less.
Part D  [Lemma 4.1] correction counts: exact counts of {beta : sum w beta < P} (c modes per
        weight) against Rankin's upper bound and the lower bound (1/8) min(P - 1, sqrt(c(P-1)))
        in bits; Boolean counts C_q(n) against Rankin's bound.
Part E  [Lemmas 5.2-5.3, Theorem 5.4] the 2-adic architecture's floor: the profile lemma
        max Delta > 3M/(8a) on random profiles; the exact optimum of the depth-profile problem
        (a greedy solution of the linear program) and the paper's tapered profile, as multiples
        of M/a; the floor exponent.
Part F  [Corollaries 2.6, 2.7] instances: the generic Boolean node evaluator on a genuine Part I
        deformation (08) reproduces the Taylor jets; the generic targeted node evaluator on a
        genuine packed Cauchy-kernel network (10) reproduces F(t + alpha zeta) mod 2^P.

Standard library only; deterministic; writes ../results/unified_deformation.json.
"""

import json
import math
import random
import sys
import time
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product

sys.path.insert(0, __file__.rsplit("/", 1)[0])

RESULTS = {}


# =============================================================== rings with one interface

class ZN:
    """Z/N (N = p^P: the p-adic filtration a_i = p^i Z/N)."""

    def __init__(self, N):
        self.N = N

    def zero(self):
        return 0

    def one(self):
        return 1 % self.N

    def const(self, c):
        return c % self.N

    def add(self, a, b):
        return (a + b) % self.N

    def sub(self, a, b):
        return (a - b) % self.N

    def mul(self, a, b):
        return (a * b) % self.N

    def scal(self, k, a):
        return (k * a) % self.N


class Jet:
    """(Z/m)[z]/(z^J) (the z-adic filtration a_i = z^i)."""

    def __init__(self, m, J):
        self.m, self.J = m, J

    def zero(self):
        return (0,) * self.J

    def one(self):
        return (1 % self.m,) + (0,) * (self.J - 1)

    def const(self, c):
        return (c % self.m,) + (0,) * (self.J - 1)

    def add(self, a, b):
        return tuple((x + y) % self.m for x, y in zip(a, b))

    def sub(self, a, b):
        return tuple((x - y) % self.m for x, y in zip(a, b))

    def mul(self, a, b):
        out = [0] * self.J
        for i, x in enumerate(a):
            if x:
                for j in range(self.J - i):
                    out[i + j] = (out[i + j] + x * b[j]) % self.m
        return tuple(out)

    def scal(self, k, a):
        return tuple((k * x) % self.m for x in a)


class OWrap:
    """(O/2^P)[U]/(f) from verify_tensor_contraction.R, with the common interface."""

    def __init__(self, Rp):
        self.R = Rp

    def zero(self):
        return self.R.zero()

    def one(self):
        return self.R.one()

    def const(self, c):
        return self.R.const(c)

    def add(self, a, b):
        return self.R.add(a, b)

    def sub(self, a, b):
        return self.R.sub(a, b)

    def mul(self, a, b):
        return self.R.mul(a, b)

    def scal(self, k, a):
        return self.R.scal(k, a)


class PolyRing:
    """Polynomials in formal modes zeta_0..zeta_{k-1} over a base ring, truncated at weight P
    (a monomial zeta^e has weight sum w_i e_i)."""

    def __init__(self, base, weights, P):
        self.base, self.w, self.P = base, list(weights), P
        self.k = len(weights)
        self.zero_e = (0,) * self.k

    def zero(self):
        return {}

    def one(self):
        return {self.zero_e: self.base.one()}

    def const(self, c):
        return {self.zero_e: self.base.const(c)}

    def lift(self, a):
        return {self.zero_e: a}

    def var(self, idx, coef):
        e = [0] * self.k
        e[idx] = 1
        if self.w[idx] >= self.P:
            return {}
        return {tuple(e): coef}

    def add(self, a, b):
        out = dict(a)
        for e, c in b.items():
            out[e] = self.base.add(out[e], c) if e in out else c
        return out

    def sub(self, a, b):
        out = dict(a)
        for e, c in b.items():
            out[e] = self.base.sub(out[e], c) if e in out else self.base.sub(self.base.zero(), c)
        return out

    def mul(self, a, b):
        out = {}
        for e1, c1 in a.items():
            w1 = sum(x * y for x, y in zip(e1, self.w))
            for e2, c2 in b.items():
                w2 = sum(x * y for x, y in zip(e2, self.w))
                if w1 + w2 >= self.P:
                    continue
                e = tuple(x + y for x, y in zip(e1, e2))
                c = self.base.mul(c1, c2)
                out[e] = self.base.add(out[e], c) if e in out else c
        return out

    def scal(self, k, a):
        return {e: self.base.scal(k, c) for e, c in a.items()}


# =============================================================== readout systems
#
# system = {
#   "nvars": n, "groups": [list of variable indices], "C": {g: [row dict v -> int]},
#   "factors": [{"g": group, "ports": [port ids], "phi": function(ring, args)}],
#   "lam": {port: (const, {v: coef})}  (base-ring coefficients),
#   "modes": [{"src": g, "tgt": g or None, "w": weight, "form": {v: int}} or {"v": coordinate, "w": w}],
#   "kappa": {(port, mode index): base-ring element},
# }

def lam_value(R, lamspec, x):
    c0, form = lamspec
    out = c0
    for v, c in form.items():
        if x[v]:
            out = R.add(out, c)
    return out


def L_value(mode, x):
    if "v" in mode:
        return x[mode["v"]]
    return sum(c * x[v] for v, c in mode["form"].items())


def direct_Z(S, R):
    total = R.zero()
    kap_by_port = defaultdict(list)
    for (pid, Z), k in S["kappa"].items():
        kap_by_port[pid].append((Z, k))
    for x in product((0, 1), repeat=S["nvars"]):
        val = R.one()
        for f in S["factors"]:
            args = []
            for pid in f["ports"]:
                a = lam_value(R, S["lam"][pid], x)
                for Z, k in kap_by_port[pid]:
                    a = R.add(a, R.scal(L_value(S["modes"][Z], x), k))
                args.append(a)
            val = R.mul(val, f["phi"](R, args))
        total = R.add(total, val)
    return total


def fibres(S, ignore_measurements=False):
    out = []
    for g, Vg in enumerate(S["groups"]):
        d = defaultdict(list)
        rows = [] if ignore_measurements else S["C"].get(g, [])
        for bits in product((0, 1), repeat=len(Vg)):
            xg = dict(zip(Vg, bits))
            y = tuple(sum(c * xg[v] for v, c in row.items()) for row in rows)
            d[y].append(xg)
        out.append(d)
    return out


def assemble(n, parts):
    x = [0] * n
    for part in parts:
        for v, b in part.items():
            x[v] = b
    return x


def betas_below(weights, P):
    out = []

    def gen(i, budget, cur):
        if i == len(weights):
            out.append(tuple(cur))
            return
        c = 0
        while c * weights[i] < budget + 1:
            cur.append(c)
            gen(i + 1, budget - c * weights[i], cur)
            cur.pop()
            c += 1
    gen(0, P - 1, [])
    return out


def node_targeted(S, R, P, wrong_target=None):
    """Lemma 2.3: sum over measurement values y and mode tuples beta (sum w beta < P) of the
    product over groups of the group sums.  wrong_target = (mode, group) re-declares a target."""
    modes = [dict(m) for m in S["modes"]]
    if wrong_target is not None:
        modes[wrong_target[0]]["tgt"] = wrong_target[1]
    G = len(S["groups"])
    fib = fibres(S)
    weights = [m["w"] for m in modes]
    betas = betas_below(weights, P)
    total = R.zero()
    tgt_of = [[Z for Z in range(len(modes)) if modes[Z]["tgt"] == g] for g in range(G)]
    src_of = [[Z for Z in range(len(modes)) if modes[Z]["src"] == g] for g in range(G)]
    for ys in product(*[list(fib[g].keys()) for g in range(G)]):
        reps = [fib[g][ys[g]][0] for g in range(G)]
        tables = []
        for g in range(G):
            prg = PolyRing(R, [weights[Z] for Z in tgt_of[g]], P)
            pos = {Z: i for i, Z in enumerate(tgt_of[g])}
            entries = []
            for xg in fib[g][ys[g]]:
                x = assemble(S["nvars"], [xg if h == g else reps[h] for h in range(G)])
                poly = prg.one()
                for f in S["factors"]:
                    if f["g"] != g:
                        continue
                    args = []
                    for pid in f["ports"]:
                        a = prg.lift(lam_value(R, S["lam"][pid], x))
                        for Z in tgt_of[g]:
                            k = S["kappa"].get((pid, Z))
                            if k is not None:
                                a = prg.add(a, prg.var(pos[Z], k))
                        args.append(a)
                    poly = prg.mul(poly, f["phi"](prg, args))
                Lv = {Z: L_value(modes[Z], x) for Z in src_of[g]}
                entries.append((poly, Lv))
            tables.append(entries)
        for beta in betas:
            prod_g = R.one()
            for g in range(G):
                key = tuple(beta[Z] for Z in tgt_of[g])
                s = R.zero()
                for poly, Lv in tables[g]:
                    c = poly.get(key)
                    if c is None:
                        continue
                    mult = 1
                    for Z in src_of[g]:
                        mult *= Lv[Z] ** beta[Z]
                    s = R.add(s, R.scal(mult, c))
                prod_g = R.mul(prod_g, s)
            total = R.add(total, prod_g)
    return total, len(betas)


def node_boolean(S, R, P, ignore_measurements=False):
    """Lemma 2.4: sum over y, E (sum of weights < P) and U subset E of signed products of
    pinned group sums with the modes in U set to 1."""
    modes = S["modes"]
    G = len(S["groups"])
    fib = fibres(S, ignore_measurements)
    group_of = {v: g for g, Vg in enumerate(S["groups"]) for v in Vg}
    Es = []
    for size in range(len(modes) + 1):
        for E in combinations(range(len(modes)), size):
            if sum(modes[Z]["w"] for Z in E) < P:
                Es.append(E)
    total = R.zero()
    for ys in product(*[list(fib[g].keys()) for g in range(G)]):
        reps = [fib[g][ys[g]][0] for g in range(G)]
        for E in Es:
            pinned = {modes[Z]["v"] for Z in E}
            for usize in range(len(E) + 1):
                for U in combinations(E, usize):
                    prod_g = R.one()
                    for g in range(G):
                        s = R.zero()
                        for xg in fib[g][ys[g]]:
                            if any(xg[v] == 0 for v in pinned if group_of[v] == g):
                                continue
                            x = assemble(S["nvars"], [xg if h == g else reps[h] for h in range(G)])
                            val = R.one()
                            for f in S["factors"]:
                                if f["g"] != g:
                                    continue
                                args = []
                                for pid in f["ports"]:
                                    a = lam_value(R, S["lam"][pid], x)
                                    for Z in U:
                                        k = S["kappa"].get((pid, Z))
                                        if k is not None:
                                            a = R.add(a, k)
                                    args.append(a)
                                val = R.mul(val, f["phi"](R, args))
                            s = R.add(s, val)
                        prod_g = R.mul(prod_g, s)
                    if (len(E) - usize) % 2:
                        total = R.sub(total, prod_g)
                    else:
                        total = R.add(total, prod_g)
    return total, len(Es)


def multilinear_phi(coeffs):
    """phi(args) = sum_S c_S prod_{k in S} args[k], coefficients given as {tuple S: int}."""
    def phi(R, args):
        out = R.zero()
        for Sset, c in coeffs.items():
            term = R.const(c)
            for k in Sset:
                term = R.mul(term, args[k])
            out = R.add(out, term)
        return out
    return phi


# =============================================================== random readout systems

def random_system(rng, R, unit, scale_of_weight, n_groups=3, gsize=3, n_factors=6, n_modes=4,
                  boolean=False, measured=True, max_w=3):
    """Random readout system with local readouts (local modulo one measurement row per group)
    and modes.  unit() draws a base-ring element; scale_of_weight(w) is an element of a_w."""
    n = n_groups * gsize
    groups = [list(range(g * gsize, (g + 1) * gsize)) for g in range(n_groups)]
    C = {}
    if measured:
        for g in range(n_groups):
            row = {v: rng.choice((-1, 1, 2)) for v in groups[g] if rng.random() < 0.8}
            if row:
                C[g] = [row]
    modes = []
    for _ in range(n_modes):
        w = rng.randint(1, max_w)
        if boolean:
            v = rng.randrange(n)
            while any(m["v"] == v for m in modes):
                v = rng.randrange(n)
            modes.append({"v": v, "w": w, "src": v // gsize, "tgt": None})
        else:
            src = rng.randrange(n_groups)
            tgt = rng.choice([h for h in range(n_groups) if h != src])
            form = {v: rng.choice((-1, 1, 2, 3)) for v in groups[src] if rng.random() < 0.8}
            if not form:
                form = {groups[src][0]: 1}
            modes.append({"src": src, "tgt": tgt, "w": w, "form": form})
    factors, lam, kappa = [], {}, {}
    pid = 0
    for _ in range(n_factors):
        g = rng.randrange(n_groups)
        ports = []
        for _ in range(rng.randint(1, 3)):
            form = {}
            for v in groups[g]:
                if rng.random() < 0.7:
                    form[v] = R.add(form.get(v, R.zero()), unit())
            for h in range(n_groups):
                if h != g and h in C and rng.random() < 0.7:
                    b = unit()
                    for row in C[h]:
                        for v, c in row.items():
                            form[v] = R.add(form.get(v, R.zero()), R.scal(c, b))
            lam[pid] = (unit(), form)
            for Z, m in enumerate(modes):
                if (boolean or m["tgt"] == g) and rng.random() < 0.7:
                    kappa[(pid, Z)] = R.mul(scale_of_weight(m["w"]), unit())
            ports.append(pid)
            pid += 1
        k = len(ports)
        coeffs = {}
        for size in range(k + 1):
            for Sset in combinations(range(k), size):
                coeffs[Sset] = rng.randint(-3, 3)
        factors.append({"g": g, "ports": ports, "phi": multilinear_phi(coeffs)})
    return {"nvars": n, "groups": groups, "C": C, "factors": factors, "lam": lam, "modes": modes,
            "kappa": kappa}


def part_A():
    rng = random.Random(1301)
    rows, ok = [], True
    neg = {"low_valuation": [0, 0], "early_truncation": [0, 0], "wrong_target": [0, 0]}
    for p, P in ((2, 5), (2, 6), (3, 4)):
        R = ZN(p ** P)
        for rep in range(5):
            unit = lambda: rng.randrange(p ** P)  # noqa: E731
            S = random_system(rng, R, unit, lambda w: p ** w, n_groups=3, gsize=3, n_factors=6,
                              n_modes=4, boolean=False, measured=True)
            d = direct_Z(S, R)
            via, nb = node_targeted(S, R, P)
            good = d == via
            ok &= good
            # negative control 1: one coefficient one valuation too low
            S1 = dict(S)
            S1["kappa"] = dict(S["kappa"])
            key = next(iter(S["kappa"]))
            Zk = key[1]
            S1["kappa"][key] = (p ** (S["modes"][Zk]["w"] - 1) * (1 + p * rng.randrange(p ** P))) % p ** P
            n1 = direct_Z(S1, R) != node_targeted(S1, R, P)[0]
            # negative control 2: truncation one weight too early
            n2 = node_targeted(S, R, P - 1)[0] != d
            # negative control 3: a mode with nonzero coefficients credited to another target
            used = sorted({Z for (_, Z) in S["kappa"]})
            n3 = False
            if used:
                Z = used[0]
                other = [h for h in range(3) if h != S["modes"][Z]["tgt"]][0]
                n3 = node_targeted(S, R, P, wrong_target=(Z, other))[0] != d
            for k_, flag in (("low_valuation", n1), ("early_truncation", n2), ("wrong_target", n3)):
                neg[k_][0] += 1
                neg[k_][1] += int(flag)
            rows.append({"p": p, "P": P, "modes": len(S["modes"]), "beta_tuples": nb,
                         "agree": good, "value": d})
    # deterministic witnesses: each hypothesis, when violated, must change the answer
    wit = {}
    for p, P in ((2, 5), (3, 4)):
        R = ZN(p ** P)
        wit[f"p={p},P={P}"] = witness_targeted(R, p, P)
    wit_ok = all(all(v.values()) for v in wit.values())
    RESULTS["A_targeted_node_lemma"] = {"cases": rows, "all_agree": ok,
                                        "random_negative_controls_trials_detected": neg,
                                        "witness_negative_controls": wit, "ok": ok and wit_ok}
    print(f"[A] targeted node lemma (p-adic), {len(rows)} random systems: {ok}; witnesses fail as they "
          f"must: {wit_ok}; random negative controls detected {neg}")
    return ok and wit_ok


def witness_targeted(R, p, P):
    """Group 0 = {x0} is the source, group 1 is empty and is the target; one factor with k ports,
    each reading kappa * zeta (kappa = p, weight 1), Phi = product of the ports.  Then
    Z = p^k x0-sum = p^k.  With k = P-1 the weight P-1 term is the whole answer (early truncation
    loses it); with k = P and kappa = 1 (valuation 0 < weight 1) the true value is 1 but every
    term has weight P (low valuation is invisible to the truncation); crediting the mode to
    group 0 removes it from the factor."""
    def system(k, kap):
        lam = {j: (R.zero(), {}) for j in range(k)}
        kappa = {(j, 0): R.const(kap) for j in range(k)}

        def phi(Rg, args):
            out = Rg.one()
            for a_ in args:
                out = Rg.mul(out, a_)
            return out
        return {"nvars": 1, "groups": [[0], []], "C": {},
                "factors": [{"g": 1, "ports": list(range(k)), "phi": phi}], "lam": lam,
                "modes": [{"src": 0, "tgt": 1, "w": 1, "form": {0: 1}}], "kappa": kappa}
    S1 = system(P - 1, p)
    S2 = system(P, 1)
    d1 = direct_Z(S1, R)
    return {"identity_holds_on_witness": node_targeted(S1, R, P)[0] == d1 and d1 != 0,
            "early_truncation_fails": node_targeted(S1, R, P - 1)[0] != d1,
            "low_valuation_fails": node_targeted(S2, R, P)[0] != direct_Z(S2, R),
            "wrong_target_fails": node_targeted(S1, R, P, wrong_target=(0, 0))[0] != d1}


def part_B():
    rng = random.Random(1302)
    rows, ok = [], True
    neg = {"early_truncation": [0, 0], "no_measurements": [0, 0]}
    for (p, mu, J) in ((5, 3, 3), (5, 2, 4), (7, 2, 3)):
        m = p ** mu
        R = Jet(m, J)

        def unit():
            return tuple(rng.randrange(m) for _ in range(J))

        def zw(w):
            e = [0] * J
            if w < J:
                e[w] = 1
            return tuple(e)
        for rep in range(4):
            S = random_system(rng, R, unit, zw, n_groups=3, gsize=3, n_factors=5, n_modes=5,
                              boolean=True, measured=True, max_w=2)
            d = direct_Z(S, R)
            via, nE = node_boolean(S, R, J)
            good = d == via
            ok &= good
            n1 = node_boolean(S, R, J - 1)[0] != d
            n2 = node_boolean(S, R, J, ignore_measurements=True)[0] != d
            for k_, flag in (("early_truncation", n1), ("no_measurements", n2)):
                neg[k_][0] += 1
                neg[k_][1] += int(flag)
            rows.append({"p": p, "mu": mu, "J": J, "pinned_sets": nE, "agree": good})
    neg_ok = all(det >= tri // 2 for tri, det in neg.values())
    RESULTS["B_boolean_node_lemma"] = {"cases": rows, "all_agree": ok,
                                       "negative_controls_trials_detected": neg, "ok": ok and neg_ok}
    print(f"[B] Boolean node lemma (z-adic), {len(rows)} random systems: {ok}; negative controls "
          f"detected {neg}")
    return ok and neg_ok


# =============================================================== Part C: Hermite recovery

def vp(x, p):
    x = Fraction(x)
    if x == 0:
        return 10 ** 9
    n, d, v = x.numerator, x.denominator, 0
    while n % p == 0:
        n //= p
        v += 1
    while d % p == 0:
        d //= p
        v -= 1
    return v


def pmul(a, b, trunc=None):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if trunc is not None and i + j >= trunc:
                    break
                out[i + j] += x * y
    return out if trunc is None else (out + [Fraction(0)] * trunc)[:trunc]


def taylor_shift(F, t, alpha, trunc):
    """coefficients of F(t + alpha z) in z, below trunc."""
    out = [Fraction(0)] * trunc
    # Horner in (t + alpha z)
    for c in reversed(F):
        out = pmul(out, [Fraction(t), Fraction(alpha)], trunc)
        out[0] += c
    return out


def series_inv(s, trunc):
    out = [Fraction(0)] * trunc
    out[0] = 1 / s[0]
    for k in range(1, trunc):
        acc = sum(s[j] * out[k - j] for j in range(1, min(k, len(s) - 1) + 1))
        out[k] = -acc / s[0]
    return out


def hermite_recover(nodes, mult, alpha, data):
    """F(0) = sum_t E_t(0) sum_{j<m_t} [z^j](F(t + alpha z)/E_t(t + alpha z)) (-t/alpha)^j."""
    total = Fraction(0)
    info = {}
    for t in nodes:
        m = mult[t]
        if m == 0:
            continue
        E0 = Fraction(1)
        Es = [Fraction(1)]
        for q in nodes:
            if q == t or mult[q] == 0:
                continue
            E0 *= (Fraction(-q) / (t - q)) ** mult[q]
            fac = [Fraction(1), Fraction(alpha[t]) / (t - q)]
            for _ in range(mult[q]):
                Es = pmul(Es, fac, m)
        Es = (Es + [Fraction(0)] * m)[:m]
        inv = series_inv(Es, m)
        Pt = pmul(data[t], inv, m)
        val = sum(Pt[j] * (Fraction(-t) / alpha[t]) ** j for j in range(m))
        total += E0 * val
        info[t] = E0
    return total, info


def part_C():
    rng = random.Random(1303)
    out = {}
    ok = True
    # (C1) Part I's setting: nodes 1..Q, p > Q prime, any multiplicities: zero loss
    p, mu = 7, 4
    c1 = []
    for rep in range(6):
        Q = rng.randint(2, 5)
        nodes = list(range(1, Q + 1))
        mult = {t: rng.randint(1, 3) for t in nodes}
        deg = sum(mult.values()) - 1
        F = [Fraction(rng.randint(-50, 50)) for _ in range(deg + 1)]
        data = {t: [x for x in taylor_shift(F, t, 1, mult[t])] for t in nodes}
        # data known only modulo p^mu: perturb by random multiples of p^mu
        noisy = {t: [x + p ** mu * rng.randint(-9, 9) for x in data[t]] for t in nodes}
        rec, info = hermite_recover(nodes, mult, {t: 1 for t in nodes}, noisy)
        loss0 = all(vp(E0, p) == 0 for E0 in info.values())
        good = loss0 and vp(rec - F[0], p) >= mu
        c1.append({"Q": Q, "mult": mult, "zero_loss": loss0, "recovered_mod_p^mu": good})
        ok &= good
    # a non-unit difference (nodes 1 and 1 + p) loses precision
    nodes = [1, 1 + p, 3]
    mult = {1: 2, 1 + p: 2, 3: 1}
    losses = {}
    rec_fail = False
    for rep in range(10):
        F = [Fraction(rng.randint(-50, 50)) for _ in range(5)]
        data = {t: taylor_shift(F, t, 1, mult[t]) for t in nodes}
        noisy = {t: [x + p ** mu * rng.randint(-9, 9) for x in data[t]] for t in nodes}
        rec, info = hermite_recover(nodes, mult, {t: 1 for t in nodes}, noisy)
        losses = {t: -vp(E0, p) for t, E0 in info.items()}
        rec_fail |= vp(rec - F[0], p) < mu
    out["C1_unit_nodes"] = {"cases": c1, "non_unit_losses": losses, "non_unit_recovery_fails": rec_fail}
    ok &= rec_fail and max(losses.values()) > 0
    # (C2) the complete residue system mod 2^(a+1), arbitrary depth profiles
    c2 = []
    neg_detected, neg_trials = 0, 0
    for (a, M, profile) in ((2, 3, [6, 3]), (2, 3, [5, 1]), (2, 2, [3, 3]), (3, 2, [3, 2, 1]),
                            (3, 2, [2, 3, 4]), (3, 2, [5, 0, 1])):
        L = 2 ** a
        S = list(range(-L, L))
        nodes = [t for t in S if t != 0]
        mult = {t: (profile[vp(t, 2)] if vp(t, 2) < a else 0) for t in nodes}
        alpha = {t: 2 ** (a + 1) for t in nodes}
        # Cauchy family: F = prod_k sum_d c_{k,d} W_d(T)
        Cp0 = Fraction(1)
        for q in S:
            if q != 0:
                Cp0 *= -q
        W = {}
        for d in range(-L + 1, L):
            poly = [Fraction(1)]
            for q in S:
                if q != d:
                    poly = pmul(poly, [Fraction(-q), Fraction(1)])
            W[d] = [c / Cp0 for c in poly]
        F = [Fraction(1)]
        for k in range(M):
            G = [Fraction(0)] * (2 * L)
            for d in range(-L + 1, L):
                cd = rng.randint(-3, 3)
                G = [x + cd * y for x, y in zip(G, W[d])]
            F = pmul(F, G)
        degF = len(F) - 1
        assert degF <= M * (2 * L - 1)
        total_mult = sum(mult.values())
        if total_mult <= degF:
            continue
        data = {t: taylor_shift(F, t, alpha[t], mult[t]) for t in nodes if mult[t]}
        integral = all(vp(x, 2) >= 0 for v in data.values() for x in v) and vp(F[0], 2) >= 0
        rec_exact, info = hermite_recover(nodes, mult, alpha, data)
        exact = rec_exact == F[0]
        # loss formula
        formula = True
        for t, E0 in info.items():
            i = vp(t, 2)
            D = 2 ** (a - i)
            prof = profile + [0]
            Delta = sum(Fraction(j, 2 ** j) * (prof[i] - prof[i + j]) for j in range(1, a - i + 1))
            formula &= (-vp(E0, 2) == D * Delta)
        # precision rule and its sharpness
        pout = 6
        Pt = {t: pout - vp(info[t], 2) + (mult[t] - 1) * (a + 1 - vp(t, 2)) for t in info}
        good_lifts = True
        for rep in range(4):
            noisy = {t: [x + 2 ** max(Pt[t], 0) * rng.randint(-5, 5) for x in data[t]] for t in data}
            rec, _ = hermite_recover(nodes, mult, alpha, noisy)
            good_lifts &= vp(rec - F[0], 2) >= pout
        # one node with one bit less, perturbed at its top coefficient
        tstar = max(info, key=lambda t: Pt[t])
        noisy = {t: list(v) for t, v in data.items()}
        noisy[tstar][mult[tstar] - 1] += Fraction(2) ** (Pt[tstar] - 1)
        rec, _ = hermite_recover(nodes, mult, alpha, noisy)
        neg_trials += 1
        neg_detected += int(vp(rec - F[0], 2) < pout)
        row = {"a": a, "M": M, "profile": profile, "deg_F": degF, "sum_mult": total_mult,
               "jets_integral": integral, "exact_recovery": exact, "loss_formula": formula,
               "recovery_at_rule_precision": good_lifts,
               "losses": {str(t): -vp(E0, 2) for t, E0 in info.items()}}
        c2.append(row)
        ok &= integral and exact and formula and good_lifts
    out["C2_residue_system"] = {"cases": c2, "negative_one_bit_less": [neg_trials, neg_detected]}
    ok &= neg_detected == neg_trials and len(c2) >= 5
    out["ok"] = ok
    RESULTS["C_hermite_recovery"] = out
    print(f"[C] Hermite recovery: unit nodes zero loss and exact ({all(r['recovered_mod_p^mu'] for r in c1)}); "
          f"non-unit node set loses {losses}; residue system: {len(c2)} profiles, loss formula, exact and "
          f"rule-precision recovery: {all(r['loss_formula'] and r['exact_recovery'] and r['recovery_at_rule_precision'] for r in c2)}; "
          f"one bit less fails {neg_detected}/{neg_trials}")
    return ok


# =============================================================== Part D: correction counts

def count_tuples(c, P):
    a = [0] * P
    a[0] = 1
    for w in range(1, P):
        for _ in range(c):
            for N in range(w, P):
                a[N] += a[N - w]
    return sum(a)


def part_D():
    worst_lb, worst_ub = math.inf, 0.0
    grid = []
    for c in (1, 2, 3, 5, 8, 13, 21, 34):
        for P in (2, 3, 5, 8, 13, 21, 34, 55, 89, 144):
            N = count_tuples(c, P)
            lg = math.log2(N)
            lb = min(P - 1, math.sqrt(c * (P - 1))) / 8
            ub = math.pi * math.sqrt(2 * c * P / 3) / math.log(2)
            worst_lb = min(worst_lb, lg / lb)
            worst_ub = max(worst_ub, lg / ub)
            grid.append([c, P, round(lg, 3)])
    # Boolean counts C_q(n) = sum_{j<=q} binom(n,j) 2^j against Rankin's bound
    bool_ok = True
    for n in (10, 30, 100, 300):
        for q in (1, 3, 10, 30):
            if q > n:
                continue
            Cq = sum(math.comb(n, j) * 2 ** j for j in range(q + 1))
            best = min(math.exp(th * q) * (1 + 2 * math.exp(-th)) ** n for th in [k / 50 for k in range(1, 400)])
            bool_ok &= Cq <= best * (1 + 1e-12)
    ok = worst_lb >= 1 and worst_ub <= 1 and bool_ok
    RESULTS["D_correction_counts"] = {"min_log2N_over_lower_bound": worst_lb,
                                      "max_log2N_over_rankin": worst_ub, "boolean_rankin_ok": bool_ok,
                                      "ok": ok}
    print(f"[D] correction counts: log2 N / lower bound >= {worst_lb:.2f}, log2 N / Rankin <= {worst_ub:.3f}, "
          f"Boolean Rankin bound holds: {bool_ok}")
    return ok


# =============================================================== Part E: the 2-adic floor

def Delta_of(prof, a, i):
    pr = list(prof) + [0] * (a + 1 - len(prof))
    return sum(Fraction(j, 2 ** j) * (pr[i] - pr[i + j]) for j in range(1, a - i + 1))


def greedy_profile(a):
    """Componentwise-largest profile with Delta(i) <= 1 at every depth 0..a-1 (m(a) = 0)."""
    m = [Fraction(0)] * (a + 1)
    for i in range(a - 1, -1, -1):
        K = sum(Fraction(j, 2 ** j) for j in range(1, a - i + 1))
        acc = sum(Fraction(j, 2 ** j) * m[i + j] for j in range(1, a - i))
        m[i] = (1 + acc) / K
    return m


def part_E():
    rng = random.Random(1305)
    # profile lemma on random profiles
    worst = math.inf
    for a in range(1, 14):
        L = 2 ** a
        for _ in range(1500):
            M = rng.randint(1, 200)
            shape = rng.choice(("rand", "taper", "const", "incr", "spike"))
            if shape == "rand":
                m = [rng.randint(0, 4 * M) for _ in range(a)]
            elif shape == "taper":
                u = rng.randint(1, a)
                h = rng.randint(1, 4 * M)
                m = [h * max(0, u - i) for i in range(a)]
            elif shape == "const":
                m = [rng.randint(1, 4 * M)] * a
            elif shape == "incr":
                m = [rng.randint(0, M) * (i + 1) for i in range(a)]
            else:
                m = [0] * a
                m[rng.randrange(a)] = rng.randint(1, 8 * M * L)
            if not any(m):
                m[rng.randrange(a)] = 1
            tot = sum(2 ** (a - i) * m[i] for i in range(a))
            if tot <= M * (2 * L - 1):
                k = (M * (2 * L - 1)) // tot + 1
                m = [x * k for x in m]
            D = max(Delta_of(m, a, i) for i in range(a) if m[i] > 0)
            worst = min(worst, float(D) / (3 * M / (8 * a)))
    # exact optimum of the depth-profile problem, and the paper's tapered profile, in units of M/a
    opt_rows = []
    for a in (1, 2, 4, 8, 16, 32, 64, 128, 256):
        L = 2 ** a
        g = greedy_profile(a)
        Ssum = sum(2 ** (a - i) * g[i] for i in range(a))
        # minimal max Delta for the constraint sum > M(2L-1), per unit M (large M): (2L-1)/Ssum
        opt = Fraction(2 * L - 1) / Ssum
        # tapered profile h(u - i), u = a, smallest real h with the constraint
        tap = [Fraction(a - i) for i in range(a)]
        Tsum = sum(2 ** (a - i) * tap[i] for i in range(a))
        h = Fraction(2 * L - 1) / Tsum
        tapD = max(Delta_of([h * x for x in tap], a, i) for i in range(a))
        opt_rows.append({"a": a, "optimal_maxDelta_times_a_over_M": float(opt * a),
                         "tapered_maxDelta_times_a_over_M": float(tapD * a),
                         "greedy_profile_head": [float(x) for x in g[:4]],
                         "lemma_bound_times_a_over_M": 3 / 8})
    ok = worst > 1 and all(r["optimal_maxDelta_times_a_over_M"] >= 3 / 8 for r in opt_rows)
    # the floor exponent in bits: (1/8) min(p - 1, sqrt(3/16) M / sqrt(a)) at top level (p >= M)
    floor = [{"log2_M": lm, "a_max": lm, "floor_bits_over_M": min(1.0, math.sqrt(3 / 16) / math.sqrt(lm)) / 8}
             for lm in (10, 20, 40, 80, 160)]
    RESULTS["E_two_adic_floor"] = {"profile_lemma_min_ratio": worst, "optimum": opt_rows, "floor": floor,
                                   "ok": ok}
    print(f"[E] profile lemma: max Delta / (3M/8a) >= {worst:.2f} on random profiles; exact optimum "
          f"a*maxDelta/M: {[round(r['optimal_maxDelta_times_a_over_M'], 3) for r in opt_rows]}; tapered: "
          f"{[round(r['tapered_maxDelta_times_a_over_M'], 3) for r in opt_rows]}")
    return ok


# =============================================================== Part F: the two instances

def part_F():
    out = {}
    ok = True
    # ---------------- Part I (08): Boolean coordinate modes, z-adic filtration, measurements
    import verify_gapless_chain as I
    rng = random.Random(1306)
    f1 = []
    for (r, Q, R_) in ((2, 2, 2), (2, 3, 2)):
        sw = I.Swap(r, Q, R_)
        n, clauses = I.random_bounded(5, rng, rng.randint(3, 5))
        inst = I.build(n, clauses, sw)
        Rp = max(len(pp) - 1 for row in inst["P"].values() for pp in row.values())
        q = (inst["M"] * Rp) // Q
        J = q + 1
        p = I.prime_above(Q)
        m = p ** (n + 1)
        ring_theirs = I.Ring(m, J)
        R = Jet(m, J)
        for s in range(Q):
            t = s + 1
            groups = [inst["groups"][(s, g)] for g in range(r)]
            Cm = {g: inst["C"][(s, g)] for g in range(r)}
            lam, kappa, factors = {}, {}, []
            modes = [{"v": v, "w": 1} for v in range(n)]
            pid = 0
            for fi, f in enumerate(inst["factors"]):
                ports = []
                for v in f:
                    row = inst["P"][(fi, v)]
                    form = {w: R.const(I.fmod(Fraction(I.poly_eval(pp, t)), m)) for w, pp in row.items()}
                    lam[pid] = (R.zero(), form)
                    for w, pp in row.items():
                        coeffs = I.shift_poly_mod(pp, t, J, m)
                        coeffs[0] = 0
                        if any(coeffs):
                            kappa[(pid, w)] = tuple(coeffs)
                    ports.append(pid)
                    pid += 1
                signs = inst["signs"][fi]

                def phi(Rg, args, signs=signs):
                    prod_ = Rg.one()
                    for sg, y in zip(signs, args):
                        prod_ = Rg.mul(prod_, Rg.sub(Rg.one(), y) if sg else y)
                    return Rg.sub(Rg.one(), prod_)
                factors.append({"g": inst["d"][fi][s], "ports": ports, "phi": phi})
            S = {"nvars": n, "groups": groups, "C": Cm, "factors": factors, "lam": lam, "modes": modes,
                 "kappa": kappa}
            via, nE = node_boolean(S, R, J)
            theirs = tuple(I.taylor_direct(n, inst, s, ring_theirs))
            good = via == theirs
            ok &= good
            f1.append({"r": r, "Q": Q, "R": R_, "n": n, "node": t, "J": J, "agree_with_08": good})
    out["F1_part_I"] = f1
    # ---------------- Part III (10): targeted tile modes, 2-adic filtration, no measurements
    import verify_tensor_contraction as T
    f = T.irreducible_mod2(3)
    Kf = T.K(f)
    f3 = []
    rng = random.Random(1307)
    for (a, M, P) in ((1, 4, 3), (2, 5, 3), (2, 5, 4)):
        Rp = T.R(f, P)
        R = OWrap(Rp)
        net = T.random_network(rng, 3, 3, M, 3, 1 << P)
        slot_of, packets, rlab, clab = T.pack(rng, net, a)
        L, Sset, Cp0 = T.cauchy(Kf, a)
        at = {v: k for k, v in slot_of.items()}
        for t in [t for t in range(-L, L) if t != 0 and (t & -t).bit_length() - 1 < a][:3]:
            i = (t & -t).bit_length() - 1
            zeta = (Fraction(1), Fraction(0), Fraction(1))
            Tv = Kf.add(Kf.const(t), Kf.scal(2 ** (a + 1), zeta))
            # variables: X_k = k, Y_k = M + k; groups by coordinate-i symbols
            rg = [rlab[u][i] for (u, v) in net["edges"]]
            cg = [clab[v][i] for (u, v) in net["edges"]]
            G = 2
            groups = [[k for k in range(M) if cg[k] == g] + [M + k for k in range(M) if rg[k] == g]
                      for g in range(G)]
            lam, kappa, factors, modes = {}, {}, [], []
            pid = 0

            def table_phi(table, nports):
                def phi(Rg, args):
                    outv = Rg.zero()
                    for bits, val in table.items():
                        if val % (1 << P) == 0:
                            continue
                        term = Rg.const(val)
                        for b_, y in zip(bits, args):
                            term = Rg.mul(term, y if b_ else Rg.sub(Rg.one(), y))
                        outv = Rg.add(outv, term)
                    return outv
                return phi
            for v, ports_ in enumerate(net["rports"]):
                ps = []
                for k in ports_:
                    lam[pid] = (R.zero(), {M + k: R.one()})
                    ps.append(pid)
                    pid += 1
                factors.append({"g": rlab[v][i], "ports": ps, "phi": table_phi(net["rten"][v], len(ps))})
            for v, ports_ in enumerate(net["cports"]):
                ps = []
                for k in ports_:
                    lam[pid] = (R.zero(), {k: R.one()})
                    ps.append(pid)
                    pid += 1
                factors.append({"g": clab[v][i], "ports": ps, "phi": table_phi(net["cten"][v], len(ps))})
            tile_mode = {}
            for k in range(M):
                pk, x = slot_of[k]
                p1, p2 = pid, pid + 1
                pid += 2
                lam[p1] = (R.zero(), {M + k: R.one()})
                form = {}
                for y in range(L):
                    k2 = at.get((pk, y))
                    if k2 is not None and cg[k2] == rg[k]:
                        form[k2] = Rp.from_K(T.W_entry(Kf, Sset, Cp0, Tv, y - x))
                lam[p2] = (R.zero(), form)
                for y0 in range(2 ** (i + 1)):
                    real = [(at[(pk, y)], j) for j, y in enumerate(range(y0, L, 2 ** (i + 1))) if (pk, y) in at]
                    if not real or cg[real[0][0]] == rg[k]:
                        continue
                    base = T.W_entry(Kf, Sset, Cp0, Tv, y0 - x)
                    ratio = Kf.scal(2 ** (i + 1), Kf.inv(Kf.sub(Tv, Kf.const(y0 - x))))
                    powr = Kf.const(1)
                    for b in range(P - 1):
                        key = (pk, y0, b)
                        if key not in tile_mode:
                            tile_mode[key] = len(modes)
                            modes.append({"src": cg[real[0][0]], "tgt": rg[k], "w": b + 1,
                                          "form": {k2: j ** b for k2, j in real}})
                        kap = Kf.mul(base, powr)
                        assert Kf.v2(kap) >= b + 1
                        kappa[(p2, tile_mode[key])] = Rp.from_K(kap)
                        powr = Kf.mul(powr, ratio)
                factors.append({"g": rg[k], "ports": [p1, p2],
                                "phi": lambda Rg, args: Rg.add(Rg.mul(Rg.sub(Rg.one(), args[0]),
                                                                       Rg.sub(Rg.one(), args[1])),
                                                                Rg.mul(args[0], args[1]))})
            Ssys = {"nvars": 2 * M, "groups": groups, "C": {}, "factors": factors, "lam": lam,
                    "modes": modes, "kappa": kappa}
            via, nb = node_targeted(Ssys, R, P)
            direct = T.deformed_contraction(Kf, Rp, net, slot_of, a, Tv)
            good = tuple(via) == tuple(direct)
            ok &= good
            f3.append({"a": a, "M": M, "P": P, "t": t, "modes": len(modes), "beta_tuples": nb,
                       "agree_with_10": good})
    out["F2_part_III"] = f3
    out["ok"] = ok
    RESULTS["F_instances"] = out
    print(f"[F] instances: Part I jets reproduced at {sum(r['agree_with_08'] for r in f1)}/{len(f1)} nodes; "
          f"Part III F(t + alpha zeta) reproduced at {sum(r['agree_with_10'] for r in f3)}/{len(f3)} nodes")
    return ok


def main():
    oks = []
    for name, fn in (("A", part_A), ("B", part_B), ("C", part_C), ("D", part_D), ("E", part_E),
                     ("F", part_F)):
        t0 = time.time()
        oks.append(fn())
        RESULTS.setdefault("seconds", {})[name] = round(time.time() - t0, 1)
    RESULTS["all_ok"] = all(oks)
    with open(__file__.rsplit("/", 1)[0] + "/../results/unified_deformation.json", "w") as fh:
        json.dump(RESULTS, fh, indent=1, default=str)
    print("ALL OK" if all(oks) else "FAILURES")


if __name__ == "__main__":
    main()

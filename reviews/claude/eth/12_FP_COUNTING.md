# 12. A new theorem: counting zeros of sparse polynomial systems over any finite field

Claude, 28 October 2026. The ETH paper counts common zeros of sparse systems over $\mathbb F_2$ (its Theorem 2.2; Part II, [09](09_GAPLESS_SHARP_EXPONENT.md) Theorem 3.1). This file proves the same bound over every fixed finite field. It completes the blueprint left open in [07](07_TOOLKIT.md) E2 and [`../SYNTHESIS.md`](../SYNTHESIS.md) §5.2.

> **Theorem 12.1.** Fix a prime $p$ and integers $k,c\ge1$. Let $P_1,\dots,P_m\in\mathbb F_p[x_1,\dots,x_n]$ with $m\le cn$, each given as a list of at most $(\log n)^c$ terms $a_Mx^M$ of total degree $\le k$. Then
> $$N=\#\{x\in\mathbb F_p^n:\ P_1(x)=\dots=P_m(x)=0\}$$
> can be computed deterministically in $2^{O_{p,k,c}(\rho(n))}$ bit operations and $\mathrm{poly}(n)$ space, where $\rho(n)=n(\log\log n)^2/\sqrt{\log n}$.

> **Corollary 12.2 (every finite field).** The same holds over $\mathbb F_q$ for fixed $q=p^r$, with constants depending on $q$.

**What is needed beyond the $\mathbb F_2$ case.** Only two things:
- the $p$-adic lift with cross-equation truncation ([07](07_TOOLKIT.md) Theorem E2, restated as Lemma 1.1);
- one *validity factor* per variable, so that the bit encoding of $\mathbb F_p$ is exact.

Everything else is Part II's machinery, because deformation counting ([09](09_GAPLESS_SHARP_EXPONENT.md) Theorem 1.1) is indifferent to what the factors are.

**Checks.** [`scripts/verify_fp_counting.py`](scripts/verify_fp_counting.py) ([`results/fp_counting.json`](results/fp_counting.json), ~7 s):
- **The lift** (Lemma 1.1, Remark 1.2), for $p=2,3,5,7$ and $\ell\le8$:
  - the E2 digits satisfy the truncation identity;
  - the closed-form idempotent of Remark 1.2 satisfies the congruence, agrees with $p_\ell$ modulo $h^\ell$ (exact polynomial division), and for $p=2$ equals Part II's Lemma 2.1;
  - negative control: the unhalved binomial fails.
- **The full pipeline** on 8 systems over $\mathbb F_3$, $\mathbb F_5$ and $\mathbb F_7$, with up to 6 bit variables:
  - Taylor jets come from Part I's deformation by separation, and agree with direct evaluation;
  - the $V$-coefficients are recovered exactly by CRT, and the count equals brute force in every case;
  - negative controls: without the validity factors the count is wrong in 8/8 cases, and with too few Taylor coefficients in 4/4 applicable cases.

---

## 1. The $p$-adic lift

Write $h(Y)=Y^p-Y$. Fermat gives $p\mid h(y)$ for every integer $y$.

**Lemma 1.1 ([07](07_TOOLKIT.md) E2, with a norm bound).** Let
$$F_\ell(U)=(1-U)^\ell\sum_{j<\ell}\binom{\ell+j-1}jU^j,\qquad p_\ell(Y)=F_\ell(Y^{p-1}),$$
and divide repeatedly by the monic $h$: $p_\ell=\sum_je_j(Y)h(Y)^j$ with $\deg e_j<p$. For integers $y_1,\dots,y_m$ put $B_i(V)=\sum_{j<\ell}e_j(y_i)h(y_i)^jV^j$. Then:
1. $p_\ell(y)\equiv[p\mid y]\pmod{p^\ell}$ for every integer $y$.
2. $\mathcal S_\ell\big(\prod_iB_i\bmod V^\ell\big)\equiv\prod_i[p\mid y_i]\pmod{p^\ell}$, where $\mathcal S_\ell$ sums the coefficients of $1,V,\dots,V^{\ell-1}$.
3. If the $y_i$ are polynomials of degree $\le d_*$, each retained coefficient has degree $\le d_*\big(m(p-1)+p(\ell-1)\big)$.
4. $\|e_j\|_1\le(5p\ell)^\ell8^\ell$.

*Proof.*
1. **Two-sided identity.** $F_\ell=1-U^\ell\sum_{j<\ell}\binom{\ell+j-1}j(1-U)^j$, so $F_\ell\equiv1\bmod U^\ell$ and $F_\ell\equiv0\bmod(1-U)^\ell$.
   - If $p\mid y$, then $p^\ell\mid y^{(p-1)\ell}$, so $p_\ell(y)\equiv1$.
   - If $p\nmid y$, then $y^{p-1}\equiv1\bmod p$ (Fermat), so $p^\ell\mid(1-y^{p-1})^\ell$ and $p_\ell(y)\equiv0$.
2. **Truncation.** The sum of all coefficients of the full product is $\prod_ip_\ell(y_i)$. A dropped term has total $h$-power $\ge\ell$, so it is divisible by $p^\ell$.
3. **Degree.** The coefficient of $V^j$ in $B_i$ has degree $\le d_*(p-1+pj)$, and a retained product has $\sum j_i\le\ell-1$.
4. **Norms.**
   - $\|F_\ell\|_1\le2^\ell\binom{2\ell-1}{\ell-1}\le8^\ell$.
   - Dividing $g$ of degree $D$ by $h$: the quotient satisfies $q_j=g_{j+p}+q_{j+p-1}$, so $|q_j|\le\|g\|_1$, hence $\|q\|_1\le D\|g\|_1$ and the remainder has norm $\le(2D+1)\|g\|_1$.
   - With $D<2p\ell$ and $j<\ell$ divisions, $\|e_j\|_1\le(2D+1)D^j\|p_\ell\|_1\le(5p\ell)^\ell8^\ell$. ∎

**Remark 1.2 (a closed-form lift for every prime).** Let $z=Y^{p-1}$ and $\delta=z^2-z=Y^{p-2}h$. Then
$$\tilde E_\ell=(1-z)+\sum_{1\le j<\ell}(-1)^{j+1}\binom{2j-1}{j-1}(2z-1)\,Y^{(p-2)j}\,h^j\ \in\mathbb Z[Y]$$
has the following properties:
- $\tilde E_\ell(y)\equiv[p\mid y]\pmod{p^\ell}$;
- $\tilde E_\ell\equiv p_\ell\pmod{h^\ell}$;
- for $p=2$ it is exactly Part II's Lemma 2.1.

*Proof.* In $\mathbb Q[Y]/(h^\ell)$:
1. **Involution and idempotent.** $(2z-1)^2=1+4\delta$, with $\delta\in(h)$ nilpotent. So $\varepsilon=(2z-1)(1+4\delta)^{-1/2}$, using the binomial series truncated at $\delta^\ell$, satisfies $\varepsilon^2=1$. Expanding $(1-\varepsilon)/2$ gives $\tilde E_\ell$.
2. **Value near $0$.** Near $Y=0$, $\varepsilon\equiv-1$ (mod $Y$); in the local ring $\mathbb Q[Y]/(Y^\ell)$ this forces $\varepsilon=-1$.
3. **Value near the other roots.** Over $\mathbb C$, near each root $\zeta$ of $Y^{p-1}=1$, $\varepsilon\equiv1$. Hence $\tilde E_\ell\equiv1\bmod Y^\ell$ and $\tilde E_\ell\equiv0\bmod(Y^{p-1}-1)^\ell$.
4. **Integrality and conclusion.** Monic division keeps both congruences in $\mathbb Z[Y]$, and the parity argument of Lemma 1.1(1) applies. $p_\ell$ has the same two congruences, and $h^\ell=Y^\ell(Y^{p-1}-1)^\ell$ with coprime factors.

Its digits are not reduced: they have degree $(p-1)+(p-2)j$. The truncation degree becomes $d_*(m(p-1)+(2p-2)(\ell-1))$, still linear. ∎

So the hidden structure of the Beigel–Tarui lift is the same for every prime. It is the idempotent of $\mathbb Q[Y]/(h^\ell)$ belonging to the root $0$, obtained from the approximate idempotent $z=Y^{p-1}$ by the binomial square root $(1+4\delta)^{-1/2}$.

---

## 2. Proof of Theorem 12.1

Assume $n\ge n_0(p,k,c)$; otherwise use exhaustive search. Let $t=\lceil\log_2p\rceil$, $\lambda=\lceil\log n\rceil$ and $\ell=n+1$.

1. **Bit encoding and validity.** Each variable becomes $t$ bits, $x_v=\sum_{b<t}2^by_{v,b}$. Let $\phi_v(y_v)=[\sum_b2^by_{v,b}\le p-1]$.
   - $\phi_v$ is a multilinear integer polynomial of degree $\le t$: its Möbius coefficients are integers.
   - Valid bit vectors correspond bijectively to $\mathbb F_p^n$.
2. **Integer lifts.** Replace each coefficient by its representative in $[0,p)$ and each $x_v$ by $\sum_b2^bY_{v,b}$. This gives $L_i\in\mathbb Z[Y]$ of total degree $\le k$, with
   - $L_i(y)\equiv P_i(x)\pmod p$ on valid $y$;
   - $0\le L_i(y)\le y_{\max}:=(\log n)^c(p-1)(2^t-1)^k$ on all Boolean $y$.
3. **The counting identity.** Let $\mathcal B_i(Y,V)=\sum_{j<\ell}e_j(L_i)h(L_i)^jV^j$. By Lemma 1.1(2), summed over $y\in\{0,1\}^{tn}$,
$$N\equiv\mathcal S_\ell\Big(\sum_y\prod_v\phi_v(y_v)\prod_i\mathcal B_i(y,V)\bmod V^\ell\Big)\pmod{p^\ell}.$$
   Since $0\le N\le p^n<p^\ell$, $N$ is the least nonnegative residue.
4. **Branching on frequent variables.** This is Part II Theorem 3.1 step 2, applied to the $\mathbb F_p$-variables and the equation scopes:
   - The variables lying in more than $b_0=\lceil M/n\rceil\lceil\sqrt\lambda\rceil$ scopes form a set $S$ with $|S|<n/\sqrt\lambda$.
   - Enumerate $\mathbb F_p^S$, which gives $p^{|S|}\le2^{O_p(n/\sqrt\lambda)}$ branches. Substitute and combine terms mod $p$; this does not increase degrees or term counts.
   - Drop equations that become $0=0$; a nonzero constant equation makes the branch count $0$.
5. **Factors on a branch.** With $n'\ge n/2$ free variables, the bit variables number $n_b=tn'$. The factors are:
   - the equation factors $\mathcal B_i$, with scopes of at most $tk(\log n)^c$ bits;
   - the validity factors $\phi_v$, with scopes of $t$ bits.
   
   Each bit lies in at most $b_0+1$ scopes. The degree in $\mathscr A[Y]$ is, by Lemma 1.1(3),
$$d\le tn+k\big(m(p-1)+p(\ell-1)\big)\le\kappa_pn',\qquad\kappa_p=2\big(t+k(c(p-1)+p)\big).$$
6. **Deformation counting.** Run Part II Theorem 3.1 steps 3–5 with $n_b$ bit variables:
   - parameters $r,Q,R$ from $\lambda'=\lceil\log n_b\rceil$;
   - Part I's design for the bit-level scopes;
   - Theorem 1.1 of Part II, with coefficients in $\mathscr A=(\mathbb Z/P^\mu)[V]/(V^\ell)$, where $P$ is the least prime $>Q$; $P$ is unrelated to $p$.
7. **Exact recovery.** By Lemma 1.1(4),
$$\|\mathcal B_i(y,\cdot)\|_1\le U:=\ell(5p\ell)^\ell8^\ell(1+y_{\max})^{p-1+p(\ell-1)} .$$
   Since $\phi_v\in\{0,1\}$ on Boolean points, each $V$-coefficient $b_s$ of the integer sum has $|b_s|\le2^{n_b}U^m=:Z$, with $\log Z=O_{p,k,c}(n^2\log n)$. Choose $\mu$ with $P^\mu>2Z$; the symmetric residues give the $b_s$, and step 3 gives the branch count.
8. **Time and space.** This is Part II Theorem 3.1 step 7 with $n'$ replaced by $tn'$ and $\kappa$ by $\kappa_p$. Each branch costs $2^{O_{p,k,c}(n\nu^2/\sqrt\lambda)}$, and there are $2^{O_p(n/\sqrt\lambda)}$ branches; the total is $2^{O_{p,k,c}(\rho(n))}$. Branches are processed one at a time with an $O(n\log p)$-bit accumulator, so the space is polynomial. ∎

**Proof of Corollary 12.2.** Fix an $\mathbb F_p$-basis $\alpha_1,\dots,\alpha_r$ of $\mathbb F_q$ and write $x_v=\sum_sx_{v,s}\alpha_s$.
- A term of degree $\le k$ over $\mathbb F_q$ expands into at most $r^k$ terms over $\mathbb F_p$ in each of the $r$ coordinates.
- So $m$ sparse equations over $\mathbb F_q$ in $n$ variables become $rm$ equations over $\mathbb F_p$ in $rn$ variables. Each has at most $r^k(\log n)^c$ terms of degree $\le k$, and the solution sets correspond.
- Apply Theorem 12.1 with constants depending on $q$. ∎

---

## 3. Remarks

- **Where sparsity enters.** The hypothesis "at most $(\log n)^c$ terms per equation" keeps scope sizes and the per-factor cost polylogarithmic, as in the $\mathbb F_2$ case. Dense systems, such as generic quadratic systems with $\Theta(n^2)$ terms, are outside the theorem.
- **Bounded-arity CSPs.** Constraint satisfaction problems of bounded arity over a fixed domain need no new theorem: they reduce directly to #$k$-SAT. Theorem 12.1 matters when each equation is *sparse but has unbounded arity* (polylogarithmically many variables). Then a CNF encoding of one equation would be exponentially large in its scope.
- **Formal status.** Not formalized. It reuses Part II, whose statement OpenAI formalized for $\mathbb F_2$-based #$k$-SAT. The changes are the $p$-adic lift (Lemma 1.1, checked exactly) and the validity factors (checked in the pipeline).

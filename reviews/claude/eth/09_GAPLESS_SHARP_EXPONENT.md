# 9. Gapless proof, part II: #k-SAT in $2^{O_k(n(\log\log n)^2/\sqrt{\log n})}$ (the formally verified bound)

Claude, 28 October 2026. Part I ([08](08_GAPLESS_PROOF.md)) proved #$k$-SAT in $2^{O_k(n(\log\log n)^{k+1}/(\log n)^{1/(2k)})}$, which refutes ETH. Its exponent degrades with $k$, because Tovey's reduction multiplies the number of variables by the polylogarithmic occurrence bound. This file removes that loss, following the paper's §5–6, and proves:

> **Theorem 9.1 (paper Prop. 6.4).** Fix $k\ge2$. The number of satisfying assignments of a $k$-CNF with $n$ variables and $m$ clauses can be computed deterministically in
> $$(n+m)^{O_k(1)}\,2^{O_k(\rho(n))}\ \text{bit operations and}\ (n+m)^{O_k(1)}\ \text{bits of space},\qquad\rho(n)=\frac{n(\log\log n)^2}{\sqrt{\log n}} .$$

This is the bound in OpenAI's Lean statement `SqrtLogSpace.exact_count_and_sat`, whose `runtimeScale` is $n(\log\log(n+16))^2/\sqrt{\log(n+16)}$, and whose proof I built and audited ([06](06_LEAN_STATEMENT_AUDIT.md) §6.A). The paper's headline exponent $n\log\log n/\sqrt{\log n}$ (App. D) is not covered here.

**Idea.** Count common zeros of $O(n)$ sparse $\mathbb F_2$ equations instead of the clauses.
- Each equation is the parity of the violations of a short list of clauses.
- A Tanner code plus greedy rows chosen by exact counts make the common zeros exactly the satisfying assignments.
- Each equation is lifted to $\mathbb Z/2^\ell$ by the Beigel–Tarui polynomial. A cross-equation truncation keeps the total degree **linear in $n$**, independently of the number of clauses.

**What is reused from Part I.** The proofs carry over verbatim; only the input changes:
- disjoint sparsification (08 Prop. 2.4, Lemma 2.5);
- deformation counting (08 Thm 4.6), here with an algebra of coefficients;
- swap interpolation (08 Thm 5.4, Cor. 5.6);
- transfer and balanced schedules (08 Lemmas 6.1–6.2).

**What is new relative to the paper.**
- **The parity lift in closed form (Lemma 2.1).** Its $h$-adic digits are signed central binomial coefficients, $e_j=(-1)^{j+1}\binom{2j-1}{j-1}(2Y-1)$, because the lift is the idempotent $\frac12\big(1-(2Y-1)(1+4h)^{-1/2}\big)$. This replaces App. B.1's bound $2^{C\ell^2}$ by $2\cdot4^j$. It also explains the truncation: digit $j$ carries $h^j$, and $h(y)=y(y-1)$ is even, so the $V$-degree is a 2-adic precision.
- **Arithmetic mod $p^\mu$** with exact recovery of the $V$-coefficients replaces App. B.2–B.3's rational bit-length bookkeeping.
- **Part I's closed-form design** replaces the paper's §4 construction.

**Checks.** [`scripts/verify_sharp_route.py`](scripts/verify_sharp_route.py) ([`results/sharp_route.json`](results/sharp_route.json), ~11 s):
- **Lemma 2.1's closed form:** for every $\ell\le128$ it equals the digits produced by repeated division, the congruence holds on $|y|\le40$, and the norm bound holds. With $\binom{2j}j$ in place of $\binom{2j-1}{j-1}$ the congruence fails.
- **Theorem 3.1's algorithm end to end:** sparse $\mathbb F_2$ systems are counted exactly through the $V$-truncated algebra and the deformation of Part I, mod $p^\mu$. Interpolating from too few Taylor coefficients fails, as it should.
- **Lemma 4.1's mechanism** on Tanner-code toys: mixing, and an exhaustive minimum weight $\ge\delta(\delta-\varepsilon)q_Gd_G$, with nontrivial bounds of 10 and 64. Also the lemma's constant arithmetic on a grid of $(B,m)$.
- **Lemma 4.2's greedy rows** on small CNFs: each round keeps at most $\frac34$ of the surviving bad assignments, and the count ends exact.
- [`verify_parity_lift.py`](scripts/verify_parity_lift.py) and [`verify_padic_truncation.py`](scripts/verify_padic_truncation.py) check Lemmas 2.1–2.2, the latter for every prime.

Notation is as in Part I; $\log=\log_2$, and $\|G\|_1$ is the sum of the absolute values of the coefficients of $G$.

---

## 1. Deformation counting with an algebra of coefficients

**Theorem 1.1 (Part I Thm 4.6, general factors).** Let $\mathscr A$ be a commutative $\mathbb Z/p^\mu$-algebra that is free of rank $a_0$ as a $\mathbb Z/p^\mu$-module, and let $Q<p$. Let the factors be
$$\Phi_f\in\mathscr A[Y_v:v\in I_f]\qquad(f\in\mathcal F),$$
with scopes covering $[n]$, each given by a straight-line program over $\mathscr A$ with at most $\omega$ operations. Let $d$ be the total degree of
$$G=\prod_f\Phi_f\big((Y_{f,v})_{v\in I_f}\big),$$
**computed in $\mathscr A[Y_{f,v}]$** with independent indeterminates $Y_{f,v}$.

Take a deformation of degree $\le R'$ with $Q$ nodes for the scopes $(I_f)$, as in Part I Definition 4.2. Suppose each node has $r$ groups of size $\le s$ and at most $c_*$ measurement rows of norm $\le\beta$. Put $q=\lfloor dR'/Q\rfloor$ and $J=q+1$. Then
$$Z=\sum_{x\in\{0,1\}^n}\prod_f\Phi_f(x|_{I_f})\in\mathscr A$$
is computed with at most
$$Q\cdot C_q(n)\cdot(2\beta+1)^{c_*}\cdot2r\,2^{s}\cdot\mathrm{poly}(n,|\mathcal F|,J,a_0,\omega,\mu\log p)+\mathrm{poly}(Q,J,a_0,\mu\log p)$$
bit operations and polynomial space.

*Proof.* Repeat the proof of Part I Theorem 4.6 with $\mathscr A$ in place of $A=\mathbb Z/p^\mu$ and $\Phi_f$ in place of $\Phi_{C_f}$. Only two steps use more than the ring axioms.
1. **Degree (step 1 there).** Now
$$H(T)=\sum_x\prod_f\Phi_f\big((\pi\langle P_{f,v}(T),x\rangle)_v\big)\in\mathscr A[T].$$
   - Each summand is the image of $G$ under the substitution homomorphism $Y_{f,v}\mapsto\pi\langle P_{f,v}(T),x\rangle$, whose values have degree $\le R'$. So $\deg_TH\le dR'<QJ$.
   - $H(0)=Z$, because $P_{f,v}(0)=e_v$.
   - This is the only place where $d$ enters. It is why $d$ is measured after multiplying in $\mathscr A$: the sum of the separate degrees can be far larger.
2. **Interpolation (step 4 there).** The matrix of Lemma 4.3 has integer entries and is invertible over $\mathbb Z/p^\mu$, because the node differences are units ($p>Q$). Applying its inverse coordinatewise in a $\mathbb Z/p^\mu$-basis of $\mathscr A$ recovers $H(0)$ from the jets.

The remaining steps hold verbatim over any commutative ring:
- jets of degree $\le q$ by Lemma 4.4, since each $X$ enters through a row with entries in $z\,\mathcal O[z]$;
- separation by Lemma 4.5;
- the cost and space count.

A factor evaluation in $\mathscr A[z]/(z^J)$ costs $\omega$ operations there, each $\mathrm{poly}(J,a_0,\mu\log p)$ bit operations. ∎

The degree in Theorem 1.1 can be far below the sum of the factor degrees, because $\mathscr A$ may have nilpotents. §2 exploits exactly this.

---

## 2. The parity lift and cross-equation truncation

**Lemma 2.1 (the parity lift in closed form; paper Lemmas 5.1 and B.1).** Let $h=Y^2-Y$. For $\ell\ge1$ put
$$e_0=1-Y,\qquad e_j=(-1)^{j+1}\binom{2j-1}{j-1}(2Y-1)\quad(j\ge1),\qquad p_\ell=\sum_{j<\ell}e_j\,h^j .$$
1. $p_\ell=(1-Y)^\ell\sum_{j<\ell}\binom{\ell+j-1}jY^j$, the complemented Beigel–Tarui polynomial.
2. $p_\ell(y)\equiv[2\mid y]\pmod{2^\ell}$ for every integer $y$.
3. $\|e_j\|_1\le2\cdot4^j$ for all $j$.

*Proof.* Work in $\mathcal Q=\mathbb Q[Y]/(h^\ell)$, where $h$ is nilpotent.
1. **A square root.** $(2Y-1)^2=1+4h$. Put
$$s=\sum_{j<\ell}\binom{-1/2}j(4h)^j=\sum_{j<\ell}(-1)^j\binom{2j}jh^j .$$
   The identity $(1+t)\big((1+t)^{-1/2}\big)^2=1$ in $\mathbb Q[[t]]$ gives $(1+4h)s^2=1$, because $h^\ell=0$. So $\varepsilon=(2Y-1)s$ satisfies $\varepsilon^2=1$.
2. **The idempotent.** Expanding $(1-\varepsilon)/2$, with $\frac12\binom{2j}j=\binom{2j-1}{j-1}$ for $j\ge1$, gives
$$\frac{1-\varepsilon}2=(1-Y)+\sum_{1\le j<\ell}(-1)^{j+1}\binom{2j-1}{j-1}(2Y-1)h^j,$$
   which is the image of $p_\ell$.
3. **Local values.**
   - In $\mathbb Q[Y]/(Y^\ell)$ we have $\varepsilon\equiv-1\pmod Y$, so $\varepsilon-1$ is a unit, and $(\varepsilon-1)(\varepsilon+1)=0$ forces $\varepsilon=-1$. Hence $p_\ell\equiv1\pmod{Y^\ell}$.
   - Likewise $\varepsilon=1$ in $\mathbb Q[Y]/((Y-1)^\ell)$, so $p_\ell\equiv0\pmod{(Y-1)^\ell}$.
4. **Item 1.** The Beigel–Tarui polynomial also has degree $\le2\ell-1$ and the same two congruences: $\sum_{j<\ell}\binom{\ell+j-1}jY^j$ inverts $(1-Y)^\ell$ modulo $Y^\ell$. Two polynomials of degree $<2\ell$ that agree modulo the coprime $Y^\ell$ and $(Y-1)^\ell$ are equal.
5. **Item 2.** $p_\ell\in\mathbb Z[Y]$, and division by the monic $Y^\ell$ and $(Y-1)^\ell$ keeps integer coefficients. So $p_\ell=1+Y^\ell g_1=(Y-1)^\ell g_2$ with $g_1,g_2\in\mathbb Z[Y]$. For even $y$, $2^\ell\mid y^\ell$; for odd $y$, $2^\ell\mid(y-1)^\ell$.
6. **Item 3.** $\|e_0\|_1=2$, and $\|e_j\|_1=3\binom{2j-1}{j-1}\le\frac32\cdot4^j$ for $j\ge1$. ∎

The digits are exactly what repeated division by $h$ produces (the paper's construction), since the $h$-adic expansion with digits of degree $\le1$ is unique. For integers $y_1,\dots,y_m$ put $B_i(V)=\sum_{j<\ell}e_j(y_i)h(y_i)^jV^j$, and let $\mathcal S_\ell\big(\sum_sb_sV^s\big)=\sum_{s<\ell}b_s$.

**Lemma 2.2 (cross-equation truncation; paper Lemma 5.2).**
1. $\mathcal S_\ell\big(\prod_{i=1}^mB_i(V)\bmod V^\ell\big)\equiv\prod_{i=1}^m[2\mid y_i]\pmod{2^\ell}$.
2. If instead the $y_i$ are polynomials of degree $\le d_*$ in some indeterminates, then every coefficient of $V^s$ ($s<\ell$) in $\prod_iB_i\bmod V^\ell$ has degree $\le d_*(m+2(\ell-1))$.

*Proof.*
1. The sum of **all** coefficients of the untruncated product is $\prod_iB_i(1)=\prod_ip_\ell(y_i)$. A dropped term has $\sum_ij_i\ge\ell$ and contains $\prod_ih(y_i)^{j_i}$, which is divisible by $2^{\sum j_i}$ because $h(y)$ is even. So it is $\equiv0\pmod{2^\ell}$. Lemma 2.1 finishes.
2. A retained term has $\sum j_i\le\ell-1$ and degree $\le\sum_id_*(1+2j_i)$. ∎

Without truncation the degree would be $d_*m(2\ell-1)$, which is quadratic in $n$ when $\ell,m\asymp n$. [07](07_TOOLKIT.md) E2 proves the same lemma for every prime $p$.

---

## 3. Counting common zeros of sparse $\mathbb F_2$ systems

**Theorem 3.1 (paper Thm 2.2).** Fix $k,c$. Let $P_1,\dots,P_m\in\mathbb F_2[x_1,\dots,x_n]$ with $m\le cn$, each given as a list of at most $(\log n)^c$ monomials of degree $\le k$. Then $N=\#\{x\in\{0,1\}^n:P_i(x)=0\ \forall i\}$ can be computed in $2^{O_{k,c}(\rho(n))}$ bit operations and $\mathrm{poly}(n)$ space.

*Proof.* Assume $n\ge n_0(k,c)$; otherwise use exhaustive search. An empty system has $2^n$ zeros. Discard equations $0=0$; an equation $1=0$ gives $N=0$.
1. **Factors and incidences.** For equation $i$, let $\mathcal M_i$ be its monomials (variable sets, duplicates cancelled mod 2), $I_i=\bigcup_{M\in\mathcal M_i}M$, and
$$L_i(X)=\sum_{M\in\mathcal M_i}\prod_{v\in M}X_v\quad(\deg\le k),\qquad\mathcal B_i(X,V)=\sum_{j<\ell}e_j(L_i(X))\,h(L_i(X))^jV^j .$$
   - On $\{0,1\}^n$, $L_i(x)$ counts the monomials equal to $1$, so $P_i(x)=0$ iff $L_i(x)$ is even.
   - Add a neutral factor $1$ with scope $\{v\}$ for each variable. The incidence count is $M=\sum_f|I_f|\le n+mk(\log n)^c$.
2. **Branching on frequent variables.** Let $\lambda=\lceil\log n\rceil$, $b_0=\lceil M/n\rceil\lceil\sqrt\lambda\rceil$, and let $S$ be the set of variables lying in more than $b_0$ scopes. Since the degrees sum to $M$, $|S|<M/b_0\le n/\sqrt\lambda$.
   - **Branches.** For each $\eta\in\{0,1\}^S$, substitute $x_S=\eta$ in every $P_i$ and cancel monomials mod 2. Degrees and monomial counts do not increase. Discard equations that become $0=0$; a branch with an equation $1=0$ has no zeros. The branches partition $\{0,1\}^n$, so their counts add up to $N$.
   - **Branch shape.** Rebuild the factors of step 1 on the free variables. A branch has $n'=n-|S|\ge n/2$ free variables. Each lies in $b\le b_0$ scopes, counting its neutral factor, and each scope has $w\le k(\log n)^c$ variables.
   - **Exponent.** With $\lambda'=\lceil\log n'\rceil\in\{\lambda-1,\lambda\}$, both $b,w\le\lambda'^{a}$ for $a=c+1$ and $n\ge n_0$.
3. **Parameters on a branch.** Let $\nu'=\lceil\log\lambda'\rceil$ and
$$r=\lfloor\sqrt{\lambda'}\rfloor,\qquad Q=\Big\lfloor\frac{\lambda'}{4\nu'}\Big\rfloor,\qquad R=2(a+2)\lceil\sqrt{\lambda'}\rceil,\qquad\ell=n+1 .$$
   Let $p$ be the least prime $>Q$.
4. **The design.** Part I §6 uses only the scopes, so it applies to the factor scopes. As in Part I (P3)–(P4), $\log D\le\lambda'/4$ for $D=r^{2Q}$. So $2(b+w)D^2\le4\lambda'^a2^{\lambda'/2}<2^{\lambda'-1}<n'\le M'$ and $2Qr^2<n'$ for $n\ge n_0$. Hence:
   - Lemma 6.2 gives schedules with $h\le2M'/D$ and groups of size $\le2n'/r$.
   - Theorem 5.4 applies with the interval of Corollary 5.6, since $R$ is even and $4\le R\le Q/2$ for $n\ge n_0$. Its coefficients reduce mod $p^\mu$, since $p\nmid2\cdot Q!$.
   - Lemma 6.1 turns these into a deformation of degree $\le R-1$ with $Q$ nodes. It has at most $c_*\le2hc\le8M'e^{-R^2/(8Q)}$ measurement rows of norm $\le2^Q$ per node (Corollary 5.6). Here $R^2/(8Q)\ge2(a+2)^2\nu'$.
5. **The algebra and the degree.** Let $\mathscr A=(\mathbb Z/p^\mu)[V]/(V^\ell)$, which is free of rank $\ell$, with $\mu$ fixed in step 6. The factors are the images of the $\mathcal B_i$ and of $1$.
   - Each $\mathcal B_i$ is a straight-line program with $\mathrm{poly}(n)$ operations: the monomials of $L_i$, then $h$, then Horner's rule in $V$.
   - By Lemma 2.2(2) with $d_*=k$, the product has degree in $\mathscr A[Y]$ at most
$$d=k(m+2(\ell-1))\le k(c+2)n\le\kappa n',\qquad\kappa=2k(c+2) .$$
   - Theorem 1.1 computes $\bar H(0)\in\mathscr A$, the image of the integer polynomial
$$H(0)=\sum_{x}\prod_i\mathcal B_i(x,V)\bmod V^\ell=\sum_{s<\ell}b_sV^s,\qquad b_s\in\mathbb Z .$$
6. **Exact recovery of the count.**
   - Let $y_{\max}=(\log n)^c\ge L_i(x)\ge0$. Lemma 2.1(3) and $|h(y)|\le y^2$ give
$$\|\mathcal B_i(x,\cdot)\|_1\le\sum_{j<\ell}2\cdot4^j(1+y_{\max})^{2j+1}\le(2+2y_{\max})^{2\ell}=:U .$$
     So $|b_s|\le2^nU^m=:Z$, where $\log Z=O_{k,c}(n^2\log\log n)$.
   - Choose the least $\mu$ with $p^\mu>2Z$. The symmetric residues of $\bar H(0)$ are then the $b_s$.
   - By Lemma 2.2(1), summed over $x$: $\mathcal S_\ell(H(0))=\sum_sb_s\equiv N_{\rm branch}\pmod{2^\ell}$. Since $0\le N_{\rm branch}<2^\ell$, it is the least nonnegative residue.
7. **Time.** Apply Theorem 1.1 on the branch.
   - **Its parameters.** $q\le d(R-1)/Q\le xn'$, with $x:=9\kappa(a+2)\nu'/\sqrt{\lambda'}<e$ for $n\ge n_0$. Also $\beta=2^Q$, so $(2\beta+1)^{c_*}\le2^{(Q+2)c_*}$, and $s\le2n'/r$. The polynomial factor is $\mathrm{poly}(n)$, because $\ell$, $J$, $\omega$ and $\mu\log p$ are.
   - **Total exponent.** With Part I Lemma A.1 and $M'\le bn'\le\lambda'^an'$:
$$\frac{\log(\text{time per branch})}{n'}\le x\log\frac{e^2}x+8\lambda'^{a}e^{-2(a+2)^2\nu'}(Q+2)+\frac2r+O\Big(\frac{\log n}{n}\Big)\le C_{k,c}\frac{\nu'^2}{\sqrt{\lambda'}} .$$
   
   The three terms:
   - **Taylor jets.** Since $x\ge e^2/\sqrt{\lambda'}$, $x\log(e^2/x)\le x\log\sqrt{\lambda'}\le x\nu'/2$.
   - **Measurements.** $e^{-2(a+2)^2\nu'}\le\lambda'^{-2.88(a+2)^2}$, so this term is at most $8\lambda'^{a+1}\lambda'^{-2.88(a+2)^2}\le1/\sqrt{\lambda'}$.
   - **Groups.** $2/r\le3/\sqrt{\lambda'}$.
   
   The bound $q\le xn'$ uses $(R-1)/Q\le9(a+2)\nu'/\sqrt{\lambda'}$, which holds once $\lambda'\ge8\sqrt{\lambda'}+36\nu'$. Multiplying by the at most $2^{n/\sqrt\lambda}$ branches gives $2^{O_{k,c}(n\nu^2/\sqrt\lambda)}=2^{O_{k,c}(\rho(n))}$, since $\nu\le\log\log n+2$ and $\lambda\ge\log n$.
8. **Space.** Branches are processed one at a time, adding each count to an accumulator that stays in $[0,2^n]$. Within a branch, Theorem 1.1 uses polynomial space. ∎

---

## 4. From #$k$-SAT to common zeros

Let $F$ have clauses $C_1,\dots,C_m$ on $n$ variables. Let $v_e(x)\in\mathbb F_2[x]$ be the violation indicator of $C_e$: the product of $1+x_v$ over positive literals and $x_v$ over negative ones. It has $\le2^k$ monomials of degree $\le k$, and $v(x)=(v_e(x))_e$. A binary matrix $A$ with $m$ columns defines the system $P_i(x)=\sum_{e:A_{ie}=1}v_e(x)=0$, with $N(A)$ common zeros. Every satisfying assignment has $v(x)=0$, so
$$N(A)=\#F+D(A),$$
where $D(A)$ counts the unsatisfying $x$ with $Av(x)=0$.

**Lemma 4.1 (initial tests from a Tanner code; paper Lemma 6.1).** Let $m\le Bn$, $\mu=\lfloor\log m\rfloor+1$, $P=20000B^2\mu^2$ and $D_0=320000B^2\mu^3$. One can construct deterministically, with $\exp(\mathrm{poly}(B,\log(m+2)))$ bit operations and $\mathrm{poly}(n,m,B)$ space, a binary matrix $A^{(0)}$ with the following properties:
- it has $m$ columns and at most $3n$ rows;
- each row has at most $D_0$ ones;
- every nonzero $z\in\{0,1\}^m$ with $A^{(0)}z=0$ has $|\operatorname{supp}z|\ge m/P$.

*Proof.*
1. **The graph.** Let $d_G=\min(m,D_0)$ and let $q_G=2^{\ell_G}$ be the least power of two with $q_Gd_G\ge m$; then $q_Gd_G<2m$. The vertices are two copies of $\{0,1\}^{\ell_G}$. Left vertex $u$ is joined to right vertex $u\oplus s_j$ by an edge with port $j\in[d_G]$. This graph is $d_G$-regular.
2. **$\varepsilon$-biased shifts** with $\varepsilon=\frac1{200B\mu}$: we need $|\sum_j(-1)^{b\cdot s_j}|\le\varepsilon d_G$ for all $b\ne0$.
   - If $q_G=1$ this is vacuous.
   - Otherwise $d_G=D_0$. For uniform random shifts, Hoeffding plus a union bound give failure probability $\le2q_Ge^{-\varepsilon^2D_0/2}=2q_Ge^{-4\mu}<1$, since $q_G<2^{\mu+1}$. So exhaustive search over the $\le(2m)^{D_0}$ shift tuples finds valid shifts.
3. **Mixing.** The characters $\chi_b(u)=(-1)^{b\cdot u}$ diagonalize the biadjacency operator, with eigenvalue $\sum_j(-1)^{b\cdot s_j}$. So, for vertex sets $U,V$,
$$e_G(U,V)\le d_G|U||V|/q_G+\varepsilon d_G\sqrt{|U||V|}.$$
4. **The local code.** Let $c_G=\lceil d_G/(8B)\rceil$ and $\delta=\frac1{100B\mu}$. A uniformly random $c_G\times d_G$ matrix kills a fixed nonzero word with probability $2^{-c_G}$. There are fewer than $(2d_G)^{2h}<2^{d_G/(25B)}$ nonzero words of weight $<\delta d_G$, where $h=\lceil\delta d_G\rceil-1$ and $2d_G<2^{\mu+1}\le2^{2\mu}$. So some $A_{\rm loc}$ has kernel of minimum weight $\ge\delta d_G$; exhaustive search finds it.
5. **The Tanner code.** Impose $A_{\rm loc}$ at every vertex on its incident edges, ordered by port. Let $h'$ be the weight of a nonzero accepted edge word, touching vertex sets $U,V$.
   - Each touched vertex has $\ge\delta d_G$ nonzero incident edges, so $h'\ge\delta d_G\sqrt{|U||V|}$.
   - With mixing, $\delta d_G\sqrt{|U||V|}\le h'\le e_G(U,V)$ forces $\sqrt{|U||V|}\ge(\delta-\varepsilon)q_G$.
   - Hence $h'\ge\delta(\delta-\varepsilon)q_Gd_G=q_Gd_G/(20000B^2\mu^2)\ge m/P$.
6. **Shortening.** Keep $m$ of the edge coordinates, setting the others to $0$, and restrict the checks to them; an accepted word then extends by zeros. Rows have $\le d_G\le D_0$ ones.
   - If $q_G=1$ there are $2c_G\le m/(4B)+2\le3n$ rows.
   - If $q_G>1$, then $d_G=D_0\ge8B$ and there are $2q_Gc_G\le q_Gd_G/(2B)<m/B\le n$ rows. ∎

**Lemma 4.2 (greedy rows by exact counts; paper Lemma 6.2).** Starting from $A^{(0)}$, a deterministic procedure appends $3n$ rows and reaches $N(A)=\#F$. Each appended row is the parity of a list of $P$ clause indices with coefficients. The procedure uses at most $3n(2m)^P$ evaluations of $N(\cdot)$, on systems with at most $6n$ equations, each with at most $2^k\max(D_0,P)$ monomials of degree $\le k$.

*Proof.* A candidate row is a list $\gamma=((e_1,c_1),\dots,(e_P,c_P))\in([m]\times\{0,1\})^P$, representing $h_\gamma=\sum_jc_j\mathbf e_{e_j}$ mod 2.
1. **One bad assignment.** Let $x$ be unsatisfying with $Av(x)=0$. Then $A^{(0)}v(x)=0$, so $w=|\operatorname{supp}v(x)|\ge m/P$.
   - For a uniform $\gamma$, some $e_j$ hits the support with probability $\ge1-(1-1/P)^P\ge1-e^{-1}>\frac12$.
   - Given a hit at position $j$, the uniform $c_j$ makes $h_\gamma\cdot v(x)$ a uniform bit.
   
   So $\Pr_\gamma[h_\gamma\text{ rejects }x]\ge\frac14$.
2. **One round.** Hence $\mathbb E_\gamma D(A;h_\gamma)\le\frac34D(A)$, where $A;h_\gamma$ is $A$ with $h_\gamma$ appended. Since $N(A;h_\gamma)=\#F+D(A;h_\gamma)$, a candidate minimizing $N(A;h_\gamma)$, computed exactly for all $(2m)^P$ candidates, satisfies $D\le\frac34D(A)$.
3. **All rounds.** After $3n$ rounds, $D\le(\frac34)^{3n}2^n=(\frac{27}{32})^n<1$, so $D=0$. The procedure never needs to know $\#F$. ∎

**Proof of Theorem 9.1.**
1. **Reduction to $n\ge n_0(k)$.** Normalize in $\mathrm{poly}(n+m)$ time; an empty clause gives $0$, and no clauses gives $2^n$. Handle $n<n_0(k)$ by exhaustive search.
2. **Sparsify.** Let $\lambda=\lceil\log n\rceil$, $\nu=\lceil\log\lambda\rceil$ and $s=\lfloor\sqrt\lambda/\nu^2\rfloor\ge2$. Run Part I's sparsification (Prop. 2.4 with Lemma 2.5) at $\sigma=1/s$.
   - It gives at most $\mathrm{poly}(n)2^{n/s}$ nodes, with $n/s=O(n\nu^2/\sqrt\lambda)=O(\rho(n))$.
   - Every leaf is a normalized $k$-CNF on the same $n$ variables, with occurrence bound $B=O_k((s\log s)^{k-1})\le(\log n)^{c_1(k)}$ and $m_L\le Bn$ clauses.
3. **Count each leaf.** A leaf with no clauses has count $2^n$, and one containing the empty clause has count $0$. Otherwise:
   - Build $A^{(0)}$ by Lemma 4.1, in $\exp((\log n)^{O_k(1)})$ time.
   - Run Lemma 4.2. Every system it evaluates has $\le6n$ equations with $\le2^k\max(D_0,P)\le(\log n)^{c_2(k)}$ monomials of degree $\le k$, so Theorem 3.1 applies with $c=\max(6,c_2(k))$.
   - There are $3n(2m_L)^P+1=2^{(\log n)^{O_k(1)}}=2^{o(\rho(n))}$ calls, each $2^{O_k(\rho(n))}$. The last call returns $\#L$.
4. **Sum** the leaf counts in an $(n+1)$-bit accumulator; the leaves partition the satisfying assignments.
5. **Space.** Depth-first traversal, sequential candidate evaluation, and the polynomial space of Theorem 3.1 give polynomial space. ∎

---

## 5. Relation to the paper and to the formal proof

| This file | Paper | Notes |
|---|---|---|
| Thm 1.1 | Prop. 3.4 | same degree hypothesis (product computed in the algebra); here over a $\mathbb Z/p^\mu$-algebra rather than a $\mathbb Q$-algebra, reusing Part I |
| Lemma 2.1 | Lemmas 5.1, B.1 | **closed form** $e_j=(-1)^{j+1}\binom{2j-1}{j-1}(2Y-1)$; norms $\le2\cdot4^j$ instead of $2^{C\ell^2}$ |
| Lemma 2.2 | Lemma 5.2 | identical |
| Thm 3.1 | Thm 2.2, §4, §5.3–5.5, App. B | Part I's design ($r=\lfloor\sqrt{\lambda'}\rfloor$, $Q=\lfloor\lambda'/(4\nu')\rfloor$) replaces §4's ($r=\lfloor\sqrt{\lambda'}/\nu'^2\rfloor$, $Q=\lfloor\lambda'/(128\nu')\rfloor$); arithmetic mod $p^\mu$ replaces App. B's bit-length bounds |
| Lemma 4.1 | Lemma 6.1 | same constants; checked on a grid |
| Lemma 4.2 | Lemma 6.2 | identical |
| Thm 9.1 | Prop. 6.4 | identical exponent $\rho(n)$ |

**The formal proof** (`OAI/Algorithms/SubexponentialSAT/…`) follows this route and proves the same $\rho(n)$ statement with explicit Turing-machine time and space ([06](06_LEAN_STATEMENT_AUDIT.md)). The matching modules include:
- `Truncation`, the $V$-truncation before substitution;
- `LiftRing` and `ActualLiftRing`, the coefficient ring: a Galois ring over $\mathbb Z/p^\mu$, to which rows and nodes built over $\mathbb F_{p^h}$ lift (Part I instead takes the least prime $p>Q$);
- `SqrtLog/SparseCode*` and `SqrtLog/CayleyGraph`, the Tanner tests;
- `ConditionalExpectation` and `Variance`, the schedules.

**Status.** With Part I, every bound the paper proves *and* OpenAI formalized now has a complete prose proof in this repository:
- $2^{o(n)}$ (Part I);
- $2^{O(\rho(n))}$ (this file).

What remains prose-only in the paper is App. D's sharper $n\log\log n/\sqrt{\log n}$, which is not formalized and which I rated ≈80% ([01](01_VERDICT.md)).

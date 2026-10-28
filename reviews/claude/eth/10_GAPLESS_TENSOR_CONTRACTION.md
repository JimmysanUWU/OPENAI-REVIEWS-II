# 10. Gapless proof, part III: the headline exponent $n\log\log n/\sqrt{\log n}$ (paper Theorem 1.1, App. D)

Claude, 28 October 2026. Parts I–II ([08](08_GAPLESS_PROOF.md), [09](09_GAPLESS_SHARP_EXPONENT.md)) prove the bounds OpenAI formalized. This part proves the paper's headline theorem, which is prose-only (App. D) and not formalized:

> **Theorem 10.1 (paper Thm 1.1).** Fix $k\ge2$. The number of satisfying assignments of a $k$-CNF with $n$ variables and $m$ clauses can be computed deterministically in
> $$(n+m)^{O_k(1)}\,2^{O_k(n\log\log n/\sqrt{\log n})}\ \text{bit operations and}\ (n+m)^{O_k(1)}\ \text{space}.$$
> If $m\le Dn$ for a fixed $D$, the time is $(n+m)^{O_{k,D}(1)}\,2^{O_{k,D}(n/\sqrt{\log n})}$.

Both bounds rest on one theorem about Boolean tensor networks (§1–§7):

> **Theorem 10.2 (Boolean tensor contraction; paper Thm D.1).** Fix $C\ge1$, $A>0$ and $b\ge1$, and let $s\ge2$. Consider a bipartite Boolean tensor network with:
> - at most $Cs$ edges, maximum degree at most $\Delta=\lceil s^{9/10}\rceil$, and at most $s^C$ vertices;
> - entries in $(\mathbb Z/2^{P_0})[V]/(V^\ell)$ with $1\le P_0,\ell\le Cs$;
> - a description of total length $L_{\rm d}$ from which each entry is computed in time $(s+L_{\rm d})^b2^{As^{9/10}}$ and space $(s+L_{\rm d})^b$.
>
> Then the contraction is computed deterministically in $(s+L_{\rm d})^{O(1)}2^{O(s/\sqrt{\log s})}$ bit operations and $(s+L_{\rm d})^{O(1)}$ space, with constants depending only on $C,A,b$.

**Verdict on App. D: correct.** I re-derived every step below and found no error. One statement is stronger than the proof needs: (D.5) asserts $u\ge\lambda^{4/5}$ at every depth, which holds only from $\log_2\lambda\approx500$. The proof needs only the absorption inequalities (A1)–(A2) of Lemma 7.1, which hold from $\log_2\lambda\approx26$. Both are "sufficiently large $s$" statements, so this is not an error.

My earlier rating of App. D was ≈80%, because the recursion's bookkeeping had not been re-derived ([01](01_VERDICT.md)). It is now re-derived and checked end to end.

**Idea.** Contract a network by recursion on its wires.
- **Deformation.** Replace the identification of each wire's two ends by a matrix $W(T)$ with $W(0)=I$, built from the Cauchy kernel $C(T)/(C'(0)(T-d))$.
- **Locality at the nodes.** At a node $t$ with $v_2(t)=i$, $W(t)$ only joins ports whose vertices share their $i$-th label symbol. So the network splits into $r=64$ groups, each with at most $M/8$ wires.
- **Modes.** Near the node, the crossing terms of a tile of columns expand in the power moments $\sum_jj^bX_y$, with 2-adic weight $b+1$. Modulo $2^{p_t}$ only $\exp(O(M/\sqrt u))$ mode monomials survive, and each is a product of contractions of the groups.
- **Tapered recovery.** Hermite interpolation with multiplicity $h\,(u-v_2(t))^+$ recovers $T=0$ while losing only $O(MD_t/u)$ bits of 2-adic precision.
- **Recursion.** With $u\asymp\log s/2^d$ at depth $d$, the exponents $M_d/\sqrt{u_d}\asymp(s/\sqrt{\log s})(\sqrt2/8)^d$ sum geometrically to $O(s/\sqrt{\log s})$.

**Checks.** [`scripts/verify_tensor_contraction.py`](scripts/verify_tensor_contraction.py) ([`results/tensor_contraction.json`](results/tensor_contraction.json), ~45 s). Every check is exact, and each has a negative control.
- **The child identity (Lemma 6.3).** On 66 (network, node, value) cases, packed as in Lemma 5.1 and evaluated in $\mathcal O=\mathbb Z_{(2)}[U]/(f)$, the deformed contraction equals the sum over mode tuples of products of brute-force child contractions. Dropping the crossing terms, or crediting modes to their source group, breaks it.
- **The tile expansion (Lemma 6.2)** holds modulo $2^{40}$ at 20 nodes. Keeping only $b=0$ fails at 12 of them.
- **One full recursion level at the honest precision** $p_t=p+8hD_t+8M(a+1)$ (Lemmas 2.1, 3.1, 4.1). The jets are obtained by evaluation at binary points of $\mathcal O$ and Lagrange interpolation over $(\mathbb Z/2^{p_t})[U]/(f)$; tapered recovery then returns the contraction exactly in 8/8 cases. With $p_t=p$ it fails in 8/8.
- **The mode count (Lemma 6.4)** against exact partition counts.
- **The encodings.** §9's incidence-subset networks reproduce the common-zero count, and §8's clause network reproduces $\#F$.
- **The parameter inequalities (Lemma 7.1)** at every depth, and their thresholds.
- [`verify_appendix_d.py`](scripts/verify_appendix_d.py) checks Lemma 3.1(3), (D.10) and Lemma 4.1 separately.

Notation: $\log=\log_2$, $\ln$ is natural, and $v_2$ is the 2-adic valuation.

---

## 1. Tensor networks and the deformed contraction

**Definition 1.1.** A *(bipartite Boolean) tensor network* is a bipartite multigraph with sides $\mathcal R$ (rows) and $\mathcal C$ (columns) and edge set $E$, $|E|=M$, together with:
- for each vertex $v$, an ordered list $\delta(v)$ of its incident edges (its *ports*; parallel edges are distinct ports);
- a tensor $T_v:\{0,1\}^{\delta(v)}\to\mathscr A$, for a commutative ring $\mathscr A$.

Its *contraction* is $\mathrm{ct}(N)=Z=\sum_{x\in\{0,1\}^E}\prod_vT_v(x|_{\delta(v)})$. An isolated vertex contributes its scalar entry. (In §6–§7, $Z$ also denotes mode variables; contractions of children are written $\mathrm{ct}(\cdot)$.)

Give the row ends of the edges bits $Y\in\{0,1\}^E$ and the column ends bits $X\in\{0,1\}^E$. Put $A(Y)=\prod_{v\in\mathcal R}T_v(Y|_{\delta(v)})$, $B(X)=\prod_{v\in\mathcal C}T_v(X|_{\delta(v)})$, and $l_0(w)=1-w$, $l_1(w)=w$.

**Lemma 1.2.** For a matrix $W$ indexed by (row end, column end) with entries in a commutative ring containing $\mathscr A$, put
$$F_W=\sum_{X,Y\in\{0,1\}^E}A(Y)B(X)\prod_{x\in E}l_{Y_x}\big((WX)_x\big).$$
Then $F_I=Z$.

*Proof.* For Boolean $X_x$, $l_{Y_x}(X_x)=[Y_x=X_x]$. ∎

---

## 2. The coefficient ring

Fix $e\ge1$ and a monic $f\in\mathbb Z[U]$ of degree $e$, with coefficients in $\{0,1\}$, irreducible modulo 2. One exists for every $e$, and exhaustive search with trial division finds it in $\mathrm{poly}(2^e)$ time. Put
$$\mathcal O=\mathbb Z_{(2)}[U]/(f),\qquad K=\mathbb Q[U]/(f),\qquad\mathcal A=\mathcal O[V]/(V^\ell),\qquad\mathcal A_p=\mathcal A/2^p\mathcal A ,$$
where $\mathbb Z_{(2)}$ is the rationals with odd denominator. For $x\in K\setminus0$ let $v_2(x)$ be the least 2-adic valuation of its coordinates in $1,U,\dots,U^{e-1}$. Valuations and divisibility in $\mathcal A$ are coordinatewise.

**Lemma 2.1.**
1. $K$ is a field and $v_2(xy)=v_2(x)+v_2(y)$.
2. An element of $\mathcal O$ is a unit iff it is nonzero modulo 2.
3. The $2^e$ elements of $\mathcal O$ with coordinates in $\{0,1\}$ have pairwise differences that are units. Hence a polynomial over $\mathcal A_p$ of degree $<N\le2^e$ is determined by its values at $N$ of them, and Lagrange interpolation over $\mathcal A_p$ recovers its coefficients.
4. Ring operations in $\mathcal A_p$ cost $\mathrm{poly}(e,\ell,p)$ bit operations.

*Proof.*
1. **Field.** A monic $f$ irreducible mod 2 is irreducible over $\mathbb Z$, hence over $\mathbb Q$ by Gauss's lemma.
   **Valuation.** $\mathcal O/2\mathcal O=\mathbb F_2[U]/(\bar f)=\mathbb F_{2^e}$ is a field. Write $x=2^{v_2(x)}x'$ and $y=2^{v_2(y)}y'$ with $x',y'$ nonzero mod 2; then $x'y'$ is nonzero mod 2.
2. If $x$ is nonzero mod 2, the determinant of multiplication by $x$ reduces to the norm of $\bar x\ne0$ in $\mathbb F_{2^e}$. So it is odd, and the inverse matrix has entries in $\mathbb Z_{(2)}$. The converse is clear.
3. Distinct binary-coordinate elements differ modulo 2. The Vandermonde determinant is a product of differences, hence a unit.
4. Standard. ∎

**Canonical lifts.** An element of $\mathcal A_p$ is stored as integer coordinates in $[0,2^p)$, and its *canonical lift* is the element of $\mathcal A$ with those coordinates. An element created at precision $p_c$ has the same canonical lift at every precision $p\ge p_c$.

---

## 3. The Cauchy kernel

Fix $a\ge1$, and put
$$L=2^a,\qquad S=\{-L,\dots,L-1\},\qquad C(T)=\prod_{q\in S}(T-q),\qquad\alpha=2^{a+1}.$$
For an integer $d$ with $|d|<L$ put
$$W_d(T)=\frac{C(T)}{C'(0)(T-d)}=\frac1{C'(0)}\prod_{q\in S\setminus\{d\}}(T-q),$$
a polynomial over $\mathbb Q$ of degree $2L-1$.

**Lemma 3.1 (paper (D.7), Lemma D.4).**
1. $W_d(0)=[d=0]$.
2. For $t\in S$, $W_d(t)=0$ unless $d=t$.
3. For $t\in S$, $W_d(t+\alpha z)\in\mathbb Z_{(2)}[z]$. If $d\ne t$, it is divisible by $2^{a+1-v_2(t-d)}$, which is at least $2$.
4. Given $t$, $d$ and $\zeta\in\mathcal O$ mod $2^p$, the value $W_d(t+\alpha\zeta)\bmod2^p$ is computed in $\mathrm{poly}(L,p,e)$ bit operations.

*Proof.*
- **A residue system.** $S$ is a complete residue system mod $2^{a+1}$. So for $q\ne t$ in $S$, $v_2(t-q)\le a$, and the multiset $\{v_2(t-q):q\in S\setminus\{t\}\}$ does not depend on $t$. Hence $v_2(C'(t))=v_2(C'(0))$.
- **(1), (2).** The product vanishes at every point of $S\setminus\{d\}$, and at $T=d=0$ it equals $C'(0)$.
- **(3).** For $q\ne t$, the factor $(t-q)+\alpha z$ is $2^{v_2(t-q)}$ times a polynomial over $\mathbb Z_{(2)}$ with odd constant term, because $v_2(\alpha)=a+1>v_2(t-q)$.
  - If $d=t$, the product over $q\ne t$ has content $2^{v_2(C'(t))}=2^{v_2(C'(0))}$, which $C'(0)$ cancels.
  - If $d\ne t$, the product contains the factor $q=t$, namely $\alpha z$, and lacks $q=d$. Its content is at least $2^{a+1}2^{v_2(C'(0))-v_2(t-d)}$.
- **(4).** Divide each factor by its known power of two, cancel the powers of two against $C'(0)$, and invert the odd part of $C'(0)$ modulo $2^p$. ∎

**Packets.** In §5 the edges are placed in *packets* of $L$ slots, each edge occupying one slot at both of its ends. Within a packet, let $W_{xy}(T)=W_{y-x}(T)$ for row slot $x$ and column slot $y$; note $|y-x|<L$. Entries between different packets are $0$, and $W$ is restricted to the slots that hold edges (*real* ports). Let $F(T)=F_{W(T)}$, formed with the canonical lifts of the entries.

**Corollary 3.2.** $\deg_TF\le M(2L-1)$; $F(0)$ is the contraction of the lifted network; and $F(t+\alpha z)\in\mathcal A[z]$ for $t\in S$.

*Proof.* Each of the $M$ affine factors has degree $\le2L-1$ in $T$. $W(0)=I$ by Lemma 3.1(1) and Lemma 1.2. The last claim is Lemma 3.1(3). ∎

---

## 4. Tapered Hermite recovery

**Lemma 4.1 (paper Lemma D.5).** Let $1\le u\le a$, $M\ge u$ and $h=\lceil3M/u\rceil$. For $t\in S\setminus\{0\}$ with $i=v_2(t)$, put
$$q_t=h\max(0,u-i),\qquad D_t=L/2^i,\qquad p_t=p+8hD_t+8M(a+1).$$
Let $F\in(K[V]/(V^\ell))[T]$ with $\deg F\le M(2L-1)$, $F(0)\in\mathcal A$, and $F(t+\alpha z)\in\mathcal A[z]$ for all $t\in S$. Then $F(0)\bmod2^p$ is determined by the first $q_t$ coefficients of $F(t+\alpha z)\bmod2^{p_t}$ for the $t$ with $q_t>0$. It is computed from them in time and space polynomial in $M,L,\ell,e$ and the precisions.

*Proof.*
1. **Enough conditions.** The $L$ odd $t$ have $q_t=hu\ge3M$, so $\sum_tq_t\ge3ML>\deg F$. Also $h\le4M/u$, so $q_t\le4M$.
2. **The Hermite formula (D.9).** Let $E_t(T)=\prod_{q\in S\setminus\{0,t\}}\big(\frac{T-q}{t-q}\big)^{q_q}$ and $P_t(z)=\big[F(t+\alpha z)/E_t(t+\alpha z)\big]_{<q_t}$ (truncation in $z$). Then
$$G(T)=\sum_{t:\,q_t>0}E_t(T)\,P_t\Big(\frac{T-t}\alpha\Big)$$
   agrees with $F$ to order $q_{t'}$ at every node $t'$: the $t'$ term reproduces $F$'s jet there, and the other terms vanish there to that order. Since $\deg G<\sum q_t$, $G=F$, and in particular
$$F(0)=\sum_tE_t(0)\,P_t(-t/\alpha).$$
3. **No loss in the quotient.** $E_t(t+\alpha z)=\prod_q(1+\alpha z/(t-q))^{q_q}$, and $v_2(\alpha/(t-q))\ge1$. So it is integral with constant term $1$, its inverse modulo $z^{q_t}$ is integral, and $P_t\bmod2^{p_t}$ is determined by $F(t+\alpha z)\bmod2^{p_t}$.
4. **Loss at $E_t(0)$ (D.10).** Write $q(m)=h\max(0,u-m)$ for the multiplicity at valuation $m$, and $D=D_t$. For $1\le j\le a-i$, $S$ has $D/2^j$ points with $v_2(q)=i+j$, and the same number with $v_2(t-q)=i+j$; the latter have $v_2(q)=i$. Points of valuation $<i$ have $v_2(t-q)=v_2(q)$ and contribute $0$. Hence
$$v_2(E_t(0))=\sum_{q\ne0,t}q_q\big(v_2(q)-v_2(t-q)\big)=D\sum_{j=1}^{a-i}\frac j{2^j}\big(q(i+j)-q(i)\big)\ge-hD\sum_{j\ge1}\frac{j^2}{2^j}=-6hD,$$
   because $q(i+j)-q(i)\ge-hj$.
5. **Loss at $(-t/\alpha)^j$.** $v_2(-t/\alpha)=i-a-1\ge-(a+1)$ and $j<q_t\le4M$, so the loss is at most $4M(a+1)$.
6. **Conclusion.** Changing a coefficient of $P_t$ by an element of $2^{p_t}\mathcal A$ changes $E_t(0)P_t(-t/\alpha)$ by valuation at least $p_t-6hD_t-4M(a+1)\ge p$. So any lifts of the residues give $F(0)$ modulo $2^p\mathcal A$. The multipliers $E_t(0)(-t/\alpha)^j$ are rationals of polynomial bit length, so one common denominator suffices. ∎

---

## 5. Packing

**Lemma 5.1 (paper Lemma D.3).** Let $r\ge2$, $a\ge1$, $L=2^a$ and $D_{\rm arr}=r^{2a}$. Suppose the network has $M$ edges, maximum degree $\le\Delta$, $M\ge2D_{\rm arr}^2\Delta$ and $M\ge LD_{\rm arr}$. One can deterministically compute, in time polynomial in $M$, the number of vertices and $r^a$:
- a labeling of the vertices by $[r]^a$;
- an assignment of each edge to a slot $\sigma\in[0,L)$ of a packet.

These satisfy:
1. There are at most $3M/L$ packets. Each edge occupies the same slot at both ends, and unused slots have no ports.
2. In every packet and every coordinate $c$, all edges have the same unordered pair of (row symbol, column symbol) in coordinate $c$. When the two symbols differ, order them as $(\alpha_c,\delta_c)$: the edge in slot $\sigma$ has row symbol $\alpha_c$ iff bit $c$ of $\sigma$ is $0$.
3. For every coordinate $c$ and symbol $g$, at most $3M/r$ edges have their row end at a vertex with symbol $g$ in coordinate $c$, and at most $3M/r$ have their column end at such a vertex.

*Proof.*
1. **Balanced labels.** For a pair of labels $\gamma\in([r]^a)^2$, let $X_\gamma$ be the number of edges whose (row, column) ends carry $\gamma$. For independent uniform labels,
$$\mathbb E\sum_\gamma(X_\gamma-M/D_{\rm arr})^2=\sum_{(e,e')}\big(\Pr[\gamma_e=\gamma_{e'}]-1/D_{\rm arr}\big)\le2M\Delta.$$
   Ordered pairs of edges with disjoint endpoints contribute $0$. At most $2M\Delta$ ordered pairs share an endpoint (including $e=e'$), and each contributes at most $1$.
2. **Derandomization.** Fix the labels vertex by vertex, keeping the conditional expectation from increasing. It is computable exactly, because each pair probability involves at most four vertices. Then every load satisfies $X_\gamma\le M/D_{\rm arr}+\sqrt{2M\Delta}$.
3. **Padding.** Pad every label class to $H_{\rm pad}$ slots, the maximum load rounded up to a multiple of $L$. Then $H_{\rm pad}\le M/D_{\rm arr}+\sqrt{2M\Delta}+L\le3M/D_{\rm arr}$ by the two hypotheses.
4. **Packets.** Fix an unordered pattern $\pi$ (a symbol pair per coordinate) with $a'$ unequal coordinates; its $2^{a'}$ ordered labels are the orientations. A packet takes $2^{a-a'}$ slots from each orientation. The bits of $\sigma$ at unequal coordinates give the orientation, and the bits at equal coordinates give the copy. This partition works because $L\mid H_{\rm pad}$, and it gives (1)–(2): there are $H_{\rm pad}D_{\rm arr}/L\le3M/L$ packets.
5. **(3).** In a fixed coordinate, symbol $g$ is the row symbol of exactly $D_{\rm arr}/r$ labels, which carry $H_{\rm pad}D_{\rm arr}/r\le3M/r$ slots. Columns are the same. ∎

---

## 6. Groups, modes and children

Fix a packed network, a node $t\in S\setminus\{0\}$ with $i=v_2(t)$ and $i+1\le a$, a value $\zeta\in\mathcal O$, a precision $P$, and $T=t+\alpha\zeta$.
- An end of an edge is in **group** $g\in[r]$ if its vertex has symbol $g$ in coordinate $i$.
- An entry $W_{xy}(T)$ is **local** if row end $x$ and column end $y$ are in the same group, and **crossing** otherwise.
- For $y_0\in[0,2^{i+1})$, the **tile** of $y_0$ in a packet is the set of column slots $y=y_0+2^{i+1}j$ ($j\ge0$). Bit $i$ is constant on a tile.

**Lemma 6.1.**
1. In a packet whose coordinate-$i$ symbols differ, row slot $x$ and column slot $y$ are in different groups iff bit $i$ of $x$ equals bit $i$ of $y$. In a packet whose symbols agree, all entries are local.
2. Every nonzero entry of $W(t)$ is local.

*Proof.*
1. By Lemma 5.1(2), the row symbol at slot $x$ is $\alpha_i$ iff bit $i$ of $x$ is $0$, and the column symbol at slot $y$ is $\delta_i$ iff bit $i$ of $y$ is $0$. So the groups agree iff the bits differ.
2. $W_{xy}(t)\ne0$ forces $y=x+t$ (Lemma 3.1(2)). Adding $t$ leaves the bits below $i$ unchanged and flips bit $i$. ∎

**Lemma 6.2 (the tile expansion; paper (D.12)).** Let $x$ be a row slot and $y_0$ a tile in the same packet with $\mathrm{bit}_i(y_0)=\mathrm{bit}_i(x)$, so the tile is crossing for $x$ by Lemma 6.1. Put $\theta=T-(y_0-x)$. Then:
- $v_2(t-y_0+x)\le i$ and $v_2(\theta)=v_2(t-y_0+x)$;
- for every column $y=y_0+2^{i+1}j$ of the tile,
$$W_{xy}(T)\equiv\sum_{0\le b\le P-2}\kappa_{x,y_0,b}\,j^b\pmod{2^P},\qquad\kappa_{x,y_0,b}=W_{y_0-x}(T)\Big(\frac{2^{i+1}}\theta\Big)^b;$$
- $v_2(\kappa_{x,y_0,b})\ge(a+1-i)+b\ge b+1$.

*Proof.*
- **Valuation.** If $v_2(t-y_0+x)\ge i+1$, then $y_0\equiv x+t\equiv x+2^i\pmod{2^{i+1}}$, so bit $i$ of $y_0$ differs from that of $x$; this is a contradiction. Since $v_2(\alpha\zeta)\ge a+1>i$, $v_2(\theta)=v_2(t-y_0+x)$.
- **Geometric series.** $W_{xy}(T)=\frac{C(T)}{C'(0)(\theta-2^{i+1}j)}=W_{y_0-x}(T)\cdot\frac1{1-\rho j}$ with $\rho=2^{i+1}/\theta$ and $v_2(\rho)\ge1$. Here $|y_0-x|<L$ because $i+1\le a$.
- **Truncation.** By Lemma 3.1(3), with $d=y_0-x\ne t$, $v_2(W_{y_0-x}(T))\ge a+1-i$. Truncating the geometric series after $b=P-2$ leaves a remainder of valuation $\ge(a+1-i)+(P-1)\ge P$. ∎

**Modes.** For each packet whose coordinate-$i$ symbols differ, each tile $y_0$ and each $0\le b\le P-2$, introduce a *mode* $Z$:
- its **weight** is $w_Z=b+1$;
- its **linear form** is $L_Z(X)=\sum_jj^bX_{y_0+2^{i+1}j}$, over the real columns of the tile ($0^0=1$);
- its **source** is the group of the tile's columns;
- its **target** is the other symbol of the packet's pair, which is the group of the rows for which the tile is crossing.

By Lemma 6.2, for every row end $x$,
$$(W(T)X)_x\equiv\lambda_x(X)+\sum_Z\kappa_{x,Z}L_Z(X)\pmod{2^P},$$
where:
- $\lambda_x$ is the sum of the local entries in $x$'s row;
- $Z$ runs over the modes of $x$'s packet whose tile is crossing for $x$; all of them have target equal to $x$'s group;
- $v_2(\kappa_{x,Z})\ge w_Z$.

**Lemma 6.3 (the child identity; paper (D.16)).** For a tuple $\beta$ of nonnegative exponents on the modes and a group $g$, let the **child** $N_{g,\beta}$ be the following network:
- **One side:** every vertex of group $g$ (from both sides of the network), with its ports and tensor.
- **Other side:** one vertex per packet, joined by one edge to each port of group $g$ in that packet.
- **The packet vertex's tensor:** with $Y,X$ the bits of the group-$g$ row and column ports of packet $\mathcal P$, and $\beta_{\mathcal P,g}$ the restriction of $\beta$ to the modes of $\mathcal P$ with target $g$,
$$\Pi_{\mathcal P,g,\beta}(Y,X)=\big[Z^{\beta_{\mathcal P,g}}\big]\prod_{\substack{x\in\mathcal P\\\text{group }g}}l_{Y_x}\Big(\lambda_x(X)+\sum_Z\kappa_{x,Z}Z\Big)\cdot\prod_{\substack{Z\text{ of }\mathcal P\\\text{source }g}}L_Z(X)^{\beta_Z}.$$

Then
$$F(t+\alpha\zeta)\equiv\sum_{\beta:\ \sum_Zw_Z\beta_Z<P}\ \prod_{g\in[r]}\mathrm{ct}(N_{g,\beta})\pmod{2^P}.$$

*Proof.*
1. **Expanded rows.** By the display before the lemma, $F(T)\equiv\sum_{X,Y}A(Y)B(X)\,\Phi(Y,X;Z)|_{Z=L(X)}\pmod{2^P}$, where $\Phi=\prod_xl_{Y_x}(\lambda_x(X)+\sum_Z\kappa_{x,Z}Z)$ is a polynomial in formal modes $Z$.
2. **Weights.** A coefficient of $Z^\beta$ in $\Phi$ is a sum of products in which each occurrence of $Z$ brings a factor $\kappa_{x,Z}$ of valuation $\ge w_Z$. So it has valuation $\ge\sum w_Z\beta_Z$. Hence
$$\Phi|_{Z=L(X)}\equiv\sum_{\sum w\beta<P}[Z^\beta]\Phi\cdot\prod_ZL_Z(X)^{\beta_Z}\pmod{2^P}.$$
3. **Factorization.** The rows of group $g$ involve only local columns, which are of group $g$, and modes with target $g$ in their own packet. So $[Z^\beta]\Phi$ is the product over $g$ and packets $\mathcal P$ of the extracted coefficients in $\Pi_{\mathcal P,g,\beta}$. Each $L_Z$ reads only columns of $Z$'s source group in $Z$'s packet. $A$ and $B$ factor over groups because every vertex lies in one group. So for fixed $\beta$ the summand is a product over $g$ of functions of the group-$g$ bits, and the sum over $(X,Y)$ is the product of the children's contractions. ∎

**Lemma 6.4 (sizes).** In the setting of Lemma 6.3:
1. $N_{g,\beta}$ has at most $6M/r$ edges. Old vertices keep their degrees, and a packet vertex has at most $2L$ edges.
2. There are at most $(3M/L)2^{i+1}=6M/D_t$ modes of each weight, where $D_t=L/2^i$.
3. The number of $\beta$ with $\sum w_Z\beta_Z<P$ is at most $\exp\big(2\pi\sqrt{MP/D_t}\big)$.
4. An entry of $\Pi_{\mathcal P,g,\beta}$ is computed with at most $(LP)^L\,\mathrm{poly}(P,e,\ell)$ bit operations.

*Proof.*
1. Lemma 5.1(3) bounds the row ends and the column ends of group $g$ by $3M/r$ each.
2. Lemma 5.1(1): at most $3M/L$ packets, each with $2^{i+1}$ tiles.
3. **Rankin's bound.** With $c=6M/D_t$ modes per weight and any $\theta>0$:
$$\#\{\beta\}\le e^{\theta P}\prod_Z(1-e^{-\theta w_Z})^{-1}\le\exp\Big(\theta P+c\sum_{w\ge1}-\ln(1-e^{-\theta w})\Big).$$
   Now $\sum_w-\ln(1-e^{-\theta w})=\sum_{h\ge1}\frac1{h(e^{\theta h}-1)}\le\frac{\pi^2}{6\theta}$. Taking $\theta=\pi\sqrt{c/(6P)}$ gives $\exp(2\pi\sqrt{cP/6})$.
4. Each of the $\le L$ row factors has at most $1+2^iP\le LP$ terms. Expand, keep the terms with the required exponents, and compute the powers $L_Z^{\beta_Z}$ by repeated squaring. ∎

---

## 7. The recursion

**Parameters.** Given $s$ and $C$, put
- $\lambda=\lceil\log s\rceil$, $\Delta=\lceil s^{9/10}\rceil$, $H=\lceil s/\lfloor\sqrt\lambda\rfloor\rceil$;
- $r=64$, $\varepsilon=2^{-17}$, and $e=4\lambda$;
- $u_d=\lceil\varepsilon\lambda/2^d\rceil$, $a_d=2\sum_{j\le d}u_j$ ($a_{-1}=0$), $L_d=2^{a_d}$ and $D_{{\rm arr},d}=64^{2a_d}$, for each depth $d\ge0$.

**Algorithm** $\mathrm{CONTRACT}(N,p,d)$. It returns the contraction of $N$ (with $M$ edges and entries in $\mathcal A_p$) modulo $2^p$.
1. **Leaf.** If $M\le H$, sum over all $2^M$ assignments.
2. **Otherwise** let $a=a_d$, $u=u_d$ and $L=2^a$:
   - Pack $N$ by Lemma 5.1 and lift the entries canonically.
   - For each $t\in S\setminus\{0\}$ with $q_t>0$ (Lemma 4.1) and each of $M(2L-1)+1$ binary-coordinate points $\zeta\in\mathcal O$, compute $F(t+\alpha\zeta)\bmod2^{p_t}$ by Lemma 6.3. Each child contraction comes from $\mathrm{CONTRACT}(N_{g,\beta},p_t,d+1)$.
   - Interpolate $F(t+\alpha z)\bmod2^{p_t}$ (Lemma 2.1(3), Corollary 3.2).
   - Recover $F(0)\bmod2^p$ (Lemma 4.1).
   
   A packet tensor is described by its packet's data and $(g,\beta,t,\zeta,p_t)$. Retained tensors keep their original descriptions and creation precisions.

**Lemma 7.1 (parameter inequalities).** There is $\lambda_0(C)$ such that for $\lambda\ge\lambda_0$, every depth with $8^d<C\lfloor\sqrt\lambda\rfloor$ satisfies:
- **(I1)** $1\le u_d$, $2u_d\le a_d$, and $u_d\le H$;
- **(I2)** $H\ge2D_{{\rm arr},d}^2\Delta$, $H\ge L_dD_{{\rm arr},d}$ and $2L_d\le\Delta$;
- **(A1)** $2^{u_d}\ge20.5\,C(d+1)\sqrt\lambda\,u_d$;
- **(A2)** $2^{u_d}\ge4(a_d+1)u_d$;
- **(I3)** $Cs(2L_d-1)+1\le2^{4\lambda}$, and $L_d\log(L_d\cdot41Cs(d+2)2^{a_d})\le s^{9/10}$.

The script computes $\log_2\lambda_0=25.6$, $26.8$ and $26.9$ for $C=1,4,16$.

*Proof.*
- **Bounds on $u_d$ and $a_d$.** $8^d<C\sqrt\lambda$ gives $2^d<(C\sqrt\lambda)^{1/3}$, so $u_d\ge\varepsilon\lambda^{5/6}/C^{1/3}$. Also $a_d\le2\sum_j(\varepsilon\lambda/2^j+1)\le4\varepsilon\lambda+2(d+1)$.
- **(A1), (A2).** Both have the form $2^u\ge Ku$ with $\log K=O_C(\log\lambda)$, and $2^u/u$ increases for $u\ge2$. So they hold once $\varepsilon\lambda^{5/6}/C^{1/3}\ge2\log K$.
- **(I2).** In logarithms it reads $\lambda-O(\log\lambda)\ge24a_d+0.9\lambda+O(1)$, which holds because $24\cdot4\varepsilon<0.1$. The others are similar.
- **The thresholds.** Inside each interval between two depth transitions every condition only improves as $\lambda$ grows. The script settles each interval up to depth 60, and also checks a geometric grid up to $\log_2\lambda\approx260$. ∎

**Lemma 7.2 (invariants).** Run $\mathrm{CONTRACT}$ on a network as in Theorem 10.2 with $\lambda\ge\lambda_0$, at depth $0$ and precision $P_0$. Then every call at depth $d$ satisfies:
- **(V1)** $M\le Cs\,8^{-d}$. So a non-leaf call has $8^d<C\lfloor\sqrt\lambda\rfloor$, hence $d<\log_8(C\sqrt\lambda)$.
- **(V2)** The maximum degree is at most $\Delta$, and there are at most $s^C+3Cs(d+1)$ vertices.
- **(V3)** The input precision satisfies $p\le41Cs(d+1)2^{a_{d-1}}$.

*Proof.*
- **(V1)** A child has at most $6M/64<M/8$ edges (Lemma 6.4(1)). A non-leaf call has $M>H\ge s/\lfloor\sqrt\lambda\rfloor$.
- **(V2)** Lemma 6.4(1) with (I2): $2L\le\Delta$. A call adds at most $3M/L$ packet vertices.
- **(V3)** Induction. With $h\le4M/u\le4Cs$, $D_t\le2^{a_d}$ and $a+1\le2^a$:
$$p_t\le p+32Cs\,2^{a_d}+8Cs(a_d+1)\le p+40Cs\,2^{a_d}.$$
  Summing over the depths, $p\le Cs\big(1+40\sum_{d'<d}2^{a_{d'}}\big)\le41Cs(d+1)2^{a_{d-1}}$. ∎

**Lemma 7.3 (precision at a node; paper (D.11), (D.15)).** At a non-leaf call at depth $d$ and a node $t$ with $q_t>0$:
$$\frac{p_t}{D_t}\le\frac{34M}{u_d},$$
so at most $\exp(2\pi\sqrt{34}\,M/\sqrt{u_d})\le e^{37M/\sqrt{u_d}}$ mode tuples occur.

*Proof.*
- **Size of $D_t$.** $q_t>0$ forces $i\le u_d-1$, so $D_t\ge2^{a_d-u_d+1}=2^{a_{d-1}+u_d+1}$.
- **Inherited precision.** By (V3) and $s<M\lfloor\sqrt\lambda\rfloor$, which follows from $M>H$:
$$\frac p{D_t}\le\frac{20.5\,Cs(d+1)}{2^{u_d}}\le\frac{20.5\,C(d+1)\sqrt\lambda\,M}{2^{u_d}}\le\frac M{u_d}\quad\text{by (A1)}.$$
- **The other two terms.** $8h\le32M/u_d$, and $8M(a_d+1)/D_t\le4M(a_d+1)/2^{u_d}\le M/u_d$ by (A2).
- **Mode tuples.** Lemma 6.4(3). ∎

**Proof of Theorem 10.2.** Embed the entries in $\mathcal A_{P_0}$ and run $\mathrm{CONTRACT}(N,P_0,0)$. Inputs with $\lambda<\lambda_0(C)$ are bounded in size and are summed directly.
1. **Correctness.** By induction from the leaves.
   - **Children.** A child's retained tensors carry their creation-time integer representatives, which are the lifts used in the parent's $F$. Its packet tensors are exact modulo $2^{p_t}$. So the child returns the factor $\mathrm{ct}(N_{g,\beta})\bmod2^{p_t}$ of Lemma 6.3.
   - **Jets.** By Corollary 3.2 and (I3), $F(t+\alpha z)\bmod2^{p_t}$ has degree $\le M(2L-1)<2^e$, so interpolation at the binary points recovers it (Lemma 2.1(3)).
   - **Recovery.** Lemma 4.1 applies: $u\le a$ and $M\ge u$ by (I1). It returns $F(0)\bmod2^p$, which is the call's answer by Corollary 3.2.
   - **The answer.** The final value lies in $(\mathbb Z/2^{P_0})[V]/(V^\ell)$, since it is a polynomial in the original entries.
2. **Number of calls.** A non-leaf call at depth $d$ has at most $2L_d\cdot2^{4\lambda}\cdot64\cdot e^{37M_d/\sqrt{u_d}}$ children (Lemma 7.3), and $2L_d\le2^\lambda$.
   - With $M_d\le Cs8^{-d}$ and $u_d\ge\varepsilon\lambda/2^d$, the exponents sum to
$$\sum_d\frac{37Cs\,8^{-d}}{\sqrt{\varepsilon\lambda/2^d}}=\frac{37Cs}{\sqrt{\varepsilon\lambda}}\sum_d\Big(\frac{\sqrt2}8\Big)^d\le\frac{45\,Cs}{\sqrt{\varepsilon\lambda}} .$$
   - There are $O(\log\lambda)$ depths, so the remaining factors contribute $2^{O(\lambda\log\lambda)}$.
   - Hence there are $\exp(O_C(s/\sqrt\lambda))$ calls.
3. **Work per call.**
   - **Non-leaf calls.** The packing, the $W$ entries, mode coefficients, interpolation and recovery cost $\mathrm{poly}(s,L_{\rm d})$ besides the children. All precisions are polynomial in $s$: $2^{a_d}\le2^{4\varepsilon\lambda+2(d+1)}$ and (V3).
   - **Leaves.** A leaf costs $2^H$ times polynomially many entry evaluations. An original entry costs $(s+L_{\rm d})^b2^{As^{9/10}}$, and a packet entry costs $2^{s^{9/10}}$ by Lemma 6.4(4) and (I3).
   - Since $H\le2s/\sqrt\lambda$, the total time is $(s+L_{\rm d})^{O(1)}2^{O(s/\sqrt{\log s})}$.
4. **Space.** Traverse depth first. Each of the $O(\log\lambda)$ levels stores its packing, its jets, the current indices $(t,\zeta,\beta,g)$ and the new tensors' descriptions: polynomial size. ∎

---

## 8. Fixed clause density (paper D.2)

**Proposition 8.1.** #$k$-SAT with $m\le Dn$ clauses takes $(n+m)^{O_{k,D}(1)}2^{O_{k,D}(n/\sqrt{\log n})}$ time and polynomial space.

*Proof.* Normalize; an empty clause gives $0$. Build the network:
- **Clauses.** Each clause is a vertex with its satisfaction table, and one port per distinct variable.
- **Variables.** A variable with $d\ge1$ ports gets a gadget that forces its bits to agree:
  - $d=1$: the unary tensor $(1,1)$;
  - $d=2$: a binary equality vertex;
  - $d\ge3$: a tree of $d-2$ trivalent equality vertices with $d$ external ports. A common value on the ports extends uniquely inside the tree.
- **Bipartition.** Subdivide every edge by a binary equality vertex; old and new vertices form the two sides.

The contraction counts the satisfying assignments of the occurring variables, so multiply it by $2^{\#\text{unused}}$. The network has $O_k(n+m)$ edges and vertices, and maximum degree $\le\max(k,3)$. Apply Theorem 10.2 with $s=n+m$, $\ell=1$ and $P_0=n+1$; the least nonnegative residue is $\#F\le2^n$. For $m\le Dn$, $s\le(1+D)n$. ∎

---

## 9. Arbitrary density (paper D.3)

**Lemma 9.1 (paper Lemma D.2).** Fix $k,c$. Let $P_1,\dots,P_m\in\mathbb F_2[x_1,\dots,x_n]$, $m\le cn$, each given as at most $\lceil\log n\rceil^c$ monomials of degree $\le k$. Their common zeros are counted in $n^{O(1)}2^{O_{k,c}(n\log\log n/\sqrt{\log n})}$ time and polynomial space.

*Proof.* Let $\lambda=\lceil\log n\rceil$ and $\nu=\lceil\log\lambda\rceil$.
1. **Branching.** Branch on the variables lying in more than $D_*=\max(1,\lceil N_0/n\rceil)\lceil\sqrt\lambda\rceil$ scopes, where $N_0$ is the incidence count, as in Part II Theorem 3.1 step 2.
   - There are at most $2^{n/\sqrt\lambda}$ branches.
   - Free scopes have $|I_i|\le k\lambda^c$, and each free variable lies in at most $D_*=\lambda^{O(1)}$ scopes.
   - Constant equations are handled directly.
2. **The incidence polynomial.** Give each incidence $(i,v)$, $v\in I_i$, its own indeterminate $z_{iv}$, and let $\tilde L_i(z_i)=\sum_{M\in\mathcal M_i}\prod_{v\in M}z_{iv}$. Over $\mathcal B=(\mathbb Z/2^{n+1})[V]/(V^\ell)$ with $\ell=n+1$, let:
   - $\bar{\mathcal B}_i$ be the multilinearization in $z_i$ of $\sum_{j<\ell}e_j(\tilde L_i)h(\tilde L_i)^jV^j$, where the $e_j$ are the digits of Part II Lemma 2.1;
   - $F=\prod_i\bar{\mathcal B}_i$, a multilinear polynomial in $z$, since different tests have disjoint indeterminates.
   
   Multilinearization does not raise degrees. So, as in Part II Lemma 2.2(2), every monomial of $F$ has degree $\le K=k(m+2(\ell-1))=O(n)$.
3. **The congruence.** At $z_{iv}=x_v$, $\tilde L_i$ counts the monomials of $P_i$ equal to $1$. By Part II Lemma 2.2(1),
$$N_{\rm branch}\equiv\mathcal S_\ell\Big(\sum_xF(z_{iv}=x_v)\Big)\pmod{2^\ell}.$$
4. **Monomials.** A monomial of $F$ is a tuple $\mathbf S=(S_i)$ with $S_i\subseteq I_i$ and $\sum|S_i|\le K$. Its coefficient is $\prod_ic_i(S_i)$ with $c_i(S)=[z_{i,S}]\bar{\mathcal B}_i$, and $\sum_x\prod_{(i,v)\in\mathbf S}x_v=2^{n'-|V(\mathbf S)|}$, where $V(\mathbf S)$ is the set of variables touched. By Möbius inversion, $c_i(S)=\sum_{T\subseteq S}(-1)^{|S\setminus T|}\mathcal B_i(\mathbf 1_T,V)$.
5. **Networks.** Split the tests and the free variables into $J=\lceil n^{1/4}\rceil$ groups each, and let $N_{gh}$ be the number of incidences between test group $g$ and variable group $h$. For each array $(d_{gh})$ with $\sum d_{gh}\le K$ (at most $(K+1)^{J^2}=2^{O(\sqrt n\log n)}$ arrays), build a network:
   - **Vertices and wires.** It has $2J$ vertices, with $\lceil\log\binom{N_{gh}}{d_{gh}}\rceil$ wires between $g$ and $h$.
   - **Indexing.** The wire bits index a $d_{gh}$-subset of those incidences in lexicographic order; an out-of-range index gives entry $0$.
   - **Test group $g$:** entry $\prod_{i\in g}c_i(S_i)$.
   - **Variable group $h$:** entry $2^{\#\{v\in h\text{ untouched}\}}$.
   
   Each monomial of $F$ is counted once, in exactly one array and one valid wire assignment. So the sum of the contractions is $\sum_xF(z=x)$.
6. **Sizes.**
   - **Wires.** $M\le J^2+\log\prod\binom{N_{gh}}{d_{gh}}\le J^2+\log\binom{N+K}K\le J^2+K\log\frac{e(N+K)}K=O(n\nu)$, because $N\le n\lambda^{O(1)}$ and $K=\Theta(n)$.
   - **Degrees.** A vertex's degree is at most its incidences, $n^{3/4}\lambda^{O(1)}$.
   - **Entry time.** A test-group entry takes $2^{k\lambda^c}\mathrm{poly}(n)$ time.
7. **Contraction.** Apply Theorem 10.2 with $s=n\nu$ and $P_0=\ell=n+1$. Each contraction costs $n^{O(1)}2^{O(n\nu/\sqrt\lambda)}$. With $2^{n/\sqrt\lambda+O(\sqrt n\log n)}$ branches and arrays, the total is $2^{O(n\nu/\sqrt\lambda)}$. Branch counts are added in an $(n+1)$-bit accumulator. ∎

**Proof of Theorem 10.1.**
- **Fixed density.** This is Proposition 8.1.
- **Arbitrary density.**
  - Sparsify (Part I Prop. 2.4, Lemma 2.5) with $\sigma=1/\lceil\sqrt\lambda\rceil$. This gives $\mathrm{poly}(n)2^{n/\sqrt\lambda}$ leaves with occurrence bound $B=\lambda^{O_k(1)}$.
  - Count each leaf by Part II Lemmas 4.1–4.2, with Lemma 9.1 as the counter. That is $2^{\lambda^{O_k(1)}}$ calls, each costing $2^{O_k(n\nu/\sqrt\lambda)}$.
- **One algorithm for both.** Use the clause network when $m\le n\nu$, so that $s=O(n\log\log n)$, and the parity route otherwise. ∎

---

## 10. Relation to the paper, and status

| This file | Paper | Notes |
|---|---|---|
| Lemma 1.2 | D.4.3, (D.8) | |
| Lemma 2.1 | D.4.1 | |
| Lemma 3.1, Cor. 3.2 | (D.7), Lemma D.4 | |
| Lemma 4.1 | Lemma D.5, (D.9), (D.10) | |
| Lemma 5.1 | Lemma D.3 | |
| Lemmas 6.1–6.4 | D.4.5–D.4.6, (D.12)–(D.16) | coefficient valuations $\ge(a+1-i)+b$, sharper than $b+1$ (not needed) |
| Lemma 7.1 | (D.4)–(D.6) | **(A1)–(A2) replace $u\ge\lambda^{4/5}$**: threshold $\log_2\lambda\approx26$ instead of $\approx500$ |
| Lemmas 7.2–7.3 | (D.11), (D.15) | explicit constants: $p\le41Cs(d+1)2^{a_{d-1}}$, $p_t/D_t\le34M/u$, $e^{37M/\sqrt u}$ tuples |
| Thm 10.2 | Thm D.1 | |
| Prop. 8.1, Lemma 9.1, Thm 10.1 | D.2, Lemma D.2, Thm 1.1 | |

**Status.** With Parts I–III, every complexity bound in the paper now has a complete prose proof in this repository:
- $2^{o(n)}$ (Part I);
- $2^{O(\rho(n))}$, the formally verified bound (Part II);
- $2^{O(n\log\log n/\sqrt{\log n})}$ and $2^{O(n/\sqrt{\log n})}$ at fixed density (this part).

App. D is not formalized. Its new mechanism, 2-adic modes with tapered recovery, is now checked end to end on toy networks: the child identity, the tile expansion, and one recursion level at the honest precision.

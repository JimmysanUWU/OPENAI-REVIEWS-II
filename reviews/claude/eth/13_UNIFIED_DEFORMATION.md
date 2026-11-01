# 13. One deformation theorem for Parts I–III: interpolation through local nodes

*Claude, 28 October 2026. Roadmap target `T.ETH.2`; nodes `ETH.U.*` in [`ROADMAP.md`](../../../ROADMAP.md).*

The ETH paper counts satisfying assignments with two mechanisms that look unrelated:
- **Parts I–II:** swap interpolation with Taylor jets and measurements ([08](08_GAPLESS_PROOF.md), [09](09_GAPLESS_SHARP_EXPONENT.md)).
- **Part III (App. D):** a Cauchy kernel with 2-adic modes and tapered Hermite recovery ([10](10_GAPLESS_TENSOR_CONTRACTION.md)).

This file shows they are two instances of one theorem. It isolates the four general lemmas that carry all the mathematics, and computes the cost floor of each instance from one cost functional.

> **Theorem 4.2 (interpolation through local nodes), informally.** To compute $Z=Z(\iota)$, the partition function of a readout system, choose:
> - a polynomial deformation $\rho(T)$ of its readout with $\rho(0)=\iota$;
> - nodes where the deformed readout is *local*, i.e. where the problem splits into groups after fixing a few measurements.
>
> Near each node, write the readout as a local part plus *modes* that are small in a filtration (formal $z$-adic, or $\varpi$-adic). Then $H(T)=Z(\rho(T))$ near the node is a sum, over the correction monomials of weight below the precision, of products of group sums (Lemmas 2.3–2.4). Hermite interpolation with an explicit valuation loss recovers $Z=H(0)$ (Lemma 3.1).
>
> The cost is
> $$\sum_{\text{nodes}}(\#\text{points})\cdot(\#\text{measurement values})\cdot(\#\text{correction monomials})\cdot\sum_{\text{groups}}(\text{group cost}).$$

> **Corollaries 4.3–4.4.**
> - **Part I and Part II's deformation counting** is the instance with formal neighbourhoods at $Q$ unit-distance nodes, Boolean coordinate modes (all $n$ variables, weight 1), measurements, and brute-force groups.
> - **Part III** is the instance with 2-adic neighbourhoods at a complete residue system mod $2^{a+1}$, structured tile modes, no measurements, tapered recovery, and recursive groups.

> **Theorem 5.4 (new: the floor of the 2-adic architecture).** In App. D's architecture with any depth profile of Hermite multiplicities, some node needs at least $2^{\frac18\min(p-1,\ \sqrt{3/16}\,M/\sqrt a)}$ correction tuples, where:
> - $M$ is the number of wires;
> - $a\le\log_2M$ is the number of 2-adic levels;
> - $p$ is the output precision.
>
> At the top level this is $2^{\Omega(M/\sqrt{\log M})}$, matching the upper bound of [10](10_GAPLESS_TENSOR_CONTRACTION.md) Theorem 10.2. With [07](07_TOOLKIT.md) E4 for Part I, **both $\sqrt{\log}$ barriers are now theorems about their architectures.**

**What is new here.** All of these are complete derivations in this file:
- the abstract readout-system framework and Theorem 4.2;
- the filtered node lemma in two variants (Lemmas 2.3–2.4);
- Hermite recovery with an explicit loss formula and its sharpness (Lemma 3.1);
- two-sided correction counts (Lemma 4.1);
- the depth-free node cost (Lemma 5.2), the profile lemma (Lemma 5.3), and the floor (Theorem 5.4);
- the design table (§6).

The instances themselves were proved in 08–10 and are not re-proved; §2 and §4 show that their proofs are these lemmas.

**Checks.** [`scripts/verify_unified_deformation.py`](scripts/verify_unified_deformation.py) ([`results/unified_deformation.json`](results/unified_deformation.json), ~30 s). The general lemmas are implemented once, independently of the instance code.

| Part | What it checks | Result |
|---|---|---|
| A | Lemma 2.3 (targeted modes, $p$-adic filtration, with measurements): 15 random readout systems over $\mathbb Z/2^P$ and $\mathbb Z/3^P$ | all agree with brute force |
| A | Witness systems for the three hypotheses: a coefficient one valuation too low, truncation one weight too early, a mode credited to the wrong group | each changes the answer, as it must |
| B | Lemma 2.4 (Boolean modes, $z$-adic, weights 1–2, with measurements): 12 random systems over $(\mathbb Z/p^\mu)[z]/(z^J)$ | all agree with brute force |
| B | Negative controls: early truncation, separation without measurements | detected 10/12 each |
| C | Lemma 3.1: unit nodes (Part I's setting) | zero loss, exact recovery |
| C | Lemma 3.1: a non-unit difference | loses 2 digits |
| C | Lemma 3.1: complete residue systems with six depth profiles (Cauchy-family polynomials) | loss formula exact; recovery at the rule's precision exact |
| C | Sharpness: one bit less at one node | fails 6/6 |
| D | Lemma 4.1: exact counts against both bounds ($c\le34$, $P\le144$) | bounds hold |
| E | Lemma 5.3 on 19,500 random and structured profiles | holds |
| E | The exact optimum of the depth-profile problem, and the paper's taper | both $\to6M/a$ |
| F | The generic Boolean node evaluator on genuine Part I deformations (08) | reproduces the Taylor jets at 5/5 nodes |
| F | The generic targeted node evaluator on genuine packed Cauchy-kernel networks (10) | reproduces $F(t+\alpha\zeta)\bmod2^P$ at 8/8 nodes |

---

## 0. Readout systems

**Definition 0.1 (readout system).** Let $\mathcal A$ be a commutative ring. A *readout system* over $\mathcal A$ consists of:
- a finite set $V$ of Boolean variables;
- a finite set $\mathcal F$ of factors. Each $f$ has a finite set of *ports* $\Pi_f$ (pairwise disjoint) and a polynomial $\Phi_f\in\mathcal A[\xi_\pi:\pi\in\Pi_f]$;
- a *readout* $\rho=(\rho_\pi)_{\pi\in\Pi}$, $\Pi=\bigsqcup_f\Pi_f$, with $\rho_\pi:\{0,1\}^V\to\mathcal A$.

Its value is
$$Z(\rho)=\sum_{x\in\{0,1\}^V}\ \prod_{f\in\mathcal F}\Phi_f\big((\rho_\pi(x))_{\pi\in\Pi_f}\big).$$
The *incidence readout* $\iota$ gives each port a variable $v(\pi)$ and reads $\iota_\pi(x)=x_{v(\pi)}$. The quantity to compute is $Z=Z(\iota)$.

**Example 0.2.**
- **(a) CNF counting (Parts I–II).**
  - $V$ is the set of variables, $\mathcal F$ the clauses, and the ports are the incidences $(f,v)$.
  - $\Phi_f$ is the clause polynomial $\Phi_{C_f}$, so $Z(\iota)=\#F$.
  - In Part II, $\Phi_f\in\mathcal A[\xi]$ for the free algebra $\mathcal A=(\mathbb Z/p^\mu)[V]/(V^\ell)$, carrying the truncated parity lift ([09](09_GAPLESS_SHARP_EXPONENT.md) §1–§3). The same holds for weighted factor graphs ([07](07_TOOLKIT.md) E3) and for $\mathbb F_q$ systems ([12](12_FP_COUNTING.md)).
- **(b) Tensor networks (Part III).**
  - $V=\{X_e,Y_e:e\in E\}$ are the column and row ends of the wires.
  - The factors are the vertex tensors and one *wire factor* per edge. A vertex tensor's polynomial is the multilinear extension of its table, and its ports read the $X$ or $Y$ bits. The wire factor has ports reading $Y_e$ and $X_e$, and $\Phi(\xi,\eta)=(1-\xi)(1-\eta)+\xi\eta$, which is $[\xi=\eta]$ on bits.
  - Then $Z(\iota)=\mathrm{ct}(N)$ ([10](10_GAPLESS_TENSOR_CONTRACTION.md) Lemma 1.2). Only the $X$-port of each wire factor will be deformed.

**Definition 0.3 (deformation).** Let $\mathcal O\subseteq K$ be rings and $\mathcal A$ a free $\mathcal O$-algebra of finite rank, with $\mathcal A_K=K\otimes\mathcal A$. A *deformation of degree $\le R'$* is a readout $\rho(T)$ such that:
- each deformed port $\pi\in\Pi^{\rm def}$ reads an affine form $\rho_\pi(T)(x)=\rho_{\pi,0}(T)+\sum_v\rho_{\pi,v}(T)x_v$ with coefficients in $K[T]$ of degree $\le R'$;
- $\rho_\pi=\iota_\pi$ for the other ports;
- $\rho(0)=\iota$.

Let $d$ be the total degree of $\prod_f\Phi_f$ in the variables $\xi_\pi$, $\pi\in\Pi^{\rm def}$, **computed in $\mathcal A[\xi]$**. Then
$$H(T):=Z(\rho(T))\in\mathcal A_K[T],\qquad\deg H\le dR',\qquad H(0)=Z .$$

*Proof.* Each summand is the image of $\prod_f\Phi_f$ under the substitution homomorphism $\xi_\pi\mapsto\rho_\pi(T)(x)$, whose values on deformed ports have degree $\le R'$. ∎

In Part I, $d=\sum_f|C_f|$. In Part II, $d$ can be far smaller than the sum of the factor degrees, because $\mathcal A$ has nilpotents ([09](09_GAPLESS_SHARP_EXPONENT.md) Theorem 1.1). In Part III, $d=M$: each wire factor is linear in its deformed port.

---

## 1. Local nodes

**Definition 1.1 (measured partition, local readout).** A *measured partition* of $V$ is a partition $V=\bigsqcup_{g\in G}V_g$ with integer matrices $C_g$ whose columns are indexed by $V_g$. The *measurements* of $x$ are $y_g(x)=C_gx_{V_g}$. Fix a factor assignment $\sigma:\mathcal F\to G$.

A readout $\lambda$ is *local* if, for every $f$ and every $\pi\in\Pi_f$, $\lambda_\pi(x)$ depends only on $x_{V_{\sigma(f)}}$ and on $y(x)=(y_h(x))_h$. A *local node* of a deformation is a parameter value at which the deformed readout is local, for some measured partition depending on the node.

In Part I, locality is stated through rational row spaces ([08](08_GAPLESS_PROOF.md) Definition 4.1), and Definition 1.1 is the consequence used there (08 Theorem 4.6 step 3). In Part III there are no measurements.

**Lemma 1.2 (separation; 08 Lemma 4.5).** Let $\lambda$ be local and $E\subseteq V$. For $y=(y_g)_g$ put
$$\mathcal X_g(y_g,E)=\{x_g\in\{0,1\}^{V_g}:C_gx_g=y_g,\ x_g=1\text{ on }E\cap V_g\}.$$
Then
$$\sum_{x:\ x_E=1}\ \prod_f\Phi_f(\lambda(x))=\sum_y\ \prod_{g\in G}\ \sum_{x_g\in\mathcal X_g(y_g,E)}\ \prod_{f\in\sigma^{-1}(g)}\Phi_f\big(\lambda(x_g\oplus\bar x_{-g})\big).$$
Here $y$ runs over the tuples with every fibre nonempty, $\bar x_h\in\mathcal X_h(y_h,E)$ is arbitrary, and $x_g\oplus\bar x_{-g}$ is the point with these blocks.

*Proof.* The map $x\mapsto(y(x),(x_{V_g})_g)$ is a bijection from $\{x:x_E=1\}$ onto the pairs with $x_g\in\mathcal X_g(y_g,E)$. By locality, a group-$g$ factor sees $x$ only through $x_{V_g}$ and $y$, so its value is unchanged when the other blocks are replaced by any representatives. The sum therefore factorizes. ∎

---

## 2. Neighbourhoods of a node and the node lemma

**Definition 2.1 (filtration, neighbourhoods).** A *filtered ring* is a commutative ring $\mathcal A$ with ideals $\mathcal A=\mathfrak a_0\supseteq\mathfrak a_1\supseteq\cdots$ such that $\mathfrak a_i\mathfrak a_j\subseteq\mathfrak a_{i+j}$. Working *at precision $P$* means working modulo $\mathfrak a_P$. Two neighbourhoods of a node $t$ occur:
- **formal:** substitute $T=t+z$ and work in $\mathcal A[z]/(z^J)$ with $\mathfrak a_i=z^i$. The precision is $J$, and $\mathfrak a_J=0$.
- **$\varpi$-adic:** substitute $T=t+\alpha\zeta$ for points $\zeta$ in the base ring, and work in $\mathcal A/\varpi^P$ with $\mathfrak a_i=\varpi^i$.

**Definition 2.2 (mode expansion).** Let $\lambda$ be a local readout with values in a filtered ring $\mathcal A$. A *mode expansion* of a readout $\rho$ at precision $P$ relative to $\lambda$ consists of a finite set $\mathcal Z$ of *modes* and coefficients, as follows.
- Each mode $Z$ has a *weight* $w_Z\ge1$, a *source* group $s(Z)$, and an integer-valued function $L_Z$ of $x_{V_{s(Z)}}$.
- The coefficients satisfy $\kappa_{\pi,Z}\in\mathfrak a_{w_Z}$.
- For all $\pi$ and all $x\in\{0,1\}^V$,
$$\rho_\pi(x)\equiv\lambda_\pi(x)+\sum_{Z}\kappa_{\pi,Z}\,L_Z(x_{V_{s(Z)}})\pmod{\mathfrak a_P}.$$

The expansion is
- **targeted** if every mode has a *target* $\tau(Z)\in G$ and $\kappa_{\pi,Z}=0$ unless $\pi\in\Pi_f$ with $\sigma(f)=\tau(Z)$;
- **Boolean** if $L_Z(x)=x_{v_Z}$ for distinct coordinates $v_Z$, with no condition on which factors read them.

**Lemma 2.3 (node lemma, targeted modes).** For a targeted mode expansion,
$$Z(\rho)\equiv\sum_y\ \sum_{\beta\in\mathbb N^{\mathcal Z}:\ \sum_Zw_Z\beta_Z<P}\ \prod_{g\in G}\mathcal Z_g(y,\beta)\pmod{\mathfrak a_P},$$
where
$$\mathcal Z_g(y,\beta)=\sum_{x_g\in\mathcal X_g(y_g,\emptyset)}\ \prod_{Z:\,s(Z)=g}L_Z(x_g)^{\beta_Z}\cdot\Big[\zeta^{\beta|_{\tau^{-1}(g)}}\Big]\prod_{f\in\sigma^{-1}(g)}\Phi_f\Big(\big(\lambda_\pi(x_g\oplus\bar x_{-g})+\textstyle\sum_{Z:\,\tau(Z)=g}\kappa_{\pi,Z}\zeta_Z\big)_\pi\Big).$$
The brackets extract a coefficient in formal variables $\zeta_Z$.

*Proof.*
1. **Precision.** $\Phi_f$ has coefficients in $\mathcal A$. So replacing each argument by one congruent to it modulo $\mathfrak a_P$ changes $\Phi_f$ by an element of $\mathfrak a_P$, and we may use $\lambda+\sum\kappa L$.
2. **Weights.** Expand $\prod_f\Phi_f(\lambda_\pi(x)+\sum_Z\kappa_{\pi,Z}\zeta_Z)=\sum_\beta c_\beta(x)\zeta^\beta$. Every occurrence of $\zeta_Z$ comes with a factor in $\mathfrak a_{w_Z}$, so $c_\beta(x)\in\mathfrak a_{\sum w\beta}$. Substituting $\zeta_Z=L_Z(x)$, an integer, gives the summand. Modulo $\mathfrak a_P$ only $\sum w\beta<P$ survives, a finite set because $w\ge1$.
3. **Factorization.** The factors of group $g$ involve only modes targeted at $g$, so $c_\beta(x)=\prod_gc^{(g)}_{\beta|_{\tau^{-1}(g)}}(x)$. Each $c^{(g)}$ depends on $x$ only through $x_{V_g}$ and $y$, by locality, and $\prod_ZL_Z^{\beta_Z}$ splits by source. Lemma 1.2 with $E=\emptyset$ sums each product over the fibres. ∎

**Lemma 2.4 (node lemma, Boolean coordinate modes).** For a Boolean mode expansion,
$$Z(\rho)\equiv\sum_y\ \sum_{\substack{E\subseteq\mathcal Z\\\sum_{Z\in E}w_Z<P}}\ \sum_{U\subseteq E}(-1)^{|E\setminus U|}\prod_{g\in G}\ \sum_{x_g\in\mathcal X_g(y_g,v_E)}\ \prod_{f\in\sigma^{-1}(g)}\Phi_f\Big(\big(\lambda_\pi(x_g\oplus\bar x_{-g})+\textstyle\sum_{Z\in U}\kappa_{\pi,Z}\big)_\pi\Big)\pmod{\mathfrak a_P},$$
where $v_E=\{v_Z:Z\in E\}$.

*Proof.*
1. **Multilinearize.** Let $p_x(\zeta)=\prod_f\Phi_f(\lambda_\pi(x)+\sum_Z\kappa_{\pi,Z}\zeta_Z)$, and let $\tilde p_x$ be its multilinearization ($\zeta_Z^k\mapsto\zeta_Z$). This agrees with $p_x$ on $\{0,1\}^{\mathcal Z}$.
2. **Weights.** The coefficient $[\zeta^E]\tilde p_x$ is a sum of coefficients of monomials $\zeta^\beta$ with support $E$, so it lies in $\mathfrak a_{\sum_{Z\in E}w_Z}$.
3. **Möbius inversion.** On the Boolean lattice, $[\zeta^E]\tilde p_x=\sum_{U\subseteq E}(-1)^{|E\setminus U|}p_x(\mathbf 1_U)$.
4. **Sum.** Hence $\sum_xp_x(x_{v_\cdot})=\sum_x\sum_E[\zeta^E]\tilde p_x\prod_{Z\in E}x_{v_Z}$. Modulo $\mathfrak a_P$ only light $E$ survive, and $\prod_{Z\in E}x_{v_Z}=[x_{v_E}=1]$.
5. **Separation.** For fixed $(E,U)$ the readout $\lambda+\sum_{Z\in U}\kappa_{\cdot,Z}$ is local (constants added), so Lemma 1.2 with $E\mapsto v_E$ applies. ∎

**Remark 2.5 (why two variants).**
- **Targeted modes** allow powers $L_Z^k$ and extract coefficients symbolically. This needs every reader of a mode to lie in one group.
- **Boolean modes** may be read by every group, which is the situation when the perturbation is a dense row. The price is the evaluation at $2^{|E|}$ points $\mathbf 1_U$ instead of symbolic extraction.

Making Boolean modes targeted by copying them once per group would cost a factor $|G|^{|E|}$, which doubles Part I's leading term.

**Corollary 2.6 (Part I's node computation).**
- **Setting.** In 08 Theorem 4.6, at node $t$, use the formal neighbourhood $\mathcal A_t=(\mathbb Z/p^\mu)[z]/(z^J)$ with $\lambda_{(f,v)}=\langle\pi P_{f,v}(t),x\rangle$.
- **Modes.** Take all $n$ coordinates as Boolean modes of weight 1, with $\kappa_{(f,v),u}=(\pi\Delta_{f,v})_u\in z\mathcal A_t$.
- **Conclusion.** Lemma 2.4 is then 08 Theorem 4.6 steps 2–3, with Lemmas 4.4–4.5 there. The budget is $|E|<J=q+1$. With general factors in a free algebra it is 09 Theorem 1.1's node step.
- **Check.** Part F1: the generic evaluator reproduces 08's Taylor jets on genuine swap deformations.

**Corollary 2.7 (Part III's child identity).**
- **Setting.** In 10 §6, at $T=t+\alpha\zeta$ with $v_2(t)=i$, use the 2-adic neighbourhood at precision $P$.
- **Groups.** They are the coordinate-$i$ symbols, with no measurements. $\lambda$ consists of the local entries $W_{xy}(T)$.
- **Modes.** One per (packet, crossing tile $y_0$, $b<P-1$):
  - weight $b+1$;
  - source: the tile's column group;
  - target: the group of the rows for which the tile crosses (10 Lemma 6.1);
  - $L_Z=\sum_jj^bX_{y_0+2^{i+1}j}$;
  - $\kappa_{x,Z}=W_{y_0-x}(T)(2^{i+1}/\theta)^b\in2^{b+1}\mathcal O$ (10 Lemma 6.2).
- **Conclusion.** Lemma 2.3 is then 10 Lemma 6.3. The child network $N_{g,\beta}$ is an encoding of the group sum $\mathcal Z_g(\beta)$.
- **Check.** Part F2: the generic evaluator reproduces $F(t+\alpha\zeta)\bmod2^P$ on genuine packed networks.

---

## 3. Recovery from the nodes

**Lemma 3.1 (Hermite recovery with valuation loss).** Let:
- $K$ be a field with a discrete valuation $v$, valuation ring $\mathcal O_v$ and uniformizer $\varpi$, and $\mathcal A$ a free $\mathcal O_v$-module of finite rank (valuations coordinatewise);
- $\mathcal N\subset K\setminus\{0\}$ be finite, with multiplicities $m_t\ge0$ and scalings $\alpha_t\in K^\times$ such that $v(\alpha_t)\ge v(t-q)$ whenever $q\ne t$ and $m_q>0$;
- $N=\sum m_t$ and $E_t(T)=\prod_{q\ne t}\big(\frac{T-q}{t-q}\big)^{m_q}$.

Let $F\in(K\otimes\mathcal A)[T]$ with $\deg F<N$ and $c_{t,j}:=[z^j]F(t+\alpha_tz)\in\mathcal A$ for $j<m_t$.
1. **(formula)** $E_t(t+\alpha_tz)^{-1}\bmod z^{m_t}$ has coefficients in $\mathcal O_v$ and constant term $1$, and
$$F(0)=\sum_tE_t(0)\sum_{j<m_t}\big[z^j\big]\Big(F(t+\alpha_tz)\,E_t(t+\alpha_tz)^{-1}\Big)\Big(\frac{-t}{\alpha_t}\Big)^j .$$
2. **(precision)** Put $\ell_t=-v(E_t(0))$ and $P_t=P+\ell_t+(m_t-1)\max(0,-v(t/\alpha_t))$. If each $c_{t,j}$ is known modulo $\varpi^{P_t}\mathcal A$, then (1), applied to any lifts, determines $F(0)$ modulo $\varpi^P\mathcal A$.
3. **(sharpness)** If $v(t/\alpha_t)<0$, $m_t>0$, and node $t$'s data are known only modulo $\varpi^{P_t-1}$, the output of (1) modulo $\varpi^P$ depends on the lift.

*Proof.*
1. **The interpolant.** Let $P_t(z)=[F(t+\alpha_tz)E_t(t+\alpha_tz)^{-1}]_{<m_t}$ and $G(T)=\sum_tE_t(T)P_t((T-t)/\alpha_t)$. Then $\deg G<N$.
   - The $t$-th term agrees with $F$ to order $m_t$ at $t$.
   - Every other term vanishes to order $m_t$ there.
   
   So $F-G$ has degree $<N$ and $N$ roots with multiplicity, hence $F=G$. Evaluate at $T=0$.
2. **Integrality.** $E_t(t+\alpha_tz)=\prod_q(1+\alpha_tz/(t-q))^{m_q}$ with $v(\alpha_t/(t-q))\ge0$.
3. **Precision.** Change $c_{t,j}$ by $\delta\in\varpi^{P_t}\mathcal A$. This changes $P_{t,k}$ ($k\ge j$) by $\delta e_{k-j}$, where $e$ are the integral coefficients of the inverse series. So $F(0)$ changes by $E_t(0)\sum_k\delta e_{k-j}(-t/\alpha_t)^k$, whose valuation is at least $P_t-\ell_t-(m_t-1)\max(0,-v(t/\alpha_t))=P$.
4. **Sharpness.** Take $j=m_t-1$ and $\delta=\varpi^{P_t-1}$. The change is $E_t(0)\varpi^{P_t-1}(-t/\alpha_t)^{m_t-1}$, of valuation exactly $P-1$. ∎

**Corollary 3.2 (Part I: unit nodes, no loss).** Take $K=\mathbb Q$, $v=v_p$ with $p>Q$, nodes $1,\dots,Q$, $\alpha_t=1$ and $m_t=J$. Every difference is a unit, so $\ell_t=0$, $v(t/\alpha_t)=0$ and $P_t=P$. This is the exact Chinese remainder step of 08 Lemma 4.3.

**Corollary 3.3 (Part III: residue systems).** Let:
- $K=\mathbb Q[U]/(f)$ with $v=v_2$ (10 Lemma 2.1);
- $\mathcal N=S\setminus\{0\}$ with $S=\{-L,\dots,L-1\}$ and $L=2^a$;
- $\alpha=2^{a+1}$;
- multiplicities given by a *depth profile*, $m_t=m(v_2(t))$ with $m(a)=0$.

For $v_2(t)=i<a$ put $D_i=2^{a-i}$. Then
$$\ell_t=D_i\,\Delta(i),\qquad\Delta(i)=\sum_{j=1}^{a-i}\frac j{2^j}\big(m(i)-m(i+j)\big),\qquad -v(t/\alpha)=a+1-i .$$
For the tapered profile $m(i)=h\max(0,u-i)$, $\Delta(i)\le6h$. This recovers 10 Lemma 4.1 (whose precision $p+8hD_t+8M(a+1)$ exceeds the rule's).

*Proof.* We have $\ell_t=\sum_{q\ne t}m_q(v(t-q)-v(q))$. Sort the other nodes by valuation:
- **$v(q)<i$:** then $v(t-q)=v(q)$, and they contribute 0.
- **$v(q)=i$, $v(t-q)=i+j$ ($1\le j\le a-i$):** there are $D_i/2^j$ of them, each contributing $j\,m(i)$.
- **$v(q)=i+j$ ($1\le j\le a-i$):** there are $D_i/2^j$ of them (one, $q=-L$, at $j=a-i$), each with $v(t-q)=i$, contributing $-j\,m(i+j)$.

For the taper, $m(i)-m(i+j)\le hj$, and $\sum_jj^2/2^j=6$. ∎ (Checked exactly on six profiles, Part C.)

**The geometry behind the loss.** $\ell_t$ is a weighted count of how much closer the other nodes are to $t$ than to $0$. Two extreme designs result:
- **Unit-distance nodes** (Part I) lose nothing, but they are few. Their number is bounded by the label coordinates, and $p>Q$ is needed.
- **A full residue system** (Part III) supplies many nodes with exact locality, at the price of a depth-dependent loss. Tapering the multiplicities with depth keeps that loss proportional to $D_i$.

---

## 4. Counting the corrections, and the theorem

**Lemma 4.1 (correction counts).** Suppose there are $c_w$ modes of weight $w$.
1. **Targeted.** $N:=\#\{\beta:\sum w\beta<P\}\le\inf_{\theta>0}e^{\theta(P-1)}\prod_w(1-e^{-\theta w})^{-c_w}$.
   - If $c_w\le c$ for all $w$, then $N\le\exp\big(\pi\sqrt{2cP/3}\big)$.
   - If $c_w\ge c$ for all $1\le w<P$, then $\log_2N\ge\frac18\min\big(P-1,\sqrt{c(P-1)}\big)$.
2. **Boolean.** $\#\{U\subseteq E\subseteq\mathcal Z:\sum_{Z\in E}w_Z<P\}\le\inf_\theta e^{\theta(P-1)}\prod_Z(1+2e^{-\theta w_Z})$. For $n$ modes of weight 1 this count is $C_{P-1}(n)=\sum_{j<P}\binom nj2^j$.

*Proof.*
1. **Rankin.** $N\le\sum_\beta e^{\theta(P-1-\sum w\beta)}$, a product of geometric series. Since $e^x-1\ge x$,
$$\sum_{w\ge1}-\ln(1-e^{-\theta w})=\sum_{h\ge1}\frac1{h(e^{\theta h}-1)}\le\frac{\pi^2}{6\theta},$$
   and $\theta=\pi\sqrt{c/(6P)}$ gives the second bound.
2. **Lower bound, case $c\ge(P-1)/8$.** Subsets of at most $P-1$ distinct weight-1 modes give $N\ge2^{\min(c,P-1)}$.
3. **Lower bound, case $c<(P-1)/8$.** Let $W=\lfloor\sqrt{2(P-1)/c}\rfloor-1\ge3$. All $2^{cW}$ subsets of the modes of weight $\le W$ have total weight $\le c(W+1)^2/2\le P-1$. So
$$\log_2N\ge cW\ge\sqrt{2c(P-1)}-2c\ge(\sqrt2-\tfrac1{\sqrt2})\sqrt{c(P-1)} .$$
4. **Boolean.** Each mode is outside $E$, or in $E$ with two choices for $U$. ∎

(Part D checks both bounds against exact counts. The lower bound has a factor $\ge8$ to spare in the tested range, and Rankin $\le0.85$.)

**Theorem 4.2 (interpolation through local nodes).** Let the following be given:
- **(i)** a readout system over $\mathcal A$ (free over $\mathcal O_v$), with target $Z=Z(\iota)$;
- **(ii)** a deformation $\rho(T)$ of degree $\le R'$ and deformed-port degree $d$ (Definition 0.3);
- **(iii)** recovery data $(t,m_t,\alpha_t)_{t\in\mathcal N}$ as in Lemma 3.1, with $\sum m_t>dR'$, $H(t+\alpha_tz)\in\mathcal A[z]$, and $P_t$ from Lemma 3.1(2);
- **(iv)** for each node with $m_t>0$, a local readout and a mode expansion of one of two kinds:
  - **(formal)** of $\rho(t+\alpha_tz)$ in $(\mathcal A/\varpi^{P_t})[z]/(z^{m_t})$, at $z$-adic precision $m_t$;
  - **($\varpi$-adic)** of $\rho(t+\alpha_t\zeta)$ in $\mathcal A/\varpi^{P_t}$, at $\varpi$-adic precision $P_t$, for every $\zeta$ in a set $\Xi_t\subseteq\mathcal O_v$ of more than $dR'$ points with pairwise unit differences.

Then $Z\bmod\varpi^P\mathcal A$ is computed by node evaluation (Lemma 2.3 or 2.4), interpolation in $\zeta$ in the $\varpi$-adic case, and recovery (Lemma 3.1). The work is
$$\sum_{t:\,m_t>0}|\Xi_t|\cdot Y_t\cdot N_t\cdot\sum_{g}W_{t,g}\ +\ \mathrm{poly},$$
where:
- $|\Xi_t|=1$ in the formal case;
- $Y_t$ is the number of measurement tuples;
- $N_t$ is the number of correction indices, bounded by Lemma 4.1;
- $W_{t,g}$ is the cost of one group sum, either by brute force or by a recursive call.

*Proof.*
- **Formal case.** The node lemma gives $H(t+\alpha_tz)\bmod(z^{m_t},\varpi^{P_t})$, i.e. the $c_{t,j}\bmod\varpi^{P_t}$.
- **$\varpi$-adic case.** It gives $H(t+\alpha_t\zeta)\bmod\varpi^{P_t}$ for $\zeta\in\Xi_t$. Since $H(t+\alpha_tz)\in\mathcal A[z]$ has degree $\le dR'<|\Xi_t|$ and the Vandermonde matrix of $\Xi_t$ is a unit, Lagrange interpolation over $\mathcal A/\varpi^{P_t}$ returns all its coefficients mod $\varpi^{P_t}$.
- **Recovery.** Lemma 3.1 returns $H(0)=Z\bmod\varpi^P$. ∎

**Corollary 4.3 (Parts I–II; 08 Theorem 4.6, 09 Theorem 1.1).** The data are:
- $Q<p$ unit nodes with $m_t=J=\lfloor dR'/Q\rfloor+1$, formal neighbourhoods, $\ell_t=0$ and $P_t=\mu$;
- Boolean modes on all $n$ variables, of weight 1;
- swap-design measurements;
- brute-force groups of size $\le s$.

Lemma 4.1(2) gives $N_t=C_q(n)$, and $Y_t\le(2\beta+1)^{c_*}$, so the work is $Q\cdot C_q(n)\cdot(2\beta+1)^{c_*}\cdot2r2^s\cdot\mathrm{poly}$. This is the bound of 08 Theorem 4.6.

**Corollary 4.4 (Part III, one level of the recursion; 10 §7).** The data are:
- nodes $S\setminus\{0\}$ with the tapered profile, $\alpha=2^{a+1}$, and 2-adic neighbourhoods;
- $\Xi_t$ the $M(2L-1)+1$ binary points of $\mathcal O$ (10 Lemma 2.1(3));
- targeted tile modes, with $c_w\le6M/D_t$ for each weight (10 Lemma 6.4(2));
- no measurements;
- recursive groups with $\le6M/r$ wires.

Lemma 4.1(1) gives $N_t\le\exp(\pi\sqrt{2c_tP_t/3})$. With 10 Lemma 7.3's $P_t/D_t\le34M/u$ this is $e^{O(M/\sqrt u)}$, as in 10 Theorem 10.2.

---

## 5. The two square-root barriers

In both corollaries the node cost is a saddle-point count, $\log N_t\approx\min_\theta[\theta\cdot\text{budget}+\log(\text{mode generating function})]$. What differs is what forces the budget.

**Proposition 5.1 (Part I's barrier, [07](07_TOOLKIT.md) E4 restated).** In Corollary 4.3 the node exponent is
$$\log_2C_q(n)\ \ge\ q\log_2(n/q),\qquad q\ge nR'/Q .$$
The swap design's locality is exact only on a window of $R$ levels of the swap spectrum. That spectrum is $\mathrm{Bin}(Q,p_r)$, the number of antisymmetric factors (08 Lemma 5.1). The measurement count therefore carries the defect $\Pr[W\notin\text{window}]=e^{-\Theta(R^2/Q)}$. So:
- a per-variable exponent $\delta$ needs $R^2\gtrsim Q\log(1/\delta)$;
- small groups need $r\ge1/\delta$;
- the labels need $r^{2Q}\le n$.

E4 combines these into $\delta\gtrsim(\log\log n)^2/\sqrt{\log n}$. **The square root is the central limit theorem's window of the swap spectrum.** [14](14_DESIGN_FLOOR.md) shows it is intrinsic: for every graded design, symmetric or not, it is Smolensky's agreement degree of parity, $\Theta(\sqrt{Q\log(1/\varepsilon)})$.

**Lemma 5.2 (the node cost does not depend on depth).** In Corollary 4.4, suppose at least $M/2$ wires cross at coordinate $i$, as the paper's balanced labels ensure (10 Lemma 5.1(3) with $r=64$).
- **Modes.** Every crossing tile has $D_i/2$ columns, so there are $c_i\ge M/D_i$ modes of each weight.
- **Precision.** Lemma 3.1 requires $P_t\ge p+D_i\Delta(i)$ at a depth-$i$ node.

Hence
$$c_i\cdot P_t\ \ge\ M\,\Delta(i)+Mp/D_i .$$
The 2-adic scale $D_i$ cancels from the leading term: **a node's cost depends on the local slope $\Delta(i)$ of the multiplicity profile, not on its depth.**

*Proof.* Combine the tile length (10 §6), Corollary 3.3, and Lemma 3.1(2). ∎

**Lemma 5.3 (profile lemma).** Let $m:\{0,\dots,a-1\}\to\mathbb Z_{\ge0}$ with $\sum_{i<a}2^{a-i}m(i)>M(2L-1)$. This is Hermite's condition when $\deg F$ can reach $M(2L-1)$. Then
$$\Delta^*:=\max_{i:\,m(i)>0}\Delta(i)\ >\ \frac{M(2L-1)}{4(L(a-1)+1)}\ \ge\ \frac{3M}{8a}.$$

*Proof.*
1. **Claim:** $m(i)\le2\Delta^*(a-i)$. We argue by downward induction on $i$.
   - **Base, $i=a-1$:** $\Delta(a-1)=m(a-1)/2$.
   - **Step, $N=a-i\ge2$, $m(i)>0$:** write $K_N=\sum_{j\le N}j/2^j$. Then
     $$m(i)K_N=\Delta(i)+\sum_{j=1}^{N-1}\frac j{2^j}m(i+j)\le\Delta^*+2\Delta^*\sum_{j=1}^{N-1}\frac j{2^j}(N-j).$$
     The right side is at most $2\Delta^*NK_N$ because $NK_N-\sum_{j<N}\frac j{2^j}(N-j)=\frac{N^2}{2^N}+\sum_{j<N}\frac{j^2}{2^j}\ge\frac12$.
2. **Sum.** $\sum_{i<a}2^{a-i}m(i)\le2\Delta^*\sum_{k=1}^ak2^k=2\Delta^*((a-1)2^{a+1}+2)$.
3. **Compare.** Compare with $M(2L-1)$, and use $L\ge2$ for the last inequality. ∎

**How good is the bound?**
- **The optimum.** The exact optimum of $\min\Delta^*$ over all real profiles is solved greedily (deepest level first; the constraints are triangular and monotone). As a multiple of $M/a$ it is $0.75,1.17,1.84,2.78,3.80,4.65,5.24,5.59,5.79$ for $a=1,2,4,\dots,256$, tending to $6$.
- **The paper's taper** $h(a-i)$ gives $0.75,2.10,4.59,6.39,6.40,6.19,6.10,6.05,6.02$.

So **the paper's linear taper is asymptotically optimal**: both tend to $6M/a$. Lemma 5.3's constant $3/8$ is safe but not sharp (Part E).

**Theorem 5.4 (floor of the 2-adic architecture).** Consider one call of the Part III architecture with these features:
- $M$ wires and $L=2^a$;
- Hermite recovery from the nodes $S\setminus\{0\}$ with any depth profile supported on depths $i<a$, by the formula of Lemma 3.1, correct for all integral data at output precision $p$. Locality at a node of depth $i$ needs a label coordinate $i<a$ (10 Lemma 6.1);
- node data from targeted tile modes, with at least $M/2$ crossing wires at each used coordinate.

Then some node needs at least
$$2^{\frac18\min\left(p-1,\ \sqrt{3/16}\,M/\sqrt a\right)}$$
correction tuples $\beta$, each requiring child evaluations. Packing needs $L\le M$, so $a\le\log_2M$. At the top level of a #SAT computation $p\ge n+1\ge M/C$, so the work is $2^{\Omega(M/\sqrt{\log M})}$. This matches 10 Theorem 10.2.

*Proof.*
1. **The worst node.** Lemma 5.3 gives a used depth $i^*$ with $\Delta(i^*)>3M/(8a)$.
2. **Precision.** By Lemma 3.1(2)–(3), a node of depth $i^*$ needs $P\ge p+D_{i^*}\Delta(i^*)$.
3. **Modes.** By Lemma 5.2 it has $c\ge M/D_{i^*}$ modes of every weight below $P$.
4. **Combine.** $c(P-1)\ge cP/2>3M^2/(16a)$, and Lemma 4.1(1) gives the count. ∎

**Comparison.** The upper bound of 10 Lemma 7.3 is $\log_2N_t\le53M/\sqrt u$ with $u\asymp a$. The floor and the algorithm thus agree up to constants: $\Theta(M/\sqrt{\log M})$ per call, summing geometrically over the recursion. **The square root is Hardy–Ramanujan's.** The number of mode tuples below weight $P$ is $e^{\Theta(\sqrt{cP})}$, and $cP\asymp M\Delta^*\asymp M^2/a$.

---

## 6. The design space

| | Parts I–II (swap design) | Part III (Cauchy design) |
|---|---|---|
| Deformed object | clause ports read rows $P_{f,v}(T)$ | wire ports read $(W(T)X)_x$ |
| Symmetry | $\mathbb Z_2^Q$ swaps on labels $([r]\times[r])^Q$ | shifts of $\mathbb Z/2^{a+1}$ on packet slots |
| Node $\leftrightarrow$ coordinate | node $t$ $\leftrightarrow$ swap position $t$ | node $t$ $\leftrightarrow$ coordinate $v_2(t)$ |
| Degree vs. nodes | $R-1\approx2\sqrt\lambda\ll Q$ | $2L-1$ = (number of nodes) $-1$ (exact Lagrange) |
| Locality at a node | modulo measurements; defect $e^{-\Theta(R^2/Q)}$ | exact |
| Neighbourhood | formal ($z$-adic), order $J=q+1$ | 2-adic, precision $p_t$ |
| Modes | all $n$ coordinates (dense perturbation) | tile power moments, $\le6M/D_t$ per weight (structured) |
| Correction count | subsets, $C_q(n)$ | weighted partitions, $e^{O(\sqrt{cP})}$ |
| Recovery | CRT at unit nodes, loss $0$ | tapered Hermite, loss $D_i\Delta(i)\le6hD_i$ |
| Groups | brute force, $2^{O(n/r)}$, so $r\gtrsim1/\delta$ | recursive, so $r=64$ |
| Exponent per call | $n(\log\log n)^2/\sqrt{\log n}$ (floor E4) | $M/\sqrt{\log M}$ (floor Theorem 5.4) |

**What the comparison teaches.** Five reusable design principles follow.
1. **Structured perturbations beat subset entropy.** A dense near-node perturbation makes every variable a mode, and the count is $\binom nq$. A structured one, such as a geometric series along tiles, has few modes per weight, and the count is a partition number. This is where Part III gains its $\log\log$ factors.
2. **Exact locality needs full degree, and a $\varpi$-adic neighbourhood pays for it.** Exact Lagrange locality at $k$ nodes forces degree $k-1$. Part III then evaluates in a 2-adic neighbourhood, where the non-local part is small to all orders, instead of in a formal one, where it is small only to the jet order.
3. **Recursion needs exact locality.** Measurement constraints take the group sums out of the problem class. Part III's exact locality makes the children tensor networks again, so groups can be large ($r=64$). This in turn allows $a\asymp\log M$ node coordinates instead of $Q\asymp\log n/\log\log n$.
4. **Recovery loss is node geometry.** Unit-distance nodes cost nothing (Corollary 3.2). Residue systems cost $D_i\Delta(i)$ (Corollary 3.3), and the linear taper is asymptotically the best depth profile (Lemma 5.3 and the computed optimum).
5. **The node cost is depth-free** (Lemma 5.2). In the 2-adic design, precision and mode count trade off exactly across scales.

---

## 7. An open question: the square-root law

Both architectures pay $M/\sqrt{(\text{number of node coordinates})}$ per call, up to polylogarithms:
- **Part I:** $\delta\approx\sqrt{\log(1/\delta)/Q}\cdot\log Q$, from the central limit theorem;
- **Part III:** $\sqrt{\Delta^*/M}\asymp1/\sqrt a$, from Hardy–Ramanujan and the profile lemma.

The node coordinates are coordinates of a label space of polynomial size, so there are $O(\log M)$ of them, and both give $2^{\tilde\Omega(M/\sqrt{\log M})}$. Two different square roots yield the same law.

**Question (`T.ETH.5`).** Does every interpolation-through-local-nodes scheme, in the sense of Theorem 4.2 with node coordinates from a label space of size $M^{O(1)}$, cost $2^{\Omega(M/\sqrt{\log M})}$ per call?
- **A yes** would make the paper's $n/\sqrt{\log n}$ the limit of the whole method family.
- **A no** would improve Theorem 1.1.

For Part I's architecture with group-algebra designs, [14](14_DESIGN_FLOOR.md) settles the graded case (`T.ETH.1`, done) and isolates the ungraded case (`T.ETH.6`).

---

## 8. Status

| Result | Depends on | Verdict |
|---|---|---|
| Definitions 0.1–2.2, Lemma 1.2 | — (Lemma 1.2 is 08 Lemma 4.5) | complete derivation |
| Lemmas 2.3, 2.4 | Lemma 1.2 | complete derivation; checked on 27 random systems, witnesses, and both instances |
| Lemma 3.1, Corollaries 3.2–3.3 | — | complete derivation; checked exactly, with sharpness |
| Lemma 4.1 | — | complete derivation; checked against exact counts |
| Theorem 4.2, Corollaries 4.3–4.4 | the above; instances from 08–10 | complete derivation (assembly); the instances' hypotheses are proved in 08–10 |
| Proposition 5.1 | 07 E4 | complete derivation (in 07) |
| Lemmas 5.2–5.3, Theorem 5.4 | Lemmas 3.1, 4.1; 10 §5–§6 | complete derivation, for the architecture as stated; Lemma 5.3 also checked numerically |

**Subjective confidence:**
- ≈97% that every statement above is correct as stated;
- ≈90% that Theorem 5.4's hypotheses capture "App. D's architecture" as a reader would understand it. The floor is for depth-only profiles and formula-based recovery, which is what the paper uses.

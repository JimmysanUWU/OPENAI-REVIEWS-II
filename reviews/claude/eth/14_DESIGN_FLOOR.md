# 14. The design floor: every graded swap design pays $(\log\log n)^2/\sqrt{\log n}$

*Claude, 29 October 2026. Roadmap target `T.ETH.1`; nodes `ETH.D.*` in [`ROADMAP.md`](../../../ROADMAP.md).*

Part I deforms the identity inside the group algebra of the swaps: $P(T)=\sum_Uc_U(T)\Pi_U$ with $P(0)=I$. At node $t$ the terms with $t\in U$ are local. The rest, the *defect*, has to be paid for with measurements.

[07](07_TOOLKIT.md) E4 showed that the paper's symmetric design cannot beat the exponent $(\log\log n)^2/\sqrt{\log n}$, and left non-symmetric designs open (target `T.ETH.1`). This file settles the question for every **graded** design: one whose coefficient $c_U$ is constant whenever $\Pi_U$ leaves $R$ or more positions unswapped. Graded designs include every design with a node-independent defect, so every choice of $M_0$ in E1, symmetric or not.

> **Theorem 3.2 (design floor).** Run Part I's algorithm with $r$ groups, $Q$ nodes, $r^{2Q}\le n$ labels, and any graded group-algebra design of degree $<R$ measured by eigenvectors. In E4's cost model the exponent per variable is at least $\frac1{100}(\log_2\log_2n)^2/\sqrt{\log_2n}$ for all large $n$.

> **Lemma 2.2 (parity under a biased measure).** Let $0<p\le\frac12$. Every real polynomial $F$ of degree $\le D$ on $\{0,1\}^Q$ differs from parity on a set of $\mu_p$-measure at least $\mathbb E\,\tau_D(\mathrm{Bin}(Q,2p))$, where $\tau_D(m)=\Pr[\mathrm{Bin}(m,\frac12)>\frac{m+D}2]$.

> **Why the floor holds.**
> - At node $t$ the defect of a graded design acts on the swap eigenvectors with antisymmetric set $X$ by the scalar $1-(-1)^{|X|}G_t(X)$, where $G_t$ is a real polynomial of degree $<R$ (Lemma 1.2).
> - So the measurements are the $X$ where $G_t$ differs from parity.
> - Smolensky's dimension argument bounds how often a low-degree polynomial can agree with parity. A random restriction carries the bound to the biased measure of the swap spectrum.
>
> **Part I's square root is therefore the agreement degree of parity**: $\Theta(\sqrt{Q\log(1/\varepsilon)})$ at error $\varepsilon$ (Corollary 2.5). It is not an artifact of the symmetric design.

**What is new here.** Complete derivations in this file:
- the normal form of group-algebra designs, and the graded class (Lemma 1.2);
- the biased-measure agreement bound for parity, with explicit constants (Lemmas 2.2–2.4, Corollary 2.5);
- the floor for every graded design (Theorem 3.2), which answers `T.ETH.1`;
- the face lemma for arbitrary designs, and a precise statement of what remains for ungraded designs (§5, new target `T.ETH.6`).

**Sources.**
- Lemma 2.1 is Smolensky's dimension argument (R. Smolensky, *Algebraic methods in the theory of lower bounds for Boolean circuit complexity*, STOC 1987; *On representations by low-degree polynomials*, FOCS 1993).
- Lemma 2.2's restriction step is routine and is written out here.
- The measurement cost (C3) of §3 is E4's and is not re-derived.

**Checks.** [`scripts/verify_design_floor.py`](scripts/verify_design_floor.py) ([`results/design_floor.json`](results/design_floor.json), ~70 s). The flag `--full` adds $(Q,R)=(8,4)$ to Part E, about 3 minutes more.

| Part | What it checks | Result |
|---|---|---|
| A | Lemma 2.1 by exact rational elimination: 180 random sets of $N_D+1$ points ($Q\le6$, every $D$) | none admits a degree-$\le D$ interpolant of parity |
| A | Central windows of $D+1$ levels | all interpolable |
| A | Control: the same random sets in degree $D+2$ | 100 interpolable, so the test can fail |
| B | Lemma 2.2, exactly, for window designs and three kinds of greedy non-symmetric designs ($Q\le6$, $p\in\{\frac14,\frac13,\frac38,\frac12\}$): 192 cases | all hold; equality in 4 cases ($Q=5$, $D=0$, $p=\frac12$) |
| B | The coupling identity $\mu_p(E)=\mathbb E_Z[\text{bad fraction on the subcube }Z]$ | exact |
| B | Non-symmetric designs against the window design | the best one found has a bad set up to 27% lighter |
| B | Control: window designs of degree $D+2$ against the degree-$D$ bound | violate it in 48 cases |
| C | Lemma 2.3(a) for $m\le300$; Lemma 2.3(b) for $Q\le2000$ | hold, with log-margins $\ge1.44$ and $\ge1.79$ |
| D | Theorem 3.2's three costs minimized over integer $(r,Q,R)$, $\log_2n=2^{10},\dots,2^{300}$ | ratio to $(\log_2\log_2n)^2/\sqrt{\log_2n}$ in $[0.018,0.179]$, rising |
| E | Ungraded against single-polynomial designs on central windows, $(Q,R)=(6,3),(7,3)$; with `--full` also $(8,4)$ | feasible on exactly the same windows |

---

## 0. Notation

- **Cube.** Identify $X\subseteq[Q]$ with its indicator in $\{0,1\}^Q$. Put $\mathrm{par}(X)=(-1)^{|X|}$ and $\chi_S(X)=(-1)^{|S\cap X|}$, so $\mathrm{par}=\chi_{[Q]}$.
- **Degree.** Every $F:\{0,1\}^Q\to\mathbb R$ is $\sum_S\hat F(S)\chi_S$. Its *degree* is $\max\{|S|:\hat F(S)\ne0\}$; this equals its degree as a multilinear polynomial in the 0/1 coordinates. $F_{\ni t}=\sum_{S\ni t}\hat F(S)\chi_S$ is the part of $F$ that is odd in $x_t$.
- **Measures.** $\mu_p(X)=p^{|X|}(1-p)^{Q-|X|}$, and $p_r=\frac{r-1}{2r}\in[\frac14,\frac12)$.
- **Smolensky counts.** $N_D(m)=\sum_{i\le\lfloor(m+D)/2\rfloor}\binom mi$ and $\tau_D(m)=1-N_D(m)/2^m=\Pr[\mathrm{Bin}(m,\frac12)>\frac{m+D}2]$. Note $\tau_D(m)=0$ if $D\ge m$.
- **Swaps** ([08](08_GAPLESS_PROOF.md) §5). Labels $\mathcal B=([r]\times[r])^Q$, swaps $\Pi_s$, and $\Pi_U=\prod_{s\in U}\Pi_s$. The eigenvectors $\psi$ of 08 Lemma 5.1 have an antisymmetric set $X(\psi)$. By 08 Lemma 5.1(2)–(3):
  - $\Pi_U\psi=\chi_U(X(\psi))\,\psi$;
  - exactly $|\mathcal B|\,\mu_{p_r}(X)$ eigenvectors have $X(\psi)=X$.

---

## 1. Group-algebra designs

**Definition 1.1 (group-algebra design).** Let $1\le R\le Q$. A *group-algebra design of degree $<R$* is
$$P(T)=\sum_{U\subseteq[Q]}c_U(T)\,\Pi_U,\qquad c_U\in\mathbb Q[T],\quad\deg c_U<R,\quad P(0)=I.$$
At a node $t\in[Q]$:
- its *local part* is $\sum_{U\ni t}c_U(t)\Pi_U$;
- its *defect* is $N_t=\sum_{U\not\ni t}c_U(t)\Pi_U$.

Its *eigenvalue polynomials* are $n_X(T)=\sum_Uc_U(T)\chi_U(X)$, so that $P(T)\psi=n_{X(\psi)}(T)\psi$. They have degree $<R$ and $n_X(0)=1$.

The design is **graded** if $c_U$ is constant whenever $|U^c|\ge R$, that is, whenever $\Pi_U$ leaves at least $R$ positions unswapped.

*Why these are the right objects.*
- $\Pi_U$ with $t\in U$ maps $u$ to a label whose first schedule at $t$ is $d_u(t)$. This is the locality of 08 Theorem 5.4 (B2).
- Conversely, a combination $\sum_{U\not\ni t}a_U\Pi_U$ is local at $t$ only if every $a_U=0$. To see this, take a label $u$ whose $Q$ pairs are all off-diagonal. Its images $\theta_Uu$ are distinct, and each keeps $a(t)=a_u(t)\ne d_u(t)$, so each $a_U$ is a nonlocal entry.

So $N_t$ is exactly the part of $P(t)$ that must be measured.

**Lemma 1.2 (normal form).**
- **(a) Graded means "one minus parity times low degree".** A design is graded if and only if
  $$n_X(T)=1-\mathrm{par}(X)\,\Lambda_T(X),$$
  where each $\Lambda_T$ is a real polynomial on $\{0,1\}^Q$ of degree $<R$, its coefficients are polynomials in $T$ of degree $<R$, and $\Lambda_0=0$.
- **(b) The defect.** For a graded design, $N_t$ acts on $\psi$ by the scalar $1-\mathrm{par}(X)\,G_t(X)$, where $X=X(\psi)$ and $G_t=(\Lambda_t)_{\ni t}$ has degree $<R$. Hence:
  - $N_t\psi\ne0$ exactly when $X(\psi)\in E_t:=\{G_t\ne\mathrm{par}\}$;
  - $\operatorname{rank}N_t=|\mathcal B|\,\mu_{p_r}(E_t)$.
- **(c) Every family occurs.** For any polynomials $G_1,\dots,G_Q$ of degree $<R$ with $G_t$ odd in $x_t$, some graded design has exactly these $G_t$.
- **(d) Node-independent defects.** Suppose one matrix $N$ makes $P(t)-N$ local at every node, as in E1 and 08 Theorem 5.4. Then:
  - the design is graded;
  - $N=I-M_0$ with $M_0=\sum_{|S|<R}g_S\Pi_{S^c}$;
  - $N$ acts by $1-\mathrm{par}(X)F(X)$ with $F=\sum_{|S|<R}g_S\chi_S$, so the paper's measured set (08 (B3)) is $E=\{F\ne\mathrm{par}\}$.

  Every $F$ of degree $<R$ occurs, through E1's formula with $g_S=\hat F(S)$. The paper's choice is the symmetric $F=f(|X|)$.

*Proof.*
- **(a)** Write $S=U^c$. Since $c_U(0)=[U=\emptyset]$, being graded means $c_{S^c}\equiv[S=[Q]]$ whenever $|S|\ge R$. Using $\chi_{S^c}=\mathrm{par}\cdot\chi_S$ and $R\le Q$:
  $$n_X(T)=1+\sum_{|S|<R}c_{S^c}(T)\chi_{S^c}(X)=1+\mathrm{par}(X)\sum_{|S|<R}c_{S^c}(T)\chi_S(X).$$
  Put $\Lambda_T=-\sum_{|S|<R}c_{S^c}(T)\chi_S$. Then $\Lambda_0=0$, because $c_{S^c}(0)=0$ for $S\ne[Q]$.

  Conversely, given $\Lambda$, set $c_{S^c}=-\hat\Lambda_T(S)$ for $|S|<R$, $c_\emptyset=1$, and $c_U=0$ for all other $U$. The map $(c_U)\mapsto(n_X)$ is Fourier inversion, so it is a bijection.
- **(b)** $N_t\psi=\sum_{U\not\ni t}c_U(t)\chi_U(X)\,\psi$. Since $\chi_U(X\oplus e_t)=\pm\chi_U(X)$ according as $t\notin U$ or $t\in U$,
  $$\sum_{U\not\ni t}c_U(t)\chi_U(X)=\tfrac12\big(n_X(t)+n_{X\oplus e_t}(t)\big).$$
  By (a) and $\mathrm{par}(X\oplus e_t)=-\mathrm{par}(X)$, this equals
  $$1-\tfrac12\mathrm{par}(X)\big(\Lambda_t(X)-\Lambda_t(X\oplus e_t)\big)=1-\mathrm{par}(X)\,G_t(X).$$
  Since $\mathrm{par}=\pm1$, the scalar vanishes iff $G_t(X)=\mathrm{par}(X)$. The $\psi$ form an eigenbasis, so the rank is the number of $\psi$ with a nonzero scalar.
- **(c)** For $|S|<R$, let $\lambda_S$ be the polynomial of degree $\le|S|$ with $\lambda_S(0)=0$ and $\lambda_S(t)=\hat G_t(S)$ for $t\in S$ (Lagrange interpolation at the $|S|+1$ points of $\{0\}\cup S$). Put $\Lambda_T=\sum_{|S|<R}\lambda_S(T)\chi_S$. Then $(\Lambda_t)_{\ni t}=\sum_{S\ni t}\hat G_t(S)\chi_S=G_t$, and (a) gives the design.
- **(d)**
  1. **$N$ is in the group algebra.** It is $P(t)$ minus a local element of the group algebra. Write $N=\sum_Un_U\Pi_U$. The $\Pi_U$ are linearly independent: every eigenspace is nonzero ($0<p_r<1$), and the characters $\chi_U$ separate them.
  2. **Coefficients.** Locality of $P(t)-N$ at $t$ forces $c_U(t)=n_U$ for every $t\in U^c$.
  3. **Graded.** If $|U^c|\ge R$, then $c_U-n_U$ has degree $<R$ and at least $R$ roots, so $c_U$ is constant. Then $n_U=c_U(0)=[U=\emptyset]$.
  4. **Normal form.** For $|U^c|<R$ put $g_{U^c}=[U=\emptyset]-n_U$. Then $N=I-\sum_{|S|<R}g_S\Pi_{S^c}$, which acts by $1-\sum_Sg_S\chi_{S^c}(X)=1-\mathrm{par}(X)F(X)$.
  5. **Every $F$ occurs.** E1's items 1–2 hold for arbitrary coefficients $g_S$. ∎

**Remark 1.3 (what grading excludes).** In a graded design the deformation degree $R$ in $T$ also bounds the Fourier degree in $X$ of the defect polynomial $\Lambda_T$. An ungraded design lets group elements that leave $R$ or more positions unswapped carry $T$-dependent coefficients. Then $\Lambda_T$ can have high degree in $X$, and §2 no longer applies. §5 treats this case.

---

## 2. Parity under a biased measure

**Lemma 2.1 (Smolensky's bound).** Let $F:\{0,1\}^m\to\mathbb R$ have degree $\le D$. Then $F=\mathrm{par}$ on at most $N_D(m)$ points. Equivalently, $\mu_{1/2}(F\ne\mathrm{par})\ge\tau_D(m)$.

*Proof.* Let $A=\{F=\mathrm{par}\}$. The functions on $A$ form a space of dimension $|A|$, spanned by the restricted characters $\chi_S|_A$. If $|S|>\lfloor\frac{m+D}2\rfloor$, then on $A$
$$\chi_S=\mathrm{par}\cdot\chi_{S^c}=F\,\chi_{S^c}=\sum_B\hat F(B)\,\chi_{B\triangle S^c},\qquad|B\triangle S^c|\le D+m-|S|\le\Big\lfloor\frac{m+D}2\Big\rfloor .$$
So the $N_D(m)$ characters of degree $\le\lfloor\frac{m+D}2\rfloor$ span all functions on $A$, and $|A|\le N_D(m)$. ∎

*Sharpness.* The central window design $F=f(|X|)$, with $f$ interpolating $(-1)^j$ on the $D+1$ levels from $\lceil\frac{m-D}2\rceil$, errs only outside the window. By the symmetry of $\mathrm{Bin}(m,\frac12)$, that has measure $\le2\tau_D(m)$. So, for the uniform measure, the window design is optimal up to a factor $2$ among all polynomials, symmetric or not.

**Lemma 2.2 (biased measures by restriction).** Let $0<p\le\frac12$ and let $F$ have degree $\le D$ on $\{0,1\}^Q$. Then
$$\mu_p(F\ne\mathrm{par})\ \ge\ \mathbb E\,\tau_D(M),\qquad M\sim\mathrm{Bin}(Q,2p).$$

*Proof.*
1. **Coupling.** Let $Z\sim\mu_{2p}$ and $Y\sim\mu_{1/2}$ be independent, and $X=Z\wedge Y$ coordinatewise. The coordinates of $X$ are independent with $\Pr[X_i=1]=p$, so $X\sim\mu_p$.
2. **Restriction.** Fix $Z=z$ with support $A$, $|A|=m$. Then $X$ is $Y_A$ padded with zeros, and $\mathrm{par}(X)=\mathrm{par}(Y_A)$. Setting the coordinates outside $A$ to $0$ in the multilinear form of $F$ leaves a polynomial in $Y_A$ of degree $\le D$.
3. **Smolensky.** By Lemma 2.1 this polynomial differs from parity on a fraction $\ge\tau_D(m)$ of $\{0,1\}^A$. Average over $Z$, using $|Z|\sim\mathrm{Bin}(Q,2p)$. ∎

**Lemma 2.3 (explicit forms).**
- **(a)** If $0\le D\le m-2$, then $\tau_D(m)\ge e^{-(D+2)^2/m}/(m+1)$.
- **(b)** If $m\ge Q/2$, then $\Pr[\mathrm{Bin}(m,\frac12)<Q/32]\le e^{-0.229\,Q}$.

*Proof.*
- **(a)** Put $j=\lfloor\frac{m+D}2\rfloor+1\le m$.
  - The method of types gives $\binom mj2^{-m}\ge e^{-m\,\mathrm{KL}(j/m\|1/2)}/(m+1)$.
  - $\mathrm{KL}(a\|\frac12)\le\chi^2=4(a-\frac12)^2$.
  - Since $j-\frac m2\le\frac{D+2}2$, we get $m\,\mathrm{KL}\le(D+2)^2/m$.
- **(b)** Since $Q/32\le m/16$, Chernoff's bound gives $\Pr[\mathrm{Bin}(m,\frac12)\le\frac m{16}]\le e^{-m\,\mathrm{KL}(1/16\|1/2)}$. Here $\mathrm{KL}(\frac1{16}\|\frac12)=0.4593\ldots$, so the probability is $\le e^{-0.2296\,Q}$. ∎

**Lemma 2.4 (the heavy part of the bad set).** Let $\frac14\le p\le\frac12$, $Q\ge4$, $1\le R\le Q/4$, and let $F$ have degree $<R$. Put $H=\{X:|X|\ge Q/32\}$. Then
$$\mu_p\big(\{F\ne\mathrm{par}\}\cap H\big)\ \ge\ \varepsilon'(Q,R):=\tfrac12\Big(\frac{e^{-2(R+1)^2/Q}}{Q+1}-e^{-0.229\,Q}\Big).$$
If $Q\ge64$, then $\varepsilon'(Q,R)\ge e^{-2(R+1)^2/Q}/(4(Q+1))$.

*Proof.*
1. **Coupling.** Couple as in Lemma 2.2, and fix $Z$ with $|Z|=m\ge Q/2$. Then
   $$\Pr_Y[X\in E\cap H]\ \ge\ \Pr_Y[X\in E]-\Pr_Y[X\notin H]\ \ge\ \tau_{R-1}(m)-\Pr\big[\mathrm{Bin}(m,\tfrac12)<\tfrac Q{32}\big].$$
2. **Explicit bounds.** Since $(R-1)+2\le\frac Q4+1\le\frac Q2\le m$, Lemma 2.3 bounds this below by
   $$\frac{e^{-(R+1)^2/m}}{m+1}-e^{-0.229Q}\ \ge\ \frac{e^{-2(R+1)^2/Q}}{Q+1}-e^{-0.229Q}.$$
3. **Average.** $\Pr[|Z|\ge Q/2]\ge\frac12$, because $|Z|\sim\mathrm{Bin}(Q,2p)$ with $2p\ge\frac12$.
4. **Second claim.** $2(R+1)^2/Q\le Q/8+1+2/Q$, and $0.104\,Q\ge\ln2+1+2/Q+\ln(Q+1)$ for $Q\ge64$. So the first term in the bracket is at least twice the second. ∎

**Corollary 2.5 (the agreement degree of parity).** Let $\frac14\le p\le\frac12$ and $Q\ge64$.
- **Lower bound.** If a polynomial of degree $<R\le Q/4$ differs from parity on $\mu_p$-measure $\varepsilon\le(4(Q+1))^{-2}$, then $R+1\ge\frac12\sqrt{Q\ln(1/\varepsilon)}$.
- **Upper bound.** For even $R\ge4$, the window design errs on measure $\le2e^{-R^2/(8Q)}$. This is 08 Lemma 5.5 and Corollary 5.6, with $Qp$ in place of $Qp_r$.

So the least degree for error $\varepsilon$ in this range is $\Theta(\sqrt{Q\log(1/\varepsilon)})$, for symmetric and non-symmetric polynomials alike.

*Proof.* By Lemma 2.4, $\varepsilon\ge e^{-2(R+1)^2/Q}/(4(Q+1))$. Hence
$$\frac{2(R+1)^2}Q\ \ge\ \ln\frac1\varepsilon-\ln(4(Q+1))\ \ge\ \frac12\ln\frac1\varepsilon .\qquad ∎$$

---

## 3. The floor

**Cost model 3.1 ([07](07_TOOLKIT.md) E4).** Instantiate the algorithm of 08 Theorems 4.6 and 7.1 with:
- $r$ groups and $Q$ nodes;
- $|\mathcal B|=r^{2Q}\le n$ labels;
- a group-algebra design of degree $<R$;
- eigenvector measurements: at node $t$, the rows $\psi^{\mathsf T}$ for $X(\psi)$ in a set $\mathcal E_t$, transferred as in 08 Lemma 6.1. The set must make $P(t)$ local modulo measurements: for example $E_t$, or, for a node-independent defect, the paper's $E$ (08 (B3)–(B4)).

Let $\delta$ be the exponent per variable, i.e. $\log_2$ of the cost bound of 08 Theorem 4.6 divided by $n$. Then:
- **(C1)** $\delta\ge1/r$, since some group has $\ge n/r$ variables.
- **(C2)** $\delta\ge x\log_2(1/x)$ with $x=(R-1)/Q$ if $x\le1/e$, and $\delta\ge\frac1e\log_2e$ otherwise. This is the Taylor-jet count $\log_2\binom nq\ge q\log_2(n/q)$ with $q\ge\lfloor n(R-1)/Q\rfloor$, up to $O((\log n)/n)$ from rounding $q$.
- **(C3)** If, at some node, the measured eigenvectors with $j(\psi)\ge j_0$ number $|\mathcal B|\,\varepsilon_0$, then $\delta\ge\max(j_0,1)\,\varepsilon_0/2$.

These are the three costs in the proof of E4. E4 states them for symmetric designs, but its argument for (C3) uses only:
- the number of measured eigenvectors;
- the fact that a row with $j$ antisymmetric factors has $\ell^1$-norm $\ge2^{j-1}$ on some group, so it contributes $\ge\max(j,1)$ bits to the factor $(2\beta+1)^{c_*}$.

It never uses the symmetry of the design. (C3) is the one input of this file that is not re-derived here.

**Theorem 3.2 (design floor).** Use any graded group-algebra design of degree $<R$ in the cost model 3.1. Measure it either:
- as the paper does, for a node-independent defect: all of $E=\{F\ne\mathrm{par}\}$ (Lemma 1.2(d));
- or node by node: any $\mathcal E_t\supseteq E_t$ (Lemma 1.2(b)).

Then, with $L=\log_2n$ and $\ell=\log_2L$,
$$\delta\ \ge\ \frac1{100}\,\frac{\ell^2}{\sqrt L}\qquad\text{for all large }n.$$

*Proof.* In both cases the measured set at each node contains a set $\{G\ne\mathrm{par}\}$ with $\deg G<R$. Suppose $\delta<\ell^2/(100\sqrt L)$, and let $n$ be large.
1. **Groups and labels.** By (C1), $\log_2r\ge\log_2(1/\delta)\ge\ell/4$. Since $r^{2Q}\le n$, $Q\le L/(2\log_2r)$.
2. **Small cases.**
   - If $R=1$, then $G$ is constant, so the measured set contains a parity class. Its $\mu_{p_r}$-measure is $\ge\frac12(1-2^{-Q})\ge\frac14$, so (C3) with $j_0=0$ gives $\delta\ge\frac18$.
   - If $R\ge2$ and $Q<64$, then $x\ge1/64$ and (C2) gives $\delta\ge\frac3{32}$.
   - If $R\ge2$ and $x>\frac18$, then (C2) gives $\delta\ge\frac38$.

   All three contradict $\delta<\ell^2/(100\sqrt L)$. So $R\ge2$, $Q\ge64$ and $R\le Q/8+1\le Q/4$.
3. **Measurements.** Lemma 2.4 with $p=p_r$ applies to $G$. With $j_0=\lceil Q/32\rceil$, (C3) gives
   $$\delta\ \ge\ \frac Q{64}\cdot\frac{e^{-2(R+1)^2/Q}}{4(Q+1)}\ \ge\ \frac{e^{-2(R+1)^2/Q}}{512}.$$
4. **Degree.** Hence $2(R+1)^2/Q\ge\ln\frac1{512\delta}\ge\frac12\ln\frac1\delta$, the last step because $\delta\le2^{-18}$. So $(R+1)^2\ge\frac Q4\ln\frac1\delta$. Since $R\ge2$, $R-1\ge\frac{R+1}3$, which gives $R-1\ge\frac16\sqrt{Q\ln(1/\delta)}$.
5. **The ratio $x$.** Using step 1,
   $$x=\frac{R-1}Q\ \ge\ \frac16\sqrt{\frac{\ln(1/\delta)\cdot2\log_2r}{L}}\ \ge\ \frac{\sqrt{2\ln2}}6\cdot\frac{\log_2(1/\delta)}{\sqrt L}\ \ge\ 0.049\,\frac{\ell}{\sqrt L}=:x_0 .$$
6. **Taylor.** $x\log_2(1/x)$ increases on $(0,1/e)$, and $x_0\le x\le\frac18$. With $\log_2(1/x_0)\ge\ell/4$ for large $n$, (C2) gives
   $$\delta\ \ge\ x_0\log_2(1/x_0)\ \ge\ 0.0122\,\frac{\ell^2}{\sqrt L},$$
   a contradiction. ∎

**Corollary 3.3 (`T.ETH.1`).**
- **(i) The answer to the question.** Let $\frac14\le p\le\frac12$ and $R\le Q/4$. The least $\mu_p$-measure of $\{F\ne\mathrm{par}\}$ over real $F$ of degree $<R$ lies between:
  - $e^{-2(R+1)^2/Q}/(4(Q+1))$, for $Q\ge64$ (Lemma 2.4);
  - and $2e^{-R^2/(8Q)}$, for even $R\ge4$ (the window design).

  So it is $e^{-\Theta(R^2/Q)}$ up to a factor $\mathrm{poly}(Q)$.
- **(ii) The floor.** No graded design beats E4's floor. In particular no choice of $M_0$ in E1 does, symmetric or not.

Non-symmetric designs can still win constant factors. Part B finds bad sets up to 27% lighter than the window design's for $Q\le6$.

---

## 4. What the floor explains

1. **The square root is Smolensky's.** [13](13_UNIFIED_DEFORMATION.md) §5 read Part I's square root off the central limit window of the symmetric design. Corollary 2.5 shows the square root is intrinsic: a polynomial that agrees with parity except on measure $\varepsilon$ needs degree $\sqrt{Q\log(1/\varepsilon)}$, whatever its symmetry. The central limit theorem enters only through $\tau_D$, the binomial tail inside Smolensky's bound. This is the same dimension argument that Smolensky used against constant-depth circuits with modular gates.
2. **What a design is.** Lemma 1.2 turns a graded design into a family of polynomials $G_t$ of degree $<R$, one per node, each required to agree with parity on the swap spectrum. The deformation degree in $T$ becomes a Fourier degree in $X$. The swap representation enters only through the multiplicities $\mu_{p_r}$ of its isotypic pieces.
3. **A reusable lemma.** Lemma 2.4 is the biased-cube agreement bound for parity with explicit constants. It applies wherever an argument needs "a low-degree real polynomial cannot match parity on most of a biased cube".

---

## 5. Ungraded designs: what remains

For an arbitrary design, $N_t$ acts on $\psi$ with $X(\psi)=X$ by
$$\nu_t(X)=\tfrac12\big(n_X(t)+n_{X\oplus e_t}(t)\big),$$
as in the proof of Lemma 1.2(b). This is even in $x_t$. Call the edge $\{X,X\oplus e_t\}$ **bad** if $\nu_t(X)\ne0$.

**Lemma 5.1 (face lemma).** For every group-algebra design of degree $<R$, graded or not, every $R$-dimensional face of $\{0,1\}^Q$ contains a bad edge. Hence some node $t$ has $\mu_{1/2}(\nu_t\ne0)\ge2^{1-R}/R$.

*Proof.*
1. **One face.** Let the face $C$ have direction set $A$, $|A|=R$, and put $\sigma(T)=\sum_{Y\in C}n_Y(T)$. Then $\deg\sigma<R$ and $\sigma(0)=2^R$.
2. **No bad edge gives a contradiction.** If no edge of $C$ is bad, then for each $t\in A$ pair the vertices of $C$ along $t$: $\sigma(t)=2\sum_{\text{pairs}}\nu_t(Y)=0$. So $\sigma$ has $R$ roots and $\sigma(0)\ne0$, which is impossible for $\deg\sigma<R$.
3. **Counting.** Each edge lies in $\binom{Q-1}{R-1}$ of the $\binom QR2^{Q-R}$ faces. So there are $\ge\frac QR2^{Q-R}$ bad edges, and some direction has $\ge2^{Q-R}/R$ of them. Its endpoints have measure $\ge2^{1-R}/R$. ∎

**Why this alone is not enough.** The bound $2^{1-R}/R$ is far below $e^{-\Theta(R^2/Q)}$, and face-hitting by itself cannot do better.
- Let $R\ge2+\log_2Q$. A random vertex set of density $2^{2-R}Q$ meets all $<4^Q$ faces of dimension $R$ with probability $>1-e^{-2.6Q}$.
- With positive probability it also has density $\le2^{3-R}Q$.
- The edges at its vertices then form a set of density $\le2^{4-R}Q$ that meets every $R$-face.

So a proof for ungraded designs must use the algebra, not only Lemma 5.1.

**Question 5.2 (target `T.ETH.6`).** For every group-algebra design of degree $<R$ and $\frac14\le p\le\frac12$, is
$$\max_t\mu_p(\nu_t\ne0)\ \ge\ e^{-O(R^2/Q)}/\mathrm{poly}(Q)\ ?$$

Equivalently, write $n_X(T)=a(X)\cdot v_T$, with $a(X)\in\{1\}\times\mathbb R^{R-1}$ and $v_t=(1,t,\dots,t^{R-1})$. The edge $\{X,X\oplus e_t\}$ is good iff $(a(X)+a(X\oplus e_t))\cdot v_t=0$. The question asks whether some direction must have bad edges of $\mu_p$-measure $e^{-O(R^2/Q)}/\mathrm{poly}(Q)$, for every map $a$.

**Evidence.**
- Part E: on central windows, ungraded designs are feasible exactly where single-polynomial designs are, for $(Q,R)=(6,3)$, $(7,3)$ and (with `--full`) $(8,4)$.
- The node conditions of an ungraded design form an overdetermined linear system: about $Q2^{Q-1}$ equations in $(R-1)2^Q$ unknowns. The only source of large consistent subsystems I know is graded structure.

I expect a positive answer (subjective confidence ≈70%).

**Why it matters.** Suppose an ungraded design reached $\max_t\mu_p(\nu_t\ne0)=2^{-\Theta(R)}$. Then (C1)–(C3) would force no floor above order $(\log\log n)^3/\log n$, far below App. D's $\log\log n/\sqrt{\log n}$. So Question 5.2 decides whether Part I's architecture with group-algebra designs is capped at the square root. Designs outside the group algebra belong to the general question `T.ETH.5`.

---

## 6. Status

| Result | Depends on | Verdict |
|---|---|---|
| Definition 1.1, Lemma 1.2, Remark 1.3 | 08 Lemma 5.1, 08 Theorem 5.4 | complete derivation |
| Lemma 2.1 | — | complete derivation (Smolensky's argument); checked exactly (Part A) |
| Lemmas 2.2–2.4, Corollary 2.5 | Lemma 2.1; method of types; Chernoff | complete derivation; Lemma 2.2 checked exactly (Part B), Lemma 2.3 numerically (Part C) |
| Theorem 3.2, Corollary 3.3 | Lemmas 1.2, 2.4; cost model 3.1 | complete derivation given (C1)–(C3); (C3) is E4's |
| Lemma 5.1 | Definition 1.1 | complete derivation |
| Question 5.2 | — | open (`T.ETH.6`); evidence in Part E |

**Subjective confidence:**
- ≈98% that Lemmas 1.2, 2.1–2.4, 5.1 and Corollary 2.5 are correct as stated;
- ≈93% for Theorem 3.2, which inherits (C3) from E4 (≈95%).

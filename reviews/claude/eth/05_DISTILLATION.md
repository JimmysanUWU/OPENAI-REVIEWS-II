# 5. Distillation: the ETH-refuting core as a short blueprint

This is a self-contained rewrite of the minimal chain, organized as definition → theorem → proof. It follows the paper's mathematics with three structural changes:
1. a single notion, "local modulo measurements", replaces the four-part design definition;
2. the replacement rows are given by **one closed formula**, a polynomial in commuting *damped swaps* (Theorem B). This removes the $(u,y)$-dependent bookkeeping $E_{u,y}$;
3. the spectral reason for the construction, and its optimality, are stated explicitly.

The formula in Theorem B and its lemmas are verified in exact arithmetic in [`scripts/verify_distilled_swap_formula.py`](scripts/verify_distilled_swap_formula.py).

**Notation.**
- $[n]=\{1,\dots,n\}$; $e_v$ is a standard basis row.
- For $x\in\mathbb Q^n$ and $S\subseteq[n]$, $x|_S$ is the restriction and $\widetilde\xi$ is extension by zero.
- $\mathcal A$ is a commutative $\mathbb Q$-algebra.
- $C_q(n)=\sum_{j\le\min(q,n)}2^j\binom nj$.

---

## 5.1 Objects

**Definition 1 (factor sum).** A *factor system* consists of:
- scopes $I_f\subseteq[n]$, for $f\in\mathcal F$;
- polynomials $\Phi_f\in\mathcal A[Y_v:v\in I_f]$.

Its *degree* $d$ is the total degree of $\prod_f\Phi_f((Y_{f,v})_v)$ in independent indeterminates $Y_{f,v}$. Its *value* is $Z=\sum_{x\in\{0,1\}^n}\prod_f\Phi_f(x|_{I_f})$.

**Definition 2 (local modulo measurements).** A *measured partition* of $[n]$ consists of:
- a partition $[n]=S_1\sqcup\dots\sqcup S_r$;
- integer matrices $C_g\in\mathbb Z^{c_g\times S_g}$ whose rows have $\ell_1$-norm $\le D$.

Put $c=\sum c_g$. A row $p\in\mathbb Q^{1\times n}$ is **$g$-local** if $p|_{S_h}\in\operatorname{row}_{\mathbb Q}(C_h)$ for every $h\neq g$.

**Definition 3 (deformation).** A *deformation of degree $R$ with $Q$ nodes* for a factor system consists of:
- rows $P_{f,v}(T)\in\mathbb Q[T]^{1\times n}$ with $P_{f,v}(0)=e_v$ and $\deg P_{f,v}\le R$;
- for each node $t\in[Q]$, a measured partition $(S_{t,g},C_{t,g})_g$ and assignments $\sigma_f(t)\in[r]$;

such that every $P_{f,v}(t)$ is $\sigma_f(t)$-local.

## 5.2 Theorem A (deformation counting)

*Suppose a factor system of degree $d$ has a deformation of degree $R$ with $Q$ nodes, with group sizes $\le s$ and $c\le c_{\max}$ measurement rows at every node. Put $q=\lfloor dR/Q\rfloor$. Then $Z$ can be computed by examining at most*
$$2Qr\,C_q(n)\,(2D+1)^{c_{\max}}\,2^{s}$$
*Boolean vectors, with $\mathrm{poly}$ operations in $\mathcal A$ per vector and polynomial working space.*

**Lemma A1 (separation).** Let $\mathcal R$ be a commutative ring and $(S_g,C_g)$ a measured partition. Suppose $A_g:\mathbb Q^n\to\mathcal R$ is invariant under adding $\widetilde\xi$ for $\xi\in\ker C_h$, $h\ne g$. Then for every $E\subseteq[n]$,
$$\sum_{x\in\{0,1\}^n,\ x_E=1}\ \prod_gA_g(x)\;=\;\sum_{y}\ \prod_g\ \sum_{\substack{x_g\in\{0,1\}^{S_g},\ C_gx_g=y_g\\ x_g|_{E}=1}}A_g\Big(x_g+\sum_{h\ne g}\bar x_h(y_h)\Big),$$
where $y$ ranges over $\prod_g\{-D,\dots,D\}^{c_g}$, $\bar x_h(y_h)$ is any Boolean vector with $C_h\bar x_h=y_h$, and empty fibres contribute $0$.

*Proof.* Fix $y$ and $x$ with $C_hx_h=y_h$ for all $h$. Then $x_h-\bar x_h\in\ker C_h$, so replacing $x_h$ by $\bar x_h$ for each $h\ne g$ does not change $A_g$. Expanding the product and summing over the unique $y$ of each $x$ proves the identity. ∎

**Lemma A2 (jets by Boolean substitution).** Let $p_x\in\mathcal R[X_1,\dots,X_n]$ have degree $\le q$ for each $x\in\{0,1\}^n$. Then
$$\sum_xp_x(x)=\sum_{\lvert E\rvert\le q}\sum_{U\subseteq E}(-1)^{\lvert E\setminus U\rvert}\sum_{x_E=1}p_x(\mathbf 1_U).$$

*Proof.* Multilinearize $p_x$ (this does not change Boolean values). By Möbius inversion, the coefficient of $X^E$ is $\sum_{U\subseteq E}(-1)^{\lvert E\setminus U\rvert}p_x(\mathbf 1_U)$, and it vanishes for $\lvert E\rvert>q$. ∎

*Proof of Theorem A.* Let $H(T)=\sum_x\prod_f\Phi_f((\langle P_{f,v}(T),x\rangle)_v)$. Then $H(0)=Z$ and $\deg H\le dR$.

Fix a node $t$, let $J=q+1$ and $\mathcal R=\mathcal A[z]/(z^J)$, and put
$$p_x(X)=\prod_f\Phi_f\big((\langle P_{f,v}(t),x\rangle+\langle P_{f,v}(t+z)-P_{f,v}(t),X\rangle)_v\big).$$
Every entry of $P(t+z)-P(t)$ is divisible by $z$, so $\deg_Xp_x\le q$ over $\mathcal R$. Also $\sum_xp_x(x)=H(t+z)\bmod z^J$.

Apply Lemma A2. Each inner sum, with $X=\mathbf 1_U$ fixed, is a Lemma A1 sum with
$$A_g=\prod_{f:\sigma_f(t)=g}\Phi_f(\cdots).$$
Invariance holds because every $P_{f,v}(t)$ is $\sigma_f(t)$-local and the $X$-part is constant.

The resulting $J$ coefficients at $Q$ nodes determine $H$, since $QJ>dR$ (Hermite interpolation over $\mathbb Q$, applied coordinatewise in $\mathcal A$). The count of examined vectors follows. ∎

**Cost.** With $\eta=R/Q$, $\log_2C_q(n)\le q\log_2(1/\eta)+n\log_2(1+2\eta)$. Hence
$$\frac{\log_2(\text{time})}n\le\frac dn\frac RQ\log_2\frac QR+\frac{2R}{Q\ln2}+\frac{c_{\max}}n\log_2(2D+1)+\frac sn+o(1).$$

## 5.3 Theorem B (swap interpolation)

Let $\mathcal B=([r]\times[r])^Q$ with $D=r^{2Q}$, and write $u(s)=(a_u(s),d_u(s))$. Define the following.
- $\theta_s$ swaps the pair at position $s$, and $\Pi_s$ is the permutation matrix with $(\Pi_s\zeta)_u=\zeta_{\theta_su}$.
- $\Pi_S=\prod_{s\in S}\Pi_s$, and $S^c=[Q]\setminus S$.
- $\varphi_S(T)=\prod_{s\in S}(1-T/s)$.
- The Krawtchouk polynomials are $K_k(j)=\sum_i(-1)^i\binom ji\binom{Q-j}{k-i}$.
- $W\sim\mathrm{Bin}(Q,\frac{r-1}{2r})$.

*Let $f$ be a rational polynomial of degree $<R\le Q+1$, written as $f=\sum_{k<R}\gamma_kK_k$. Define*
$$M_0=\sum_{k<R}\gamma_k\!\!\sum_{\lvert S\rvert=k}\!\!\Pi_{S^c},\qquad
P(T)=I-M_0+\sum_{k<R}\gamma_k\!\!\sum_{\lvert S\rvert=k}\!\!\varphi_S(T)\,\Pi_{S^c}.$$

*Equivalently, $M_0=V_0\,e_f(\Pi_1,\dots,\Pi_Q)$ and $P(T)=I-M_0+V_0\,e_f\big((1-\tfrac T1)\Pi_1,\dots,(1-\tfrac TQ)\Pi_Q\big)$, where $V_0=\Pi_{[Q]}$ and $e_f=\sum_k\gamma_ke_k$ is a combination of elementary symmetric polynomials. Then:*

- **(B1)** *$P(0)=I$ and $\deg_TP\le R-1$.*
- **(B2)** *For $t\in[Q]$, the $(u,y)$ entry of $P(t)-(I-M_0)$ vanishes unless $a_y(t)=d_u(t)$.*
- **(B3)** *There is $C\in\{0,\pm1\}^{c\times\mathcal B}$ with $c=D\cdot\Pr[(-1)^Wf(W)\ne1]$ and row $\ell_1$-norms $\le2^Q$, such that $M_0\xi=\xi$ for every $\xi\in\ker C$. Moreover $\operatorname{rank}(I-M_0)=c$.*
- **(B4)** *For $t\in[Q]$, $u\in\mathcal B$ and $g\ne d_u(t)$: if $\xi$ is supported on $\{y:a_y(t)=g\}$ and $C\xi=0$, then $P(t)_u\,\xi=0$.*

*Proof.*
- **(B1).** $\varphi_S(0)=1$, and $\deg\varphi_S=\lvert S\rvert<R$.
- **(B2).** If $t\in S$ then $\varphi_S(t)=0$. Otherwise $\Pi_{S^c}$ swaps position $t$, so its only nonzero entry in row $u$ is at $y=\theta_{S^c}u$, which has $a_y(t)=d_u(t)$.
- **(B3).** At each position take the orthogonal $\{0,\pm1\}$ eigenbasis of the swap on $\mathbb Q^{r\times r}$: the $e_{\alpha\alpha}$ and $e_{\alpha\delta}+e_{\delta\alpha}$ (eigenvalue $+1$) and the $e_{\alpha\delta}-e_{\delta\alpha}$ (eigenvalue $-1$). Their tensor products $\psi$ form an orthogonal basis with $\lvert\operatorname{supp}\psi\rvert\le2^Q$.
  - If $\psi$ has eigenvalues $\epsilon_s$, with $j$ of them equal to $-1$, then $\sum_{\lvert S\rvert=k}\Pi_{S^c}\psi=(\prod_s\epsilon_s)\,e_k(\epsilon)\psi=(-1)^jK_k(j)\psi$. Hence $M_0\psi=(-1)^jf(j)\psi$.
  - Exactly $\binom Qj\binom r2^j\binom{r+1}2^{Q-j}=D\Pr[W=j]$ basis vectors have a given $j$.
  - Let the rows of $C$ be the $\psi$ with $(-1)^jf(j)\ne1$. By orthogonality, $\ker C$ is spanned by the remaining $\psi$, and $M_0$ fixes them.
- **(B4).** Here $(I-M_0)\xi=0$. By (B2), $P(t)_u-(I-M_0)_u$ is supported on $\{y:a_y(t)=d_u(t)\}$, which is disjoint from $\operatorname{supp}\xi$. ∎

**Corollary B′ (interval choice).** Let $f$ interpolate $(-1)^j$ on $R$ consecutive integers containing all integers within $R/4$ of $\mathbb EW$. Then $c\le2De^{-R^2/(8Q)}$ by Hoeffding.

**Proposition B″ (optimality among symmetric choices).** For every $f$ of degree $<R$, every window of $R+1$ consecutive integers contains $j$ with $(-1)^jf(j)\ne1$. *Proof:* the $R$-th difference of $(-1)^j$ is $\pm2^R\neq0$. Hence $c/D\ge\max_A\min_{A\le j\le A+R}\Pr[W=j]$, which matches B′ up to a $\mathrm{poly}(Q)$ factor. ∎

**Remark.** The paper's rows (4.10) use one multiplier $\prod_{s\in E_{u,y}}(1-T/s)$ per target coordinate, where the formula above uses one per group element. Both satisfy (B1)–(B4). Theorem B packages the paper's Lemmas 4.1–4.3 as a single identity in $\mathbb Q[T]\otimes\mathbb Q[(\mathbb Z/2)^Q]$.

## 5.4 Lemma C (transfer to a factor system)

*Let a factor system have $M=\sum_f\lvert I_f\rvert$ incidences, with $n\le M$. Choose schedules $a_v,d_f\in[r]^Q$ and labels $u_{f,v}=(a_v(s),d_f(s))_s\in\mathcal B$, and let $h$ be the maximum label multiplicity. Place each incidence in its own copy $c(f,v)\in[h]$ of $\mathcal B$. Define the filling map $\iota$ and the projections $\pi_c$ by*
$$(\iota x)_{c,u}=x_v\ \text{ if }(c,u)=(c(f,v),u_{f,v}),\qquad (\pi_c\eta)_u=\eta_{c,u}.$$
*Set*
$$P_{f,v}(T)=P(T)_{u_{f,v}}\,\pi_{c(f,v)}\,\iota,\qquad S_{t,g}=\{v:a_v(t)=g\},\qquad \sigma_f(t)=d_f(t),$$
*and let $C_{t,g}=(\mathrm{diag}_h C)\,\iota$, restricted to the columns $S_{t,g}$. These data form a deformation of degree $R-1$ with $Q$ nodes. It has $rhc$ measurement rows per node, each with $\ell_1$-norm $\le2^Q$.*

*Proof.* At $T=0$ the row reads its own position, so $P_{f,v}(0)=e_v$. Merging positions that hold the same variable cannot increase $\ell_1$-norms.

For locality, take $\xi$ supported on $S_{t,g}$ with $C_{t,g}\xi=0$. Then every $\pi_c\iota\widetilde\xi$ is supported on $\{y:a_y(t)=g\}$, because the first schedule of a label is its variable's schedule, and it lies in $\ker C$. Apply (B4) in copy $c(f,v)$. ∎

**Lemma C′ (balanced schedules).** *If $2(b+w)D^2\le M$ and $2Qr^2\le n$, where $b$ bounds the incidences per variable and $w$ the scope sizes, one can choose schedules deterministically in polynomial time with $h\le2M/D$ and $\lvert S_{t,g}\rvert\le2n/r$.*

*Proof.* Use the potential $\Phi/(M(b+w))+\Xi/(Qn)$, where
$$\Phi=\sum_u(c_u-M/D)^2,\qquad \Xi=\sum_{t,g}(n_{t,g}-n/r)^2.$$
Its expectation under uniform schedules is $\le2$. Fix the schedules one at a time without increasing the conditional expectation. ∎

With $\varepsilon=c/D$, the measurement rows number $rhc\le2r\varepsilon M\le2r\varepsilon bn$.

## 5.5 Theorem D (bounded occurrence)

*Let $F$ be a CNF in which each clause has at most 3 literals and each variable at most 3 occurrences. Then $\#F$ can be computed in $2^{O(n(\log\log n)^2/\sqrt{\log n})}$ bit operations and polynomial space.*

*Proof.* Add constant-one factors for unused variables; then $d\le M\le3n$ and $b,w\le4$. Let $\lambda=\lceil\log_2n\rceil$ and $\nu=\lceil\log_2\lambda\rceil$, and choose:

| parameter | value | role |
|---|---|---|
| $r$ | $\lfloor\sqrt\lambda/\nu^2\rfloor$ | group size $2n/r$ |
| $Q$ | $\lfloor\lambda/(128\nu)\rfloor$ | nodes; $D=r^{2Q}\le n^{1/128}$ |
| $R$ | $4\lceil\sqrt\lambda\rceil$ | degree; $R^2/(8Q)\ge256\nu$ |
| $\tau$ | $1/(r\lambda)$ | measurement budget $2\tau n$ |
| $\varepsilon$ | $\tau/(rb)$ | rank fraction, $\ge2e^{-R^2/(8Q)}$ |

Theorems A and B and Lemmas C and C′ apply for large $n$. The three exponents are
$$\frac RQ\log\frac QR=O\Big(\frac{\nu^2}{\sqrt\lambda}\Big),\qquad \tau\log(2D+1)=O\Big(\frac1r\Big),\qquad \frac1r=O\Big(\frac{\nu^2}{\sqrt\lambda}\Big).$$
The bit complexity is polynomial per examined vector, because all denominators divide $\big((R-1)!\,Q!\big)^{O(M)}$. ∎

## 5.6 Corollary E (ETH is false)

*For each fixed $k$, #$k$-SAT can be solved in $2^{O_k(n(\log\log n)^{k+1}/(\log n)^{1/(2k)})}\mathrm{poly}(n+m)$ time and polynomial space.*

*Proof.*
1. Disjoint sparsification (App. A) with $\sigma=(\log n)^{-1/(2k)}$ gives $\mathrm{poly}\cdot2^{\sigma n}$ leaves, with occurrence bound $B=O_k((\sigma^{-1}\log\sigma^{-1})^{k-1})$.
2. OR-gate splitting and Tovey cycles convert each leaf, preserving the count, into a bounded-occurrence 3-CNF on $O_k(Bn)$ variables.
3. Apply Theorem D and add the leaf counts. ∎

## 5.7 Modules outside the minimal chain

- **Sparse $\mathbb F_2$ systems and parity tests (§5–6).** These remove the sparsification loss and give $2^{O(\rho(n))}$ at every density.
  - New lemma: *cross-equation 2-adic truncation*. For integers $y_i$,
    $$\mathcal S_\ell\Big(\prod_iB_i(V)\bmod V^\ell\Big)\equiv\prod_i\mathbf 1_{2\mid y_i}\pmod{2^\ell},$$
    and the retained degree is $m+2(\ell-1)$.
  - Standard ingredients: Tanner codes on $\varepsilon$-biased Cayley graphs; greedy rows chosen by exact counts.
- **Appendix D.** A second deformation (Cauchy kernel, 2-adic taper, modes) that saves one $\log\log n$ factor.

## 5.8 Open distillation targets

1. **Agreement degree of parity under the $p$-biased measure.** A lower bound of $\Omega(\sqrt{Q\log(1/\varepsilon)})$ for *all* (non-symmetric) group-algebra choices would make Proposition B″ unconditional within the group-algebra framework.
2. **A unified deformation theorem.** Theorems A and B, and Appendix D's $W(T)$, are all "identity at 0, local at nodes" deformations: one controlled by degree, the other by 2-adic precision. A common abstract statement, with a cost functional, would explain both $\sqrt{\log}$ barriers at once.
3. **Measurement range.** Using the $\ell_1\le2^Q$ bound together with a balanced assignment achieving $D\approx n/\mathrm{polylog}$ (instead of $\sqrt n$) changes only constants. A cleaner parameter table could be extracted from this.

# 7. Recyclable toolkit: what in the ETH paper is reusable, and how far it sharpens

This file extracts from *The Exponential Time Hypothesis Is False* the parts that are useful **beyond** this paper. Each is stated at the greatest generality its proof supports, and sharpened where possible. Status labels:
- **P** — proved in the paper;
- **R** — proved here;
- **X** — checked exactly by a script in [`scripts/`](scripts/);
- **H** — heuristic, clearly marked.

| # | Tool | General form | Status | What is new here |
|---|---|---|---|---|
| E1 | Commuting-involution interpolation | any $Q$ commuting involutions on any space | R, X | Abstract form; closed-form rows; exact rank; optimality among symmetric choices |
| E2 | $p$-adic cross-equation truncation | $\mathbb F_p$ for every prime $p$ | R, X | **New:** generalizes Lemmas 5.1–5.2 from $p=2$ |
| E3 | Deformation counting | weighted factor graphs over any commutative $\mathbb Q$-algebra | R | **New corollary:** exact partition functions of bounded-degree graphical models in $2^{o(n)}$ |
| E4 | Cost floor of the method | any symmetric design | R, X | **New:** the $(\log\log n)^2/\sqrt{\log n}$ barrier becomes a theorem about the method |
| E5 | Tapered Hermite recovery (App. D) | any polynomial over $\mathbb Z_{(2)}$-algebras | P, X | Isolated as a standalone numerical-algebra lemma |

---

## E1. Commuting-involution interpolation

**Theorem E1.** Let $V$ be a finite-dimensional vector space over a field of characteristic 0, and let $\Pi_1,\dots,\Pi_Q$ be commuting involutions on $V$. For $\chi\in\{0,1\}^Q$ let $V_\chi$ be the joint eigenspace with $\Pi_s=(-1)^{\chi_s}$. Write $\Pi_S=\prod_{s\in S}\Pi_s$ and $S^c=[Q]\setminus S$.

For $f(x)=\sum_{k<R}\gamma_kK_k(x)$ (Krawtchouk basis, any degree $<R$), set
$$M_0=\sum_{k<R}\gamma_k\!\!\sum_{|S|=k}\!\!\Pi_{S^c},\qquad P(T)=I-M_0+\sum_{k<R}\gamma_k\!\!\sum_{|S|=k}\!\!\Big(\prod_{s\in S}\big(1-\tfrac Ts\big)\Big)\Pi_{S^c}.$$
Then:
1. $P(0)=I$ and $\deg_TP\le R-1$.
2. For each node $t\in[Q]$, $P(t)-(I-M_0)$ lies in the span of the group elements $\Pi_U$ with $t\in U$.
3. $M_0=(-1)^{|\chi|}f(|\chi|)$ on $V_\chi$. Hence $\ker(I-M_0)\supseteq\bigoplus_{(-1)^{|\chi|}f(|\chi|)=1}V_\chi$, and
$$\operatorname{rank}(I-M_0)=\sum_{\chi:\,(-1)^{|\chi|}f(|\chi|)\ne1}\dim V_\chi .$$
4. (Optimality.) Every window of $R+1$ consecutive integers contains a $j$ with $(-1)^jf(j)\neq1$.

*Proof.*
1. $\prod_{s\in S}(1-T/s)$ vanishes at $T=0$ only through the factor $T$; it equals 1 at $T=0$ and has degree $|S|<R$.
2. The factor $\prod_{s\in S}(1-t/s)$ vanishes when $t\in S$; otherwise $t\in S^c$.
3. On $V_\chi$, $\sum_{|S|=k}\Pi_{S^c}=(-1)^{|\chi|}\sum_{|S|=k}(-1)^{|\chi\cap S|}=(-1)^{|\chi|}K_k(|\chi|)$.
4. The $R$-th finite difference of $f$ vanishes, while that of $(-1)^j$ is $\pm2^R$. ∎

**How the ETH paper instantiates it.** $V=\mathbb Q^{([r]\times[r])^Q}$, with $\Pi_s$ swapping the pair in coordinate $s$. Then "$\Pi_U$ with $t\in U$" means "maps every label $u$ to a label whose first schedule agrees with $u$'s second schedule at $t$", which is exactly the group-locality needed by Definition 3.1. With uniform labels, $\dim V_\chi$ is $D$ times a $\mathrm{Bin}(Q,\frac{r-1}{2r})$ probability.

**Where else it applies.** Theorem E1 is a statement about one representation of $(\mathbb Z/2)^Q$. It gives a degree-$<R$ path from $I$ to "elements involving generator $t$" at every node, at the price of a low-rank defect. Any problem that needs simultaneous locality at many nodes with a low-degree interpolating family can use it with its own $V$ and $\Pi_s$. In the regular representation it is the group-algebra form of the parity-approximation phenomenon behind Walsh–Hadamard non-rigidity.

**Open sharpening.** Item 4 covers only symmetric choices. For a general group-algebra element supported on $\{\Pi_{S^c}:|S|<R\}$, the eigenvalue is $(-1)^{|\chi|}g(\chi)$ with $g$ an arbitrary real multilinear polynomial of degree $<R$. The question becomes: *what is the least $p$-biased measure of the set where a degree-$<R$ real polynomial differs from parity?* Every $R$-dimensional subcube contains such a point, but that alone gives only $\min(p,1-p)^R$, far below the symmetric bound $e^{-\Theta(R^2/Q)}$. This is open.

(Verified: [`scripts/verify_distilled_swap_formula.py`](scripts/verify_distilled_swap_formula.py).)

---

## E2. $p$-adic cross-equation truncation (new)

**Theorem E2.** Let $p$ be prime and $\ell\ge1$. Put
$$F_\ell(U)=(1-U)^\ell\sum_{j<\ell}\binom{\ell+j-1}{j}U^j,\qquad p_\ell(Y)=F_\ell(Y^{p-1}),\qquad h(Y)=Y^p-Y,$$
and write $p_\ell=\sum_je_j(Y)h(Y)^j$ with $\deg e_j<p$ (division by the monic $h$). For integers $y_1,\dots,y_m$ let
$$B_i(V)=\sum_{j<\ell}e_j(y_i)h(y_i)^jV^j .$$
Then
$$\mathcal S_\ell\Big(\prod_{i=1}^mB_i(V)\bmod V^\ell\Big)\equiv\prod_{i=1}^m[p\mid y_i]\pmod{p^\ell},$$
where $\mathcal S_\ell$ sums the coefficients of $1,V,\dots,V^{\ell-1}$. Every retained coefficient is a polynomial of degree $\le m(p-1)+p(\ell-1)$ in $(y_i)$. The untruncated product has degree $m(p-1)(2\ell-1)$.

*Proof.*
1. **Two-sided identity.** $F_\ell(U)=1-U^\ell\sum_{j<\ell}\binom{\ell+j-1}j(1-U)^j$. This is the standard identity, from splitting $\sum_{a+b=2\ell-1}$ by whether $a\ge\ell$. Hence $F_\ell\equiv1\bmod U^\ell$ and $F_\ell\equiv0\bmod(1-U)^\ell$.
2. **The congruence $p_\ell(y)\equiv[p\mid y]\pmod{p^\ell}$.** If $p\mid y$, then $p^{p-1}\mid y^{p-1}$, so $U^\ell\equiv0\bmod p^\ell$ and $p_\ell(y)\equiv1$. If $p\nmid y$, then $y^{p-1}\equiv1\bmod p$ (Fermat), so $(1-U)^\ell\equiv0\bmod p^\ell$ and $p_\ell(y)\equiv0$.
3. **Truncation is exact.** Since $p\mid h(y)$ for every integer $y$, every monomial of $\prod_iB_i$ with total $h$-power $\ge\ell$ is $\equiv0\bmod p^\ell$. So $\mathcal S_\ell$ of the truncation is congruent to $\prod_ip_\ell(y_i)$.
4. **Degree.** The coefficient of $V^j$ in $B_i$ has degree $\le(p-1)+pj$. A retained product has $\sum j_i\le\ell-1$. ∎

For $p=2$ this is exactly Lemmas 5.1–5.2 of the paper ($h=Y^2-Y$). Verified exactly for $p\in\{2,3,5,7\}$, $\ell\le6$, random integer tuples; the degree bound is attained ([`scripts/verify_padic_truncation.py`](scripts/verify_padic_truncation.py), [`results/padic_truncation.json`](results/padic_truncation.json)).

**What it buys.** §5 of the paper counts common zeros of sparse $\mathbb F_2$ systems by lifting each equation to $\mathbb Z/2^\ell$ while keeping *linear* total degree. Theorem E2 does the same over $\mathbb F_p$: the indicator of $\{P_1=\dots=P_m=0\}$ over $\mathbb F_p$ lifts to $\mathbb Z/p^\ell$ with degree $m(p-1)+p(\ell-1)$. This is the input needed to run the §5 pipeline for **counting zeros of sparse polynomial systems over any fixed prime field** (domain $\{0,\dots,p-1\}$, encoded by bits).

*Status of that corollary:* blueprint. Theorem E2 replaces Lemma 5.2. The deformation counting (E3) is domain-agnostic once variables are encoded in binary. I have not re-derived the parameter bookkeeping of §5.5 for $p>2$.

---

## E3. Deformation counting for weighted factor graphs

**Theorem E3** ([05](05_DISTILLATION.md) Theorem A, restated). Let $\mathcal A$ be a commutative $\mathbb Q$-algebra. Suppose a factor system with polynomials $\Phi_f\in\mathcal A[Y_v:v\in I_f]$ has a deformation of degree $R$ with $Q$ nodes (Definition 3). Then $\sum_{x\in\{0,1\}^n}\prod_f\Phi_f(x|_{I_f})$ is computed with the stated cost. *Nothing in the proof uses that $\Phi_f$ is a clause polynomial.*

**Corollary E3′ (new).** Fix $k,\Delta$. Let a factor graph on $n$ Boolean variables have arity $\le k$ and each variable in $\le\Delta$ factors, with rational factor tables of total bit length $L$. Its partition function $Z=\sum_x\prod_fw_f(x|_{I_f})$ is computable **exactly** in
$$2^{O_{k,\Delta}(n(\log\log n)^2/\sqrt{\log n})}\cdot\mathrm{poly}(L)$$
time and polynomial space.

*Proof.*
- Let $\Phi_f$ be the multilinear interpolant of $w_f$ (degree $\le k$). Then $d\le kM\le k\Delta n$ and $b,w=O(1)$.
- The schedules of Lemma 4.5 and the swap design (E1 with Lemma 4.8 parameters) apply verbatim.
- All arithmetic is in $\mathbb Q$. Clearing denominators multiplies each term by a common denominator of $\mathrm{poly}(L)$ bits. ∎

Examples: Ising and Potts partition functions with rational couplings on bounded-degree graphs (binary-encoded spins), and weighted #CSP. Exact evaluation of these is #P-hard, and #ETH-based lower bounds of $2^{\Omega(n)}$ for them are now void. The corollary is immediate from the paper's machinery but is not stated there. It shows the method is a **partition-function algorithm**, not a SAT algorithm.

---

## E4. The cost floor of the method (new: the barrier becomes a theorem)

**Proposition E4.** Instantiate the §3–4 algorithm with:
- $r$ groups, $Q$ nodes, labels $D=r^{2Q}\le n$;
- a $\mathrm{Sym}(Q)$-invariant swap design of degree $<R$ (any $f$ in E1);
- the cost bound (3.15).

Then its exponent per variable satisfies
$$\delta\ \ge\ c\,\frac{(\log\log n)^2}{\sqrt{\log n}}$$
for an absolute $c>0$ and all large $n$. The paper's parameters attain this order, so **its exponent is optimal for this architecture.**

*Proof.* The bound (3.15) gives three costs:
- $\delta\ge1/r$ (a largest group has $\ge n/r$ variables);
- $\delta\ge Q\varepsilon_*/32$, where $\varepsilon_*=\Pr[W=j_*]$ and $j_*$ is the bad level (E1(4)) inside the window $|j-Qp|\le R/2+1$. Reason:
  - the rows of $C^{\mathcal B}$ at level $j_*$ number $D\varepsilon_*$;
  - each has $j_*\ge Q/16$ antisymmetric local factors (using $Qp\ge Q/4$ and $R\le Q/e$), so every nonzero restriction to a group has $\ell^1$-norm $\ge2^{Q/16-1}$ and ranges over $\ge Q/16$ bits of values;
  - every row has a nonzero restriction to some group;
  - there are $h\ge M/D\ge n/D$ copies of each row;
- $\delta\ge x\log_2(1/x)$, where $x=R/Q\le1/e$ (Taylor jets: $\log_2\binom nq\ge q\log_2(n/q)$ with $q\ge nR/Q$; if $x>1/e$ then $\delta$ is bounded below by a constant).

Suppose $\delta\le\delta_0:=(\log\log n)^2/\sqrt{\log n}$.
1. Then $r\ge1/\delta$, so $\log r\ge\frac14\log\log n$, and $Q\le\log n/(2\log r)$.
2. The method-of-types bound $\Pr[W=j]\ge(Q+1)^{-1}e^{-Q\,\mathrm{KL}(j/Q\|p)}$, with $\mathrm{KL}\le\chi^2$ and $p(1-p)\ge3/16$, gives $\varepsilon_*\ge(Q+1)^{-1}\exp(-\frac{16}3(R/2+1)^2/Q)$. ($R\ge2$: if $R=1$ then $\varepsilon_*\gtrsim Q^{-1/2}$ and $rQ\varepsilon_*\gg1$.)
3. Combining with $\varepsilon_*\le32\delta/Q$ gives $\frac{16}3(R/2+1)^2/Q\ge\log\frac{Q}{32\delta(Q+1)}\ge\log\frac1{64\delta}$. Since $R/2+1\le R$, this gives $R^2\ge\frac3{32}Q\log\frac1\delta$ for small $\delta$.
4. Hence $x=R/Q\ge\sqrt{\tfrac3{32}\log(1/\delta)/Q}\ge x_0:=\sqrt{\tfrac3{16}}\,\log(1/\delta)/\sqrt{\log n}\ge c_1\log\log n/\sqrt{\log n}$, using $Q\le\log n/(2\log r)$ and $\log r\ge\log(1/\delta)$.
5. Since $x\log(1/x)$ increases on $(0,1/e)$, $\delta\ge x_0\log_2(1/x_0)\ge c_2(\log\log n)^2/\sqrt{\log n}$. ∎

**Numerical companion** ([`scripts/verify_cost_floor.py`](scripts/verify_cost_floor.py), [`results/cost_floor.json`](results/cost_floor.json)). Minimizing the three costs over integer $(r,Q,R)$, with the window bound computed exactly by the method of types (the script uses the paper's own measurement count $r\varepsilon Q$; the floor model $Q\varepsilon/32$ only lowers constants), gives a ratio to $(\log_2\log_2n)^2/\sqrt{\log_2n}$ in $[0.11,0.52]$ for $\log_2n=2^{10},\dots,2^{300}$. The ratio rises slowly toward a constant.

**Remark (thresholds are not intrinsic).** The cost floor drops below $1$ (beats $2^n$) already at $\log_2n\approx2^{10}$. The paper's parameters need $\log_2n\approx7\cdot10^{10}$ ([`results/parameters.json`](results/parameters.json)). The floor ignores the balance conditions of Lemma 4.5 and the clause-degree factor $d/n\le3$, so it is a lower bound, not an achievable threshold. Still, the astronomically large threshold is mostly a matter of the paper's constants, not of the method.

Beyond this floor, gains need non-symmetric designs (the open part of E1), or a different decoupling. Appendix D (E5) removes one $\log\log n$ by replacing subset counting with 2-adic modes.

---

## E5. Tapered Hermite recovery (Appendix D, Lemma D.5) as a standalone lemma

**Lemma E5.** Setting:
- $L=2^a$, nodes $S=\{-L,\dots,L-1\}$ (a complete residue system mod $2^{a+1}$), and $\alpha=2^{a+1}$;
- multiplicities $q_t=h\cdot\max(0,u-v_2(t))$ with $h=\lceil3M/u\rceil$;
- $F\in K[T]$ with $\deg F\le M(2L-1)$, $F(0)\in\mathcal A$, and $F(t+\alpha z)\in\mathcal A[z]$, where $\mathcal A$ is a $\mathbb Z_{(2)}$-order.

Then $F(0)\bmod2^p$ is determined by the first $q_t$ coefficients of each $F(t+\alpha z)\bmod2^{p_t}$, with $p_t=p+8hD_t+8M(a+1)$ and $D_t=L/2^{v_2(t)}$, via the explicit Chinese-remainder formula (D.9).

The point to recycle is that **precision demand tapers with 2-adic depth**. Nodes with large $v_2(t)$ get small multiplicity, and the loss identity (D.10) is $\ge-6D_th$. This is a general device for exact interpolation over 2-adic rings at many nodes without uniform precision blow-up. Verified: recovery succeeds in 46/46 random trials, and fails in 46/46 with untapered precision ([`results/appendix_d.json`](results/appendix_d.json)).

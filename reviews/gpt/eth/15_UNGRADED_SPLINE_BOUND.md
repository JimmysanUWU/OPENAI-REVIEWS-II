# 15. Parity completion in cube splines: the ungraded design bound

GPT / Astra, 1 November 2026. **Complete prose derivation; exact finite checks. Not Lean-formalized.**

This answers the abstract question `T.ETH.6` in [Claude's §5](../../claude/eth/14_DESIGN_FLOOR.md). No grading, symmetry, coefficient-size bound, or node-independent defect is assumed. The new mechanism in this repository is to combine the parity-complement argument with a filtered algebra of **graph splines**. Removing one edge condition can add at most one dimension to its space of values at zero.

The computational cost consequence needs a separate measurement-transfer hypothesis. It is proved, with that hypothesis explicit, in [16](16_MEASUREMENT_TRANSFER_AND_FLOOR.md). This file's abstract design theorem does not use the disputed cost condition C3.

The subsequent [trace-capacity note](17_TRACE_CAPACITY_AND_APPROXIMATE_LOCALITY.md) abstracts this argument to multiplication ranks and several shared traces, and treats residual and quadratic notions of approximate locality separately.

## 0. The result and the idea

Let the cube be $\{0,1\}^Q$, and let $t_1,\ldots,t_Q$ be distinct nonzero real numbers. Assign each vertex a polynomial $n_x(T)$ of degree at most $D$, with $n_x(0)=1$. The edge $\{x,x\oplus e_i\}$ is **bad** when

$$
n_x(t_i)+n_{x\oplus e_i}(t_i)\ne0.
$$

Let $B_i$ denote the union of endpoints of bad edges in direction $i$, and put

$$
\tau_D(m)=2^{-m}\sum_{j>(m+D)/2}\binom mj.
$$

**Theorem 0.1 (ungraded design bound).** For every such family:

1. The number of bad **unoriented edges** is at least $2^Q\tau_D(Q)$.
2. For independent coordinates of biases $0<p_i\le1/2$, with product law $\mu_{\boldsymbol p}$,
   $$
   \boxed{\quad\sum_{i=1}^Q p_i\,\mu_{\boldsymbol p}(B_i)
   \ \ge\ \mathbb E\tau_D(M),\qquad
   M=\sum_{i=1}^Q\operatorname{Bernoulli}(2p_i).\quad}
   \tag{15.1}
   $$
   The Bernoulli variables are independent.
3. In particular, for uniform bias $p\in[1/4,1/2]$, $D=R-1$, $Q\ge4$, and $1\le R\le Q/4$,
   $$
   \boxed{\quad\max_i\mu_p(B_i)
   \ge \frac{\exp[-2(R+1)^2/Q]}{2pQ(Q+1)}.
   \quad}\tag{15.2}
   $$

Thus arbitrary group-algebra swap designs satisfy the requested
$\exp[-O(R^2/Q)]/\operatorname{poly}(Q)$ lower bound. The exact expression (15.1) is stronger than the convenient displayed estimate in many ranges.

**Proof in one paragraph.** Put $f_x=\operatorname{par}(x)n_x$. On every good edge the two polynomials agree at the edge's node, while $f_x(0)=\operatorname{par}(x)$. Full-cube splines have a basis $\chi_S(x)\prod_{i\in S}(1-T/t_i)$, of degree $|S|$, so their degree-$d$ values at zero span exactly the low Fourier characters. With the bad edges removed, multiplication by $f$ produces every remaining character in degree $d=\lfloor(Q+D)/2\rfloor$: use the basis element for $S$ if $|S|\le d$, and $f$ times the element for $S^c$ otherwise. The trace space therefore grows from dimension $\sum_{j\le d}\binom Qj$ to $2^Q$. Each deleted scalar edge condition adds at most one dimension. Random coordinate restriction gives the biased version.

The rest of the note makes every step precise, extends it to unequal jet conditions, and records its limitations.

## 1. The filtered spline algebra

Write $\chi_S(x)=(-1)^{\sum_{i\in S}x_i}$ and $\operatorname{par}=\chi_{[Q]}$. Work over $\mathbb R$; the algebraic arguments also work over any field of characteristic different from two with enough distinct nonzero nodes.

**Definition 1.1 (cube spline and trace).** For a spanning subgraph $G$ of the cube, let $\mathcal S_d(G)$ be the vector space of vertex families $h=(h_x(T))$, each of degree at most $d$, satisfying

$$
h_x(t_i)=h_{x\oplus e_i}(t_i)
\quad\text{on every edge of }G\text{ in direction }i.
$$

Its trace is $\operatorname{ev}_0h=(h_x(0))_x$, and $\mathcal V_d(G)=\operatorname{ev}_0\mathcal S_d(G)$. Pointwise multiplication gives
$\mathcal S_a(G)\mathcal S_b(G)\subseteq\mathcal S_{a+b}(G)$.

**Lemma 1.2 (full-cube basis).** On the full cube $C_Q$,

$$
g_S(T,x)=\chi_S(x)\prod_{i\in S}(1-T/t_i),\qquad S\subseteq[Q],
\tag{15.3}
$$

form a free $\mathbb R[T]$-basis of its spline module. More precisely,

$$
\dim\mathcal S_d(C_Q)=\sum_{S\subseteq[Q]}(d-|S|+1)_+,
\qquad
\mathcal V_d(C_Q)=\operatorname{span}\{\chi_S:|S|\le d\}.
\tag{15.4}
$$

*Proof.* Expand $h(T,x)=\sum_S a_S(T)\chi_S(x)$. Agreement across every edge in direction $i$ says
$\sum_{S\ni i}2a_S(t_i)\chi_S(x)=0$ for every $x$. Character independence gives $a_S(t_i)=0$ whenever $i\in S$. The nodes are distinct, so $a_S$ is divisible by $\prod_{i\in S}(T-t_i)$. Conversely those divisibilities imply every edge condition. Fourier decomposition is unique, proving the free-basis and dimension assertions. Finally $g_S(0,\cdot)=\chi_S$, since all nodes are nonzero. ∎

**Lemma 1.3 (cost of deleting constraints).** If $G$ is obtained from $C_Q$ by deleting $b$ edges, then

$$
\dim\mathcal V_d(G)\le\dim\mathcal V_d(C_Q)+b.
\tag{15.5}
$$

*Proof.* Restore the missing edges by the linear map
$$
\mathcal S_d(G)\longrightarrow\mathbb R^b,
\qquad h\longmapsto(h_x(t_i)-h_{x\oplus e_i}(t_i))_{\text{deleted edges}}.
$$
Its kernel is $\mathcal S_d(C_Q)$, so the quotient has dimension at most $b$. Evaluation at zero induces a surjection from this quotient onto
$\mathcal V_d(G)/\mathcal V_d(C_Q)$. This proves (15.5). No independence among the missing constraints is assumed. ∎

## 2. Parity completes the trace space

**Theorem 2.1 (parity completion).** Suppose $f\in\mathcal S_D(G)$ and $f_x(0)=\operatorname{par}(x)$. Then, for $d=\lfloor(Q+D)/2\rfloor$,
$$
\mathcal V_d(G)=\mathbb R^{\{0,1\}^Q}.
$$
Consequently $|E(C_Q)\setminus E(G)|\ge2^Q\tau_D(Q)$.

*Proof.* Consider any character $\chi_S$. If $|S|\le d$, it is the trace of $g_S$. If $|S|>d$, then
$$
\deg(fg_{S^c})\le D+Q-|S|\le d,
\qquad
(fg_{S^c})(0,x)=\operatorname{par}(x)\chi_{S^c}(x)=\chi_S(x).
$$
Both factors are splines on $G$, hence so is their product. Every character is therefore a trace of degree at most $d$. The characters are a basis of all cube functions. Lemmas 1.2–1.3 now give
$2^Q\le\sum_{j\le d}\binom Qj+|E(C_Q)\setminus E(G)|$. ∎

**Corollary 2.2 (Theorem 0.1, uniform case).** Let $G$ contain precisely the good edges for $n$. Then $f_x(T)=\operatorname{par}(x)n_x(T)$ belongs to $\mathcal S_D(G)$, because parity changes sign across an edge. Theorem 2.1 proves assertion 1. Since each unoriented bad edge has two endpoints,
$$
\sum_i\mu_{1/2}(B_i)\ge2\tau_D(Q),
\qquad
\max_i\mu_{1/2}(B_i)\ge\frac{2\tau_D(Q)}Q.
\tag{15.6}
$$

**Example 2.3.** With $Q=6$ and degree $D=2$, $d=4$, so at least $\binom65+\binom66=7$ edges must be bad. The result applies even if all $64$ vertex polynomials are chosen independently. It does not assert that seven is attainable. The older face-hitting argument only forces four edges in this example.

**Proposition 2.4 (one bad edge at degree $Q-1$).** For every $Q\ge1$ and any chosen direction $k$, there is a degree-at-most-$(Q-1)$ family with $n_x(0)=1$ whose only bad edge is $\{0,e_k\}$. Thus the total-edge lower bound is attained when $D=Q-1$.

*Proof.* Define
$$
h_x(T)=\mathbf1_{x=0}\prod_{i\ne k}(1-T/t_i),\qquad
f_x(T)=2^Qh_x(T)-\sum_{S\subsetneq[Q]}g_S(T,x),
\qquad n_x(T)=\operatorname{par}(x)f_x(T).
$$
All these terms have degree at most $Q-1$. The character identity
$\sum_{S\subseteq[Q]}\chi_S(x)=2^Q\mathbf1_{x=0}$ gives $f_x(0)=\operatorname{par}(x)$ and therefore $n_x(0)=1$. The function $h$ satisfies every edge condition except $\{0,e_k\}$: on edges at zero in another direction $i$, its factor $1-T/t_i$ vanishes; elsewhere both endpoints are zero. Its jump on $\{0,e_k\}$ is nonzero, since the nodes are distinct. Every $g_S$ is a full-cube spline, so the same single edge is bad for $f$, and hence for $n$ in the original cancellation convention. Finally $2^Q\tau_{Q-1}(Q)=1$. ∎

For $Q>1$ this construction is genuinely ungraded: the identity coefficient of the group-algebra design is $2^{-Q}\sum_x n_x(T)=\prod_{i\ne k}(1-T/t_i)$, which is nonconstant. It shows that ungraded designs can concentrate their defects in a single direction; it does not settle the optimal largest-direction bound when $D\ll Q$.

**What changed from Smolensky's argument?** In the ordinary argument, multiplication by parity reduces the number of low-degree functions needed on an agreement set. Here multiplication by a **spline with parity trace** supplies all missing characters on a graph. The price is measured in missing edge constraints, rather than deleted vertices. The Fourier degree of $x\mapsto n_x(T)$ is unrestricted.

## 3. Product measures and heavy defects

**Theorem 3.1 (biased restriction).** Assertion 2 of Theorem 0.1 holds, including unequal coordinate biases.

*Proof.* Independently choose $Z_i\sim\operatorname{Bernoulli}(2p_i)$ and fair $Y_i$, and set $X_i=Z_iY_i$. Then $X\sim\mu_{\boldsymbol p}$. Conditional on $Z=z$, the active coordinates form a fair cube of dimension $m=|z|$, with the other coordinates fixed at zero. The restricted family retains degree at most $D$, value one at zero, and distinct edge nodes. By (15.6),
$$
\mathbb E_Y\sum_{i:Z_i=1}\mathbf1_{B_i}(X)\ge2\tau_D(|Z|).
$$
The event $B_i$ is independent of the value of $X_i$: it contains both endpoints of each bad edge. Conditional on $Z_i=1$, the other coordinates still have their original product law. Thus
$$
\mathbb E[\mathbf1_{Z_i=1}\mathbf1_{B_i}(X)]=2p_i\mu_{\boldsymbol p}(B_i).
$$
Average and divide by two. ∎

**Corollary 3.2 (explicit answer to `T.ETH.6`).** Assertion 3 of Theorem 0.1 holds.

*Proof.* For $D\le m-2$, put $j=\lfloor(m+D)/2\rfloor+1$. The binomial mass bound
$2^{-m}\binom mj\ge e^{-m\mathrm{KL}(j/m\|1/2)}/(m+1)$, together with
$\mathrm{KL}(a\|1/2)\le4(a-1/2)^2$, gives
$$
\tau_D(m)\ge\frac{e^{-(D+2)^2/m}}{m+1}.
$$
For $M\sim\mathrm{Bin}(Q,2p)$, $\Pr(M\ge Q/2)\ge1/2$ when $p\ge1/4$. On that event, with $D=R-1$,
$\tau_D(M)\ge e^{-2(R+1)^2/Q}/(Q+1)$. Use (15.1) and
$p\sum_i\mu_p(B_i)\le pQ\max_i\mu_p(B_i)$. ∎

For completeness, if $Q\ge4$ and $Q/4<R\le Q$, keeping just the event $M=Q$ in (15.1) yields
$\max_i\mu_p(B_i)\ge p^{Q-1}/Q\ge4^{-(Q-1)}/Q$. Since $R^2/Q>Q/16$, this also has the requested $e^{-O(R^2/Q)}/\operatorname{poly}(Q)$ form. The small dimensions are covered by (15.1) itself.

The polynomial prefactor in (15.2) is harmless for the abstract question but too wasteful for the desired algorithmic floor. The next two estimates retain the mass of a whole binomial tail, rather than a single level.

**Lemma 3.3 (a binomial tail without a polynomial prefactor).** For $m\ge256$ and $0\le D\le m/8$,
$$
\tau_D(m)\ge\frac34\exp(-25-2D^2/m).
\tag{15.7}
$$

*Proof.* Set $\theta=1/2+s$, where
$s=D/(2m)+m^{-1/2}+m^{-1}$, and let $Y\sim\mathrm{Bin}(m,\theta)$. We have $0\le s\le1/4$. Chebyshev's inequality gives
$\Pr_\theta(|Y-m\theta|<\sqrt m)\ge3/4$. Every integer in this interval exceeds $(m+D)/2$. On the interval the ratio of the fair mass to the tilted mass is at least
$$
\exp\left[-m\mathrm{KL}(\theta\|1/2)-\sqrt m\log\frac\theta{1-\theta}\right]
\ge\exp[-4ms^2-6\sqrt m\,s].
$$
Here $\mathrm{KL}(1/2+s\|1/2)\le4s^2$, and
$\log[(1+2s)/(1-2s)]\le6s$ for $0\le s\le1/4$, by differentiation. Write $a=D/\sqrt m$. Since $m\ge256$, $\sqrt m\,s\le a/2+17/16$, so
$$
4ms^2+6\sqrt m\,s
\le a^2+\frac{29}4a+\frac{697}{64}
\le2a^2+25.
$$
Summing the likelihood-ratio bound over the interval proves (15.7). ∎

**Theorem 3.4 (heavy defects).** Suppose $Q\ge512$, $0\le D\le Q/16$, and $1/4\le p\le1/2$. For $H=\{x:|x|\ge Q/32\}$,
$$
\sum_i\mu_p(B_i\cap H)\ge\frac38e^{-25-4D^2/Q},
\qquad
\max_i\mu_p(B_i\cap H)\ge\frac{3}{8Q}e^{-25-4D^2/Q}.
\tag{15.8}
$$

*Proof.* Use the coupling of Theorem 3.1 and restrict to $M\ge Q/2$. On a restricted fair cube,
$$
\mathbb E_Y\sum_{i:Z_i=1}\mathbf1_{B_i\cap H}(X)
\ge2\tau_D(M)-M\Pr(\mathrm{Bin}(M,1/2)<Q/32).
$$
For $M\ge Q/2$, the last threshold is at most $M/16$. The multiplicative Chernoff bound gives a probability at most
$e^{-49M/256}\le e^{-Q/12}$. Lemma 3.3 and $\Pr(M\ge Q/2)\ge1/2$ therefore give the unconditional lower bound
$$
\frac34e^{-25-4D^2/Q}-Qe^{-Q/12}.
$$
Pointwise $\sum_{i:Z_i=1}\mathbf1_{B_i\cap H}(X)\le\sum_i\mathbf1_{B_i\cap H}(X)$. Also, for $Q\ge512$,
$$
Qe^{-Q/12}\le\frac38e^{-25-Q/64}
\le\frac38e^{-25-4D^2/Q}.
$$
The first inequality is equivalent to
$\log Q+25+\log(8/3)\le13Q/192$; it holds at 512 and its right side grows faster thereafter. This proves the sum bound. Divide by $Q$ for the maximum. ∎

## 4. Application to every group-algebra swap design

Let $\Pi_1,\ldots,\Pi_Q$ be the commuting coordinate swaps and let
$$
P(T)=\sum_{U\subseteq[Q]}c_U(T)\Pi_U,
\qquad P(0)=I,\qquad\deg c_U\le D.
$$
Its eigenvalue on the joint eigenspace indexed by the antisymmetric set $x$ is
$n_x(T)=\sum_Uc_U(T)\chi_U(x)$. Thus $n_x(0)=1$.

**Corollary 4.1 (defect-rank bound, without grading).** At node $t_i$, let
$N_i=\sum_{U\not\ni i}c_U(t_i)\Pi_U$ be the nonlocal defect. Then
$$
N_i|_x=\nu_i(x)I,
\qquad
\nu_i(x)=\tfrac12(n_x(t_i)+n_{x\oplus e_i}(t_i)),
$$
so $B_i=\{\nu_i\ne0\}$. In the label representation
$\mathcal B=([r]\times[r])^Q$, of dimension $L_{\mathcal B}=r^{2Q}$,
$$
\operatorname{rank}N_i=L_{\mathcal B}\mu_{p_r}(B_i),
\qquad p_r=(r-1)/(2r).
$$
Theorems 0.1 and 3.4 therefore give lower bounds on the defect ranks and on their high-antisymmetry part for **all** such designs.

*Proof.* Averaging the two eigenvalues cancels exactly the characters with $i\in U$. The multiplicity formula is the tensor-product swap spectrum, [Claude 08, Lemma 5.1](../../claude/eth/08_GAPLESS_PROOF.md) (`ETH.I.SWAP`). Apply the previous theorems. ∎

No implication about arbitrary faster SAT algorithms follows. Even within deformation counting, this statement addresses group-algebra designs, not deformations outside that algebra. A defect-rank statement also does not itself bound measurements after filling; [16](16_MEASUREMENT_TRANSFER_AND_FLOOR.md) supplies that missing interface under explicit assumptions.

## 5. Extension: unequal orders of edge agreement

The proof works with jets, which connects it to Hermite interpolation without identifying the two constructions.

**Definition 5.1 (jet splines).** Choose positive integers $w_i$, and require on an edge in direction $i$ that
$h_x-h_{x\oplus e_i}$ be divisible by $(T-t_i)^{w_i}$. A relaxed system instead imposes order $a_e\in\{0,\ldots,w_i\}$ on edge $e$. Its **number of erased scalar conditions** is
$b=\sum_e(w_{i(e)}-a_e)$. Put $W=\sum_iw_i$ and $w(S)=\sum_{i\in S}w_i$.

**Theorem 5.2 (weighted jet obstruction).** If a relaxed system has a spline $f$ of degree at most $D$ with parity trace at zero, then
$$
\boxed{\quad b\ge\#\{S\subseteq[Q]:w(S)>\lfloor(W+D)/2\rfloor\}.\quad}
\tag{15.9}
$$

*Proof.* The full system has basis
$g_S=\chi_S\prod_{i\in S}(1-T/t_i)^{w_i}$, with degree $w(S)$ and trace $\chi_S$. This follows by the same Fourier divisibility argument as Lemma 1.2, now using derivatives through order $w_i-1$. Restoring erased jets imposes $b$ scalar linear conditions, so Lemma 1.3 holds with this $b$. The relaxed conditions are still preserved by multiplication, because each is divisibility by an ideal. If $d=\lfloor(W+D)/2\rfloor$, use $g_S$ when $w(S)\le d$, and $fg_{S^c}$ otherwise; its degree is at most $D+W-w(S)\le d$. The trace is always $\chi_S$. The trace-space dimension deficit is exactly the right side of (15.9). ∎

This counts erased **jets**, not just edges. An edge of order $w_i$ costs $w_i$ conditions to remove. The square-root scale in the equal-weight case comes from the distribution of subset sizes; unequal weights replace it by the corresponding subset-sum distribution.

## 6. Attribution, checks, and scope

Graph splines and polynomial matching across colored edges are established tools:

- Gilbert, Tymoczko, and Viel, [*Generalized splines on arbitrary graphs*](https://msp.org/pjm/2016/281-2/pjm-v281-n2-p03-p.pdf), Pacific Journal of Mathematics 281 (2016): vertex labels differ by elements of prescribed edge ideals.
- Hambardzumyan, Hatami, and Qian, [*Lower bounds for graph bootstrap percolation via properties of polynomials*](https://arxiv.org/abs/1708.04640), JCTA 174 (2020): low-degree vertex polynomials agree at edge colors, yielding dimension obstructions for percolation.
- Morrison and Ogden, [*Using polynomials to find lower bounds for r-bond bootstrap percolation*](https://arxiv.org/abs/2410.24130): further developments of that polynomial method.
- Smolensky's parity-complement dimension argument is reproduced in [Claude 14, Lemma 2.1](../../claude/eth/14_DESIGN_FLOOR.md).

The contribution here is the particular parity-trace completion argument, its deleted-edge and biased-measure consequences, and its application to the repository's ungraded-design question. **It is new relative to the repository. A limited literature check did not establish worldwide priority.** The claims above are proved here rather than inferred from an unsuccessful search.

[`scripts/verify_spline_bound.py`](scripts/verify_spline_bound.py) uses exact rational elimination to check the full-module dimensions, trace dimensions, deletion inequality, character completion, biased coupling, genuinely ungraded feasible designs, the explicit one-edge construction, and unequal jet orders. It also checks the numerical tail bound on a finite grid. Negative controls remove distinctness, the parity trace, or the degree restriction. Results: [`results/spline_bound.json`](results/spline_bound.json). The proof, rather than those finite checks, establishes the arbitrary-dimensional assertions.

The remaining sharp question is the factor $Q$ in (15.6): can it be removed for the largest directional defect, or can ungraded designs attain a genuine polynomial improvement? The exponential scale requested in `T.ETH.6` is settled by this note; exact constants and extremizers are not.

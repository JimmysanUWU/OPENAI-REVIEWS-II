# 17. Trace capacity and approximate locality

GPT / Astra, 1 November 2026. **Complete prose proofs with exact finite checks; not Lean-formalized.** The results below extend the mechanism of [15](15_UNGRADED_SPLINE_BOUND.md). They address shared constraints, reduction modulo a prime, and the size of real defects. None asserts an improved SAT algorithm. The orthogonal-polynomial calculation is a classical ingredient, derived here rather than claimed as a new special-function identity.

## 0. Three questions with different answers

The spline theorem counts nonzero failures of locality. It leaves three natural questions.

1. **Can several low-degree constructions share the same exceptional edges?** Their traces multiply. A joint bound sees characters produced by these products, beyond the bounds obtained from the constructions separately.
2. **Can the defects all become small modulo an odd prime?** Not while the interpolation nodes remain distinct and nonzero modulo that prime: the same support obstruction survives in the residue field.
3. **Can every real defect be tiny?** Under the uniform cube measure there is a sharp quadratic bound. Its exact value is a discrete Christoffel function, with a transition at degree of order $\sqrt Q$.

These are different notions of capacity. The first two concern rank and exact vanishing. The third uses a norm and the spacing of the nodes. They must not be substituted for one another in a complexity argument.

## 1. A general trace-multiplication bound

Use the cube $C_Q=\{0,1\}^Q$, distinct nonzero nodes $t_i$ in a field $K$ of characteristic different from two, and characters $\chi_A(x)=(-1)^{\sum_{i\in A}x_i}$. Let

$$
H_d=\operatorname{span}_K\{\chi_A:|A|\le d\},
$$

with $H_d=0$ for $d<0$. For a spanning subgraph $G$, a degree-$d$ spline has vertex polynomials of degree at most $d$ agreeing at $t_i$ along every retained edge in direction $i$. Denote its space of traces at $T=0$ by $V_d(G)$. Let $b$ be the number of deleted edges.

The two facts proved in [15, §1](15_UNGRADED_SPLINE_BOUND.md#1-the-filtered-spline-algebra) are

$$
V_d(C_Q)=H_d,\qquad \dim V_d(G)-\dim H_d\le b.
\tag{17.1}
$$

**Theorem 1.1 (multiplication-rank certificate).** If a degree-$D$ spline on $G$ has trace $h$, then, for every integer $d\ge D$,

$$
\boxed{\quad b\ge
\operatorname{rank}\big(\pi_d M_h\big|_{H_{d-D}}\big).
\quad}\tag{17.2}
$$

Here $M_h$ is pointwise multiplication by $h$ and $\pi_d:K^{C_Q}\to K^{C_Q}/H_d$ is the quotient map. Equivalently, the right side is $\dim(H_d+hH_{d-D})-\dim H_d$.

*Proof.* Lift every element of $H_{d-D}$ to a full-cube spline of degree at most $d-D$. Multiplying that lift by the given spline produces a spline on $G$ of degree at most $d$. Its trace is the product of the two traces. Therefore $H_d+hH_{d-D}\subseteq V_d(G)$, and (17.1) applies. ∎

For $h=\chi_{[Q]}$, multiplication replaces a character by its complement. Choosing $d=\lfloor(Q+D)/2\rfloor$ for $D\le Q$ recovers the binomial-tail bound of [15]. The statement also makes sense when $h$ is not a character: compute a finite multiplication matrix and take its rank after discarding the low-degree rows.

**Corollary 1.2 (low-degree trace corrections).** Suppose $D\le Q$ and a degree-$D$ spline on $G$ has trace $\chi_{[Q]}+\ell$, where $\ell\in H_D$. Then

$$
b\ge\sum_{j>(Q+D)/2}\binom Qj.
$$

*Proof.* Since $\ell H_{d-D}\subseteq H_d$, the correction vanishes under $\pi_d$. Alternatively, subtract a full-cube degree-$D$ spline with trace $\ell$ and apply [15]. ∎

The hypothesis on the correction is exact; this is not a stability theorem for arbitrary small perturbations of the trace.

## 2. Several traces sharing one exceptional set

**Theorem 2.1 (joint character completion).** Suppose the same graph $G$ supports splines $f_1,\ldots,f_m$ of respective degrees at most $D_1,\ldots,D_m$, with traces $\chi_{A_1},\ldots,\chi_{A_m}$. For $J\subseteq[m]$, write

$$
A_J=\mathop{\triangle}_{j\in J} A_j,\qquad c_J=\sum_{j\in J}D_j,
$$

where $\triangle$ is symmetric difference. For every integer $d\ge0$, define a union of Hamming balls in the set of subsets of $[Q]$:

$$
U_d=\bigcup_{\substack{J\subseteq[m]\\c_J\le d}}
\{A_J\triangle S: |S|\le d-c_J\}.
$$

Then

$$
\boxed{\quad b\ge |U_d|-\sum_{k=0}^{\min(d,Q)}\binom Qk.\quad}
\tag{17.3}
$$

*Proof.* The product $\prod_{j\in J}f_j$ has degree at most $c_J$ and trace $\chi_{A_J}$. Multiply it by the full-cube character spline for $S$. When $|S|+c_J\le d$, this proves $\chi_{A_J\triangle S}\in V_d(G)$. Distinct characters are independent over $K$. The empty $J$ already contributes $H_d$. Apply (17.1). ∎

Repeated use of a generator is unnecessary: its trace squares to one, so an even repetition adds degree without changing the trace. If several $J$ produce the same $A_J$, keep the least cost. This expresses the theorem as a weighted covering problem in a binary linear code. The theorem proves a lower bound for a shared exceptional set, not the existence of an optimal graph attaining it.

**Example 2.2 (a product supplies a genuinely additional constraint).** Take $Q=4$, $A_1=\{1,2\}$, $A_2=\{3,4\}$, and $D_1=D_2=1$. At $d=2$, $H_2$ contains the eleven characters of weight at most two. Multiplying $f_1$ by the singleton splines in directions 3 and 4 produces two triple characters. Multiplying $f_2$ by the singleton splines in directions 1 and 2 produces the other two. The product $f_1f_2$ supplies the weight-four character. Thus $U_2$ is the whole set of sixteen characters and

$$
b\ge16-11=5.
$$

The two separate multiplication-rank certificates at $d=2$ each give only two; their joint span without the product gives four. The fifth dimension comes from the product. This is a precise mechanism for detecting an interaction between constructions which share a failure set.

For more general traces $h_j$, the same proof uses the span of $H_d$ and the spaces $(\prod_{j\in J}h_j)H_{d-c_J}$; unlike characters, arbitrary traces can benefit from repeated powers. A useful next question is which actual algorithmic constructions force several such traces on a common graph. That application has not been established here.

## 3. What survives reduction modulo a prime

Let $\mathcal O$ be a discrete valuation ring, with maximal ideal $\mathfrak m$, residue field $k$, and residue characteristic different from two. Suppose $t_i\in\mathcal O$ have pairwise distinct nonzero residues, and $n_x\in\mathcal O[T]$ have degree at most $D$ and $n_x(0)=1$. Define

$$
\nu_i(x)=\frac{n_x(t_i)+n_{x\oplus e_i}(t_i)}2,
\qquad E_i=\{x:\nu_i(x)\notin\mathfrak m\}.
$$

Thus $E_i$ records **unit defects**, a stronger condition than merely being nonzero in the fraction field.

**Theorem 3.1 (residual obstruction).** There are at least $\sum_{j>(Q+D)/2}\binom Qj$ unoriented edges whose defect is a unit. Moreover, for real probabilities $0<p_i\le1/2$,

$$
\sum_i p_i\,\mu_{\boldsymbol p}(E_i)
\ge \mathbb E\left[2^{-M}\sum_{j>(M+D)/2}\binom Mj\right],
\qquad M=\sum_i\operatorname{Bernoulli}(2p_i).
\tag{17.4}
$$

*Proof.* Reduce all polynomials and nodes modulo $\mathfrak m$. The hypotheses of the characteristic-not-two version of the spline theorem remain true. An edge fails its reduced equation exactly when its original defect is a unit. Apply the algebraic edge count, then the same random restriction as in [15, §3]. Probabilities count finite sets and remain ordinary real numbers; they are not elements of the residue field. ∎

In particular, if $D<Q$, it is impossible for all defects to lie in $\mathfrak m$, much less in a higher power of $\mathfrak m$. The heavy-defect conclusion of [15, Theorem 3.4] likewise applies to $E_i$ under its numerical hypotheses, since its proof only restricts coordinates and counts edges.

For integer nodes $1,\ldots,Q$ reduced modulo an ordinary prime $p$, the separation hypothesis requires $p>Q$. It would be wrong to apply (17.4) to these nodes at a fixed small odd prime as $Q$ grows. One may instead use other nodes over a larger residue field, but that changes the arithmetic setting.

**Example 3.2 (colliding residues).** Take nodes $t_i=1+ip$ in $\mathbb Z_{(p)}$ for an odd prime $p$, and set every $n_x(T)=1-T$. These distinct characteristic-zero nodes all reduce to one, and every $\nu_i=-ip$ lies in the maximal ideal. Degree one now suffices. This demonstrates the need for residue separation; it does not give locality over characteristic zero.

**Example 3.3 (characteristic two).** Over $\mathbb F_2$, $n_x=1$ satisfies $n_x+n_{x\oplus e_i}=0$ on every edge. Parity has become constant, division by two is unavailable, and the character decomposition used above has collapsed. This is a failure of this proof mechanism, not a proof that the original matrix deformation is local.

In fact, for the two-point swap $\Pi$, the matrix $N=\Pi-I$ over $\mathbb F_2$ satisfies $N^2=0$ but $N(1,0)^\mathsf T=(1,1)^\mathsf T\ne0$. Nilpotence of the swap algebra's augmentation ideal does not make arbitrary scalar readouts vanish. A proposed characteristic-two algorithm must still control those readouts and the counting modulus.

## 4. A sharp bound for the size of real defects

Return to real coefficients, use nodes $t_i=i$, and assume $0\le D<Q$. For a family $n_x$ of degree at most $D$ with $n_x(0)=1$, define $\nu_i$ as in §3 and its total quadratic defect

$$
\mathcal E(n)=\sum_{i=1}^Q 2^{-Q}\sum_{x\in C_Q}\nu_i(x)^2.
\tag{17.5}
$$

This section uses the **uniform cube measure**. The biased spectral measure of the counting algorithm is different.

**Theorem 4.1 (exact defect capacity).** Put

$$
\rho_{Q,D}=\prod_{j=0}^{D}\frac{Q-j}{Q+j+1}
=\frac{(Q!)^2}{(Q-D-1)!(Q+D+1)!}.
$$

Then

$$
\boxed{\quad\min_n\mathcal E(n)=\frac{\rho_{Q,D}}{1-\rho_{Q,D}}.\quad}
\tag{17.6}
$$

The minimum is attained by a family independent of $x$, $n_x(T)=q_*(T)$, where $q_*$ is specified in §5. Therefore if all $|\nu_i(x)|\le\varepsilon$, necessarily

$$
\varepsilon\ge\sqrt{\frac{\rho_{Q,D}}{Q(1-\rho_{Q,D})}}.
\tag{17.7}
$$

The energy optimum is sharp; no claim is made that the consequence (17.7) solves the maximum-norm minimization problem.

*Reduction to one variable.* Expand

$$
n_x(T)=\sum_{S\subseteq[Q]}a_S(T)\chi_S(x).
$$

The normalization gives $a_\varnothing(0)=1$ and $a_S(0)=0$ for nonempty $S$. Averaging across direction $i$ retains exactly those characters with $i\notin S$. Parseval's identity yields

$$
\mathcal E(n)=\sum_{S\subseteq[Q]}\sum_{i\notin S}a_S(i)^2
\ge\sum_{i=1}^Qa_\varnothing(i)^2.
\tag{17.8}
$$

Conversely, setting every nonconstant Fourier coefficient to zero attains the right side for any polynomial $q=a_\varnothing$ with $q(0)=1$. Hence

$$
\min_n\mathcal E(n)
=\min_{\deg q\le D,\ q(0)=1}\sum_{i=1}^Qq(i)^2.
\tag{17.9}
$$

The remaining calculation is given in full below. There may be other minimizing vertex families because some nonconstant character coefficients have zero energy; the scalar minimizer $q_*$ itself is unique for $D<Q$. ∎

The same Fourier argument works for arbitrary distinct nonzero real nodes and positive direction weights $w_i$. The optimum is then $1/(e_0^\mathsf T A^{-1}e_0)$, where

$$
A_{rs}=\sum_i w_i t_i^{r+s},\quad 0\le r,s\le D,
\qquad e_0=(1,0,\ldots,0)^\mathsf T.
\tag{17.10}
$$

Indeed $A$ is positive definite for $Q>D$, and Cauchy–Schwarz in the $A$-inner product gives the constrained minimum. Node spacing therefore matters for energy although it did not matter for the exact edge count.

## 5. The discrete orthogonal-polynomial calculation

For clarity, $x^{\underline j}=x(x-1)\cdots(x-j+1)$ is a falling factorial, and $\Delta p(x)=p(x+1)-p(x)$. Define, for $0\le j\le Q$,

$$
R_j(x)=x^{\underline j}(x-Q-1)^{\underline j},\qquad
P_j(x)=\Delta^jR_j(x),\qquad
H_j(x)=P_j(x)/P_j(0).
$$

**Lemma 5.1 (discrete Rodrigues calculation).** On the points $0,1,\ldots,Q$ with weight one, $H_j$ is orthogonal to all polynomials of degree less than $j$, $H_j(0)=1$, and

$$
h_j:=\sum_{x=0}^QH_j(x)^2
=\frac{(Q+j+1)!(Q-j)!}{(2j+1)(Q!)^2}.
\tag{17.11}
$$

*Proof.* Expand the finite difference and substitute $y=x+a$ in its $a$-th term. $R_j$ vanishes at $0,\ldots,j-1$ and at $Q+1,\ldots,Q+j$, so only $j\le y\le Q$ survives. With $\nabla p(y)=p(y)-p(y-1)$ this gives

$$
\sum_{x=0}^QP_j(x)p(x)
=(-1)^j\sum_{y=j}^QR_j(y)\nabla^jp(y).
\tag{17.12}
$$

It is zero if $\deg p<j$. Since $R_j$ is monic of degree $2j$, the leading coefficient of $P_j$ is $(2j)!/j!$, and $\nabla^jP_j=(2j)!$. Also

$$
P_j(0)=R_j(j)=(-1)^j j!\,Q^{\underline j},\qquad
(-1)^jR_j(y)=(j!)^2\binom yj\binom{Q+j-y}j.
$$

Apply (17.12) with $p=P_j$. The convolution

$$
\sum_{y=j}^Q\binom yj\binom{Q+j-y}j
=\binom{Q+j+1}{2j+1}
$$

counts subsets of size $2j+1$ by their middle element. It follows that

$$
\sum_{x=0}^QP_j(x)^2
=(2j)!(j!)^2\binom{Q+j+1}{2j+1}.
$$

Divide by $P_j(0)^2$ to obtain (17.11). ∎

**Lemma 5.2 (endpoint kernel).** Let

$$
K_D(x,0)=\sum_{j=0}^D\frac{H_j(x)}{h_j}.
$$

Then $K_D(0,0)=1-\rho_{Q,D}$, and the minimizer in (17.9) is

$$
q_*(x)=\frac{K_D(x,0)}{K_D(0,0)}.
\tag{17.13}
$$

*Proof.* Every polynomial $q$ of degree at most $D$ satisfies

$$
q(0)=\sum_{x=0}^Qq(x)K_D(x,0).
$$

Thus Cauchy–Schwarz shows that the minimum of $\sum_{x=0}^Qq(x)^2$ with $q(0)=1$ is $1/K_D(0,0)$, attained precisely by (17.13). Set $\rho_{Q,-1}=1$. Direct subtraction gives

$$
\rho_{Q,j-1}-\rho_{Q,j}
=\frac{(2j+1)(Q!)^2}{(Q-j)!(Q+j+1)!}
=\frac1{h_j}.
$$

The kernel diagonal telescopes to $1-\rho_{Q,D}$. Finally subtract the compulsory value $q(0)^2=1$ from the minimum on $0,\ldots,Q$, giving $\rho/(1-\rho)$ on $1,\ldots,Q$. ∎

These $H_j$ are the Hahn polynomials $Q_j(x;0,0,Q)$, also called discrete Chebyshev or Gram polynomials. The specialization of the orthogonality and norm table in [NIST DLMF §18.19](https://dlmf.nist.gov/18.19#T1) agrees with (17.11). The proof above makes the normalization and endpoint subtraction explicit. The special-function formula and the reproducing-kernel minimization are classical; their use here separates a precise amplitude question from the earlier support question.

## 6. The square-root transition and its limitations

**Corollary 6.1 (quantitative scales).** For $0\le D<Q$,

$$
\exp\left(-\frac{(D+1)^2}{Q-D}\right)
\le\rho_{Q,D}\le
\exp\left(-\frac{(D+1)^2}{Q+D+1}\right).
\tag{17.14}
$$

Consequently, if $Q\to\infty$ and $(D+1)^2/Q\to\lambda\in(0,\infty)$, then

$$
\min_n\mathcal E(n)\longrightarrow\frac1{e^\lambda-1}.
\tag{17.15}
$$

If $(D+1)^2/Q\to0$, the minimum is asymptotic to $Q/(D+1)^2$.

*Proof.* For the upper bound in (17.14), write each factor as $1-(2j+1)/(Q+j+1)$ and use $\log(1-u)\le-u$. For the lower bound use $\log(1-u)\ge-u/(1-u)$ and $\log(1+v)\le v$ on $(Q-j)/(Q+j+1)=(1-j/Q)/(1+(j+1)/Q)$. Summing gives the displayed bounds because $\sum_{j=0}^D(2j+1)=(D+1)^2$. In both limiting regimes $D/Q\to0$, so these bounds squeeze $-\log\rho$ to a quantity asymptotic to $(D+1)^2/Q$. Apply (17.6). ∎

Checks at the extremes are useful:

$$
\min\mathcal E=Q\quad(D=0),\qquad
\min\mathcal E=\frac{Q(Q-1)}{2(2Q+1)}\quad(D=1),\qquad
\min\mathcal E=\frac1{\binom{2Q}Q-1}\quad(D=Q-1).
$$

At $D=Q$, the polynomial $\prod_{i=1}^Q(1-T/i)$ gives zero energy; that case lies outside (17.6).

The support and amplitude obstructions share a degree-squared-over-dimension scale, but have different proofs and hypotheses. Small energy can coexist with many nonzero defects. Large energy can be concentrated on few defects. The actual measurement distribution is biased, and the original counting argument is arithmetic rather than an approximate real computation. A transfer theorem for rounding, error tolerance, or image compression would therefore be additional mathematics, not a consequence silently supplied by (17.6).

## 7. Verification and open interfaces

Run `python reviews/gpt/eth/scripts/verify_defect_capacity.py`. It checks the multi-trace character count and spline feasibility, finite-field trace ranks and feasible reduced designs, the residue-collision and characteristic-two controls, all Rodrigues inner products and norms in a finite range, and the energy optimum against an independent rational Gram-matrix solve. It also checks Parseval's energy identity on nonsymmetric polynomial families. Results are written to [defect_capacity.json](results/defect_capacity.json).

The worthwhile open interfaces are specific: identify an algorithmic source of several traces with a shared defect set; analyze the exact energy optimum under biased product measures; and determine whether arithmetic structure permits useful compression of the measurement image. None is claimed solved in this note. The [research ledger](../RESEARCH_LEDGER.md) records explored routes and the reasons several shortcuts fail.

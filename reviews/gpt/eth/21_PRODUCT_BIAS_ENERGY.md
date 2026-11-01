# 21. Product bias is spectral thinning: the exact defect-energy optimum

GPT, 1 November 2026. **Complete independent prose derivation with exact rational checks; not Lean-formalized.** Answers target **T.ETH.8**. The ingredients are tensor diagonalization of squared readouts, scalar least squares, and Cauchy–Binet. No worldwide novelty is asserted.

## 0. Definitions

Fix distinct nonzero real nodes $t_1,\ldots,t_Q$, biases $0<p_i<1$, positive direction weights $w_i$, and $0\le D<Q$. For the product law $\mu_{\boldsymbol p}$ on $\{0,1\}^Q$, define
$$
E_{\boldsymbol p,w}=\min_{\substack{\deg n_x\le D\\n_x(0)=1}}
\sum_{i=1}^Qw_i\,
\mathbb E_{\mu_{\boldsymbol p}}\left[
\left(\frac{n_x(t_i)+n_{x\oplus e_i}(t_i)}2\right)^2\right].
\tag{21.1}
$$
The polynomials have real coefficients and no restriction on their dependence on $x$. For $S\subseteq[Q]$, put
$$
C_D(S;v)=\min_{\deg q\le D,\ q(0)=1}\sum_{i\in S}v_iq(t_i)^2.
\tag{21.2}
$$
Set $\theta_i=4p_i(1-p_i)$. Let $\mathsf S$ include $i$ independently with probability $\theta_i$.

## 1. Exact reduction and optimizer

**Theorem 21.1 (spectral thinning).**
$$
\boxed{E_{\boldsymbol p,w}=\mathbb E_{\mathsf S}C_D(\mathsf S;w/\theta).}
\tag{21.3}
$$
For each $S$, choose a scalar minimizer $q_S$ in the right side. Define the one-coordinate functions
$$
a_i(0)=2p_i,\quad a_i(1)=2(1-p_i),\qquad b_i=1-a_i.
$$
An explicit optimal family is
$$
\boxed{n_x(T)=\sum_{S\subseteq[Q]}q_S(T)
\prod_{i\in S}a_i(x_i)\prod_{i\notin S}b_i(x_i).}
\tag{21.4}
$$

*Proof.* Work in $H=L^2(\mu_{\boldsymbol p})$, and let $P_ih=(h+h\circ\mathrm{flip}_i)/2$. The swaps are not self-adjoint under bias. Their squared readouts $K_i=P_i^*P_i$, however, act on separate coordinates and commute. In one coordinate,
$$
K=\frac14\begin{pmatrix}
(1-p)^{-1}&(1-p)^{-1}\\ p^{-1}&p^{-1}
\end{pmatrix}
$$
in value coordinates with inner product $\operatorname{diag}(1-p,p)$. This positive rank-one operator has eigenvalues $1/\theta,0$. Thus $A=\theta K$ and $B=I-A$ are orthogonal rank-one projections. Directly,
$$
A1=a,\quad B1=b,\quad a+b=1,\quad
\langle a,b\rangle_p=0,\quad
\|a\|_p^2=\theta,\quad\|b\|_p^2=1-\theta .
\tag{21.5}
$$
These formulas include $p=1/2$, where $b=0$, although the projection $B$ is still nonzero.

The tensor projections
$$
\Pi_S=\prod_{i\in S}A_i\prod_{i\notin S}B_i
$$
are orthogonal, rank one, sum to $I$, and satisfy $K_i\Pi_S=\theta_i^{-1}\Pi_S$ for $i\in S$, zero otherwise. Decompose the $H$-valued polynomial $n(T)$ into $n_S=\Pi_Sn$. Its energy is
$$
\sum_S\sum_{i\in S}(w_i/\theta_i)\|n_S(t_i)\|_H^2,
\qquad n_S(0)=\Pi_S1 .
$$
If $\Pi_S1\ne0$, rank one forces $n_S(T)=q_S(T)\Pi_S1$, $q_S(0)=1$; its optimal cost is $\|\Pi_S1\|^2 C_D(S;w/\theta)$. If $\Pi_S1=0$, the zero polynomial is optimal. Finally,
$$
\|\Pi_S1\|^2=\prod_{i\in S}\theta_i\prod_{i\notin S}(1-\theta_i)
=\Pr(\mathsf S=S).
$$
Minimize the sectors separately. Expanding $\Pi_S1$ gives (21.4), whose normalization is $\prod_i(a_i+b_i)=1$. ∎

**Explicit scalar minimizers.** If $|S|\le D$, take $q_S(T)=\prod_{i\in S}(1-T/t_i)$, giving capacity zero. Otherwise, with
$$
u(t)=(1,t,\ldots,t^D)^{\mathsf T},\quad
M_S(v)=\sum_{i\in S}v_i u(t_i)u(t_i)^{\mathsf T},\quad
e_0=(1,0,\ldots,0)^{\mathsf T},
$$
distinct nodes make $M_S$ positive definite, and
$$
C_D(S;v)=\frac1{e_0^{\mathsf T}M_S^{-1}e_0},\qquad
q_S(T)=\frac{u(T)^{\mathsf T}M_S^{-1}e_0}
{e_0^{\mathsf T}M_S^{-1}e_0}.
\tag{21.6}
$$
Indeed a coefficient vector $c$ has $e_0^{\mathsf T}c=1$ and cost $c^{\mathsf T}M_Sc$; Cauchy–Schwarz gives the formula and equality case. This specifies (21.4) by finite matrix inverses. It has $2^Q$ sectors; no efficient evaluation algorithm is claimed.

## 2. The degree controls sensitivity to bias

Write
$$
\Delta_{D+1}(S;v)=\sum_{\substack{J\subseteq S\\|J|=D+1}}
\prod_{i\in J}v_i\prod_{\substack{i<j\\i,j\in J}}(t_j-t_i)^2,
$$
$$
\Gamma_D(S;v)=\sum_{\substack{J\subseteq S\\|J|=D}}
\prod_{i\in J}(v_it_i^2)\prod_{\substack{i<j\\i,j\in J}}(t_j-t_i)^2 .
\tag{21.7}
$$
Empty products give $\Gamma_0=1$.

**Lemma 21.2 (Cauchy–Binet capacity).** For $|S|>D$,
$$
C_D(S;v)=\Delta_{D+1}(S;v)/\Gamma_D(S;v).
\tag{21.8}
$$
For every $S$, this is at least $\Delta_{D+1}(S;v)/\Gamma_D([Q];v)$.

*Proof.* Cauchy–Binet gives $\det M_S=\Delta_{D+1}(S;v)$. Its minor deleting row and column zero is $\Gamma_D(S;v)$: extract a factor $t_i$ from each Vandermonde column on powers $1,\ldots,D$. Now apply the inverse cofactor formula in (21.6). All denominator summands are nonnegative, giving the inequality. When $|S|\le D$, both its numerator and left side are zero. ∎

**Theorem 21.3 (sharp bias comparison).**
$$
\boxed{\frac{\Delta_{D+1}([Q];w)}
{\Gamma_D([Q];w/\theta)}
\le E_{\boldsymbol p,w}\le C_D([Q];w).}
\tag{21.9}
$$
In particular, with $\theta_{\min}=\min_i\theta_i$,
$$
\boxed{\theta_{\min}^{D}C_D([Q];w)
\le E_{\boldsymbol p,w}\le C_D([Q];w).}
\tag{21.10}
$$
The lower bound in (21.9) is exact at $D=Q-1$.

*Proof.* Apply Lemma 21.2 to the random set in (21.3). Every size-$(D+1)$ monomial satisfies
$$
\mathbb E\left[\mathbf1_{J\subseteq\mathsf S}
\prod_{i\in J}(w_i/\theta_i)\right]=\prod_{i\in J}w_i.
$$
Thus its expected numerator is $\Delta_{D+1}([Q];w)$, proving the lower bound. For the upper bound choose the scalar optimal family $n_x=q$ at every vertex.

Every size-$D$ monomial of $\Gamma_D([Q];w/\theta)$ is at most $\theta_{\min}^{-D}$ times the corresponding monomial of $\Gamma_D([Q];w)$. Equation (21.8) proves (21.10). If $D=Q-1$, only $\mathsf S=[Q]$ can have positive capacity; on that event the denominator bound is equality. ∎

The product of the $D$ smallest $\theta_i$ improves $\theta_{\min}^D$, by the same monomial comparison. The degree, rather than the number of nodes, controls possible relative gain.

**Corollary 21.4 (equally spaced nodes).** For $t_i=i,w_i=1$ and common bias $p$, put $\theta=4p(1-p)$. With the uniform optimum from [17, Theorem 4.1](17_TRACE_CAPACITY_AND_APPROXIMATE_LOCALITY.md),
$$
U_{Q,D}=\frac{\rho}{1-\rho},\qquad
\rho=\prod_{j=0}^{D}\frac{Q-j}{Q+j+1},
$$
one has $\theta^D U_{Q,D}\le E_p\le U_{Q,D}$. At $D=Q-1$,
$$
E_p=\frac{\theta^{Q-1}}{\binom{2Q}{Q}-1}.
\tag{21.11}
$$
For $p_r=(r-1)/(2r)$, $r\ge2$,
$$
\boxed{0\le 1-E_{p_r}/U_{Q,D}\le D/r^2.}
\tag{21.12}
$$
Hence $E_{p_r}=(1+o(1))U_{Q,D}$ if $D=o(r^2)$. Part I's $r=\lfloor\sqrt\lambda\rfloor,\ D=2\lceil\sqrt\lambda\rceil-1$ give relative error $O(\lambda^{-1/2})$.

*Proof.* Apply Theorem 21.3 and the uniform formula. At $D=Q-1$, $\rho=1/\binom{2Q}{Q}$. For the actual bias, $\theta=1-r^{-2}$ and $1-(1-r^{-2})^D\le D/r^2$. ∎

## 3. Monotonicity, example, and scope

**Proposition 21.5 (balance is hardest in each coordinate).** The optimum is nondecreasing and concave in each $\theta_i\in(0,1]$ separately, holding the others fixed. Reflecting $p_i$ to $1-p_i$ leaves it unchanged.

*Proof.* Condition on $\mathsf S\setminus\{i\}=S$. Set $f_S(q)=\sum_{j\in S}(w_j/\theta_j)q(t_j)^2$ and $c_S=\min_{q(0)=1,\deg q\le D}f_S(q)$. The conditional expectation in (21.3) is
$$
(1-\theta_i)c_S+\theta_i\min_q[f_S(q)+(w_i/\theta_i)q(t_i)^2]
=c_S+\min_q[\theta_i(f_S(q)-c_S)+w_iq(t_i)^2].
$$
This is an infimum of affine functions with nonnegative slopes. It is nondecreasing and concave. Average over $S$. Reflection follows from dependence on $4p_i(1-p_i)$ alone. Joint concavity is not claimed. ∎

**Example 21.6 (scalar optimality fails under bias).** At $Q=2,D=1,t=(1,2),p=1/4$, an optimizer has constant coefficients one and slopes, in order $(0,0),(1,0),(0,1),(1,1)$,
$$
-21/40,\quad-43/40,\quad-23/40,\quad-9/40.
$$
Its energy is $3/20$, below the vertex-independent optimum $1/5$. The sharp factor in (21.10) is $\theta=3/4$.

**Scope.** T.ETH.8 is settled as an exact finite optimization problem, with an explicit optimizer and sharp relative asymptotics for the actual bias. Energy is a real sum of squares: it does not control exact defect support, odd-prime unit defects, or compressed measurement-image sizes without further hypotheses. No faster SAT algorithm or lower bound for all counting algorithms follows here.

## 4. Reproduction

~~~text
python reviews/gpt/eth/scripts/verify_biased_energy.py
~~~

The deterministic standard-library script uses exact fractions. It forms the original normal equations in all $D2^Q$ nonconstant vertex coefficients and independently compares their optimum with (21.3) and direct substitution of (21.4). Its 24 QPs include unequal biases, a balanced coordinate, reflected biases, unequal weights, and negative nodes. It also checks the local spectra, determinant identities and bounds, 24 sharp high-degree cases, coordinate monotonicity, concavity and reflection, and 21 exact instances of the actual-swap relative bound. Controls reject scalar optimality, the wrong thinning law, a stronger false bound, and omission of $1/\theta_i$. Results: [biased_energy.json](results/biased_energy.json).

These finite checks support the separate proofs; they do not establish them by enumeration.

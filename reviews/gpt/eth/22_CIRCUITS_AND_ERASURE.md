# 22. Minimal contradictions control deletion and refitting

GPT, 1 November 2026. **Independent prose proofs with exact checks; not Lean-formalized.** This abstracts the energy comparison in ETH 21 and the frustrated-cycle route in round 9. The basic objects are linear equations and their smallest inconsistent subsets. The decomposition into sign-compatible circuits is classical; it is proved below and explicitly credited in §7. No worldwide novelty is asserted.

## 0. Definitions

Let $A\in\mathbb R^{M\times d}$ have rows $a_i^{\mathsf T}$, let $b\in\mathbb R^M$, and fix $w_i>0$. Redundant parameter directions can be removed, so take $\operatorname{rank}A=d$. For $s>1$ and nonnegative weights $v$, put
$$
\mathcal L_s(v)=\min_{\beta\in\mathbb R^d}
\sum_{i=1}^M v_i|b_i-a_i^{\mathsf T}\beta|^s .
\tag{22.1}
$$
The minimum is attained: this is distance to a finite-dimensional linear subspace in the weighted $\ell_s$ space on the positive-weight coordinates. It is nondecreasing in each weight, homogeneous of degree one, and superadditive:
$$
\mathcal L_s(u+v)\ge\mathcal L_s(u)+\mathcal L_s(v).
\tag{22.2}
$$
The last inequality follows by minimizing a sum of two objectives.

A **row circuit** is the support $J$ of a nonzero vector $g\in\ker A^{\mathsf T}$ with inclusion-minimal support. Its dependence is unique up to scale, every $g_i$ on $J$ is nonzero, and $|J|\le d+1$. Call it **inconsistent** if $g^{\mathsf T}b\ne0$. These are exactly the inclusion-minimal inconsistent sets of equations $a_i^{\mathsf T}\beta=b_i$.

Let $\mathsf S$ be a random set of retained equations with marginals $\Pr(i\in\mathsf S)=\theta_i\in(0,1]$. Reweight the retained equations by inverse retention probability:
$$
V_i=w_i\mathbf1_{i\in\mathsf S}/\theta_i.
\tag{22.3}
$$
Thus $\mathbb E V=w$, but generally $\mathbb E\mathcal L_s(V)\ne\mathcal L_s(w)$: the parameters are refitted after seeing $\mathsf S$.

## 1. Exact support is a hitting problem

**Lemma 22.1 (minimal contradictions).** Every inconsistent subsystem contains an inconsistent row circuit. Consequently, deleting a set $D$ leaves a consistent system if and only if $D$ meets every inconsistent circuit.

*Proof.* Choose independent rows spanning the subsystem's row space and fit their equations. If all equations cannot be fitted, some additional row disagrees with the forced value. Its unique dependence on the chosen basis has nonzero pairing with $b$; discarding zero coefficients gives a row circuit. Conversely, the equations on an inconsistent circuit cannot all hold, since their dependence would give $0=g^{\mathsf T}b$. These two observations prove the deletion criterion. ∎

For any nonnegative cost $F(D)$ nondecreasing under inclusion,
$$
\boxed{\min_\beta F(\{i:a_i^{\mathsf T}\beta\ne b_i\})
=\min_{\substack{D:\ D\cap J\ne\varnothing\\
\text{for every inconsistent circuit }J}} F(D).}
\tag{22.4}
$$
Indeed, every actual failure set hits all contradictions; a hitting set permits a fit whose actual failures are a subset of it. This includes weighted support and the largest load among prescribed direction classes. It is an exact reformulation, not an efficient algorithm: the circuit family may be exponentially large.

## 2. The loss splits into minimal contradictions

**Lemma 22.2 (sign-compatible circuit splitting).** There are nonnegative weight vectors $w^{(1)},\ldots,w^{(h)}$, each supported on one inconsistent circuit $J_k$, such that
$$
\sum_k w^{(k)}\le w,\qquad
\boxed{\mathcal L_s(w)=\sum_k\mathcal L_s(w^{(k)})}.
\tag{22.5}
$$
The same full-data minimizer minimizes every summand. If the full system is consistent, take the empty sum.

*Proof.* Let $\beta_*$ minimize (22.1), and write $r=b-A\beta_*$. Its first-order condition is
$$
\lambda_i=w_i\operatorname{sign}(r_i)|r_i|^{s-1},
\qquad A^{\mathsf T}\lambda=0.
\tag{22.6}
$$
We first decompose $\lambda$ without coordinate cancellation. Among nonzero kernel vectors having the signs of $\lambda$, choose one of smallest support. It is a circuit vector: otherwise a kernel vector on a proper subset lets us move to the boundary of that orthant, losing a coordinate without reaching zero. Scale this circuit until subtracting it from $\lambda$ first makes a nonzero coordinate vanish. The remainder still has the signs of $\lambda$. Repeat, at most $M$ times, to obtain
$$
\lambda=\sum_k g^{(k)},\qquad
g_i^{(k)}\lambda_i>0\ \text{whenever }g_i^{(k)}\ne0 .
\tag{22.7}
$$
Each circuit is inconsistent because
$$
(g^{(k)})^{\mathsf T}b=(g^{(k)})^{\mathsf T}r
=\sum_i g_i^{(k)}r_i>0.
$$
Define
$$
w_i^{(k)}=
\begin{cases}
w_i g_i^{(k)}/\lambda_i,&\lambda_i\ne0,\\
0,&\lambda_i=0.
\end{cases}
\tag{22.8}
$$
These weights sum to $w_i$ wherever $r_i\ne0$. The omitted rows have zero residual. At $\beta_*$ the $k$th objective has gradient $-sA^{\mathsf T}g^{(k)}=0$, so convexity makes $\beta_*$ its minimizer. Evaluating all objectives there gives (22.5). ∎

The signs are essential: an arbitrary decomposition of a kernel vector need not produce nonnegative split weights.

## 3. Only the survival of contradictions matters

**Theorem 22.3 (erasure bound from circuit size).** Suppose the full system is inconsistent, and let $c$ be the largest size of an inconsistent circuit. Assume
$$
\Pr(J\subseteq\mathsf S)\ge\prod_{i\in J}\theta_i
\quad\text{for every inconsistent circuit }J .
\tag{22.9}
$$
Writing $\theta_{(1)}\le\cdots\le\theta_{(M)}$, one has, for every $s>1$,
$$
\boxed{\left(\prod_{j=1}^{c-1}\theta_{(j)}\right)\mathcal L_s(w)
\le\mathbb E\mathcal L_s(V)\le\mathcal L_s(w).}
\tag{22.10}
$$
An empty product is one. In particular, common retention $\theta$ gives the factor $\theta^{c-1}$; since $c\le d+1$, it gives $\theta^d$. Condition (22.9) holds for $c$-wise independent Bernoulli indicators; independence beyond that order is unnecessary.

*Proof.* The upper bound evaluates each random objective at $\beta_*$ and uses only $\mathbb E V=w$.

For the lower bound, use Lemma 22.2 and put $V_i^{(k)}=w_i^{(k)}\mathbf1_{i\in\mathsf S}/\theta_i$. Superadditivity and monotonicity give
$$
\mathcal L_s(V)\ge\sum_k\mathcal L_s(V^{(k)}).
$$
Every proper subset of a circuit consists of independent rows, so its equations can be fitted exactly. Therefore
$$
\boxed{\mathbb E\mathcal L_s(V)\ge
\sum_k\Pr(J_k\subseteq\mathsf S)\,
\mathcal L_s(w^{(k)}/\theta).}
\tag{22.11}
$$
This inequality holds for any deletion law. Under (22.9), each summand is at least
$$
\frac{\prod_{i\in J_k}\theta_i}{\max_{i\in J_k}\theta_i}
\mathcal L_s(w^{(k)})
\ge\left(\prod_{j=1}^{c-1}\theta_{(j)}\right)
\mathcal L_s(w^{(k)}).
$$
Here monotonicity and homogeneity handle the first inequality; the second uses $|J_k|\le c$ and $0<\theta_i\le1$. Sum and use (22.5). ∎

One may use only the circuits occurring in a chosen splitting, with their individual survival probabilities in (22.11). Thus $c$ is a convenient uniform bound, not always the sharpest description of a particular data set.

The power loss is not essential either. The same proof works for differentiable convex row losses $\rho_i$, provided the minimum is attained, $\rho_i(0)=0$, and $\rho_i'(r)r>0$ for $r\ne0$: use $\lambda_i=w_i\rho_i'(r_i)$ in the splitting.

**Proposition 22.4 (both the exponent and independence order are sharp).** For every $d\ge1$, there is a system of $d+1$ equations in $d$ parameters with positive full loss, such that:

1. Independent common retention gives $\mathbb E\mathcal L_s(V)=\theta^d\mathcal L_s(w)$.
2. A $d$-wise independent fair retention law gives $\mathbb E\mathcal L_s(V)=0$.

*Proof.* Take distinct real $t_1,\ldots,t_{d+1}$,
$$
a_i=(1,t_i,\ldots,t_i^{d-1}),\qquad b_i=t_i^d.
$$
Every proper subset interpolates, while the full set does not. Only retaining all $d+1$ equations has positive loss. Independent retention has probability $\theta^{d+1}$ for that event, and inverse weighting multiplies its loss by $\theta^{-1}$.

For the second claim, use the uniform law on binary masks whose parity is opposite to the all-ones mask. Any $d$ coordinates are uniform: they uniquely determine the remaining coordinate. The full mask never occurs, so every sample can be fitted exactly. ∎

This also makes the $c$-wise threshold sharp by setting $d=c-1$. A universal guarantee based on the marginals alone is impossible once $c>1$.

**Example 22.5 (rank can overstate the obstruction).** Take $d$ unrelated contradictory pairs $\beta_j=1,\ \beta_j=-1$, with unit weights. Every inconsistent circuit has size two. For squared loss, $\mathcal L_2(w)=2d$ and independent common retention gives $\mathbb E\mathcal L_2(V)=2d\theta$. The circuit factor $\theta$ is exact; the rank-only factor $\theta^d$ can be much weaker. The same equality holds for every $s>1$.

## 4. The quadratic determinant refinement

**Theorem 22.6 (determinant bound).** Assume that inclusion moments match the independent Bernoulli law for sets of size at most $d+1$. Put $W=\operatorname{diag}w$. Then
$$
\boxed{\kappa\,\mathcal L_2(w)\le\mathbb E\mathcal L_2(V)
\le\mathcal L_2(w),\qquad
\kappa=\frac{\det(A^{\mathsf T}WA)}
{\det(A^{\mathsf T}\operatorname{diag}(w/\theta)A)}.}
\tag{22.12}
$$
Moreover $\kappa\ge\prod_{j=1}^d\theta_{(j)}$. For common $\theta$, $\kappa=\theta^d$.

*Proof.* For nonnegative $v$, define
$$
\Gamma(v)=\det(A^{\mathsf T}\operatorname{diag}(v)A),\qquad
\Delta(v)=\det([b\ A]^{\mathsf T}\operatorname{diag}(v)[b\ A]).
$$
The Schur complement gives
$$
\Delta(v)=\Gamma(v)\mathcal L_2(v).
\tag{22.13}
$$
This identity also holds when the retained design has rank below $d$: then both determinants are zero. **The loss itself can still be positive**; no division by zero is permitted.

Cauchy–Binet expresses $\Gamma$ as a nonnegative homogeneous multiaffine polynomial of degree $d$, and $\Delta$ as one of degree $d+1$. Since $V\le w/\theta$ coordinatewise,
$$
\Gamma(w/\theta)\mathcal L_2(V)\ge\Delta(V),
\qquad \mathbb E\Delta(V)=\Delta(w).
$$
This proves the lower bound. The upper bound was proved in Theorem 22.3 using only marginals. Comparing the degree-$d$ monomials of $\Gamma(w/\theta)$ with those of $\Gamma(w)$ gives the stated simpler factor. ∎

The circuit and determinant bounds supply different refinements: short contradictions improve the exponent, while the determinant retains the full geometry and unequal weights. No ordering between those two refinements is asserted.

For squared loss on a single circuit with dependence $g$ and positive weights $u_i$ on its support,
$$
\mathcal L_2(u)=\frac{(g^{\mathsf T}b)^2}{\sum_{i\in J}g_i^2/u_i}.
\tag{22.14}
$$
Indeed its optimal weighted residual lies in $\ker A_J^{\mathsf T}=\mathbb R g$, and its pairing with $g$ is $g^{\mathsf T}b$. Formula (22.14) makes every summand of (22.11) explicit.

There is also a classical volume-sampling interpretation. Give each size-$d$ row basis probability
$$
\pi(J)=\frac{\det(A_J)^2\prod_{i\in J}w_i}{\Gamma(w)}.
$$
Then $\kappa^{-1}=\mathbb E_\pi\prod_{i\in J}\theta_i^{-1}$. This auxiliary basis law is **not** the deletion law.

**Corollary 22.7 (fixed ridge penalty).** For $\lambda>0$, let
$$
\mathcal L_{2,\lambda}(v)=\min_\beta
\left[\sum_i v_i(b_i-a_i^{\mathsf T}\beta)^2+\lambda\|\beta\|^2\right].
$$
With the same moment assumption,
$$
\frac{\det(R+\lambda I)}
{\det(A^{\mathsf T}\operatorname{diag}(w/\theta)A+\lambda I)}
\mathcal L_{2,\lambda}(w)
\le\mathbb E\mathcal L_{2,\lambda}(V)
\le\mathcal L_{2,\lambda}(w),
\quad R=A^{\mathsf T}WA.
\tag{22.15}
$$
*Proof.* Add the $d$ equations $\beta_j=0$ with weight $\lambda$ and retention probability one, then apply Theorem 22.6. These rows make the augmented design full rank even if the original $A$ is not. ∎

For common random-row retention $\theta$ and eigenvalues $\sigma_j$ of $R$, the lower factor is
$$
\prod_j\left(1+(\theta^{-1}-1)\frac{\sigma_j}{\sigma_j+\lambda}\right)^{-1}.
$$
The penalty is fixed, not erased or rescaled.

## 5. What returns to the ETH problem

**Energy.** In each scalar sector of ETH 21, take
$$
a_i=(t_i,\ldots,t_i^D),\qquad b_i=1.
$$
The sign of the fitted coefficients is immaterial. Distinct nonzero nodes make every $D$ rows independent and every $D+1$ rows inconsistent. Thus $c=D+1$. Theorem 22.3 recovers the factor $\theta^D$ with no Vandermonde formula, even for scalar $\ell_s$ loss. Theorem 22.6 recovers ETH 21's determinant comparison. Only the squared original vertex energy has the tensor reduction of ETH 21; a nonquadratic scalar theorem does not establish a nonquadratic vertex reduction.

**Exact locality.** Write
$$
n_x(T)=1+\sum_{k=1}^D c_{xk}T^k.
$$
A good edge $\{x,y\}$ in direction $i$ is the linear equation
$$
\sum_{k=1}^D t_i^k(c_{xk}+c_{yk})=-2.
\tag{22.16}
$$
Lemma 22.1 turns the minimum defect support, or the smallest possible maximum directional load, exactly into hitting its minimal inconsistent circuits. Here there are $D2^Q$ coefficient variables before removing redundancy: **the scalar parameter count $D$ must not be substituted for this much larger rank.**

For $D=1$, $h_x=\operatorname{par}(x)c_{x1}$ turns (22.16) into prescribed differences. The row circuits are simple graph cycles; an inconsistent circuit has nonzero total prescribed circulation. With nodes $t_i=i$, exhaustive exact enumeration gives:

| Cube | Inconsistent circuits | Circuit lengths | Minimum total failed edges | Minimum largest directional count |
|---|---:|---|---:|---:|
| $Q=2$ | 1 | 4 | 1 | 1 |
| $Q=3$ | 24 | 4, 6, 8 | 4 | 2 |

Checking only surviving square faces is insufficient after deletion: a retained inconsistent six-cycle can have every square face broken. The circuit reformulation therefore needs all minimal contradictions, not just the original plaquettes.

These are small exact optima, not a sharper asymptotic answer to **T.ETH.7**. A useful next step is a circuit packing or hitting certificate that respects direction classes and improves a proved bound.

## 6. Checks and limits

~~~text
python reviews/gpt/eth/scripts/verify_circuit_erasure.py
~~~

The standard-library verifier uses exact rational normal equations, independent Gram determinants and minor sums, enumerated sign-compatible splittings, and all masks of 29 small full-rank designs. It checks unequal weights and probabilities, zero-rank models, rank-deficient retained models, $(d+1)$-wise parity laws, sharpness and the $d$-wise failure for $d=1,\ldots,4$, short-circuit and fourth-power fixtures, three fixed-penalty models, and all deletion sets in the two cubes above. Results: [circuit_erasure.json](results/circuit_erasure.json).

Moment matching guarantees bounds, not the exact optimized expectation. At $A=(1,2,3)^{\mathsf T}, b=(1,1,1)$ and $\theta_i=1/2$, independent retention gives $503/1820$, whereas the pairwise independent even-parity law gives $22/65$. Both satisfy the rank-one theorem's bounds.

**Example 22.8 (nonlinear models can fail).** Let the allowed residual vectors be the six permutations of $(0,1/4,1)$, and minimize $\sum_i v_i r_i^2$ over those states. Each singleton can be fitted exactly, while every pair is inconsistent, so the largest minimal contradiction has size two. Nevertheless, the full loss is $17/16$, and independent retention with $\theta=1/2$ gives expected inverse-weighted loss $5/16<17/32$. The affine model and its sign-compatible kernel decomposition are substantive hypotheses; a hypergraph of contradictions alone does not imply the quantitative loss bound.

The general $s>1$ assertion has a separate convexity proof; the script's nonquadratic checks are the stated fourth-power fixtures, not a generic nonlinear solver. No Lean proof, exact-counting consequence, efficient circuit enumeration, or new SAT exponent is claimed.

## 7. Attribution

Sign-compatible decomposition into support-minimal kernel vectors is classical. See Stefan Müller and Georg Regensburger, [*Elementary Vectors and Conformal Sums in Polyhedral Geometry and their Relevance for Metabolic Pathway Analysis*](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2016.00090/full), Theorem 3 (2016), which credits Rockafellar's 1969 elementary-vector theorem and corrects an overstrong auxiliary support claim. We use only the decomposition and give its elementary support-reduction proof.

Cauchy–Binet, the Schur complement and least-squares projection are standard. For volume sampling and regression identities, see Michał Dereziński and Manfred K. Warmuth, [*Reverse Iterative Volume Sampling for Linear Regression*](https://jmlr.csail.mit.edu/papers/v19/17-781.html), JMLR 19(23), 2018. Their prediction-loss problem and sampling law differ from the inverse-probability training residual after deletion studied here.

The repo contribution is this explicit passage from circuit splitting to erasure bounds and a sharp independence threshold, together with its two distinct ETH interfaces. Novelty beyond this project has not been established.

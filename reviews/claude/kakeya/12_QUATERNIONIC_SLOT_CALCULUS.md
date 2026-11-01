# 12. The quaternionic slot calculus: a standalone toolkit

*Claude, 28 October 2026. Roadmap target `T.KAK.3`; nodes `KAK.Q.CALC`, `KAK.FRAME`, `KAK.MAENERGY` in [`ROADMAP.md`](../../../ROADMAP.md).*

The paper's analysis runs on a small algebraic calculus:
- quaternionic gradients and Hessians of real functions;
- the Moore determinant;
- the fact that its cofactor is divergence-free, a quaternionic null Lagrangian.

[08](08_GAPLESS_PROOF.md) §2 proves what the paper needs, in the paper's setting of 3 or 4 real coordinates per slot. This file states the calculus so it can be used outside the paper, in its natural generality: each slot carries an **arbitrary finite frame of quaternions**.

> **Theorem 12.1 (cofactor divergence for arbitrary frames).** Give each slot $a\le m$ a finite list of quaternions $f_{a,1},\dots,f_{a,k_a}$, of any length and not necessarily independent. For a real $C^3$ function $u$ of the variables $y_{a,i}$, put
> $$H_f(u)_{ab}=\sum_{i,j}f_{a,i}\,\overline{f_{b,j}}\,\partial_{a,i}\partial_{b,j}u .$$
> Then the Moore cofactor of $H_f(u)$ is divergence-free in the frame sense:
> $$\sum_{a,i}\overline{f_{a,i}}\,\partial_{a,i}\,C(H_f(u))_{ab}=0\qquad(b\le m).$$

> **Lemma 12.2 (frame reduction).** Every frame calculus is a restriction of the 4-slot calculus. If $A_{a,i\alpha}$ are the real components of $f_{a,i}$, then $H_f(u)\circ A=H^{(4)}(u\circ A)$, and the frame divergence pulls back to the 4-slot divergence.

> **Theorem 12.3 (the Monge–Ampère operator is the gradient of a concave energy).** Let the frames be independent. On $\mathcal P=\{u:u|_{\partial\Omega}=0,\ H_f(B_0+u)>0\}$, the energy $E(u)=\int_\Omega u\int_0^1\det H_f(B_0+\theta u)\,d\theta$ has
> $$dE(u)[v]=\int_\Omega v\,\det H_f(B_0+u),\qquad d^2E(u)[v,v]=-\int_\Omega\operatorname{Re}\big(G(v)^*C(H_f(B_0+u))G(v)\big)<0\quad(v\ne0).$$
> So $\det H_f(B_0+u)=g$ is the Euler–Lagrange equation of the strictly concave functional $E(u)-\int gu$.

> **Proposition 12.4 (ellipticity ⇔ independent frames).** For $Q>0$, the linearization $w\mapsto Q^{-1}:H_f(w)$ is elliptic iff every slot frame is $\mathbb R$-linearly independent. In particular each slot has at most $4$ real coordinates.

> **Proposition 12.5 (the 3-slot Hessian is real ⊕ complex).** In a 3-slot $X+P\mathrm i$, $P=p+\mathrm kq$, write $\bar\partial=\frac12(\partial_p+\mathrm k\partial_q)$ and $\partial=\frac12(\partial_p-\mathrm k\partial_q)$. Then
> $$H_{ab}(u)=\partial_{X_a}\partial_{X_b}u+4\,\bar\partial_a\partial_bu+2\big(\bar\partial_a\partial_{X_b}u-\partial_{X_a}\bar\partial_bu\big)\,\mathrm i .$$
> That is: the real Hessian in the real parts, plus four times the complex Hessian in the imaginary parts, plus an $\mathrm i$-twisted coupling in $\operatorname{span}(\mathrm i,\mathrm j)$.

**What is new.**
- Theorem 12.1 for arbitrary frames: non-orthonormal, dependent, and of any length. The paper and 08 state it for the frames $(1,\mathrm i,\mathrm j)$ and $(1,\mathrm i,\mathrm j,\mathrm k)$, and OpenAI's Lean development for the slot systems it uses ([07](07_TOOLKIT.md) K6). The proof is the same; the general statement shows that only the rank-one structure and associativity matter. I have not searched the literature for it.
- Lemma 12.2, which explains Theorem 12.1.
- Theorem 12.3, the variational structure behind the energy of 08 §6.
- Propositions 12.4 and 12.5.

**Checks.** [`scripts/verify_slot_calculus.py`](scripts/verify_slot_calculus.py) ([`results/slot_calculus.json`](results/slot_calculus.json), ~12 s; exact sympy unless stated):
- **Part F:** Theorem 12.1 for five frames, with $m=2,3$ slots and random quartic $u$. The frames include mixed sizes $(1,2,4)$ and a dependent 5-vector frame. Two negative controls are nonzero: the divergence taken with $f$ instead of $\bar f$, and a Hermitian field that is not a frame Hessian.
- **Part R:** Lemma 12.2.
- **Part K:** the real and complex cases.
- **Part D:** Proposition 12.5.
- **Part E:** Proposition 12.4, numerically: positive symbol for an independent $(3,4,2)$ frame, a null covector for a 5-vector frame.
- **Part U:** Euler's identity and $C(A)>0$ for $A>0$, numerically.

The paper's cases are also checked by [`scripts/verify_quaternion_calculus.py`](scripts/verify_quaternion_calculus.py).

---

## 1. Frames, gradients and Hessians

**Conventions.** These are as in 08 §1:
- vectors are columns with scalars acting on the right, and $^*$ is the quaternionic conjugate transpose;
- $A:B=\operatorname{Re}\operatorname{tr}(AB)$;
- $\det$ is the Moore determinant of a quaternionic Hermitian matrix.

**Definition 1.1 (frame calculus).** A *frame* on $m$ slots is a family $f=(f_{a,i})_{a\le m,\ i\le k_a}$ of quaternions. The variables are real numbers $y_{a,i}$. For real $u$:
$$G_a(u)=\sum_if_{a,i}\,\partial_{a,i}u,\qquad H_f(u)_{ab}=\sum_{i,j}f_{a,i}\overline{f_{b,j}}\,\partial_{a,i}\partial_{b,j}u,\qquad A^{\mathbb R}_{(a,i),(b,j)}=\operatorname{Re}\big(\overline{f_{a,i}}\,A_{ab}\,f_{b,j}\big).$$
A Hermitian field $A$ is **conserved** if $\sum_{a,i}\overline{f_{a,i}}\partial_{a,i}A_{ab}=0$ for every $b$. Then every real column of $A^{\mathbb R}$ is divergence-free.

**Examples.**

| Frame per slot | $H_f(u)$ | $\det H_f(u)$ |
|---|---|---|
| $(1)$ | real Hessian $D^2u$ | real Monge–Ampère |
| $(1,\mathrm i)$ | $4\,\partial_{\bar z_a}\partial_{z_b}u$ | complex Monge–Ampère |
| $(1,\mathrm i,\mathrm j,\mathrm k)$ | Alesker's quaternionic Hessian | quaternionic Monge–Ampère |
| $(1,\mathrm i,\mathrm j)$ | the paper's 3-slot Hessian, real ⊕ complex (Prop. 12.5) | the paper's fibre equation |

**Lemma 1.2 (Hessian calculus for frames).**
1. $H_f(u)$ is Hermitian. $H_f(\xi\xi^{\mathsf T})=gg^*$ with $g_a=\sum_if_{a,i}\xi_{a,i}$. Hence $D^2u\ge0$ implies $H_f(u)\ge0$.
2. **Product rule.** $H_f(uv)=uH_f(v)+vH_f(u)+G(u)G(v)^*+G(v)G(u)^*$.
3. **Green's form.** $A:H_f(u)=A^{\mathbb R}:D^2u$, and $\operatorname{Re}\big(G(v)^*AG(u)\big)=\nabla v^{\mathsf T}A^{\mathbb R}\nabla u$. If $A$ is conserved, $A:H_f(u)=\operatorname{div}(A^{\mathbb R}\nabla u)$.
4. $\operatorname{tr}H_f(u)=\sum_a\sum_{i,j}\operatorname{Re}(f_{a,i}\overline{f_{a,j}})\,\partial_{a,i}\partial_{a,j}u$. This is the Laplacian for orthonormal frames.

*Proof.*
1. **Hermitian.** $\overline{H_{ba}}=\sum\overline{f_{b,i}\bar f_{a,j}}\partial\partial u=\sum f_{a,j}\bar f_{b,i}\partial_{a,j}\partial_{b,i}u=H_{ab}$.
2. **Rank one.** Insert $D^2u=\xi\xi^{\mathsf T}$; the map $D^2u\mapsto H_f(u)$ is real-linear.
3. **Product rule.** $D^2(uv)=uD^2v+vD^2u+\nabla u\otimes\nabla v+\nabla v\otimes\nabla u$. The last two terms are $\frac12[(\nabla u+\nabla v)^{\otimes2}-(\nabla u-\nabla v)^{\otimes2}]$; apply 1.
4. **Green's form.** $\operatorname{Re}(A_{ba}f_{a,i}\bar f_{b,j})=\operatorname{Re}(\bar f_{b,j}A_{ba}f_{a,i})$. Conservation lets the derivative pass through $A^{\mathbb R}$. ∎

## 2. The Moore determinant toolkit

**Lemma 2.1.** For Hermitian $A\in\mathbb H^{m\times m}$, with cofactor $C(A)$ defined by $C(A):B=\frac d{dt}\det(A+tB)|_{t=0}$:
1. $C(A)=(\det A)A^{-1}$ for $A>0$. $\log\det$ is concave and increasing on $\{A>0\}$, and $C(A)>0$ there.
2. If $A\ge0$ is singular, then $C(A)\ge0$ vanishes on the image of $A$.
3. **Rank-one expansion.** $v^*C(A)v=\det(A|_{v^\perp})$ and $\det(A+tvv^*)=\det A+t\,v^*C(A)v$ for unit $v$. Along a segment from a rank-one $B$:
$$\det((1-\theta)B+\theta D)=\theta^m\det D+(1-\theta)\theta^{m-1}C(D):B .$$
4. **Schur.** For $M=\begin{pmatrix}\ell&c^*\\c&D\end{pmatrix}$ with $D>0$: $\det M=\det D\,(\ell-c^*D^{-1}c)$, and $M-(\ell-c^*D^{-1}c)e_0e_0^*\ge0$.
5. **Euler.** $C(A):A=m\det A$.

*Proof.*
- **Items 1–4.** These are 08 Lemma 2.2, which does not involve slots. The segment formula is 3 applied to $D+\frac{1-\theta}\theta B$, using homogeneity (08 Lemma 6.3).
- **Item 5.** $\det$ is homogeneous of degree $m$, so $\frac d{dt}\det((1+t)A)|_0=m\det A$. ∎

## 3. Frame reduction

**Lemma 12.2.** Write $f_{a,i}=\sum_{\alpha=0}^3A_{a,i\alpha}e_\alpha$ with $(e_\alpha)=(1,\mathrm i,\mathrm j,\mathrm k)$. Let $A$ act slotwise as a linear map $\mathbb R^{4m}\to\mathbb R^{\sum k_a}$, $(Ax)_{a,i}=\sum_\alpha A_{a,i\alpha}x_{a,\alpha}$. Then for every real $u$ and every Hermitian field $M(y)$:
$$H^{(4)}(u\circ A)=H_f(u)\circ A,\qquad\sum_{a,\alpha}\bar e_\alpha\partial_{a\alpha}(M\circ A)_{ab}=\Big(\sum_{a,i}\overline{f_{a,i}}\,\partial_{a,i}M_{ab}\Big)\circ A .$$

*Proof.* The chain rule gives $\partial_{a\alpha}(u\circ A)=\sum_iA_{a,i\alpha}(\partial_{a,i}u)\circ A$. Also $\sum_\alpha A_{a,i\alpha}e_\alpha=f_{a,i}$ and $\sum_\alpha A_{a,i\alpha}\bar e_\alpha=\overline{f_{a,i}}$. ∎

So every frame identity that holds in 4 slots holds for every frame. This needs one remark: if $A$ is not surjective (dependent frames), an identity proved at the points $Ax$ is a polynomial identity in the jets of $u$ at those points. It therefore holds at $0=A0$ for every $u$, and everywhere by translation.

## 4. Cofactor divergence

**Proof of Theorem 12.1.** This is the proof of 08 Lemma 2.3. It uses only Lemma 1.2.1 and Lemma 2.1.2, so it applies to any frame. For $m=1$, $C\equiv1$. Otherwise let $j_0=m-1$ and let $Q(M_1,\dots,M_{j_0})$ be the symmetric multilinear form with $Q(M,\dots,M)=j_0!\,C(M)$.
1. **Rank-one inputs.** For $M_j=g_jg_j^*$, the matrix $\sum t_jg_jg_j^*\ge0$ ($t_j>0$) has rank $<m$. By Lemma 2.1.2 its cofactor kills every $g_i$; by polynomiality in $t$, $Q(g_1g_1^*,\dots)g_i=0$.
2. **Exponentials.** For $u_j=e^{\xi_j\cdot y}$, Lemma 1.2.1 gives $H_f(u_j)=u_jg_jg_j^*$. The frame divergence of $e^{(\sum\xi_j)\cdot y}Q(g_1g_1^*,\dots)$ is $e^{(\sum\xi_j)\cdot y}\big(\sum_jg_j\big)^*Q(\dots)=0$, since $\sum_{a,i}\overline{f_{a,i}}(\sum_j\xi_j)_{a,i}=(\sum_jg_j)^*_a$.
3. **General inputs.** The divergence is a multilinear constant-coefficient expression in $(u_1,\dots,u_{j_0})$. On exponentials it is $e^{(\cdot)}$ times a polynomial in the $\xi_j$, which vanishes identically. Polarize to $u_j=u$. ∎

**Remarks.**
1. **What the proof needs.** It needs only two facts:
   - the Hessian maps rank-one real forms to rank-one quaternionic forms;
   - the Moore cofactor of a positive semidefinite matrix of rank $<m$ kills its image.
   
   The second needs the quaternionic spectral theorem, hence **associativity**. Among the real division algebras, frames in $\mathbb H$ cover all the associative ones ($\mathbb R\subset\mathbb C\subset\mathbb H$); the octonions are not associative.
2. **Lemma 12.2 gives a second proof.** Apply it to the 4-slot case.

**Corollary 4.1 (Green's identity and compression).** These are 08 Lemma 2.4 for any frame, with the same proofs; they use only Lemma 1.2.3:
- for conserved $A$ on $\Omega=\{r<0\}$, $\int(vA:H_f(u)-uA:H_f(v))=\int_{\partial\Omega}\operatorname{Re}(v\nu^*AG(u)-u\nu^*AG(v))\,dS$;
- if $u=0$ on $\partial\Omega$, then $H_f(u)|_{\nu^\perp}=\rho H_f(r)|_{\nu^\perp}$ there.

## 5. Ellipticity and comparison

**Proof of Proposition 12.4.** Let $L_Qw=Q^{-1}:H_f(w)$. By Lemma 1.2, its principal symbol at a covector $\xi$ is $\sigma(\xi)=\hat g(\xi)^*Q^{-1}\hat g(\xi)$ with $\hat g_a=\sum_if_{a,i}\xi_{a,i}$. Since $Q^{-1}>0$, $\sigma(\xi)=0$ iff $\hat g(\xi)=0$. The map $\xi\mapsto\hat g(\xi)$ is injective iff each slot frame is linearly independent in $\mathbb H\cong\mathbb R^4$, which forces $k_a\le4$. ∎

**Corollary 5.1 (comparison).** For independent frames, 08 Lemma 5.1 holds verbatim. If $H_f(\psi_i)>0$, $\det H_f(\psi_1)\ge\det H_f(\psi_2)$ in $\Omega'$, and $\psi_1\le\psi_2$ on $\partial\Omega'$, then $\psi_1\le\psi_2$.

**Why the paper uses 3 coordinates per slot.**
- **Upper bound.** Proposition 12.4 gives at most 4 coordinates per slot.
- **The time slot.** It must have exactly 3: $T$ plus a two-dimensional auxiliary time. This is the capacity-dimension principle ([07](07_TOOLKIT.md) K4).
- **The spatial slots.** These have 3 because the lifted current needs one real coordinate plus two auxiliary flows. Two flows combine into one $\mathbb C_{\mathrm k}$-coordinate ([08](08_GAPLESS_PROOF.md) §3).
- **The fourth unit.** It enters only in the proof of the boundary $C^2$ estimate (08 Thm 5.2 step 5), through Lemma 12.2, so that quaternionic rigid motions are available.

## 6. The Monge–Ampère energy

**Proof of Theorem 12.3.** Let the frames be independent and $\Omega$ bounded and smooth, possibly periodic in some directions. Let $B_0$ be smooth with $H_f(B_0)\ge0$; then $\mathcal P$ is convex. Let $u\in\mathcal P$ and $v\in C^2(\bar\Omega)$ with $v|_{\partial\Omega}=0$. Put $M_\theta=C(H_f(B_0+\theta u))$, conserved by Theorem 12.1.
1. **First variation.** Differentiate under the integral:
$$\frac d{d\epsilon}E(u+\epsilon v)\Big|_0=\int v\int_0^1\det H_f(B_0+\theta u)\,d\theta+\int u\int_0^1\theta\,M_\theta:H_f(v)\,d\theta .$$
   Green's identity (Corollary 4.1, with $u=v=0$ on $\partial\Omega$) turns $\int uM_\theta:H_f(v)$ into $\int vM_\theta:H_f(u)=\int v\,\partial_\theta\det H_f(B_0+\theta u)$. The two terms combine to
$$\int v\int_0^1\partial_\theta\big(\theta\det H_f(B_0+\theta u)\big)\,d\theta=\int v\det H_f(B_0+u).$$
2. **Second variation.** Differentiate $\int v\det H_f(B_0+u+\epsilon v)$:
$$\int vM_1:H_f(v)=\int v\operatorname{div}(M_1^{\mathbb R}\nabla v)=-\int\nabla v^{\mathsf T}M_1^{\mathbb R}\nabla v=-\int\operatorname{Re}(G(v)^*M_1G(v)).$$
   Here $M_1=C(H_f(B_0+u))>0$ (Lemma 2.1.1). Independent frames make $M_1^{\mathbb R}$ positive definite, so this is $<0$ unless $\nabla v\equiv0$, that is $v=0$.
3. **Conclusion.** $J(u)=E(u)-\int gu$ is strictly concave on the convex set $\mathcal P$, with $dJ(u)=\det H_f(B_0+u)-g$. So a solution of the Dirichlet problem is the unique critical point and the unique maximizer of $J$. ∎

**Relation to 08 §6.** The energy there is $\mathcal F=X_0+E_z(u)$, where $E_z(u)=\lambda^{-1}E(u)$ and $X_0=\int(-u)\mu$. Since $\det D_u=\lambda\mu$, $\mathcal F=-\big(\lambda^{-1}\int u\det H(B_0+u)-E_z(u)\big)$, the Legendre-type transform of the concave $E$ evaluated at the solution. Proposition 6.2 there computes how this transform varies with the parameters $z$, and Lemma 6.3 compares it with $X_0$. In complex geometry this is the Aubin–Mabuchi energy: its gradient is the Monge–Ampère operator, and the energy is concave along the space of potentials.

## 7. The 3-slot structure

**Proof of Proposition 12.5.**
- **Rewriting the operator.** Write the frame operator as $\sum_\alpha e_\alpha\partial_{a\alpha}=\partial_{X_a}+2\bar\partial_a\,\mathrm i$, using $\mathrm j=\mathrm k\mathrm i$. Its conjugate is $\partial_{X_b}-\mathrm i\,2\partial_b$.
- **Expanding the product.** $H_{ab}=(\partial_{X_a}+2\bar\partial_a\mathrm i)(\partial_{X_b}-\mathrm i2\partial_b)u$. Here $\mathrm i\,\mathrm i=-1$, and $\mathrm ic=\bar c\,\mathrm i$ for $c\in\mathbb C_{\mathrm k}$. ∎

**Consequences.**
1. A function of the real parts only has $H=D^2_Xu$, the real Hessian.
2. A function of the imaginary parts only has $H=4\bar\partial\partial u$, the complex Hessian with Wirtinger derivatives for $\mathrm k$. In particular $H(\sum|F_\nu|^2)=4J^*J$ for holomorphic $F$ (08 Lemma 2.1(4)).
3. In the paper's fibres the background $B_0=\lambda|\Phi|^2$ is a function of the imaginary parts, so $H(B_0)$ is $\mathbb C_{\mathrm k}$-Hermitian of rank one. The solution $u$ depends on both: the fibre equation $\det H_{\rm fib}(B_0+u)=\lambda\mu$ is a mixed real–complex Monge–Ampère equation.

## 8. Dictionary

| Here | 08 | Paper |
|---|---|---|
| Lemma 1.2 | Lemma 2.1 | Lemma 3.1 |
| Lemma 2.1 | Lemma 2.2, Lemma 6.3 (segment) | Lemma 3.2, Lemma 4.3 |
| Theorem 12.1 | Lemma 2.3 (frames $(1,\mathrm i,\mathrm j)$, $(1,\mathrm i,\mathrm j,\mathrm k)$) | Lemma 3.3 |
| Corollary 4.1 | Lemma 2.4 | §3 |
| Corollary 5.1 | Lemma 5.1 | App. A |
| Lemma 12.2 | Theorem 5.2 step 5 (the 4th slot) | App. A |
| Theorem 12.3 | §6 (energy) | §4 |
| Propositions 12.4, 12.5 | — | — |

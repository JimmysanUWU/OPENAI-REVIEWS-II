# 22. Curvature forces an orthogonal design; every multiple-sheet full fibre diverges

GPT, 1 November 2026. **Independent prose derivation, supported by exact algebra and Jacobian checks; not Lean-formalized.** This extends the obstruction in [Claude 20](../../claude/kakeya/20_GAUGE_ORBIT_OBSTRUCTION.md) beyond Hermitian gauge groups. In particular, it settles the octonionic full-fibre convergence question left open in [Claude 21, §2](../../claude/kakeya/21_GAUGE_QUOTIENT.md).

The mechanism is the second variation of a constant-multiplicity eigenvalue. At the dimension required for tangency, it becomes a square orthogonal design. Its determinant computes the conditional measure on the complexified polarization sphere. No classification of pencils or new dispersive estimate is claimed.

## 0. Setting

Use the setting of [Claude 17](../../claude/kakeya/17_STATE_DIMENSION_RIGIDITY.md), nodes **KAK.RIGID** and **KAK.C.LIFT**. After compression to the span of the extremal eigenvectors:

- $C^0=I_N,\ C^1,\ldots,C^n$ are real symmetric;
- $\varphi(\xi)=\lambda_{\max}(\sum_j\xi_jC^j)$ has constant multiplicity $m$ on an open cone;
- its Hessian is positive definite transverse to the radial direction;
- the tangency dimension is $N=n+m-1$, and $n\ge2$.

Normalize $\xi_0=e_n$. Set $k=n-1$, $E=\ker(\varphi_0I-C^n)$, $F=E^\perp$, so $\dim E=m,\dim F=k$. In orthonormal coordinates $z=(a,b)\in E_{\mathbb C}\oplus F_{\mathbb C}$, put
$$
S=(\varphi_0I-C^n)|_F>0,\quad
d_j=\partial_j\varphi(\xi_0),\quad
H=(\partial_j\partial_l\varphi(\xi_0))_{j,l\le k}>0.
$$
All complexification uses the auxiliary unit $\mathrm k$; dots below are bilinear unless a norm is displayed.

For $a\in E_{\mathbb C}$, define the square matrix
$$
R_a:F_{\mathbb C}\longrightarrow\mathbb C^k,\qquad
(R_ab)_j=b^{\mathsf T}(C^j-d_jI)a.
\tag{22.1}
$$
Hellmann–Feynman gives $P_EC^jP_E=d_jP_E$, so each vector $(C^j-d_jI)a$ lies in $F_{\mathbb C}$.

## 1. The second variation is an orthogonal design

**Lemma 22.1 (curvature identity).**
$$
\boxed{R_aS^{-1}R_a^{\mathsf T}=\tfrac12(a\cdot a)H
\quad\text{for every }a\in E_{\mathbb C}.}
\tag{22.2}
$$
Consequently
$$
(\det R_a)^2=\det S\,\det(H/2)\,(a\cdot a)^k.
\tag{22.3}
$$

*Proof.* First let $a\in E$ be real. Along the perturbation $\xi_0+s e_l$, choose a differentiable eigensection $a(s)$ with $a(0)=a$, constant norm, and $P_Ea'(0)=0$. This follows by choosing an orthonormal frame for the constant-rank spectral projection and removing its internal rotation at zero. Differentiating the eigenvector equation gives
$$
a'(0)=S^{-1}(C^l-d_lI)a.
$$
Differentiate $a(s)^{\mathsf T}C^ja(s)=\partial_j\varphi(\xi_0+s e_l)\,|a|^2$. The result is
$$
H_{jl}|a|^2
=2\,((C^j-d_jI)a)^{\mathsf T}S^{-1}((C^l-d_lI)a).
$$
This is (22.2) for real $a$. Both sides are quadratic polynomials in $a$, so the identity extends to every complex $a$. Take determinants. ∎

**Corollary 22.2 (a necessary Clifford structure).** After fixed invertible scalings, $a\mapsto R_a$ is a linear orthogonal design on $\mathbb R^k$. More precisely, there are real $k\times k$ matrices $J_2,\ldots,J_m$ with
$$
J_\alpha^{\mathsf T}=-J_\alpha,\qquad
J_\alpha J_\beta+J_\beta J_\alpha=-2\delta_{\alpha\beta}I.
\tag{22.4}
$$
In particular, $k=n-1$ is even if $m>1$.

*Proof.* Set $L(a)=\sqrt2\,H^{-1/2}R_aS^{-1/2}$. Then $L(a)L(a)^{\mathsf T}=|a|^2I$. Right-multiply every $L(a)$ by the orthogonal matrix $L(e_1)^{\mathsf T}$, making $L(e_1)=I$. Expand this identity at $a=e_1+t e_\alpha$ and at $a=e_\alpha+t e_\beta$. The linear and quadratic coefficients give (22.4). An invertible skew-symmetric real matrix has even size. ∎

This is a necessary local structure, not a sufficiency or classification theorem. No choice of the lower blocks of the pencil has been reconstructed from these relations.

## 2. Conditional measure on the central complex sphere

Use the single lift from **KAK.C.LIFT**, with
$$
Q_\mu(z)=\tfrac12z^{\mathsf T}C^\mu z,\qquad
\Theta(z,\eta)=(Q_0,\ldots,Q_{n-1},Q_n+\eta).
$$
The passive direction is transverse at the centre: its pairing with the normal $(-\varphi_0,\xi_0)$ equals one. The determinant-one target change
$$
y_0=Q_0,\quad y_j=Q_j-d_jQ_0\ (1\le j\le k),\quad
y_n=Q_n-\varphi_0Q_0+\eta
\tag{22.5}
$$
sends the centre to $(1/2,0,\ldots,0)$. Its fibre contains
$$
\mathcal Q=\{(a,0,0):a\in E_{\mathbb C},\ a\cdot a=1\}.
\tag{22.6}
$$

**Lemma 22.3 (quadric coarea).** Every point of $\mathcal Q$ is a regular point of $\Theta$, and $\mathcal Q$ is open in that fibre. On the chart $a_1\ne0$, write $a'=(a_2,\ldots,a_m)$. Disintegration of ordinary complex Lebesgue measure in $(z,\eta)$, relative to target Lebesgue measure, gives
$$
\boxed{d\nu=\frac{\gamma}{|a_1|^2}\,d^{\,2m-2}a',
\qquad \gamma=\frac1{\det S\,\det(H/2)}.}
\tag{22.7}
$$

*Proof.* Adjoin $a'$ to the target coordinates (22.5). At $(a,0,0)$, the complex Jacobian from $(a_1,a',b,\eta)$ to $(y,a')$ has determinant, up to sign,
$$
a_1\det R_a.
\tag{22.8}
$$
Indeed the $y_0$ derivative in $a_1$ is $a_1$; the $y_1,\ldots,y_k$ derivatives in $b$ form $R_a$; and the $y_n$ derivative in $\eta$ is one. The other blocks do not change this determinant. On $a\cdot a=1$, (22.3) gives $|\det R_a|^2=1/\gamma$. A complex Jacobian contributes the squared absolute determinant to real volume. Ordinary change of variables now gives (22.7), including its coarea factor.

Some coordinate of $a$ is nonzero at every point of $\mathcal Q$; use it in place of $a_1$. Hence the target rank is $n+1$ everywhere on $\mathcal Q$. Its dimension $m-1$ equals the fibre dimension, so it is an open fibre submanifold. ∎

## 3. Universal divergence, without a gauge group

**Theorem 22.4 (every multiple-sheet full fibre diverges).** Under the setting of §0, if $m>1$, the unweighted full-fibre single lift has infinite physical-time density at every polarization centre for the plane wave polarized at that centre. More precisely, choose a unit polarization $e_1\in E$, frequency $\kappa>0$, and
$$
f(T,X)=e_1e^{\mathrm i\kappa(X\cdot\xi_0-T\varphi_0)}.
$$
On the part of $\mathcal Q$ with $a_1\ne0$, $|a'|<R$, its contribution is exactly
$$
\boxed{\frac{2\gamma\pi^{m-1}}{(m-1)!}\,R^{2m-2}.}
\tag{22.9}
$$
Thus the obstruction applies to all rigid constant-multiplicity pencils, not only complex or quaternionic Hermitian pencils. For the octonionic Hopf pencil, $m=8,\gamma=1$, and the contribution is $2\pi^7R^{14}/7!$.

*Proof.* Write $z=w_1+\mathrm k w_2$. On $\mathcal Q$, both $w_l$ lie in $E$, so the real displacement $q^\mu=\frac12\sum_lw_l^{\mathsf T}C^\mu w_l$ satisfies
$$
\varphi_0q^0-\xi_0\cdot q'
=\tfrac12\sum_l w_l^{\mathsf T}(\varphi_0I-C^n)w_l=0.
$$
The analytically continued plane wave is therefore undamped. Its physical-time entry before push-forward is $|a_1|^2$, by the lifted-current formula in Claude 08 Lemma 3.1 (also GPT 18, equation (18.8)). The map $\Theta$ leaves physical time unchanged. Multiplying by (22.7) cancels $|a_1|^2$, leaving $\gamma$ times flat volume in $a'$.

Except on the zero-measure hypersurface $1-\sum_{j\ge2}a_j^2=0$, there are exactly two roots
$$
a_1=\pm\sqrt{1-\sum_{j\ge2}a_j^2}.
$$
Summing the two sheets over the complex $(m-1)$-ball gives (22.9). It tends to infinity for $m>1$. Any other fibre components add nonnegative mass. Nonunit centres follow by rescaling the state variables; divergence persists. ∎

For $m=1$, the corresponding quadric consists of two points and this contribution is $2\gamma$; the argument does not assert divergence.

**Corollary 22.5 (scope of the unweighted single-lift route).** Under the standing extremality and curvature assumptions, a tangent unweighted full-fibre single lift can have finite central density only if $m=1$ and $N=n$, after compression to the eigenvector span.

*Proof.* Tangency forces $N=n+m-1$ by **KAK.RIGID**. Theorem 22.4 excludes $m>1$. ∎

Thus this particular route cannot extend App. E to larger effective state dimensions merely by integrating the same lifted current. A current built directly on a quotient, a modified tensor, or a reduction with compensating flux remains open (**T.KAK.7**). Neither all Hermitian engines nor all possible regularizations are ruled out.

## 4. Exact checks and limitations

~~~text
python reviews/gpt/kakeya/scripts/verify_curvature_fibre.py
~~~

The standard-library checker uses exact rational and Gaussian-rational arithmetic:

- the quadratic coefficients in (22.2), not just sample values;
- the complete ambient complex Jacobian, computed from the original pencil, against (22.8);
- the conditional density, damping, and the cancellation in (22.9);
- Hopf pencils over $\mathbb R,\mathbb C,\mathbb H,\mathbb O$, including their Clifford generators;
- a non-Clifford realified complex-Hermitian $3\times3$ pencil in five spatial variables, with nonzero lower blocks and unequal curvatures;
- controls for a missing $a_1$ Jacobian factor, positive central damping, confusing $q$ with $Q$, a perturbed orthogonal design, and removal of the square-dimension hypothesis.

Results: [curvature_fibre.json](results/curvature_fibre.json). The exact finite checks support the general second-variation and coarea proofs; no asymptotic claim is inferred from them. Unlike the Weyl complex-target profile in Claude 20, this note proves **infinite central density**, not non-integrability throughout a full target neighbourhood for every pencil.

# 18. A noncompact-fibre obstruction for the Hermitian single lift

GPT / Astra, 1 November 2026. **Complete derivation for the explicit Weyl example below; exact algebraic checks.**

Claude's [state-dimension analysis](../../claude/kakeya/17_STATE_DIMENSION_RIGIDITY.md) (`KAK.RIGID`, `T.KAK.7`) isolates an important distinction: a holomorphic passive coordinate can exist even when the push-forward requires integration over positive-dimensional fibres. Tangency alone does not establish convergence of that integral.

This note proves that the natural **unweighted integral over the full complex fibre** already diverges for the smallest complex Hermitian Hopf pencil, on a single admissible plane wave at its polarization centre. Thus the remaining issue is not just proving Gaussian estimates for that integral: a different reduction, measure, or current is needed. The Weyl sheet is the ordinary light cone, whose smoothing estimates are already available. No counterexample to local smoothing or to the original wave current is asserted.

## 1. The realified Weyl pencil

On $\mathbb R^4$, in coordinates corresponding to $(x_1+\mathrm i x_2,x_3+\mathrm i x_4)\in\mathbb C^2$, take
$$
C^0=I_4,\quad
C^1=\begin{pmatrix}0&0&1&0\\0&0&0&1\\1&0&0&0\\0&1&0&0\end{pmatrix},\quad
C^2=\begin{pmatrix}0&0&0&1\\0&0&-1&0\\0&-1&0&0\\1&0&0&0\end{pmatrix},\quad
C^3=\operatorname{diag}(1,1,-1,-1).
$$
They satisfy $C^iC^j+C^jC^i=2\delta_{ij}I$. The top eigenvalue of $C(\xi)$ is $|\xi|$, of real multiplicity two. Here $N=4$, $n=3$, $m=2$, so $N=n+m-1$: Claude's tangency criterion is met.

For auxiliary real states $w_1,w_2\in\mathbb R^4$, let $z=w_1+\mathrm k w_2$. All complex coordinates below use $\mathrm k$; the distinct Fourier unit is $\mathrm i$. Define
$$
Q_\mu(z)=\tfrac12z^{\mathsf T}C^\mu z,
\qquad
q^\mu(z)=\tfrac12\sum_{\ell=1}^2w_\ell^{\mathsf T}C^\mu w_\ell.
$$
**The bilinear map $Q$ and the real displacement $q$ are different.** In particular, $q$ need not be constant along a fibre of $Q$. This distinction is essential for the argument below.

Set
$$
u=(z_1+\mathrm k z_2,z_3+\mathrm k z_4),\qquad
v=(z_1-\mathrm k z_2,z_3-\mathrm k z_4).
$$
Then
$$
Q_0=\tfrac12(u_1v_1+u_2v_2),\quad
Q_3=\tfrac12(u_1v_1-u_2v_2),\quad
Q_1=\tfrac12(u_1v_2+u_2v_1),\quad
Q_2=\frac{u_2v_1-u_1v_2}{2\mathrm k}.
\tag{18.1}
$$
Consequently $Q_0^2=Q_1^2+Q_2^2+Q_3^2$, and $Q$ is invariant under
$$
(u,v)\longmapsto(\lambda u,\lambda^{-1}v),\qquad\lambda\in\mathbb C_{\mathrm k}^{\times}.
\tag{18.2}
$$
This is the complexification of the phase rotation of a complex state. Its noncompact radial part will cause the divergence.

## 2. The fibre measure is multiplicative Haar measure

Add the passive variable $\eta\in\mathbb C_{\mathrm k}$ and use the proposed single-lift map
$$
\Theta(z,\eta)=(Q_0,Q_1,Q_2,Q_3+\eta).
\tag{18.3}
$$
Physical spacetime coordinates are carried along unchanged. Near the north-point target $\omega_*=(1/2,0,0,1/2)$, choose the branch with $\eta$ near zero.

**Lemma 2.1 (explicit disintegration).** Put
$a=u_1v_1$, $b=u_1v_2$, $c=u_2v_1$, and $\lambda=u_1$. For $a\ne0$,
$$
u_1=\lambda,\quad u_2=c\lambda/a,\quad v_1=a/\lambda,\quad v_2=b/\lambda,
\qquad u_2v_2=bc/a.
\tag{18.4}
$$
At $(a,b,c,\eta)=(1,0,0,0)$, disintegrating Lebesgue measure $dV_z\,d^2\eta$ with respect to target Lebesgue measure under $\Theta$ gives on this fibre exactly
$$
\frac{d^2\lambda}{|\lambda|^2}.
\tag{18.5}
$$

*Proof.* The complex Jacobian of (18.4), from $(\lambda,a,b,c)$ to $(u_1,u_2,v_1,v_2)$, has determinant $1/(a\lambda)$, up to an irrelevant sign. The change from $z$ to $(u_1,u_2,v_1,v_2)$ has determinant of absolute value four. Thus
$$
dV_z\,d^2\eta=\frac{d^2\lambda\,d^2a\,d^2b\,d^2c\,d^2\eta}{16|a|^2|\lambda|^2}.
$$
By (18.1), the map from $(a,b,c,\eta)$ to $\omega$ has complex Jacobian determinant of absolute value $1/4$ at $(1,0,0,0)$, hence real Jacobian $1/16$. Divide by that Jacobian to obtain (18.5). This calculation uses ordinary change of variables, so it includes the coarea weight rather than only surface area on the fibre. ∎

The entire branch of the fibre through that point is
$$
z_1=\tfrac12(\lambda+\lambda^{-1}),\quad
z_2=\tfrac1{2\mathrm k}(\lambda-\lambda^{-1}),\quad
z_3=z_4=0,\quad\eta=0.
\tag{18.6}
$$
It is a copy of $\mathbb C^\times$, not a compact circle.

## 3. An admissible plane wave has infinite pushed time density

Take $\kappa>0$ and the top-sheet plane wave
$$
f(T,X)=e_1\exp(\mathrm i\kappa(X_3-T)).
$$
It solves $(\partial_T+\sum C^j\partial_{X_j})f=0$. In the lifted-current construction of [Claude 08, Lemma 3.1](../../claude/kakeya/08_GAPLESS_PROOF.md),
$G=\overline{f(\mathsf X-\mathrm i q)}$. Since
$$
q^0-q^3=|z_3|^2+|z_4|^2,
\qquad |G|^2=e^{-2\kappa(q^0-q^3)},
\tag{18.7}
$$
the Gaussian damping is identically one on (18.6). This is proved from $q$, not from the invariance of $Q$.

The physical-time entry of the original real positive tensor is
$$
B^{TT}=\sum_{\ell=1}^2|(C^0w_\ell)\cdot G|^2
=|z_1|^2e^{-2\kappa(q^0-q^3)}.
\tag{18.8}
$$
Because $\Theta$ leaves physical time unchanged, its pushed time entry is the integral of (18.8) with the measure in Lemma 2.1.

**Theorem 3.1 (divergence of the full-fibre single lift).** At $\omega_*$, the unweighted full-fibre push-forward of this positive current has infinite time density. More precisely, the contribution from $R^{-1}\le|\lambda|\le R$, for $R\ge1$, is
$$
\int_{R^{-1}\le|\lambda|\le R}
\left|\frac{\lambda+\lambda^{-1}}2\right|^2
\frac{d^2\lambda}{|\lambda|^2}
=\frac\pi2(R^2-R^{-2}).
\tag{18.9}
$$

*Proof.* Write $\lambda=r e^{\mathrm k\theta}$. The angular integral of the numerator is
$(\pi/2)(r^2+r^{-2})$; its cross term integrates to zero. Haar measure is $dr\,d\theta/r$. Integrating $(\pi/2)(r+r^{-3})$ from $1/R$ to $R$ gives (18.9). It tends to infinity. All contributions are nonnegative, so additional fibre branches cannot cancel the divergence. ∎

This excludes the natural full-fibre integral from the smooth finite engine data of `KAK.ENGINE`, even before checking (A3)–(A4). The example uses one plane wave in an arbitrarily small cap, exactly at the centre where a Gaussian estimate should be strongest.

**Remark 3.2 (what happens off the centre).** In coordinates (18.4), the same wave has
$$
|G|^2=\exp\left[-\kappa\left(\frac{|c|^2}{|a|^2}r^2+\frac{|b|^2}{r^2}\right)\right].
$$
When both $b,c\ne0$, this damps both ends of the radial integral. The damping disappears at the centre. It is therefore invalid to infer a uniform smooth Gaussian row estimate from generic-fibre convergence.

## 4. What a different reduction must supply

**Proposition 4.1 (cutoff identity).** If a symmetric tensor $B$ has divergence-free columns and $w$ is a smooth scalar weight, then
$$
\operatorname{div}(wB)=B\nabla w.
$$
Thus weighting or truncating the state fibres does not automatically preserve conservation.

*Proof.* Apply the product rule in each column. ∎

For example, at $w_1=e_1,w_2=0$, the weight $e^{-\epsilon(|w_1|^2+|w_2|^2)}$ gives a nonzero $w_{11}$-column divergence $-2\epsilon e^{-\epsilon}$ for the lifted plane wave above. This is an ambient calculation; a proposed push-forward with a weight still needs its own flux calculation. The note does not assert that every possible regularization or quotient is impossible.

A compact phase quotient, a gauge section, a modified current, or a compensating boundary flux could be useful. Each requires a proof of positivity, conservation, and quantitative row estimates after the change. Dividing an infinite positive integral by a supposed gauge volume is not itself such a construction: the time integrand varies along the gauge orbit.

The reusable lesson is precise: **tangency is algebraic; a finite positive push-forward also requires control of the noncompact symmetry directions of each fibre.** Here those directions are invisible to the plane-wave damping at the polarization centre.

## 5. Status and checks

This resolves one proposed implementation of `T.KAK.7` negatively, not the whole target. No classification of all Hermitian engines or improved logarithmic estimate is claimed. The original equal-state-dimension construction has discrete local fibres and is unaffected.

[`scripts/verify_weyl_fibre.py`](scripts/verify_weyl_fibre.py) checks the real matrices, Clifford relations, the bilinear identities, both Jacobians, the real displacement on the fibre, and the coefficients of the exact annular integral. It also records the Gaussian-cutoff source term. All algebra is rational or Gaussian-rational; no external packages are required. Results are in [`results/weyl_fibre.json`](results/weyl_fibre.json). The divergence follows from the displayed exact integral, not from numerical growth alone.

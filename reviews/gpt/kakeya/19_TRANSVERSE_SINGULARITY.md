# 19. The Weyl fibre current is not locally integrable

GPT / Astra, 1 November 2026. **Complete prose derivation for the current specified below; exact algebraic and finite integral checks. Not Lean-formalized.**

[18](18_NONCOMPACT_FIBRE_OBSTRUCTION.md) showed that the full complex-fibre time density is infinite at a plane wave's polarization centre. An infinite value on a set of measure zero would not by itself rule out a locally integrable current. This note closes that loophole for the same unweighted construction: on the real target slice, its nearby time density behaves as an inverse square of the transverse distance. The transverse dimension is two, so the integral diverges logarithmically in every neighbourhood of the centre.

The conclusion concerns the proposed Hermitian single lift. It neither refutes the original wave current nor rules out every gauge reduction or altered positive current.

## 1. The full coarea factor near the centre

Retain all definitions of [18, §§1–3](18_NONCOMPACT_FIBRE_OBSTRUCTION.md): $z\in\mathbb C_{\mathrm k}^4$, the realified Weyl matrices, $u,v$, and

$$
\Theta(z,\eta)=(Q_0,Q_1,Q_2,Q_3+\eta).
$$

Set $a=u_1v_1$, $b=u_1v_2$, $c=u_2v_1$, $\lambda=u_1$, and $d=bc/a$. Work in the chart $a\ne0$, $a-d\ne0$ and on one local branch of the target map. Then

$$
(u_1,u_2,v_1,v_2)=(\lambda,c\lambda/a,a/\lambda,b/\lambda),
$$

$$
\omega_0=(a+d)/2,\quad \omega_1=(b+c)/2,\quad
\omega_2=(c-b)/(2\mathrm k),\quad \omega_3=(a-d)/2+\eta.
$$

**Lemma 1.1 (nearby fibre measure).** Disintegration of $dV_z\,d^2\eta$ with respect to complex target Lebesgue measure gives the conditional measure

$$
\boxed{\quad \frac1{|a-d|^2}\frac{d^2\lambda}{|\lambda|^2}.\quad}
\tag{19.1}
$$

*Proof.* The first change of variables, calculated in [18], gives the factor $1/(16|a|^2|\lambda|^2)$. The complex Jacobian of $(a,b,c,\eta)\mapsto\omega$ has determinant

$$
\frac{1-bc/a^2}{4\mathrm k}
$$

up to sign. Its real Jacobian is $|1-bc/a^2|^2/16$. Divide by this factor and use $a(1-bc/a^2)=a-d$. ∎

Equation (19.1) specifies the fibre density before taking real target values. Below, “the density on the real slice” means the resulting explicit fibre integral evaluated at $\omega\in\mathbb R^4$. It does not mean that an eight-dimensional measure is restricted to a four-dimensional set by ordinary set restriction.

## 2. Exact reduction of the nearby time density

Write

$$
s^2=\omega_1^2+\omega_2^2,
\qquad \omega_0>s,\qquad
a=\omega_0+\sqrt{\omega_0^2-s^2},\quad
d=\omega_0-\sqrt{\omega_0^2-s^2}.
\tag{19.2}
$$

Choose the north branch. For real $\omega$,

$$
b=\omega_1-\mathrm k\omega_2,\quad c=\omega_1+\mathrm k\omega_2,
\quad |b|=|c|=s,\quad ad=s^2,
\quad \eta=\omega_3-\sqrt{\omega_0^2-s^2}.
$$

The plane wave is still $f(T,X)=e_1\exp(\mathrm i\kappa(X_3-T))$, $\kappa>0$. As computed from the real displacement in [18], its lifted time entry on this fibre is

$$
\left|\frac{\lambda+a/\lambda}{2}\right|^2
\exp\left[-\kappa s^2\left(\frac{|\lambda|^2}{a^2}+\frac1{|\lambda|^2}\right)\right].
\tag{19.3}
$$

Let $F_\kappa(\omega)$ be the integral of (19.3) against (19.1) on this one entire $\mathbb C^\times$ branch. The full positive push-forward is at least this contribution.

**Theorem 2.1 (exact profile).** For $s>0$ in the region (19.2), the branch integral is finite and equals

$$
\boxed{\quad F_\kappa(\omega)
=\frac{\pi a^2}{2\kappa s^2(a-d)^2}
J\left(\frac{\kappa^2s^4}{a^2}\right),\qquad
J(v)=\int_0^\infty e^{-u-v/u}\,du.\quad}
\tag{19.4}
$$

*Proof.* Set $\lambda=re^{\mathrm k\theta}$. The cross term in the squared modulus integrates to zero, leaving

$$
F_\kappa=\frac{\pi}{2(a-d)^2}
\int_0^\infty(r^2+a^2r^{-2})
e^{-Ar^2-B/r^2}\frac{dr}{r},
\quad A=\frac{\kappa s^2}{a^2},\quad B=\kappa s^2.
$$

Both $A,B$ are positive, so both ends converge. In the $r^2$ term substitute $u=Ar^2$; its integral is $J(AB)/(2A)$. In the $a^2r^{-2}$ term substitute $u=B/r^2$; its integral is $a^2J(AB)/(2B)$. These are equal. Multiplication by $\pi/[2(a-d)^2]$ gives (19.4). ∎

This derivation needs no special-function asymptotics. The elementary integral above directly proves every estimate used here.

## 3. A positive current cannot be recovered merely by ignoring the axis

**Lemma 3.1 (elementary profile bounds).** For $v\ge0$,

$$
e^{-2\sqrt v}\le J(v)\le1.
\tag{19.7}
$$

In particular $J(v)\to1$ as $v\downarrow0$, and $J(v)\ge e^{-3}$ for $0\le v\le1$.

*Proof.* The upper bound is the integral of $e^{-u}$. For the lower bound integrate over $u\ge\sqrt v$, where $v/u\le\sqrt v$. The resulting integral is $e^{-\sqrt v}\int_{\sqrt v}^\infty e^{-u}\,du=e^{-2\sqrt v}$. Dominated convergence also gives the limit directly. ∎

**Theorem 3.2 (inverse-square singularity and failure of local finiteness).** Uniformly for $\omega_0$ in a compact subinterval of $(0,\infty)$ and for arbitrary $\omega_3$,

$$
\lim_{s\downarrow0}s^2F_\kappa(\omega)=\frac{\pi}{2\kappa}.
\tag{19.5}
$$

In particular, near $\omega_*=(1/2,0,0,1/2)$, the density on the real slice is not locally integrable with respect to $d\omega_0\,d\omega_1\,d\omega_2\,d\omega_3$. The full unweighted positive current therefore does not define a locally finite positive measure with this time density.

*Proof.* In (19.2), $a/(a-d)\to1$ and $\kappa^2s^4/a^2\to0$, uniformly when $\omega_0$ stays bounded away from zero. Equations (19.4) and Lemma 3.1 give (19.5).

For an explicit lower bound, let $\omega_0\ge m>0$, and choose $s$ small enough that $s<m/2$ and $\kappa s^2\le m$. Then $a\ge m$, the argument of $J$ is at most one, and $a/(a-d)\ge1$. Thus

$$
F_\kappa(\omega)\ge\frac{\pi e^{-3}}{2\kappa s^2}.
\tag{19.6}
$$

Any neighbourhood of $\omega_*$ contains a product of nontrivial intervals in $\omega_0,\omega_3$ and a small transverse disk. Polar integration in $(\omega_1,\omega_2)$ gives $\int_0^\delta s^{-2}s\,ds=\infty$. Positivity prevents cancellation by additional branches. ∎

**Corollary 3.3 (universal logarithmic coefficient).** For fixed $\omega_0\ge m>0$, any $\omega_3$, and $0<\epsilon<\delta\le m/2$,

$$
\int_{\epsilon<\sqrt{\omega_1^2+\omega_2^2}<\delta}
F_\kappa(\omega)\,d\omega_1\,d\omega_2
=\frac{\pi^2}{\kappa}\log\frac\delta\epsilon+O_{m,\kappa}(\delta^2).
\tag{19.8}
$$

The error is uniform in $\epsilon,\omega_0,\omega_3$ in the stated range.

*Proof.* Since $d=s^2/a$, $a\ge m$, and $a-d\ge\sqrt3m$,

$$
1\le\frac{a^2}{(a-d)^2}
=\left(1+\frac d{a-d}\right)^2
\le1+2s^2/m^2\qquad(s\le m/2).
$$

Use (19.7) and $e^{-x}\ge1-x$ in (19.4) to obtain

$$
\frac{\pi}{2\kappa s^2}-\frac\pi m
\le F_\kappa(\omega)
\le\frac{\pi}{2\kappa s^2}+\frac{\pi}{\kappa m^2}.
$$

Integration against $2\pi s\,ds$ proves (19.8). ∎

Changing the density only on the axis $s=0$ cannot repair the divergence. A distributional finite part would not be the same positive measure: a positive distribution must be a locally finite measure.

## 4. What this adds, and what it leaves open

The obstruction is now stronger than failure of a smooth pointwise estimate. The naive real-slice current has infinite mass on neighbourhoods of the polarization axis. An alternative construction must change the fibre integral or the current, rather than merely overlook exceptional target values.

This example uses a single plane wave, as in the Gaussian testing of the proposed engine. It does not by itself prove a statement for every finite-energy wave packet. Nor does it rule out a gauge section, a compact quotient with a different density, or a compensating current. Those alternatives must establish positivity and conservation afresh; the cutoff source in [18, Proposition 4.1](18_NONCOMPACT_FIBRE_OBSTRUCTION.md#4-what-a-different-reduction-must-supply) remains relevant.

The structural mechanism is that a noncompact fibre symmetry has two radial ends, while the damping coefficient vanishes quadratically in transverse target displacement. The radial integral turns that vanishing into an inverse-square singularity. It is a reusable diagnostic for proposed positive push-forwards: calculate the transverse order of damping and compare the resulting singularity with transverse dimension.

Run `python reviews/gpt/kakeya/scripts/verify_weyl_singularity.py`. The script checks the general coarea determinant using Gaussian rationals, the real north-branch relations on rational samples, the two radial substitution factors, the exact prefactor limit, and rigorous elementary bounds for the one-dimensional profile. Its output is [weyl_singularity.json](results/weyl_singularity.json). Numerical profile enclosures illustrate convergence; the proof of (19.5) is dominated convergence, not a fitted numerical asymptotic.

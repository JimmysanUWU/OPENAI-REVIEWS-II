# 20. The gauge-orbit obstruction: every Hermitian single lift diverges

*Claude, 1 November 2026. This note verifies GPT's [18](../../gpt/kakeya/18_NONCOMPACT_FIBRE_OBSTRUCTION.md) and [19](../../gpt/kakeya/19_TRANSVERSE_SINGULARITY.md) and extends them. Roadmap nodes `KAK.WEYL.PROFILE` and `KAK.GAUGE.DIV`; target `T.KAK.7`.*

[17](17_STATE_DIMENSION_RIGIDITY.md) §6 restated `T.KAK.7` for the Hermitian pencils that pass the tangency test, and flagged the open point: "whether the fibre integrals converge … is part of the problem: the fibres of $\Theta$ are complexified spheres, non-compact in their imaginary directions." GPT settled that point for the Weyl pencil:
- **18:** the unweighted full-fibre time density is infinite at the polarization centre;
- **19:** on real target points near the centre it behaves like $\pi/(2\kappa s^2)$.

I checked both notes line by line and agree with every computation (§0). This note adds two things.

> **Proposition 20.1 (the complex profile).** Near the centre, at every complex target point, the Weyl time density of the north branch is
> $$F_\kappa=\frac{\pi|a|^2}{4\kappa|a-d|^2}\,J\!\Big(\frac{\kappa^2|b|^2|c|^2}{|a|^2}\Big)\Big(\frac1{|b|^2}+\frac1{|c|^2}\Big),\qquad J(v)=\int_0^\infty e^{-u-v/u}\,du .$$
> So $F_\kappa=+\infty$ on the two complex hypersurfaces $\{b=0\}$ and $\{c=0\}$, and $F_\kappa$ is not locally integrable on the complex target. This closes a step that 19 left implicit (§1).

> **Theorem 20.2 (gauge-orbit divergence).** Let the pencil be a realified complex or quaternionic Hermitian pencil. Then:
> - its complexified unitary gauge group $K_{\mathbb C}$ ($\mathbb C^\times$, or $SL(2,\mathbb C)$) preserves the flux map $Q$ and Lebesgue measure;
> - the $K_{\mathbb C}$-orbit of any real extremal eigenvector stays inside the extremal eigenspace, where the plane wave with that polarization has no Gaussian damping.
>
> Hence the unweighted full-fibre single lift has **infinite time density at every polarization centre**. In the tangency class $N'=n+m-1$ the fibre is locally this orbit, and over $|\operatorname{Im}\theta|\le T$ it contributes exactly $\rho\,\pi\sinh 2T$ in the complex case (§2).

> **Update (1 Nov, later).** GPT's [22](../../gpt/kakeya/22_CURVATURE_AND_QUADRIC_FIBRES.md) Theorem 22.4 proves the divergence of Theorem 20.2 for **every** pencil in the tangency class with $m>1$, without a gauge group. That covers the octonionic case too. Its mechanism: the curvature identity $R_aS^{-1}R_a^{\mathsf T}=\frac12(a\cdot a)H$ makes the coarea density $\gamma/|a_1|^2$, which cancels the time entry. Theorem 20.2 below is its gauge-group special case; Proposition 20.1 and §3 are not covered there. Verified in [`../GPT_ROUND2_REVIEW.md`](../GPT_ROUND2_REVIEW.md).

**What this changes.** The fibre-integrated route to `T.KAK.7` is closed for the whole Hermitian class, not just for Weyl. The cause is structural: the bilinear flux map is invariant under the complexification of a compact gauge group, and the plane wave sees none of the imaginary gauge directions. Any Hermitian engine must quotient or fix that gauge, and then prove positivity and conservation afresh (§3).

**Checks.** [`scripts/verify_gauge_orbit.py`](scripts/verify_gauge_orbit.py) writes [`results/gauge_orbit.json`](results/gauge_orbit.json) and takes about 10 s. Its methods are independent of GPT's scripts: finite-difference real Jacobians, quadrature, and random pencils.

| Part | What it checks |
|---|---|
| W | Weyl:<br>• the chart reproduces $Q$<br>• the coarea density $1/(|a-d|^2|\lambda|^2)$ by $10\times10$ and $8\times8$ finite-difference real Jacobians<br>• Proposition 20.1 against 2D quadrature at random complex $(a,b,c)$, to $10^{-14}$ relative<br>• the reduction to 19 (19.4) on real targets<br>• $|b|^2F_\kappa\to\pi/(4\kappa)$ as $b\to0$<br>• growth with the cutoff when $c=0$ |
| H | Complex Hermitian pencils $m'=2,3,4$ ($n=2m'-1$; $m'=2$ is GPT's Weyl form):<br>• $J$ commutes with the pencil<br>• $\det g_\theta=1$, $Q$ invariant, the orbit inside $E(\xi_0)_{\mathbb C}$, no damping<br>• the time entry along the orbit<br>• the coarea density constant along the orbit, and equal to 1 for Weyl, matching 18's $d^2\lambda/|\lambda|^2$<br>• orbit integral $\rho\,\pi\sinh 2T$ |
| Q | Quaternionic Hermitian pencils $2\times2$, $3\times3$ ($n=5,9$):<br>• the right multiplications $I_a$ commute with the pencil<br>• $\exp(\sum\theta_aI_a)$ preserves $Q$, has determinant 1 and keeps the orbit undamped<br>• $\langle gz_0,z_0\rangle=\cos s$<br>• constant coarea density on the 3-dimensional orbit<br>• $\mathbb E_u\langle uz_0,z_0\rangle^2=\frac14$ |
| N | Controls:<br>• a random antisymmetric "gauge" neither commutes nor preserves $Q$<br>• off-orbit points are damped<br>• the arrow pencil has no fibre directions |

---

## 0. Verdict on GPT 18–19

| Claim | Verdict | How checked |
|---|---|---|
| 18 §1: realified Weyl matrices, Clifford relations, the $(u,v)$ form (18.1) of $Q$, light-cone identity, $\mathbb C^\times$ invariance (18.2) | agree | by hand; Part W |
| 18 Lemma 2.1: coarea density $d^2\lambda/|\lambda|^2$ at $(1,0,0,0)$ | agree | Jacobians $1/(a\lambda)$ and $1/(4\mathrm k)$ by hand; finite differences in Part W |
| 18 (18.7)–(18.8): $q^0-q^3=|z_3|^2+|z_4|^2$, time entry $|z_1|^2e^{-2\kappa(q^0-q^3)}$ | agree | by hand; Part W |
| 18 Theorem 3.1: $\frac\pi2(R^2-R^{-2})$ | agree | by hand. In the orbit coordinate $\lambda=e^{\mathrm k\theta}$ this is $\pi\sinh 2T$ with $R=e^T$ (Theorem 20.2) |
| 18 Proposition 4.1 (cutoffs create sources) | agree | product rule |
| 19 Lemma 1.1: general coarea factor $1/(|a-d|^2|\lambda|^2)$ | agree | Jacobian $(1-bc/a^2)/(4\mathrm k)$ by hand; Part W |
| 19 Theorem 2.1, Lemma 3.1, Theorem 3.2, Corollary 3.3 | agree | by hand, including the constants $\sqrt3m$, $1+2s^2/m^2$, $e^{-3}$; Part W checks (19.4) |
| 19's conclusion "the full current does not define a locally finite positive measure" | **agree, but the proof needs one more step** | 19 proves non-integrability of the density on the **real** slice $\omega\in\mathbb R^4$. The current lives on the complex target $\mathbb C_{\mathrm k}^4$. Non-integrability on a slice does not in general pass to the ambient space (Remark 1.3). Proposition 20.1 proves it directly |

Neither note claims more than it proves about Hermitian engines in general. Both say explicitly that other reductions remain open. Theorem 20.2 shows that the naive route fails for every Hermitian pencil, not just for Weyl.

---

## 1. The complex profile of the Weyl density

Use GPT 18's setting:
- the realified Weyl pencil and the map $\Theta(z,\eta)=(Q_0,Q_1,Q_2,Q_3+\eta)$;
- the chart $(u_1,u_2,v_1,v_2)=(\lambda,\,c\lambda/a,\,a/\lambda,\,b/\lambda)$, with $d=bc/a$, on the north branch near $(a,b,c,\eta)=(1,0,0,0)$;
- the plane wave $f=e_1e^{\mathrm i\kappa(X_3-T)}$.

The target $\omega=\big(\frac{a+d}2,\frac{b+c}2,\frac{c-b}{2\mathrm k},\frac{a-d}2+\eta\big)$ is now an arbitrary complex point.

**Proposition 20.1.** Let $a\ne0$, $a\ne d$ and $b,c\ne0$. The north-branch fibre integral of the time entry is the $F_\kappa$ displayed above.

*Proof.*
1. **Fibre measure.** By 19 Lemma 1.1, the fibre measure is $d^2\lambda/(|a-d|^2|\lambda|^2)$.
2. **Damping.** $q^0-q^3=|z_3|^2+|z_4|^2=\frac12(|u_2|^2+|v_2|^2)$. Here $z_3=\frac{u_2+v_2}2$, $z_4=\frac{u_2-v_2}{2\mathrm k}$, and the parallelogram law holds in $\mathbb C_{\mathrm k}$. With $\lambda=re^{\mathrm k\vartheta}$ this gives
$$|G|^2=e^{-2\kappa(q^0-q^3)}=\exp\!\big(-Ar^2-B/r^2\big),\qquad A=\kappa|c|^2/|a|^2,\quad B=\kappa|b|^2 .$$
3. **Angular average.** $z_1=\frac12(\lambda+a/\lambda)$, and the cross term of $|z_1|^2$ averages to zero over $\vartheta$. So
$$\int_0^{2\pi}|z_1|^2\,d\vartheta=\frac\pi2\big(r^2+|a|^2r^{-2}\big).$$
4. **Radial integral.** Since $d^2\lambda/|\lambda|^2=d\vartheta\,dr/r$,
$$F_\kappa=\frac{\pi}{2|a-d|^2}\int_0^\infty\big(r^2+|a|^2r^{-2}\big)e^{-Ar^2-B/r^2}\frac{dr}r .$$
The substitutions $u=Ar^2$ and $u=B/r^2$ give $J(AB)/(2A)$ and $|a|^2J(AB)/(2B)$, as in 19 Theorem 2.1.
5. **Assemble.** $AB=\kappa^2|b|^2|c|^2/|a|^2$, which gives the formula. ∎

On real targets $b=\bar c$ in $\mathbb C_{\mathrm k}$, $|b|=|c|=s$ and $a$ is real, so the formula reduces to 19 (19.4). Part W checks this.

**Corollary 20.1a (singular hypersurfaces).**
1. If exactly one of $b$, $c$ vanishes, one end of the radial integral is undamped and $F_\kappa=+\infty$.
2. Near $(1,0,0,0)$, $F_\kappa\ge\frac{\pi}{8\kappa}\big(|b|^{-2}+|c|^{-2}\big)$.
3. So $\int F_\kappa\,d^8\omega=\infty$ over every neighbourhood of the centre. The pushed time density is not locally integrable on the complex target.

*Proof.*
1. With $c=0$, $A=0$, and the $r^2$ term integrates $r\,dr$ to infinity. The case $b=0$ is symmetric, through the term $|a|^2r^{-3}\,dr$ at $r\to0$.
2. Near the centre $|a|^2/|a-d|^2\ge\frac12$ and $J\ge e^{-2\sqrt v}\ge\frac12$, by 19 Lemma 3.1.
3. The map $(a,b,c,\eta)\mapsto\omega$ is a local biholomorphism with real Jacobian $\approx\frac1{16}$. And $\int_{|b|<\delta}d^2b/|b|^2=2\pi\int_0^\delta dr/r=\infty$. ∎

**Remark 1.3 (why 19's slice argument needs this).** A function on $\mathbb R^2_s\times\mathbb R^4_y$ such as $(|s|^2+|y|^2)^{-1}$ is not integrable on the slice $y=0$, yet it is locally integrable on $\mathbb R^6$. So the non-integrability of the real-slice density in 19 does not by itself prove that the 8-dimensional density fails to be locally integrable. Proposition 20.1 proves it directly: the singular set is a complex hypersurface, of real codimension 2, and the singularity is $|b|^{-2}$ in a transverse complex coordinate. 19's **conclusion** is right; it needed this.

---

## 2. The mechanism: complexified gauge orbits

**Setting.** Let $C^0=I_N,C^1,\dots,C^n$ be real symmetric. Suppose a compact connected group $K\subset O(N)$ commutes with every $C^\mu$ and acts freely and transitively on the unit sphere of each extremal eigenspace $E(\xi)$, $\xi\in\Gamma$. Let $\mathfrak k\subset\mathfrak{so}(N)$ be its Lie algebra and $K_{\mathbb C}=\exp(\mathfrak k\otimes\mathbb C_{\mathrm k})$, acting $\mathbb C_{\mathrm k}$-linearly on $\mathbb C_{\mathrm k}^N$. The two cases are:

| Pencil | $N$ | $m$ | $K$ | $\mathfrak k$ | $K_{\mathbb C}$ |
|---|---|---|---|---|---|
| complex Hermitian $m'\times m'$, realified | $2m'$ | 2 | $U(1)=\{e^{\alpha J}\}$, $J$ = multiplication by $i$ | $\mathbb RJ$ | $\mathbb C^\times=\{\cos\theta+J\sin\theta:\theta\in\mathbb C_{\mathrm k}\}$ |
| quaternionic Hermitian $m''\times m''$, realified | $4m''$ | 4 | $Sp(1)$, right multiplication by unit quaternions | $\operatorname{span}(I_1,I_2,I_3)$ | $SL(2,\mathbb C)$ |

In each case the extremal eigenspace is a complex or quaternionic line $z_0\cdot\mathbb C$ or $z_0\cdot\mathbb H$. Hence $m=\dim_{\mathbb R}K+1$.

**Theorem 20.2.** In this setting, use the single lift of 17 Definition 17.1 with auxiliary space $U=V$. Fix $\xi_0\in\Gamma$, $z_0\in E(\xi_0)$ real with $|z_0|=1$, $\kappa>0$, and the plane wave $f=z_0e^{\mathrm i\kappa(X\cdot\xi_0-T\varphi(\xi_0))}$.
1. **(Symmetry.)** Every $g\in K_{\mathbb C}$ satisfies $Q(gz)=Q(z)$ and $\det g=1$. So $g$ preserves Lebesgue measure on $\mathbb C_{\mathrm k}^N$ and commutes with $\Theta$. It does **not** preserve the real displacement $q$; only $K$ does.
2. **(No damping on the orbit.)** For $z\in K_{\mathbb C}z_0$, the real and imaginary parts $w_1,w_2$ lie in $E(\xi_0)$. Then $|G|^2\equiv1$ there, and the time entry is $|\langle gz_0,z_0\rangle|^2$, with $\langle\cdot,\cdot\rangle$ the $\mathbb C_{\mathrm k}$-bilinear extension of the dot product.
3. **(Fibre = orbit, Haar measure.)** Suppose the pencil is in the tangency class $N'=n+m-1$ (lift with $U=V$; $V$ is $K$-invariant because each $E(\xi)$ is), so that 17 Theorem 17.2 applies. Then near the centre, the fibre of $\Theta$ through $(z_0,0)$ contains $K_{\mathbb C}z_0$ as an open subset, injectively parametrized by $K_{\mathbb C}$. On it, the fibre measure is $\rho\cdot dg$ for Haar measure $dg$ and a constant $\rho>0$.
4. **(Divergence.)** $\int_{K_{\mathbb C}}|\langle gz_0,z_0\rangle|^2\,dg=\infty$. In the complex case:
$$\rho\int_{|\operatorname{Im}\theta|\le T}|\cos\theta|^2\,d^2\theta=\rho\,\pi\sinh 2T .$$

So the time density of the unweighted full-fibre push-forward is $+\infty$ at every centre $\Theta(z_0,0)$. Other parts of the fibre only add nonnegative mass.

*Proof.*
1. **Symmetry.**
   - Let $A\in\mathfrak k\otimes\mathbb C_{\mathrm k}$ and $g=e^A$. $A$ is a $\mathbb C_{\mathrm k}$-combination of real antisymmetric matrices that commute with every $C^\mu$.
   - So $g^{\mathsf T}=e^{-A}=g^{-1}$, and $g^{\mathsf T}C^\mu g=g^{-1}C^\mu g=C^\mu$. Hence $Q(gz)=\frac12z^{\mathsf T}g^{\mathsf T}C^\mu gz=Q(z)$.
   - $\det g=e^{\operatorname{tr}A}=1$, and the real Jacobian of a $\mathbb C_{\mathrm k}$-linear map is $|\det|^2=1$.
2. **No damping.**
   - $E(\xi_0)$ is invariant under $\mathfrak k$, which commutes with $C(\xi_0)$. So $E(\xi_0)\otimes\mathbb C_{\mathrm k}$ is invariant under $K_{\mathbb C}$.
   - The real and imaginary parts of a vector in $E\otimes\mathbb C_{\mathrm k}$ lie in $E$.
   - By 17 Lemma 17.3(4), $\varphi q^0-\xi_0\cdot q'=\frac12\sum_\ell w_\ell^{\mathsf T}S w_\ell$ with $S=\varphi I-C(\xi_0)$. This vanishes, since $Sw_\ell=0$.
   - $G=\overline{f(\mathsf X-\mathrm iq)}=z_0\,e^{-\mathrm i\kappa(\cdots)}e^{-\kappa(\varphi q^0-\xi_0\cdot q')}$ has modulus one.
   - The time entry is $\sum_\ell|w_\ell\cdot G|^2=|w_1\cdot z_0|^2+|w_2\cdot z_0|^2=|\langle z,z_0\rangle|^2$.
3. **Fibre = orbit.**
   - By the proof of 17 Theorem 17.2, $\operatorname{rank}DQ_V=n$ near $z_0$ and the passive direction is transverse. So $\Theta$ is a submersion near $(z_0,0)$, with fibres of complex dimension $N'+1-(n+1)=m-1=\dim_{\mathbb C}K_{\mathbb C}$.
   - By part 1, $\Theta$ is $K_{\mathbb C}$-invariant, so it is a submersion along the whole orbit.
   - The orbit map $g\mapsto gz_0$ is injective. In the complex case $gz_0=z_0$ forces $\cos\theta=1,\sin\theta=0$, since $z_0,Jz_0$ are independent. In the quaternionic case $gz_0=z_0q=z_0$ with $q\in\mathbb H\otimes\mathbb C_{\mathrm k}$ forces $q=1$, because some quaternionic coordinate of $z_0$ is a nonzero real quaternion, hence invertible.
   - Its differential $A\mapsto Az_0$ is injective: $z_0\mathrm i,z_0\mathrm j,z_0\mathrm k$, respectively $Jz_0$, are real-independent, hence $\mathbb C_{\mathrm k}$-independent.
   - So the orbit is an injectively immersed submanifold of the same dimension as the fibre, hence open in it.
   - Finally, Lebesgue measure and $\Theta$ are $K_{\mathbb C}$-invariant. In coordinates $(g,y)\mapsto g\,\sigma(y)$ for a local transversal $\sigma$, left translation preserves both the volume and a left-invariant frame, so the conditional measure is $\rho(y)\,dg$. Part H checks that the density is constant along the orbit, and Part Q does so on the 3-dimensional orbit.
4. **Divergence.**
   - **Complex case.** $\langle g_\theta z_0,z_0\rangle=\cos\theta$ because $Jz_0\perp z_0$. Haar measure in $\theta$ is $d^2\theta$. With $\theta=\alpha+\mathrm kt$,
$$|\cos\theta|^2=\cos^2\alpha\cosh^2t+\sin^2\alpha\sinh^2t,$$
     whose $\alpha$-integral is $\pi\cosh2t$. Integrating over $|t|\le T$ gives $\pi\sinh 2T$.
     For Weyl, $\lambda=e^{\mathrm k\theta}$, $\rho=1$ (Part H), and $R=e^T$ recovers 18 (18.9): $\frac\pi2(R^2-R^{-2})=\pi\sinh2T$.
   - **Quaternionic case.**
     - On $E(\xi_0)\otimes\mathbb C_{\mathrm k}$, $g$ acts as right multiplication by some $q\in\mathbb H\otimes\mathbb C_{\mathrm k}$ of norm one, and $\langle z_0q,z_0\rangle=\operatorname{Re}_{\mathbb H}(q)$.
     - Use the Cartan decomposition $SL(2,\mathbb C)=Sp(1)\cdot\exp(\mathrm k\,\operatorname{Im}\mathbb H)$: $q=u\,e^{\mathrm kY}$ with $u$ a real unit quaternion and $Y$ a real pure quaternion.
     - Haar measure is $du\,j(Y)\,dY$ with $j(Y)=\big(\sinh 2|Y|/2|Y|\big)^2\ge1$. This is the volume growth of $SL(2,\mathbb C)/SU(2)=\mathbb H^3$, up to normalization.
     - Since $e^{\mathrm kY}=\cosh|Y|+\mathrm k\sinh|Y|\,\hat Y$,
$$|\operatorname{Re}_{\mathbb H}(ue^{\mathrm kY})|^2=\cosh^2|Y|\,(\operatorname{Re}u)^2+\sinh^2|Y|\,\operatorname{Re}(u\hat Y)^2\ge\cosh^2|Y|\,(\operatorname{Re}u)^2 .$$
     - The $u$-average of $(\operatorname{Re}u)^2$ is $\frac14$ (Part Q). So the integral is at least $\frac14\operatorname{vol}(Sp(1))\int\cosh^2|Y|\,dY=\infty$. ∎

**Remarks.**
1. **Where the obstruction sits.**
   - Real $K$-orbits are compact, and the integrand is bounded on them by $|z|^2$.
   - The divergence comes entirely from the imaginary directions $\exp(\mathrm k\mathfrak k)$. They move $(w_1,w_2)$ within $E(\xi_0)$, so they are invisible to the damping, and they grow $|z|$.
   - In GPT's coordinates these are the radial directions $|\lambda|\to0,\infty$.
2. **Why the tangency class has $N'=n+m-1$.**
   - $V_{\mathbb C}/\!\!/K_{\mathbb C}$ has dimension $N'-(m-1)$. Tangency ($N'=n+m-1$) says exactly that this quotient has dimension $n$, the dimension of the image hypersurface of $Q$.
   - So the flux map descends to a local biholomorphism on the gauge quotient. This is the $N=n$ situation of the original wave current, but with the gauge integrated out, which is where the divergence enters.
   - For complex pencils the $\mathbb C^\times$-invariants are the bilinears $u_av_b$, so the quotient is the rank-one determinantal variety of dimension $2m'-1=n$.
3. **What does not hold.** Theorem 20.2 is a statement about this push-forward. It does not show that local smoothing fails for these pencils, nor that every positive current fails. For isotropic members such as Weyl the sheet is the round cone, where the estimate is already known.

---

## 3. What this means for `T.KAK.7`

**The route that is now closed.** The tangency class of 17 was supposed to carry engine data by integrating the lifted current over the fibres of $\Theta$. By Theorem 20.2:
- that integral is infinite at every polarization centre, for every complex and quaternionic Hermitian pencil;
- (A4) needs the time row exactly there;
- by Proposition 20.1, for Weyl it is not even locally finite nearby.

So no estimate of the form (A3)–(A4) can hold for this current.

**What a repair must do.** It must remove the imaginary gauge directions $\exp(\mathrm k\mathfrak k)$ while keeping a positive conserved current. Two natural attempts fail as they stand:
- **Dividing by the gauge volume** fails because the integrand is not gauge-invariant (GPT 18 §4).
- **The Kempf–Ness slice** $\{\mu=0\}$ fails too. Here $\mu$ is the moment map of the real $K$-action; for $U(1)$, $\mu(w_1,w_2)=w_1^{\mathsf T}Jw_2$. The slice meets each closed $K_{\mathbb C}$-orbit in one compact $K$-orbit; on the central orbit, exactly the circle $t=0$. But the lifted current of 08 Lemma 3.1 is not tangent to it. At $w_1=e_\alpha$, $w_2=0$ (in $E(\xi_0)$), $G=z_0$:
  - $\mathcal Z_2(\mu)=-z_0^{\mathsf T}Je_\alpha=\sin\alpha$;
  - so $B\nabla\mu=\sin\alpha\,z_0\ne0$ in the $w_2$-block.

  By GPT 18 Proposition 4.1, restricting to the slice creates a source.

A Hermitian engine therefore needs a current that is built on the gauge quotient from the start, or a flux correction that compensates the slice. That is the refined `T.KAK.7`.

---

## 4. Status

| Result | Status | Uses |
|---|---|---|
| §0 verdicts on GPT 18–19 | checked by hand and by Part W | GPT 18 Lemma 2.1, 19 Lemma 1.1 |
| Proposition 20.1, Corollary 20.1a | proved here; Part W (quadrature to $10^{-14}$) | 19 Lemma 1.1, Lemma 3.1 |
| Theorem 20.2 | proved here; Parts H, Q, N | 17 Theorem 17.2, Lemma 17.3(4); Cartan decomposition of $SL(2,\mathbb C)$ (standard) |
| §3 Kempf–Ness remark | one-line computation | 08 Lemma 3.1 |

**Subjective confidence:**
- ≈97% for Proposition 20.1 and the complex case of Theorem 20.2;
- ≈92% for the quaternionic case, which cites the Haar-measure form of the Cartan decomposition.

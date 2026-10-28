# 11. The positivity engine as a theorem

*Claude, 28 October 2026. Roadmap target `T.KAK.1`; nodes `KAK.ENGINE`, `KAK.PENCIL`, `KAK.APP.E`, `KAK.ARROW` in [`ROADMAP.md`](../../../ROADMAP.md).*

The paper's critical estimate (Theorem 1.2) and its App. E reuse one mechanism: a positive conserved current, Monge–Ampère potentials on designed fibres, and propagation in a two-dimensional auxiliary time. [07](07_TOOLKIT.md) K7 described this mechanism as a blueprint. This file turns it into a theorem.

> **Theorem 11.1 (engine theorem).** Suppose a cap carries *engine data* (Definition 1.4):
> - a holomorphic fibre chart with the scales of (G1)–(G2);
> - a positive, conserved current tangent to a passive coordinate, whose rows are Fourier multipliers with Gaussian symbols (A1)–(A3);
> - an invertible time row (A4).
>
> Then, for some aperture $\gamma>0$ and all $0<h\le1$, every $f$ on $\mathbb T^n$ with Fourier support in the $\gamma$-cap at scale $h$ satisfies
> $$\int_1^2\|e^{-\mathrm iT\varphi(D)}f\|_{L^{p}(\mathbb T^n)}^{p}\,dT\le C\,[\log(2/h)]^2\,\|f\|_{L^{p}(\mathbb T^n)}^{p},\qquad p=p_c=\frac{2n}{n-1}.$$

> **Theorem 11.2 (pencils carry engine data).** Let $C_1,\dots,C_n$ be real symmetric $n\times n$ matrices, and let $\varphi(\xi)=\lambda_{\max}(\sum_j\xi_jC_j)$. Let $\omega\in S^{n-1}$ be a direction where this eigenvalue is simple and $D^2\varphi(\omega)>0$ on $\omega^\perp$. Then a cap around $\omega$ carries engine data. The current is App. E's push-forward of the lifted current.

> **Corollary 11.3 (paper Theorem E.1, complete).** Suppose $\lambda_{\max}$ is simple with positive definite transverse Hessian on an open cone $\Gamma$. Then for every $b\in C_c^\infty(\Gamma)$,
> $$\|e^{-\mathrm iT\varphi(D)}b(hD)f\|_{L^{p_c}(\mathbb R^n\times[1,2])}\le C_b[\log(2/h)]^{2/p_c}\|f\|_{L^{p_c}(\mathbb R^n)}.$$
> The same holds for a simple $\lambda_{\min}$ with negative definite transverse Hessian.

> **Corollary 11.4 (the wave is the arrow pencil).** The pencil
> $$A(\xi)=\begin{pmatrix}\xi_n&\xi'^{\mathsf T}\\\xi'&-\xi_nI_{n-1}\end{pmatrix}$$
> has $\lambda_{\max}=|\xi|$, simple off the ray $\{\xi'=0,\ \xi_n<0\}$. Its polarization is the stereographic projection $v(\xi)=\xi'/(|\xi|+\xi_n)$, and its group-velocity map is the inverse stereographic projection. Theorem 1.2 is Corollary 11.3 for two rotated copies of $A$. The paper's wave system (08 §3.1) is $A$ after a change of frame and a rescaling of the states.

**What is new here.**
- The axioms are stated precisely.
- An abstract proof is written, line by line in terms of the axioms.
- App. E's quantitative chart and current claims, which [10](10_GAPLESS_REMAINING.md) §5 had checked only at the level of exponents and identities, are now proved. App. E's statement therefore moves from PARTIAL to PROVED.
- The identification of the wave with the arrow pencil is new. It gives a second derivation of the cap estimate behind Theorem 1.2, with a spatial (not null) passive coordinate.
- A hypothesis-usage matrix (§7) explains what each hypothesis does and what breaks without it.

**Checks.** [`scripts/verify_positivity_engine.py`](scripts/verify_positivity_engine.py) ([`results/positivity_engine.json`](results/positivity_engine.json)), in 60-digit arithmetic:
- **Part G: pencil geometry**, for the arrow pencils ($n=2,3,4$), the cubic pencil of Example E.2, and a perturbed arrow pencil.
  - Hellmann–Feynman, $Dv=S_\perp^{-1}L$, $dB(0)=2L^{\mathsf T}$, stationarity and the characteristic conormal hold to $10^{-40}$. The curvature formula, checked by second differences, holds to $10^{-25}$.
  - The cubic pencil's transverse Hessian is $\mathrm{diag}(1,\frac23)$.
- **Part C: the current**, for the same five pencils, computed directly in the target coordinates.
  - Conservation holds to $10^{-38}$ relative; tangency, $J=|R|^{2n}j(W)$, the Gaussian identity and the row identity all hold to $10^{-39}$ or better.
  - Two negative controls (wrong conjugation; no volume factor) fail at order $10^{-4}$ to $1$.
- **Part M: a middle eigenvalue.** The current is still conserved, but its "window" $e^{-w^*Sw/2}$ exceeds $1$.
- **Part W: the arrow pencil, exactly (sympy), for $n=2,3,4$.** The paper's null-coordinate system is the arrow pencil, and the characteristic polynomial factors as $(t+\xi_n)^{n-2}(t^2-|\xi|^2)$. The stereographic formulas for $v$ and $B$ hold, and $B(v(\xi))=\xi/|\xi|$.

---

## 0. Conventions

All conventions are those of [08](08_GAPLESS_PROOF.md) §0–§2.
- **Dimensions.** $n\ge2$, $d=n-1$, $p=p_c=2n/d$ and $q=n/d$.
- **Slots.** A slot has type $\mathcal D=\operatorname{span}_{\mathbb R}(1,\mathrm i,\mathrm j)$. Its coordinate is $\mathsf X+\mathsf P\mathrm i$ with $\mathsf X\in\mathbb R$ and $\mathsf P\in\mathbb C_{\mathrm k}$; the real coordinates $(\mathsf X,\operatorname{Re}\mathsf P,\operatorname{Im}_{\mathrm k}\mathsf P)$ carry the units $(1,\mathrm i,\mathrm j)$.
- **Calculus.** $G$, $H$, the Moore determinant and cofactor $C(\cdot)$, and conserved Hermitian fields are as in 08 §2. "Holomorphic" means holomorphic in the $\mathbb C_{\mathrm k}$-structure of the imaginary parts, with Wirtinger derivatives $\partial$ as in 08 Lemma 2.1(4).
- **The $n+1$ slots.** Slot $0$ is $T+\tau\mathrm i$ with $\tau\in\mathbb C_{\mathrm k}$. Slots $1,\dots,n$ are $X_j+P_j\mathrm i$. The parameter is $z=(T,\tau)\in\mathbb R^3$, and $s=|\tau|$.
- **Torus and multipliers.** $\mathbb T^n=\mathbb R^n/(100\mathbb Z)^n$. For a symbol $m$, $m(hD)g=\sum_\xi m(h\xi)\hat g(\xi)e^{\mathrm iX\cdot\xi}$. The kernels $J_\ell$ (with $M=2n+4$) and the cutoff $g(Q)$ are those of 08 §6.1.
- **Scales.** $\sigma=s/h$.

## 1. The axioms

**Definition 1.1 (cap).** A *cap* $(\varphi,K',v)$ consists of:
- a compact $K'\subset\mathbb R^n\setminus0$;
- a real $\varphi$, smooth and homogeneous of degree one on a conic neighbourhood of $K'$;
- a smooth map $v:K'\to\mathbb R^d$, the *polarization*.

For a compact $K\subset\operatorname{int}K'$ and $\gamma>0$, put $K_\gamma=\{\eta\in K:|v(\eta)|\le\gamma\}$. A finite Fourier sum $f$ on $\mathbb T^n$ is **$\gamma$-cap data at scale $h$** if $h\xi\in K_\gamma$ for every frequency $\xi$ of $f$. Then $g_T=e^{-\mathrm iT\varphi(D)}f$.

**Definition 1.2 (engine geometry).** Fix constants $w_0\in(0,\frac1{10}]$ and $C_0\ge1$, and let $0<a\le\frac14$. An engine geometry consists of holomorphic functions $W:\mathcal V\to\mathbb C_{\mathrm k}^d$ and $\Phi:\mathcal V\to\mathbb C_{\mathrm k}$, independent of $(T,X)$. Here $\mathcal V\subset\{a<|\tau|<4\}\times\mathbb C_{\mathrm k}^n$ is open. Put $Q=\Phi/a$. For every $\tau$, with $\mathcal V_\tau$ the slice:
- **(G1) Chart.** $P\mapsto(W,Q)$ maps $\mathcal V_\tau$ biholomorphically onto an open set containing $\bar E=\{|W|^2/w_0^2+100|Q|^2\le1\}$. It depends smoothly on $\tau$.
- **(G2) Scales.** Let $\mathcal W=\{(\tau,P)\in\mathcal V:(W,Q)\in\bar E\}$. On $\mathcal W$:
  - $|\partial W|\le C_0/s$ and $|\partial\Phi|\le C_0$, with $\partial$ the holomorphic differential in all of $(\tau,P)$;
  - $|\partial_P\Phi|\ge1/C_0$ (the **background** direction);
  - the inverse chart satisfies $|\partial P/\partial W|\le C_0s$, $|\partial P/\partial Q|\le C_0a$ and $|\det_{\mathbb C}\partial P/\partial(W,Q)|\ge as^d/C_0$.

Put $\mathrm r=|W|^2/w_0^2+100|Q|^2-1$. The fibres are $\Omega_z=\mathbb T^n\times\{P\in\mathcal V_\tau:\mathrm r<0\}$, and the total domain is $\Omega=\{(z,\zeta):a<|\tau|<4,\ \zeta\in\Omega_z\}$.

**Definition 1.3 (engine current).** Fix constants $c_0>0$, $C_1,C_2\ge1$, $m_1\ge0$ and $L\in\mathbb N$, and a cap with $|v|\le w_0/3$ on $K'$. Let $0<h\le a$. An engine current assigns to every $\gamma$-cap data $f$ at scale $h$ a smooth quaternionic column $\chi=(\chi_0,\dots,\chi_n)$ on $\mathbb R_T\times\mathbb T^n\times\mathcal W'$, where $\mathcal W'\supset\mathcal W$ is open. It must satisfy:
- **(A1) Positive conserved current.** $A^0=\chi\chi^*$ is conserved. It is positive semidefinite by construction.
- **(A2) Tangency.** $\sum_\mu(\partial_\mu\Phi)\,\chi_\mu=0$, with left multiplication by $\partial_\mu\Phi\in\mathbb C_{\mathrm k}$.
- **(A3) Gaussian rows.** At every $(\tau,P)\in\mathcal W$ and for every $\mu$ there are real $m_{\mu\ell}=m^{\tau,P}_{\mu\ell}\in C^\infty(K')$, $\ell\le L$, independent of $(T,X,f)$, with
$$|\chi_\mu(T,\tau,X,P)|^2=\sum_\ell|m_{\mu\ell}(hD)g_T(X)|^2 .$$
  For $|\beta|\le2n+4$, $\eta\in K'$ and $W=W(\tau,P)$ they satisfy
$$|\partial_\eta^\beta m_{\mu\ell}(\eta)|\le C_1(hs)^{-d/2}\sigma^{m_1}\exp\big(-c_0\sigma|W-v(\eta)|^2\big).$$
- **(A4) Invertible time row.** At every $(\tau,P)\in\mathcal W$ there is a real $m_0\in C^\infty(K')$ with $|\chi_0|^2\ge|m_0(hD)g_T|^2$. Moreover, for every $\gamma\in(0,w_0/8]$: if $|W(\tau,P)|\le2\gamma$, then on $\{\eta\in K':|v(\eta)|\le2\gamma\}$ we have $m_0\ne0$ and
$$|\partial_\eta^\beta(1/m_0)|\le C_1(hs)^{d/2}\sigma^{m_1}\exp(C_2\gamma^2\sigma)\qquad(|\beta|\le2n+4).$$
  The constant $C_2$ does not depend on $\gamma$.

**Definition 1.4 (engine).** An *engine on the cap $(\varphi,K',v)$* is a choice, for all $0<h\le a\le\frac14$, of an engine geometry and an engine current. The constants $(w_0,C_0,c_0,C_1,C_2,m_1,L)$ must not depend on $h$, $a$ or $\gamma$.

**Remarks.**
1. **What is not assumed.** The axioms never mention the equation. They do not assume homogeneity of $\varphi$ beyond the cap, a specific chart, or a specific current. The equation enters only through (A1)–(A4).
2. **The Gaussian is centred at the polarization.** On the wall of the fibre $|W|\ge\frac{\sqrt3}2w_0$, while the data have $|v|\le w_0/3$. So (A3) forces exponential decay $e^{-c_0'\sigma}$ there, with $c_0'=c_0w_0^2/4$. This is the "Gaussian decay on the wall" of App. E.
3. **Two numbers set the aperture.** The decay rate $c_0$ and the inversion cost $C_2$ compete, and (4.0) below chooses $\gamma$ so that the decay wins.

---

## 2. Fibre geometry from the chart

**Lemma 2.1 (fibre geometry; 08 Lemma 4.1 under (G1)–(G2)).** Let $F=(W/w_0,10Q)$, with $n\times(n+1)$ complex Jacobian $J=\partial F$ and spatial block $J_{\rm fib}=\partial_PF$. Constants below depend only on $n,w_0,C_0$.
1. $H_{\rm fib}(\mathrm r)=4J_{\rm fib}^*J_{\rm fib}>0$. The full $H(\mathrm r)=4J^*J\ge0$ has a one-dimensional quaternionic kernel $N$, and $N\subset G(\mathrm r)^\perp\cap\ker H(|\Phi|^2)$.
2. On $\mathcal W$:
   - $H_{\rm fib}(\mathrm r)\ge cs^{-2}I$ and $G_{\rm fib}(\mathrm r)G_{\rm fib}(\mathrm r)^*\le H_{\rm fib}(\mathrm r)$;
   - $dP\asymp a^2s^{2d}\,dW\,dQ$ (real volume elements) and $\det H_{\rm fib}(\mathrm r)\asymp a^{-2}s^{-2d}$.
3. $H_{\rm fib}(|\Phi|^2)=4(\partial_P\Phi)^*\partial_P\Phi$ has rank one, with eigenvalue in $[4C_0^{-2},4C_0^2]$. Also $|G(|W|^2)|\le Cs^{-1}$ on $\mathcal W$.
4. On each fibre, $\int\delta(\mathrm r)\,dP\le Ca^2s^{2d}$ and $\int_{\mathrm r<0}dP\le Ca^2s^{2d}$.

*Proof.* $F$ is holomorphic in the imaginary parts and independent of $(T,X)$. So 08 Lemma 2.1(4) gives $H(\mathrm r)=H(|F|^2)=4J^*J$ and $G(\mathrm r)=2J^*F\mathrm i$.
1. **Kernel.**
   - $J_{\rm fib}=\partial_P(W/w_0,10Q)$ is invertible by (G1), so $H_{\rm fib}(\mathrm r)>0$ and $J$ has rank $n$. Let $\kappa$ span its complex kernel.
   - Write a quaternionic column as $x=x_1+x_2\mathrm j$ with $x_1,x_2\in\mathbb C_{\mathrm k}^{n+1}$. Then $Jx=Jx_1+(Jx_2)\mathrm j$. So $\ker J=\kappa\mathbb H=:N$, which is $\ker H(\mathrm r)$.
   - $G(\mathrm r)$ lies in the image of $J^*$, hence is orthogonal to $N$.
   - $\partial\Phi=\frac a{10}\times$ (last row of $J$). Hence $(\partial\Phi)\kappa=0$, and $H(|\Phi|^2)=4(\partial\Phi)^*\partial\Phi$ vanishes on $N$.
2. **Bounds.**
   - $J_{\rm fib}^{-1}=\partial P/\partial F=[w_0\,\partial P/\partial W,\ \frac1{10}\partial P/\partial Q]$ has norm $\le C_0(w_0s+a/10)\le Cs$. Hence $H_{\rm fib}(\mathrm r)\ge4\|J_{\rm fib}^{-1}\|^{-2}\ge cs^{-2}$.
   - $G_{\rm fib}G_{\rm fib}^*=4J_{\rm fib}^*FF^*J_{\rm fib}\le4|F|^2J_{\rm fib}^*J_{\rm fib}\le H_{\rm fib}(\mathrm r)$, since $|F|\le1$ on $\bar E$.
   - The real Jacobian of a holomorphic map is $|\det_{\mathbb C}|^2$. By (G2) and Hadamard, $as^d/C_0\le|\det_{\mathbb C}\partial P/\partial(W,Q)|\le C(as^d)$.
   - $\det H_{\rm fib}(\mathrm r)=4^n|\det_{\mathbb C}J_{\rm fib}|^2$ (08 Lemma 4.1's proof), and $|\det J_{\rm fib}|=10w_0^{-d}|\det\partial P/\partial(W,Q)|^{-1}$.
3. The first claim is 08 Lemma 2.1(4) with (G2). For the second, $G(|W|^2)=2\sum_j(\partial W_j)^*W_j\mathrm i$, with $|\partial W_j|\le C_0/s$ and $|W_j|\le w_0$.
4. Pull back to $(W,Q)$. There $\mathrm r$ is a fixed quadratic with $|\nabla\mathrm r|\asymp1$ on $\partial E$, and the Jacobian is $\le C(as^d)^2$. ∎

---

## 3. The regularized current, the potential and the energy

Fix the exponents:
- $\Lambda_0=2d+2m_1+n$;
- $\lambda=h^{-\Lambda_0-1}$;
- $\Lambda_r=50n$.

As in 08 (6.1), put
$$A_{\rm w}=g(Q)\,(J_h*_XA^0),\qquad A=A_{\rm w}+h^{\Lambda_r}C(H(\mathrm r)),\qquad\mu=A_{TT},\qquad B_0=\lambda|\Phi|^2,\qquad B=H(B_0).$$

**Lemma 3.1 (the data; 08 Lemma 6.1 under the axioms).** Let $\|f\|_2=1$, and let $h\le h_1$ with $h_1$ depending on the constants.
1. $A\ge0$ is conserved on $\Omega$, and $\mu\ge h^{\Lambda_r}\det H_{\rm fib}(\mathrm r)>0$.
2. $A:B=0$ and $\det B=0$. $B_{\rm fib}=\lambda H_{\rm fib}(|\Phi|^2)$ has rank one, with eigenvalue $\beta\in[4C_0^{-2}\lambda,4C_0^2\lambda]$.
3. $\nu^*C(H(\mathrm r))\nu=0$ on the wall, where $\nu=G(\mathrm r)/|\nabla\mathrm r|$.
4. $\mu\le Ch^{-\Lambda_0}$, hence $\mu^{1/d}\le\lambda$.

*Proof.*
1. **Conservation.**
   - $X$-convolution commutes with the constant-coefficient divergence, and $\Phi$ does not depend on $X$. So $J_h*_XA^0$ is conserved by (A1) and annihilates $(\partial\Phi)^*$ by (A2): indeed $A^0(\partial\Phi)^*=\chi\big((\partial\Phi)\chi\big)^*=0$.
   - $g(Q)$ is a function of the holomorphic $\Phi$, so $G(g(Q))=2(\partial\Phi)^*(\partial_{\bar\Phi}g)\mathrm i$ (08 Lemma 2.1(4)). The divergence of $g\cdot M$ is $g\operatorname{div}M+G(g)^*M$, so $g(Q)(J_h*A^0)$ is conserved.
   - $C(H(\mathrm r))\ge0$ is conserved (08 Lemmas 2.2–2.3). Its $TT$ entry is $\det H_{\rm fib}(\mathrm r)$ (08 Lemma 2.2(3)).
2. **The background.**
   - $H(|\Phi|^2)=4(\partial\Phi)^*\partial\Phi$ has rank one, so $\det B=0$.
   - $A_{\rm w}(\partial\Phi)^*=0$ by tangency.
   - The image of $C(H(\mathrm r))$ is $N\subset\ker H(|\Phi|^2)$ (Lemma 2.1.1). Hence $A:B=0$.
   - The value of $\beta$ is Lemma 2.1.3.
3. The image of $C(H(\mathrm r))$ is $N\perp G(\mathrm r)$.
4. **Size of $\mu$.**
   - By (A3) with $\beta=0$, $|\chi|^2\le LC_1^2(hs)^{-d}\sigma^{2m_1}\big(\sum_\xi|\hat f(\xi)|\big)^2$.
   - At most $Ch^{-n}$ lattice frequencies have $h\xi\in K$, so $\sum|\hat f|\le Ch^{-n/2}$.
   - With $s\ge h$ and $\sigma\le4/h$, this gives $\|A_{\rm w}\|\le Ch^{-\Lambda_0}$.
   - The regularizer adds $\le Ch^{\Lambda_r-2-2d}\le C$ (Lemma 2.1.2).
   
   Finally $\mu^{1/d}\le(Ch^{-\Lambda_0})^{1/d}\le h^{-\Lambda_0-1}$ for small $h$. ∎

**The potential.** Lemma 2.1 and Lemma 3.1 verify the hypotheses of 08 Theorem 5.2 on each fibre, with:
- $m=n$, the lattice of $\mathbb T^n$, and $V_z$ the chart image of the open ellipsoid;
- the defining function $\mathrm r$, the background $B_0$, and $g=\lambda\mu>0=\det B_{\rm fib}$;
- smooth dependence on $z$, since the chart depends holomorphically on $\tau$ and $A$ smoothly on $(T,\tau)$.

So there is $u$, smooth on $\bar\Omega$, with
$$D_u:=H_{\rm fib}(B_0+u)>0,\qquad\det D_u=\lambda\mu,\qquad u=0\text{ on the wall},\qquad u\le0 .$$
Put $\rho=u/\mathrm r\ge0$.

**Proposition 3.2 (energy balance and comparison).** Let
$$X_0(z)=\int_{\Omega_z}(-u)\mu,\qquad\mathcal F=X_0+\frac1\lambda\int_{\Omega_z}u\int_0^1\det D_{\theta u}\,d\theta .$$
Then
$$\text{(a)}\ \ \Delta_z\mathcal F\le\mathcal B(z):=\int\rho\,G(\mathrm r)^*A_{\rm w}G(\mathrm r)\,\delta(\mathrm r)\,d\zeta,\qquad\text{(b)}\ \ \tfrac{n-1}{n+1}X_0\le\mathcal F\le\tfrac n{n+1}X_0 .$$

*Proof.* The proofs of 08 Proposition 6.2 and Lemma 6.3 use only:
- the conservation of $A$ and $\mu>0$ (Lemma 3.1.1);
- $A:B=0$, $\det B=0$, and $B_{\rm fib}$ of rank one (3.1.2);
- the vanishing of $\nu^*C(H(\mathrm r))\nu$ on the wall (3.1.3);
- the common null line $N\subset\ker B\cap G(\mathrm r)^\perp$ of $B$ and $H(\mathrm r)$ (Lemma 2.1.1);
- the potential;
- the cofactor divergence, Green's identity and the Schur split (08 Lemmas 2.2–2.4).

All of these hold here, so those proofs apply verbatim. ∎

---

## 4. Estimates on one fibre

Two regions of a fibre matter:
- **the wall region:** $\mathrm r=0$ and $|Q|\le\frac1{20}$. There $|W|\ge\frac{\sqrt3}2w_0$, so $|W-v(\eta)|\ge w_0/2$ for $\eta\in K'$.
- **the centre region:** $|W|\le\gamma$ and $|Q|\le\frac1{400}$.

Fix $\psi\in C_c^\infty(\operatorname{int}K')$ with $\psi=1$ on $K$. Every symbol below is multiplied by $\psi$, which does not change its action on cap data.

**Proposition 4.1 (fixed fibre; 08 Proposition 7.1 under the axioms).** Let $c_0'=c_0w_0^2/4$ and fix
$$\gamma=\min\Big(\frac{w_0}8,\ \Big(\frac{c_0'}{4qC_2}\Big)^{1/2}\Big).\tag{4.0}$$
There are $K_0$, $h_0$, $c$ and $C$, depending only on $n$, the engine constants, $K$ and $K'$, with the following property. Let $a=K_0h$ with $0<h<h_0$, let $f$ be $\gamma$-cap data with $\|f\|_2=1$, and let $T\in\mathbb R$. Then
$$\mathcal F(0,\tau)\le Ca^2\|f\|_p^p\ \ (|\tau|=4a),\qquad\mathcal F(T,\tau)\ge ca^2\|g_T\|_p^p\ \ (|\tau|=2a),\qquad\Delta_z\mathcal F\le Cs^{-2}e^{-c_0's/h}\mathcal F\ \ (2a\le s\le3).$$

*Proof.* Each step names the axioms it uses.

**Step 1 (spatial average).** For a constant $L_*\ge1$ fixed in Step 3, put $P_*=(J_{L_*s}*|g_T|^2)^{1/d}$. As in 08 (7.2),
$$|\nabla_X^jP_*|\le C(L_*s)^{-j}P_*\quad(j=1,2),\qquad\int P_*^{d+1}\,dX\le\|g_T\|_p^p,\qquad P_*^d\ge c\,h^{n+4} .$$
The last bound is the nearest lattice translate, with $\int|g_T|^2=1$.

**Step 2 (upper bounds for the current; (A3), cap).** At a fixed point $(\tau,P)$, apply 08 Lemma 7.2 to each $\psi m_{\mu\ell}(hD)$. By (A3), its derivative bound is $\mathsf B=C(hs)^{-d/2}\sigma^{m_1}$, times $e^{-c_0'\sigma}$ on the wall. Summing the squares and using $J_h*J_h\le CJ_h$ and $J_h\le(L_*s/h)^nJ_{L_*s}$:
$$\|A_{\rm w}\|\le Cs^{-2d}\sigma^{m}P_*^d\ (\mathrm r\le0),\qquad\|A_{\rm w}\|\le Cs^{-2d}\sigma^me^{-2c_0'\sigma}P_*^d\ (\text{wall}),\qquad\mu\le Cs^{-2d}\sigma^mP_*^d\ (\mathrm r\le0),\tag{4.1}$$
with $m=2m_1+n+d$. Here $(hs)^{-d}(s/h)^n=s^{-2d}\sigma^{n+d}$. The regularizer obeys $h^{\Lambda_r}\det H_{\rm fib}(\mathrm r)\le Ch^{\Lambda_r}a^{-2}s^{-2d}\le s^{-2d}P_*^d$ by Step 1, since $\Lambda_r\ge n+7$.

**Step 3 (upper comparison; Lemma 2.1, Lemma 3.1).** Let $R=H_{\rm fib}(\mathrm r)$, so $R\ge cs^{-2}I$ and $G_{\rm fib}G_{\rm fib}^*\le R$. For quaternionic $e$, the product rule and Step 1 give
$$e^*H_{\rm fib}(\mathrm rP_*)e\ge(1-C/L_*)P_*\,e^*Re .$$
- The cross terms are $\le C(e^*Re)^{1/2}(L_*s)^{-1}P_*|e|$, and $|e|\le Cs(e^*Re)^{1/2}$.
- The term $\mathrm rH(P_*)$ is $\le CL_*^{-2}P_*e^*Re$, since $|\mathrm r|\le1$.

Fix $L_*$ with $1-C/L_*\ge\frac12$. Let $v_*=C_3\sigma^{m_3}\mathrm rP_*$. Then $H_{\rm fib}(v_*)\ge\kappa I$ with $\kappa=cC_3\sigma^{m_3}P_*s^{-2}$. Since $B_{\rm fib}\ge0$ has rank one with eigenvalue $\beta\ge c\lambda$,
$$\det(B_{\rm fib}+H_{\rm fib}(v_*))\ge(\beta+\kappa)\kappa^d\ge c\lambda\kappa^d\ge\lambda\mu ,$$
by (4.1), once $C_3$ and $m_3$ are large. Both potentials vanish on the wall, so 08 Lemma 5.1 gives $v_*\le u\le0$. Hence
$$0\le\rho=u/\mathrm r\le C_3\sigma^{m_3}P_*\qquad(\mathrm r\le0).\tag{4.2}$$
With Proposition 3.2(b), (4.1), the volume bound of Lemma 2.1.4 and Step 1:
$$\mathcal F(T,\tau)\le X_0\le Ca^2\sigma^{m_4}\|g_T\|_p^p .\tag{4.3}$$

**Step 4 (the flux; (A2), Lemma 2.1).**
- $G(\mathrm r)=w_0^{-2}G(|W|^2)+100a^{-2}G(|\Phi|^2)$, and $A_{\rm w}G(|\Phi|^2)=2A_{\rm w}(\partial\Phi)^*\Phi\mathrm i=0$ by tangency.
- So $G(\mathrm r)^*A_{\rm w}G(\mathrm r)\le Cs^{-2}\|A_{\rm w}\|$ by Lemma 2.1.3.

With the wall bound in (4.1), (4.2) and Lemma 2.1.4, Proposition 3.2(a) gives
$$\Delta_z\mathcal F\le\mathcal B\le C(a/s)^2\sigma^{m_5}e^{-2c_0'\sigma}\|g_T\|_p^p .\tag{4.4}$$

**Step 5 (recovering the wave at the centre; (A4), (G2), Lemma 3.1.4).** Take a centre $(X,P)$ of the centre region, and the ball of radius $h$ about it in the $3n$ real spatial coordinates.
1. **The ball stays inside.**
   - Along a segment in $P$ of length $\le h$, (G2) moves $W$ by $\le C_0h/s\le C_0/(2K_0)$ and $Q$ by $\le C_0/K_0$. Moving $X$ changes neither.
   - Take $K_0\ge\max(C_0/\gamma,200C_0)$. A continuity argument along the segment shows the ball lies in $\{|W|\le2\gamma,\ |Q|\le\frac1{100}\}$.
   - There $g(Q)=1$, and $\mathrm r\le4\gamma^2/w_0^2+\frac1{100}-1<0$.
2. **Inverting the time row.** Fix $\psi_\gamma\in C_c^\infty(\operatorname{int}K'\cap\{|v|<2\gamma\})$ with $\psi_\gamma=1$ on $K_\gamma$. At every auxiliary point of the ball, $|W|\le2\gamma$. So by (A4), $\psi_\gamma/m_0$ has derivatives through order $2n+4$ bounded by
$$\mathsf B_{\rm inv}=C_\gamma(hs)^{d/2}\sigma^{m_1}e^{C_2\gamma^2\sigma}.$$
   Since $g_T=(\psi_\gamma/m_0)(hD)\,m_0(hD)g_T$, 08 Lemma 7.2 and (A4) give
$$|g_T(X)|^2\le C\mathsf B_{\rm inv}^2\,(J_h*A^0_{TT}(\cdot,P'))(X)$$
   at every auxiliary point $P'$ of the ball. Moving $X$ by $\le h$ changes $J_h$ by a bounded factor. So on the whole ball
$$\mu\ge m_X:=c_\gamma(hs)^{-d}\sigma^{-2m_1}e^{-2C_2\gamma^2\sigma}|g_T(X)|^2 .\tag{4.5}$$
3. **Depth at the centre.**
   - Let $l$ be the distance to the centre, and $v_b=\epsilon m_X^{1/d}(l^2-h^2)$. Its fibre Hessian is $6\epsilon m_X^{1/d}I$.
   - Since $m_X^{1/d}\le\mu^{1/d}\le\lambda$ (Lemma 3.1.4), $\det(B_{\rm fib}+6\epsilon m_X^{1/d}I)\le(C\lambda+6\epsilon\lambda)(6\epsilon)^dm_X\le\lambda\mu$ for small $\epsilon=\epsilon(n,C_0)$.
   - $u\le0=v_b$ on the sphere, so 08 Lemma 5.1 on the ball gives $u\le v_b$.
   - Hence $(-u)\mu\ge\epsilon h^2m_X^{1+1/d}=\epsilon h^2m_X^q$ at the centre, and $m_X^q\propto|g_T(X)|^p$.
4. **Integration over the centres.**
   - By (G2) the centre region has volume $\ge ca^2s^{2d}\gamma^{2d}$.
   - $h^2s^{2d}(hs)^{-n}=\sigma^{n-2}\ge1$.
   
   With Proposition 3.2(b):
$$\mathcal F(T,\tau)\ge\tfrac{n-1}{n+1}X_0\ge c_\gamma a^2\sigma^{-2m_1q}e^{-2qC_2\gamma^2\sigma}\|g_T\|_p^p .\tag{4.6}$$

**Step 6 (constants and conclusion).**
- By (4.0), $2qC_2\gamma^2\le c_0'/2$. Divide (4.4) by (4.6):
$$\frac{\Delta_z\mathcal F}{\mathcal F}\le C_\gamma s^{-2}\sigma^{m_6}e^{-(2c_0'-c_0'/2)\sigma}\le Cs^{-2}e^{-c_0'\sigma}.$$
- On $|\tau|=4a$ we have $\sigma=4K_0$, and (4.3) at $T=0$ gives the first bound.
- On $|\tau|=2a$ we have $\sigma=2K_0$, and (4.6) gives the second.
- $h_0$ is chosen for Lemma 3.1 and Step 1. ∎

---

## 5. Proof of Theorem 11.1

By homogeneity, let $\|f\|_2=1$.

**Small $h$.** Let $h<h_0$, $a=K_0h$ and $b=2a<\frac14$.
- $\mathcal F$ is smooth on $a<|\tau|<4$ (Proposition 3.2) and nonnegative (3.2(b)).
- With $V(s)=Cs^{-2}e^{-c_0's/h}$ and $s=bt$, $\int_b^3s\log(s/b)V\le C\int_1^\infty\frac{\log t}te^{-2c_0'K_0t}dt$, uniformly in $h$. Also $V$ is bounded on $[\frac12,\frac32]$.
- 08 Lemma 8.1 with $r=2b$ gives
$$H(0,2b)\ge c\,\frac{\log2}{D^2}\int_1^2H(T,b)\,dT,\qquad D=\log(3/b)\le C\log(2/h).$$
- By Proposition 4.1, $H(0,2b)\le Ca^2\|f\|_p^p$ and $H(T,b)\ge ca^2\|g_T\|_p^p$. Cancel $a^2$.

**Large $h$.** For $h\ge h_0$, at most $C$ lattice frequencies have $h\xi\in K$. So $\|g_T\|_p\le\sum|\hat f(\xi)|\le C\|f\|_1\le C\|f\|_p$. ∎

**What the proof used.** The proof reached the equation only through (A1)–(A4). It used the geometry only through (G1)–(G2), and the auxiliary parameter only through its dimension: $z\in\mathbb R^3$, so $\tau\in\mathbb R^2$. §7 records exactly where each one enters.

---

## 6. Symmetric pencils

Throughout this section:
- $C_1,\dots,C_n$ are real symmetric $n\times n$ matrices;
- $C(\xi)=\sum\xi_jC_j$ and $\varphi=\lambda_{\max}\circ C$;
- the system is $\mathscr P=\partial_T+\sum_jC_j\partial_{X_j}$, that is $C^0=I$ and $C^j=C_j$;
- $\omega\in S^{n-1}$ is a direction where $\varphi$ is simple and $D^2\varphi(\omega)|_{\omega^\perp}>0$.

All products of $\mathbb C_{\mathrm k}$-vectors written with a dot are bilinear.

**Lemma 6.1 (polarization).** Let $e_0$ be the unit top eigenvector of $C(\omega)$, and use an orthonormal basis $(e_0,e_\perp)$ of the states. Put $L:\mathbb R^n\to e_0^\perp$, $Lw=(C(w)e_0)_\perp$. Then on a conic neighbourhood $\Gamma_\omega$ of $\omega$:
1. $\varphi$ is smooth (indeed real-analytic). There is a top eigenvector $E(\xi)=(1,v(\xi))$, with $v$ smooth, homogeneous of degree $0$, and $v(\omega)=0$.
2. $S(\xi):=\varphi(\xi)I-C(\xi)\ge0$ on all of $\mathbb R^n$, $S(\xi)E(\xi)=0$, and $c|\xi|I\le S_\perp(\xi)\le C|\xi|I$ on $\Gamma_\omega$. Here $S_\perp$ is the $e_0^\perp$ block.
3. $D^2\varphi(\omega)[w,w]=2(Lw)^{\mathsf T}S_\perp(\omega)^{-1}Lw$.
4. $\ker L=\mathbb R\omega$, $L(\omega^\perp)=e_0^\perp$, and $Dv(\omega)=S_\perp(\omega)^{-1}L$. So $v$ is a local angular coordinate.

*Proof.*
1. The eigenprojection $\Pi(\xi)=(2\pi\mathrm i)^{-1}\oint(\zeta-C(\xi))^{-1}d\zeta$ around the simple eigenvalue is real-analytic (Kato, Ch. II). Put $E=\Pi e_0/(e_0^{\mathsf T}\Pi e_0)$.
2. $\varphi=\lambda_{\max}$ gives $S\ge0$. At $\omega$, $S_\perp=(\varphi I-C(\omega))|_{e_0^\perp}$ is at least the spectral gap. Continuity and homogeneity extend this to a cone.
3. **Curvature.** Let $e(\sigma)$ be the unit top eigenvector of $C(\omega+\sigma w)$.
   - Hellmann–Feynman gives $\varphi'=e^{\mathsf T}C(w)e$, hence $\varphi''=2e'^{\mathsf T}C(w)e$.
   - Differentiating $(C-\varphi)e=0$ at $\sigma=0$, with $e'\perp e_0$, gives $S_\perp e'=Lw$.
4. **The map $L$.** $L\omega=(\varphi(\omega)e_0)_\perp=0$. By 3 and definiteness, $L$ is injective on $\omega^\perp$, and $\dim\omega^\perp=d=\dim e_0^\perp$. Finally $E'=e'$ at $\sigma=0$, so $Dv(\omega)w=S_\perp^{-1}Lw$. ∎

**Lemma 6.2 (group velocity and chart).** For $W\in\mathbb C_{\mathrm k}^d$ with $|W|<\frac12$, let $e(W)=(1,W)$ and
$$B_j(W)=\frac{e(W)^{\mathsf T}C_je(W)}{1+W\cdot W}.$$
1. $B$ is holomorphic, and $B(v(\xi))=\nabla\varphi(\xi)$ for $\xi\in\Gamma_\omega$.
2. $dB(0)=2L^{\mathsf T}$, which maps $e_0^\perp$ isomorphically onto $\omega^\perp$.
3. For real $\xi\in\Gamma_\omega$, the holomorphic function $W\mapsto\xi\cdot B(W)$ is stationary at $W=v(\xi)$.
4. There is $w_0\in(0,\frac1{10}]$ such that $\pi_\perp\circ B$ is a biholomorphism from $\{|W|<2w_0\}$ onto its image, and $|\det[dB(W),\omega]|\ge c>0$ there. Here $\pi_\perp$ is the orthogonal projection onto $\omega^\perp$.

*Proof.*
1. For real $W$, $\xi\cdot B(W)$ is the Rayleigh quotient of $C(\xi)$ at $e(W)$. Hellmann–Feynman gives $\partial_j\varphi=E^{\mathsf T}C_jE/E^{\mathsf T}E$ at $E=e(v)$.
2. At $W=0$, $d(e^{\mathsf T}C_je)=2(C_je_0)_\perp\cdot dW$ and $d(1+W\cdot W)=0$. The image is $(\ker L)^\perp=\omega^\perp$.
3. A Rayleigh quotient is stationary at eigenvectors. The complex derivative of a holomorphic function at a real point equals its real derivative.
4. By the inverse function theorem, $\pi_\perp dB(0)$ is invertible. Choose $w_0$ so that $\|\pi_\perp(dB(W)-dB(0))\|\le\frac12\|(\pi_\perp dB(0))^{-1}\|^{-1}$ on the convex ball $|W|<2w_0$. This gives injectivity. $[dB(0),\omega]=[2L^{\mathsf T},\omega]$ is invertible, and continuity gives the determinant bound. ∎

**Lemma 6.3 (the pencil chart satisfies (G1)–(G2)).** For $0<|\tau|<4$, $|W|<2w_0$ and $\Phi\in\mathbb C_{\mathrm k}$, put
$$P=\tau B(W)+\Phi\,\omega,\qquad R^2=\frac{2\tau}{1+W\cdot W},\qquad y=|R|^2 .$$
This is a biholomorphism onto its image $\mathcal V_\tau$, holomorphic in $\tau$ too, with inverse
$$W=(\pi_\perp B)^{-1}(\pi_\perp P/\tau),\qquad\Phi=\omega\cdot P-\tau\,\omega\cdot B(W).$$
For every $a\le\frac14$ it satisfies (G1)–(G2), with $C_0$ depending only on the pencil, $\omega$ and $w_0$. Also $y\asymp s$, and at $W=v(\xi)$ the characteristic conormal is
$$d\Phi=\frac{\xi\cdot dP-\varphi(\xi)\,d\tau}{\xi\cdot\omega}.$$

*Proof.*
- **Bijectivity.** This is Lemma 6.2.4. Given $P$, $\pi_\perp P/\tau=\pi_\perp B(W)$ determines $W$, and $\omega\cdot P$ then determines $\Phi$. The $(W,Q)$-image is $\{|W|<2w_0\}\times\mathbb C_{\mathrm k}\supset\bar E$.
- **Derivatives of the inverse chart.**
  - $\partial W/\partial P=(\pi_\perp dB)^{-1}\pi_\perp/\tau$ and $\partial W/\partial\tau=-(\pi_\perp dB)^{-1}\pi_\perp B(W)/\tau$. Both are $O(1/s)$.
  - $\partial_P\Phi=\omega^{\mathsf T}-\omega\cdot dB\,(\pi_\perp dB)^{-1}\pi_\perp$ and $\partial_\tau\Phi=-\omega\cdot B+\omega\cdot dB\,(\pi_\perp dB)^{-1}\pi_\perp B$. Both are bounded.
  - $(\partial_P\Phi)\,\omega=1$ because $\pi_\perp\omega=0$. So $|\partial_P\Phi|\ge1$.
- **Derivatives of the chart.** $\partial P/\partial W=\tau dB(W)=O(s)$ and $\partial P/\partial Q=a\omega$. Also $\det_{\mathbb C}\partial P/\partial(W,Q)=\tau^da\det[dB(W),\omega]$, which is $\ge cas^d$ in modulus.
- **Size of $y$.** $y=2s/|1+W\cdot W|\in[2s/(1+4w_0^2),2s/(1-4w_0^2)]$.
- **Conormal.** At $W=v(\xi)$, Lemma 6.2.3 kills the $dW$ term of $\xi\cdot dP=\xi\cdot B\,d\tau+\tau\,\xi\cdot dB\,dW+(\xi\cdot\omega)d\Phi$, and $\xi\cdot B(v(\xi))=\varphi(\xi)$. ∎

**Lemma 6.4 (the pencil current satisfies (A1)–(A4)).** Let $f$ be cap data. Put:
- the vector wave $f_{\rm sys}=\sum_\xi|\xi|^{d/2}E(\xi)u_\xi$, with $u_\xi=\hat f(\xi)e^{\mathrm i(X\cdot\xi-T\varphi(\xi))}$, a finite sum of exponential solutions of $\mathscr Pf_{\rm sys}=0$;
- $z=Re(W)\in\mathbb C_{\mathrm k}^n$, $q^\mu=\frac12z^*C^\mu z$ (real), and $G=\overline{f_{\rm sys}(T-\mathrm iq^0,X-\mathrm iq')}$;
- $M_{\mu a}=(C^\mu z)_a$ and $J=|\det_{\mathbb C}[M\,|\,(0,\omega)]|^2$;
- $\chi_\mu=J^{-1/2}\sum_aM_{\mu a}\,\mathrm i\,G_a$.

Then:
1. **(A1)** $A^0=\chi\chi^*$ is single-valued and conserved on $\{R\ne0,\ |W|<2w_0\}$.
2. **(A2)** $(\partial\Phi)\chi=0$.
3. **(Explicit form)** $J=y^nj(W)$ with $j(W)=|\det[(C^\mu e(W))_a\,|\,(0,\omega)]|^2\asymp1$. With
$$\Psi_\xi=\frac y2(W-v)^*S_\perp(\xi)(W-v)\qquad\text{and}\qquad c_\mu(W,\xi)=j(W)^{-1/2}e(W)^{\mathsf T}C^\mu E(\xi)\in\mathbb C_{\mathrm k},$$
   we have $\xi\cdot q'-\varphi q^0=-\Psi_\xi$, and
$$\chi_\mu=\frac R{|R|}\sum_\xi c_\mu(W,\xi)\,|\xi|^{d/2}y^{-d/2}e^{-\Psi_\xi}\,\mathrm i\,\overline{u_\xi}.$$
4. **(A3)** with $L=2$. The symbols are
$$m'_\mu,\ m''_\mu(\eta)=\big(\operatorname{Re}_{\mathrm k},\operatorname{Im}_{\mathrm k}\big)\,c_\mu(W,\eta)\cdot|\eta|^{d/2}(hy)^{-d/2}e^{-\Psi},\qquad\Psi=\frac y{2h}(W-v(\eta))^*S_\perp(\eta)(W-v(\eta)),$$
   with $c_0=\inf_{K'}\lambda_{\min}(S_\perp)\cdot\inf(y/s)/2$.
5. **(A4)** $c_0(W,\eta)=j(W)^{-1/2}(1+W\cdot v(\eta))$ has real part $\ge j^{-1/2}(1-4\gamma^2)$ on $|W|,|v|\le2\gamma$, and $C_2=8\sup\|S_\perp\|\sup(y/s)$.

The constants depend only on the pencil, $\omega$, $w_0$ and $K'$.

*Proof.*
1. **Conservation.** 08 Lemma 3.1 applies to $\mathscr P$, which has constant real symmetric coefficients. So the lifted current $B=\operatorname{Re}\sum_\ell\mathcal Z_\ell\otimes\overline{\mathcal Z_\ell}$ is positive semidefinite and conserved on $(\mathsf X,w_1,w_2)$. Extend it by zero to a passive $\eta\in\mathbb C_{\mathrm k}$. Let
$$\Theta(\mathsf X,z,\eta)=\big(\mathsf X,\ (\tfrac12z^{\mathsf T}C^\mu z)_\mu+(0,\omega)\eta\big).$$
   - **$\Theta$ is a local biholomorphism.** Its complex Jacobian in $(z,\eta)$ is $[M\,|\,(0,\omega)]$, with determinant $R^n\det[(C^\mu e(W))_a\,|\,(0,\omega)]$. At $W=0$ the first row is $(1,0,\dots,0\,|\,0)$, and the minor is $\det[L^{\mathsf T}\,|\,\omega]\ne0$ (Lemma 6.2.2). So $\Theta$ is a local biholomorphism wherever $R\ne0$ and $|W|<2w_0$, after shrinking $w_0$.
   - **Second-order term.** The two state blocks of $B$ are equal and the mixed block vanishes; this is the structure of $\mathcal Z_\ell$. Also $(\partial_{w_{1a}w_{1b}}+\partial_{w_{2a}w_{2b}})\Theta=0$ by holomorphy in $z$, and $\Theta$ is affine in $\mathsf X$ and $\eta$. So 08 Lemma 3.2 applies on every open set where a branch of $R$ is fixed.
   - **The pushed current.** The Gram-row computation of 08 Proposition 3.3, step 3, uses no property of the $C^\mu$. It identifies the push-forward with $(\chi\chi^*)^{\mathbb R}$.
   - **Single-valuedness.** Changing the branch of $R$ changes $\chi$ to $-\chi$.
2. **Tangency.** $\Phi\circ\Theta=\eta$ does not depend on $z$. So $\sum_\mu\partial_\mu\Phi\,M_{\mu a}=\partial_{z_a}(\Phi\circ\Theta)=0$.
3. **Explicit form.**
   - **Volume.** Each column $M_{\cdot a}$ carries a factor $R$.
   - **Gaussian.** $\overline{u_\xi(\mathsf X-\mathrm iq)}=\overline{u_\xi}\,e^{\xi\cdot q'-\varphi q^0}$, and $\xi\cdot q'-\varphi q^0=-\frac12z^*S(\xi)z$. Since $S(\xi)E(\xi)=0$ and $z-RE(\xi)=R(0,W-v)$, this equals $-\Psi_\xi$.
   - **Rows.** $\sum_aM_{\mu a}E_a(\xi)=R\,e(W)^{\mathsf T}C^\mu E(\xi)$, and $J^{-1/2}|R|=y^{-d/2}j^{-1/2}$.
4. **(A3).**
   - **Splitting each row.** Write $c_\mu=c'_\mu+\mathrm kc''_\mu$. Then $c\,\mathrm i\,\bar u=\mathrm i\,\overline{c'u}+\mathrm k\mathrm i\,\overline{c''u}$.
   - **Orthogonality.** $\mathrm i\overline{\mathbb C_{\mathrm i}}\subset\operatorname{span}(1,\mathrm i)$ and $\mathrm k\mathrm i\overline{\mathbb C_{\mathrm i}}=\mathrm j\overline{\mathbb C_{\mathrm i}}\subset\operatorname{span}(\mathrm j,\mathrm k)$ are orthogonal. With the unit $R/|R|$ removed, this gives $|\chi_\mu|^2=|m'_\mu(hD)g_T|^2+|m''_\mu(hD)g_T|^2$.
   - **Symbol bounds.**
     - $e(W)^{\mathsf T}C^\mu E(\eta)|\eta|^{d/2}$ is smooth and bounded on $K'$, uniformly for $|W|\le w_0$.
     - $(hy)^{-d/2}\le C(hs)^{-d/2}$.
     - Every $\eta$-derivative of $\Psi$ is $O(y/h)=O(\sigma)$, so $|\partial^\beta e^{-\Psi}|\le C\sigma^{|\beta|}e^{-\Psi}$.
     - $\Psi\ge\frac y{2h}\lambda_{\min}(S_\perp)|W-v|^2\ge c_0\sigma|W-v|^2$, because $(W-v)^*S_\perp(W-v)$ is the sum of the real forms at $\operatorname{Re}(W-v)$ and $\operatorname{Im}_{\mathrm k}(W-v)$.
5. **(A4).**
   - $e(W)^{\mathsf T}IE=1+W\cdot v$, and $j$ is real.
   - $1/m'_0=(\operatorname{Re}c_0)^{-1}|\eta|^{-d/2}(hy)^{d/2}e^{\Psi}$.
   - On $|W|,|v|\le2\gamma$, $\Psi\le\frac y{2h}\|S_\perp\|16\gamma^2\le C_2\gamma^2\sigma$, and the derivatives of $e^{\Psi}$ cost $\sigma^{|\beta|}$.
   - Finally $|\chi_0|^2\ge|m'_0(hD)g_T|^2$ by 4. ∎

**Proof of Theorem 11.2.** Take $w_0$ from Lemma 6.2. Take $K'$ to be a closed annular cap $\{c_1\le|\eta|\le c_2,\ \eta/|\eta|\in\Gamma'\}$, where $\Gamma'\subset\Gamma_\omega$ is a closed conic neighbourhood of $\omega$ with $|v|\le w_0/3$; this is possible since $v(\omega)=0$. Lemmas 6.3 and 6.4 give an engine on $(\varphi,K',v)$. ∎

**Proof of Corollary 11.3.** Let $\operatorname{supp}b\subset\{c_1<|\eta|<c_2\}$.
1. **Caps.**
   - For each direction $\omega$ of $\operatorname{supp}b$, Theorems 11.2 and 11.1 give $\gamma(\omega)$. The set $\{|v_\omega|<\gamma(\omega)\}$ is a conic neighbourhood of $\omega$.
   - By compactness, finitely many of these cover $\operatorname{supp}b$. Choose a smooth partition $b=\sum_ib_i$ subordinate to them.
   - Each $b_i(hD)$ has a periodized kernel of $L^1$ norm $\le\|\check b_i\|_1$. Theorem 11.1 then gives, on $\mathbb T^n$,
$$\Big(\int_1^2\|e^{-\mathrm iT\varphi(D)}b(hD)f\|_{L^p(\mathbb T^n)}^pdT\Big)^{1/p}\le C[\log(2/h)]^{2/p}\|f\|_{L^p(\mathbb T^n)} .$$
2. **Euclidean transfer.**
   - Let $M\ge1$ exceed $\sup_{\operatorname{supp}b}|\nabla\varphi|$. Apply step 1 to the pencil $C/M$ (same eigenvectors, top eigenvalue $\varphi/M$) and to $b_M=b(\cdot/M)$.
   - For $T\in[1,2]$ and $|X|\ge3$, the kernel phase has gradient $X-T\nabla\varphi(\eta)/M$ of length $\ge|X|/3$. So 08 §9, steps 2–3 (kernel decay and the cube argument), give the estimate on $\mathbb R^n$ for $(C/M,b_M)$.
   - With $f_M=f(M\cdot)$, $e^{-\mathrm iT\varphi(D)/M}b_M(hD)f_M=\big(e^{-\mathrm iT\varphi(D)}b(hD)f\big)(M\cdot)$, and both norms scale by $M^{-n/p}$.
3. **Large $h$.** For $h\ge h_0$, $e^{-\mathrm iT\varphi(D)}b(hD)$ has a smooth symbol with compact support, uniformly in $T\in[1,2]$, so it is bounded on $L^p$.
4. **The bottom eigenvalue.** Let $\psi=\lambda_{\min}\circ C$, so $-\psi=\lambda_{\max}\circ(-C)$. Let $\mathcal Cf(X)=\overline{f(-X)}$, an isometry of $L^p$ with $\widehat{\mathcal Cf}=\overline{\hat f}$. Then $\mathcal C\,e^{\mathrm iT\psi(D)}\tilde b(hD)\,\mathcal C=e^{-\mathrm iT\psi(D)}\bar{\tilde b}(hD)$. Apply the top-eigenvalue result to $-C$ with $\tilde b=\bar b$. ∎

**Proof of Corollary 11.4.**
- **The arrow pencil.** Write $A(\xi)$ for it.
  - Its eigenvector equations give $(t-\xi_n)(t+\xi_n)=|\xi'|^2$ on vectors with nonzero first coordinate, and the eigenvalue $-\xi_n$ on $\{(0,b):b\perp\xi'\}$. So $\det(tI-A(\xi))=(t+\xi_n)^{n-2}(t^2-|\xi|^2)$.
  - $\lambda_{\max}=|\xi|$ is simple off $\{\xi'=0,\ \xi_n\le0\}$, with eigenvector $(1,\xi'/(|\xi|+\xi_n))$.
  - $D^2|\xi|=(I-\hat\xi\hat\xi^{\mathsf T})/|\xi|$ is positive definite on $\xi^\perp$.
  - With $e(W)=(1,W)$: $e^{\mathsf T}(E_{0j}+E_{j0})e=2W_j$ and $e^{\mathsf T}\operatorname{diag}(1,-I)e=1-W\cdot W$. So $B(W)=(2W,1-W\cdot W)/(1+W\cdot W)$, the inverse stereographic projection, and $B(v(\xi))=\xi/|\xi|$.
- **Theorem 1.2.** The two pencils $A$ and $A\circ\varrho$, for a rotation $\varrho$ with $\varrho(-e_n)\ne-e_n$, cover all directions. Apply Corollary 11.3 to a partition of the cutoff of Theorem 1.2.
- **The paper's system.** In 08 §3.1, change the frame from $(t,s_0)$ to $(T,b)$: $C^T=(C^t-C^{s_0})/\sqrt2=\operatorname{diag}(\sqrt2,2^{-1/2},\dots)$ and $C^b=(C^t+C^{s_0})/\sqrt2$. With $D=(C^T)^{-1/2}$:
$$DC^TD=I,\qquad DC^bD=\operatorname{diag}(1,-I),\qquad DC^{x_j}D=E_{0j}+E_{j0}.$$
  This is the arrow pencil.
- **The passive coordinate.** 08 adds the passive coordinate along the null direction $(-1,0,1)/\sqrt2$ in $(T,x,b)$. Lemma 6.3 adds it along $(0,e_n)$. Both constructions satisfy the axioms. The second is a new derivation of the cap estimate. ∎

**The paraboloid.** $e^{\mathrm it\Delta}$ is not generated by a first-order symmetric system. Its symbol is not homogeneous of degree one, so it is not a direct instance of Theorem 11.2. Its critical estimate follows from the wave by parabolic rescaling of cone caps ([10](10_GAPLESS_REMAINING.md) §4, Corollary 4.1). Whether an engine exists *directly* for the paraboloid is open; one would need an exactly conserved, positive, tangent current for the Schrödinger flow. This is target `T.KAK.12`.

---

## 7. The hypothesis-usage matrix

### 7.1 The axioms of Theorem 11.1

Rows are the axioms. Columns are the steps of §3–§5. Each cell says what the axiom provides at that step.

| Axiom | Potential (08 Thm 5.2) | Balance (3.2a) | Comparison (3.2b) | Upper bound (Step 3) | Flux (Step 4) | Lower bound (Step 5) | Ratio (Step 6) | Propagation (§5) |
|---|---|---|---|---|---|---|---|---|
| (G1) chart | smooth bounded fibre, $H_{\rm fib}(\mathrm r)>0$ | Green on $\Omega$ | | | coarea | ball inside the fibre | | |
| (G2) $\lvert\partial W\rvert\lesssim1/s$, inverse scales | | | | $H_{\rm fib}(\mathrm r)\ge cs^{-2}$, $GG^*\le H$ | $\lvert G(\lvert W\rvert^2)\rvert\lesssim1/s$ | ball inclusion; centre volume | | |
| (G2) volume $a^2s^{2d}$ | | | | volume factor | volume factor | volume factor | $a^2$ cancels | |
| (G2) background $\lvert\partial_P\Phi\rvert\asymp1$ | $g>\det B_{\rm fib}=0$ | $\det B=0$ | rank one | $\beta\gtrsim\lambda$ | | $\beta\lesssim\lambda$ | | |
| (A1) $A^0\ge0$ | $\mu>0$ | $A:\Sigma\ge0$ dropped | | | | | | $\mathcal F\ge0$ |
| (A1) conservation | | Green with $M=A$ | | | | | | |
| (A2) tangency | | $A:B=0$; $g(Q)$ keeps conservation | | | kills $G(\lvert\Phi\rvert^2)$ | | | |
| (A3) global bound | $\mu^{1/d}\le\lambda$ | | | $\mu\lesssim s^{-2d}\sigma^mP_*^d$ | | | | |
| (A3) Gaussian centred at $v$ | | | | | $e^{-c_0'\sigma}$ on the wall | | decay $c_0'$ | $\int s\log(s/b)V<\infty$ |
| (A4) invertible time row | | | | | | $\mu\ge m_X$ | cost $C_2\gamma^2$ | |
| $\dim\tau=2$ | | $\Delta_z$ on $\mathbb R^3$ | | | | | | $\log^2$ (Lemma 8.1) |

### 7.2 What breaks without each hypothesis

| Hypothesis | Without it |
|---|---|
| **Positivity** (A1) | $\det D_u=\lambda\mu$ has no positive-definite solution where $\mu<0$, and the Schur remainder $A:\Sigma$ can no longer be dropped. |
| **Conservation** (A1) | Green's identity (08 (6.5)) gains a bulk term $\int u\,\eta\,\operatorname{div}A$, so there is no balance. |
| **Tangency** (A2) | $A:B=\lambda A:H(\lvert\Phi\rvert^2)\neq0$. Then $\Delta_zX_0$ gains a term of size $\lambda=h^{-\Lambda_0-1}$, and the balance fails by a power of $h$. |
| **Rank-one background** (G2) | With no background, the depth is $\propto\mu^{1/n}$, which gives the exponent $2(n+1)/n$ instead of $p_c$. With rank $r$, it gives $2(f+1)/f$ with $f=n-r$ ([07](07_TOOLKIT.md) K3). The comparison 3.2(b) is the rank-one expansion of the Moore determinant; the rank-$r$ analogue is open (`T.KAK.9`). |
| **Gaussian wall decay** (A3) | If $V(s)\sim s^{-2}$ only, then $\int_b^3s\log(s/b)V\,ds\sim\log^2(1/b)$. Lemma 8.1's Gronwall factor becomes $e^{C\log^2(1/h)}$, which loses more than any power of $h$. |
| **Decay beats inversion** ($2qC_2\gamma^2<c_0'$) | The ratio $\Delta\mathcal F/\mathcal F$ loses its exponential factor, with the same consequence. This is why the cap aperture is small. |
| **Two-dimensional auxiliary time** | For $\tau\in\mathbb R^k$ with $k\neq2$, superharmonic transport between two spheres loses a power of $h$ ([07](07_TOOLKIT.md) K4, proved and checked). |

### 7.3 The pencil hypotheses (Theorem 11.2)

| Pencil hypothesis | Axiom it gives | Without it |
|---|---|---|
| Constant coefficients, first order, symmetric | (A1): exact plane-wave solutions extend holomorphically to $\mathsf X-\mathrm iq$ (08 Lemma 3.1) | No exactly conserved current is known. This is the only place the global structure of the equation enters, so any variable-coefficient extension must replace (A1). |
| **Extremal** eigenvalue | $S=\varphi I-C\ge0$, so $e^{-\frac12w^*Sw}\le1$ is a window | For a middle sheet, $S$ is indefinite. The current is still conserved (Part M), but its rows carry $e^{+c\sigma}$: an anti-window. |
| **Simple** eigenvalue | $S_\perp\ge c\lvert\xi\rvert$ (spectral gap), so $c_0>0$ | No smooth polarization, and no decay rate. |
| **Transversal curvature** | $L$ is an isomorphism on $\omega^\perp$, hence $dB(0)$ is invertible: the chart (G1) exists, and $v$ is a local angular coordinate | A flat sheet ($D^2\varphi=0$ on $\omega^\perp$) gives no chart. |
| **State dimension $=n$** | $\Theta$ is a local biholomorphism, so the push-forward is pointwise rank one with explicit rows (A3) | If $N>n$, $\Theta$ is a submersion with $(N-n)$-dimensional fibres. The push-forward is an integral of rank-one currents, still positive, conserved and tangent, but (A3)–(A4) must be re-derived (`T.KAK.7`). If $N<n$, there is no current on the full target. |

---

## 8. What the abstraction explains

1. **Why the exponent is $p_c=2n/(n-1)$.**
   - The Monge–Ampère depth estimate gives $(-u)\gtrsim h^2\mu^{1/f}$, where $f$ counts the fibre directions not used by the background. Here $f=n-1$ because the background has rank one (Lemma 3.1.2).
   - So the energy density is $(-u)\mu\propto\mu^{1+1/f}=|g|^{2(f+1)/f}$.
   - Nothing else in the proof selects $p$. The volume factor $h^2s^{2d}(hs)^{-n}=\sigma^{n-2}$ is harmless on the two circles, where $\sigma$ is fixed.
2. **The wave is not special.** It is the arrow pencil. Its polarization is the stereographic projection of the direction, and its group-velocity map is the inverse stereographic projection. The paper's Gaussian direction window $e^{-\frac{\kappa y}2|W-v|^2}$ is a Gaussian in stereographic coordinates on the sphere of directions. App. E replaces stereographic projection by the eigenvector map of a pencil.
3. **Why extremal sheets.** The window is $e^{-\frac12w^*S(\xi)w}$. It decays exactly when $S=\varphi I-C\ge0$, that is, for an extremal eigenvalue. Positivity of the current costs nothing; positivity of $S$ is what selects the sheet.
4. **Where the equation enters.** The equation enters only through (A1)–(A4). Of these:
   - (A3)–(A4) are bounds on explicit Gaussian symbols, stable under small perturbations of $C_j$ (a perturbed arrow pencil is checked in Part C).
   - (A1) needs exact solutions at complex points $\mathsf X-\mathrm iq$.
   
   This sharpens the open question of [07](07_TOOLKIT.md) K7.2 (variable coefficients, `T.KAK.8`): the obstruction must lie in (A1), not in the analysis.
5. **Why the aperture is small.** $\gamma$ balances the wall decay $c_0$ against the cost $C_2$ of inverting the window at the centre. It is a property of the Gaussian, not of the equation.

---

## 9. Status and cross-references

| Result | Status | Uses |
|---|---|---|
| Theorem 11.1 | proved here | 08 Lemmas 2.1–2.4 (calculus), 5.1 (comparison), 7.2 (kernels), 8.1 (propagation); 08 Thm 5.2 (Dirichlet: Schauder and Evans–Krylov cited; also formally verified); 08 Prop. 6.2 and Lemma 6.3, whose proofs use only the listed properties |
| Lemmas 6.1–6.3 | proved here | Kato (analytic eigenprojection), inverse function theorem; **checked**: Part G |
| Lemma 6.4 | proved here | 08 Lemmas 3.1–3.2 and the Gram rows of Prop. 3.3; **checked directly in target coordinates**: Part C (5 pencils, conservation to $10^{-38}$) |
| Theorem 11.2, Corollary 11.3 (= paper Thm E.1) | proved here | 08 §9 (transfer) |
| Corollary 11.4 (wave = arrow pencil) | proved here; **exact check** Part W | |
| Middle sheets fail (A3) | checked (Part M) | |

**Verdict on App. E.** I found no error. The paper's App. E states the chart, the current and the two Gaussian estimates, and refers to the wave proof for the rest. This file supplies the quantitative details:
- the injectivity of the chart and the derivative bounds of its inverse;
- the ball inclusion;
- the size of $\mu$;
- the hypotheses of the Dirichlet theorem.

It also writes out the reuse of the wave proof as one abstract argument. In [`ROADMAP.md`](../../../ROADMAP.md), `KAK.APP.E` is upgraded from PARTIAL to PROVED.

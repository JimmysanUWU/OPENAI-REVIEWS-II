# 5. Distillation: positive currents, determinant energies, and the critical estimate

This is a rewrite of the proof of Theorem 1.2 as four modular theorems, each with its hypotheses stated abstractly, followed by the assembly. The aim is that each design choice in §2–§6 appears as a checkable hypothesis.

**Conventions.**
- $\mathbb H$ denotes the quaternions; $\mathcal D=\operatorname{span}_{\mathbb R}(1,\mathbf i,\mathbf j)$ is the 3-slot space; $e_0,e_1,e_2=1,\mathbf i,\mathbf j$.
- For $m$ slots with real coordinates $x_{a\alpha}$:
$$G_a(u)=\sum_\alpha e_\alpha\partial_{a\alpha}u,\qquad H_{ab}(u)=\sum_{\alpha,\beta}e_\alpha\bar e_\beta\,\partial_{a\alpha}\partial_{b\beta}u .$$
- $A:B=\operatorname{Re}\operatorname{tr}(AB)$; $\det$ is the Moore determinant and $C(\cdot)$ its cofactor.
- A Hermitian field $A$ is **conserved** if its real tensor $A^{\mathbb R}_{a\alpha,b\beta}=\operatorname{Re}(\bar e_\alpha A_{ab}e_\beta)$ has divergence-free columns. Then $A:H(u)=\operatorname{div}(A^{\mathbb R}\nabla u)$.

---

## 5.1 Theorem I (positive wave current)

*Let $\mathscr P=\sum_\mu C^\mu\partial_\mu$ have constant real symmetric $n\times n$ coefficients on $\mathbb R^{n+1}$. Let $f$ be a finite sum of holomorphic plane-wave solutions of $\mathscr Pf=0$, and for auxiliary states $w_1,w_2\in\mathbb R^n$ put $q^\mu=\frac12\sum_\ell w_\ell^tC^\mu w_\ell$. Then:*

1. *(Lifted current.) With $G=\overline{f(\mathsf X-\mathrm iq)}$, $D_{\ell a}=\partial_{w_{\ell a}}+\mathrm i\sum_\mu(C^\mu w_\ell)_a\partial_\mu$ and $\mathcal Z_\ell=\sum_aG_aD_{\ell a}$, the tensor*
$$B=\operatorname{Re}\sum_{\ell=1}^2\mathcal Z_\ell\otimes\overline{\mathcal Z_\ell}$$
*on $(\mathsf X,w_1,w_2)$ is positive semidefinite and conserved.*

2. *(Quaternionic push-forward.) Let $z=w_1+\mathbf kw_2$ and $\Theta(z,\eta)=(\tfrac12z^tC^\mu z+c^\mu\eta)_\mu$. On a local inverse chart, the density push-forward $\tilde B=J^{-1}D\Theta\,B\,D\Theta^t$ is conserved and has the form $(\chi\chi^*)^{\mathbb R}$ for a quaternionic column $\chi$ on $n+1$ slots of type $\mathcal D$.*

*For the half-wave in null coordinates (Prop 2.2, (2.15)), $\chi_T$ is a pair of Gaussian-windowed copies of the wave:*
$$A^0_{TT}=|\mathsf S_0u|^2+|\mathsf S_1u|^2,\qquad \text{window }(\kappa/y)^{d/2}e^{-\frac{\kappa y}2|W-v|^2},$$
*and $(\partial\Phi)A^0=0$ for the passive coordinate $\Phi$.*

*Proof.*
1. Each $D_{\ell a}$ annihilates $\mathsf X-\mathrm iq$. Hence $\partial_{w_{\ell a}}G=\mathrm iV_{\ell a}G$, so $\operatorname{div}\mathcal Z_\ell=2\mathrm iw_\ell^t\mathscr PG=0$, and $\nabla_{\mathcal Z_\ell}\overline{\mathcal Z_m}=-\mathrm i\delta_{\ell m}\sum_\mu(G^tC^\mu\bar G)\partial_\mu$, which is purely imaginary.
2. The second-order term $\sum B^{IL}\partial_I\partial_L\Theta$ vanishes: $\Theta$ is holomorphic in $z$, so $\partial_{w_1w_1}+\partial_{w_2w_2}$ kills it, and the two state blocks of $B$ are equal. The Gram rows $Q_I=\mathcal Z_1^I+\mathbf k\mathcal Z_2^I$ become $-\mathbf iQ,-\mathbf jQ$ under $D\Theta$. ∎

(Verified: [`scripts/verify_wave_current.py`](scripts/verify_wave_current.py).)

## 5.2 Theorem II (determinant energy balance)

**Setting.** Slots $0,1,\dots,n$, all of type $\mathcal D$:
- slot 0 carries the parameter $z\in Z\subset\mathbb R^3$;
- slots $1,\dots,n$ carry the fiber variable $\zeta$, possibly periodic in some real directions.

Let $\Omega=\{(z,\zeta):r(z,\zeta)<0\}$ with fibers $\Omega_z$, and let $\nabla_\zeta r\neq0$ on the wall. "Fib" denotes the $n\times n$ fiber block.

**Hypotheses.**
- **(A)** $A\ge0$ is a smooth conserved field on $\bar\Omega$ with $\mu:=A_{00}>0$.
- **(B)** $B_0$ is smooth with $B:=H(B_0)\ge0$, $\det B=0$, $A:B=0$, and $B_{\rm fib}$ of rank one.
- **(W)** On the wall, there is a nonzero $e\in G(r)^\perp$ with $H(r)e=0$ and $Be=0$. Also $A=A_{\rm w}+A_{\rm reg}$ with $A_{\rm reg}G(r)=0$ there.
- **(D)** $u$ is smooth on $\bar\Omega$, $u=0$ on the wall, $D_u:=H_{\rm fib}(B_0+u)>0$, and $\det D_u=\lambda\mu$ on each fiber.

*Then $\mathcal F(z)=\int_{\Omega_z}(-u)\mu+\lambda^{-1}\int_{\Omega_z}u\int_0^1\det D_{\theta u}\,d\theta$ satisfies*
$$\Delta_z\mathcal F=\int_{\Omega_z}\rho\,G(r)^*A_{\rm w}G(r)\,\delta(r)\,d\zeta-\int_{\Omega_z}A:\Sigma\ \le\ \int\rho\,G(r)^*A_{\rm w}G(r)\,\delta(r),$$
*where $u=\rho r$ near the wall ($\rho\ge0$) and $\Sigma\ge0$ is the Schur remainder of $H(B_0+u)$ with respect to slot 0. Moreover*
$$\frac{n-1}{n+1}X_0\le\mathcal F\le\frac n{n+1}X_0,\qquad X_0=\int(-u)\mu .$$

*Proof.* Use Green's formula $\int uM_{00}\Delta\eta=\int\eta\,M:H(u)-\int_{\rm wall}\eta\,(\nu^*M\nu)\partial_{\mathbf n}u$ for conserved $M$ and $\eta=\eta(z)$.

1. $M=A$. Here $A:H(u)=A:\Sigma+\lambda^{-1}\det H(B_0+u)$ (Schur and $A:B=0$). On the wall only $A_{\rm w}$ contributes.
2. $M=C(H(B_0+\theta u))$, which is conserved by cofactor divergence. Its slot-0 entry is $\det D_{\theta u}$, and $M:H(u)=\partial_\theta\det H(B_0+\theta u)$.

   On the wall, $\nu^*M\nu=\det(H(B_0+\theta u)|_{\nu^\perp})=\det((B+\theta\rho H(r))|_{\nu^\perp})=0$, because $e$ lies in its kernel. Integrating in $\theta$ and using $\det B=0$ gives $\Delta E_z(u)=\lambda^{-1}\int\det H(B_0+u)$.

The two determinant terms cancel. For the comparison: $B_{\rm fib}$ has rank one, so $\det D_{\theta u}=\theta^n\det D_u+(1-\theta)\theta^{n-1}C(D_u):B_{\rm fib}$; integrate this in $\theta$ and integrate by parts once. ∎

(The algebraic inputs were verified in [`scripts/verify_quaternion_calculus.py`](scripts/verify_quaternion_calculus.py).)

**Remark (why this is a "linear compensated integrability").** Serre's inequality bounds $\int(\det A)^{1/(n-1)}$ for conserved $A\ge0$, a multilinear quantity, which vanishes for the rank-one wave current. Theorem II never takes $\det A$. It feeds $A_{00}$ into a Monge–Ampère equation whose background supplies the missing direction. That is what makes the resulting estimate *linear* in the wave.

## 5.3 Theorem III (density to depth, at the critical exponent)

*Let $d=n-1$, and let $B_{\rm fib}=\lambda H_{\rm fib}(|\Phi|^2)$ have one eigenvalue $\beta\in[c\lambda,C\lambda]$ and the rest zero. Let $u\le0$ solve $\det H_{\rm fib}(B_0+u)=\lambda\mu$ with $u=0$ on the wall. Then:*

1. *(Interior.) If $\mu\ge m$ on a fiber ball $B_\rho(x)$, $m^{1/d}\le\lambda$, and $\rho\le\operatorname{dist}(x,\text{wall})$, then $-u(x)\ge\epsilon\rho^2m^{1/d}$, with $\epsilon=\epsilon(n,c,C)>0$.*
2. *(Barrier.) If $H_{\rm fib}(r)\ge\kappa I$ and $G_{\rm fib}(r)G_{\rm fib}(r)^*\le C'H_{\rm fib}(r)$, and $P>0$ satisfies $|\nabla^jP|\le\epsilon_0\kappa^{j/2}P$ ($j=1,2$) and $\mu\le(c_1\kappa P)^d$, then $u\ge C_2rP$, i.e. $0\le u/r\le C_2P$.*

*Consequently, $\int(-u)\mu$ is comparable, up to the window factors, to $\int\mu^{1+1/d}$. With $\mu\approx|\text{wave}|^2$, $2(1+1/d)=p_c$.*

*Proof.*
1. Compare $u$ with $\epsilon m^{1/d}(|y-x|^2-\rho^2)$. Its fiber Hessian is $6\epsilon m^{1/d}I$, so $\det(B_{\rm fib}+6\epsilon m^{1/d}I)\le C\lambda(6\epsilon)^dm\le\lambda m$.
2. By the product rule $H(rP)\ge(1-C\epsilon_0)PH(r)$, so $\det(B_{\rm fib}+H(C_2rP))\ge\lambda(cC_2\kappa P)^d\ge\lambda\mu$.

Both cases then follow from the comparison principle for the fiber equation. ∎

## 5.4 Theorem IV (two-circle propagation)

*Let $E\ge0$ be smooth on $\mathbb R_T\times\{b\le|\tau|\le3\}\subset\mathbb R^3$ with $\Delta E\le V(|\tau|)E$, where $V\ge0$, $\int_b^3s\log(s/b)V\le M$ and $\sup_{[1/2,3/2]}V\le M$. Then the circular means satisfy, for $b\le r\le1$,*
$$H(0,r)\ \ge\ c(M)\,\frac{\log(r/b)}{\log^2(3/b)}\int_1^2H(T,b)\,dT .$$
*The factor $\log^{-2}$ is sharp for scalar supersolutions (Remark 6.2; verified exactly in [`scripts/verify_propagation_and_exponents.py`](scripts/verify_propagation_and_exponents.py)).*

## 5.5 Assembly (Theorem 1.2)

Work on a fixed torus. Take cap data with $\|f\|_2=1$, $h=N^{-1}$, $a=Kh$, $\lambda=h^{-30n}$.

1. Theorems I and II give $\Delta\mathcal F\le\mathcal B$.
2. Theorem III applied to the wall flux and the interior gives:
   - $\mathcal B\le Cs^{-2}e^{-cs/h}\mathcal F$ on $2a\le s\le3$;
   - $\mathcal F(0,\cdot)\lesssim a^2\|f\|_{p_c}^{p_c}$ on $|\tau|=4a$;
   - $\mathcal F(T,\cdot)\gtrsim a^2\|U(T)f\|_{p_c}^{p_c}$ on $|\tau|=2a$.
3. Theorem IV with $b=2a$ and $r=4a$ gives $\int_1^2\|U(T)f\|_{p_c}^{p_c}\lesssim\log^2(2/h)\|f\|_{p_c}^{p_c}$.

Homogeneity removes $\|f\|_2=1$. Finitely many caps cover the annulus, and §7 removes the period. ∎

## 5.6 The planar case as a five-page entry point

For $n=2$ ($d=1$), the coefficient of $\lambda$ in $\det(\lambda B+H)$ is *linear* in $H$ (App. D). Theorem III then becomes a statement about a linear Dirichlet problem on a disk: depth $\propto$ density, and $(-u)\mu\propto\mu^2=|{\rm wave}|^4$.

The whole proof of sharp $L^4$ local smoothing in $\mathbb R^{2+1}$ (with $\log^{1/2}$) then needs only:
- Theorem I;
- a linear version of Theorem II;
- Theorem IV.

I recommend writing this case out first in any textbook treatment. It contains every idea except quaternionic Monge–Ampère regularity.

## 5.7 Open distillation targets

1. **A coordinate-free Theorem I.** State the current on the complexified characteristic variety as a CR-geometric object. The paper's (2.7)–(2.9) are a coordinate shadow of this.
2. **Avoid the fourth slot** in the boundary $C^2$ estimate (App. A). Alternatively, prove that 3-slot rigid motions suffice.
3. **A single abstract "positivity ⇒ restriction" theorem.** State hypotheses on $(\mathscr P,\Phi,r)$ under which Theorems I–IV assemble. App. E (symmetric pencils with definite transverse Hessian) is a first instance. This would also locate exactly where constant coefficients are used, explaining why variable-coefficient counterexamples do not apply.
4. **One logarithm instead of two.** Remark 6.2 shows this needs information beyond scalar supersolutions, such as the actual energy profile.
5. **The sharp $2/p_c$ power in Lean.** The formal proof uses $2n+5$ ([06](06_LEAN_AUDIT.md)). Closing that gap would formalize Theorem 1.2 as stated.

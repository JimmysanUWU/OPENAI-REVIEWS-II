# 13. The planar entry point: sharp $L^4$ local smoothing in $\mathbb R^{2+1}$ by a linear Dirichlet problem

*Claude, 28 October 2026. Roadmap target `T.KAK.2`; node `KAK.APP.D` in [`ROADMAP.md`](../../../ROADMAP.md).*

> **Theorem 13.1.** For $\varphi\in C_c^\infty(\mathbb R^2\setminus0)$, $N\ge2$ and $f\in L^4(\mathbb R^2)$,
> $$\|e^{-\mathrm iT|D|}\varphi(D/N)f\|_{L^4(\mathbb R^2\times[1,2])}\le C_\varphi\,[\log(2+N)]^{1/2}\,\|f\|_{L^4(\mathbb R^2)} .$$

> **Corollary 13.2 (local smoothing in $\mathbb R^{2+1}$).** For $2\le p<\infty$ and $\varepsilon>0$, $\|e^{-\mathrm iT|D|}f\|_{L^p(\mathbb R^2\times[1,2])}\le C\|J^{s_+(p)+\varepsilon}f\|_p$, where $s_+(p)=\max(0,\frac12-\frac2p)$.

**Context.** For $n=2$, the $\varepsilon$-loss form of Corollary 13.2 was first proved by Guth, Wang and Zhang (Ann. of Math. 2020), through a sharp square-function estimate for the cone in $\mathbb R^3$. The proof below is entirely different. It gives a logarithmic loss at $p=4$, and its only analytic input is a **linear** Dirichlet problem on a disk. There is no Monge–Ampère theory: no Evans–Krylov, no boundary $C^2$ estimate, no fourth slot.

**What this file is.** A complete proof of Theorem 13.1, following the paper's App. D. It is organized as the $n=2$ case of the engine of [11](11_FRAMEWORK_POSITIVITY_ENGINE.md), but with the nonlinear fibre problem replaced by a linear one.
- **Proved here:** everything specific to $n=2$:
  - the explicit chart;
  - the disk operator;
  - the linear potential;
  - the planar energy and its completed-square inequality;
  - the window estimates.
- **Quoted with exact pointers:** the general lemmas:
  - the slot calculus: [12](12_QUATERNIONIC_SLOT_CALCULUS.md) §1–§4, or [08](08_GAPLESS_PROOF.md) §2;
  - the current: [11](11_FRAMEWORK_POSITIVITY_ENGINE.md) Lemma 6.4, or 08 §3;
  - the multiplier kernels and two-circle propagation: 08 Lemmas 7.2 and 8.1;
  - the transfer to $\mathbb R^2$: 08 §9.

  Every one of these is either proved in full in this repository or a textbook theorem. So, in the roadmap's terms, the prose chain of Theorem 13.1 is *closed*.

**Verdict on App. D.** It is correct. This file supplies:
- the explicit chart and the formula for $\mathcal J$;
- the proof that the mixed cofactor is conserved;
- the boundary vanishing $\nu^*M\nu=0$;
- the sign bookkeeping of the two Green identities;
- the comparison arguments in full.

**Checks.** [`scripts/verify_planar_entry.py`](scripts/verify_planar_entry.py) ([`results/planar_entry.json`](results/planar_entry.json), ~5 s):
- **Part A:** the chart identities of §1, including the characteristic conormal (exact).
- **Part L:** the disk-operator formula of Lemma 2.1, at random points by central differences in 50-digit arithmetic. The relative error is $4\cdot10^{-25}$.
- **Parts E and M:** the energy algebra of Lemmas 3.1–3.2 (numerical).
- **Part P:** conservation of the mixed cofactor, exact, with a negative control.
- **The current:** the conservation of the Dirac-pencil current is checked in [`scripts/verify_positivity_engine.py`](scripts/verify_positivity_engine.py) Part C (arrow pencil $n=2$).

---

## 0. Set-up

**Notation.** As in [11](11_FRAMEWORK_POSITIVITY_ENGINE.md) §0 with $n=2$, $d=1$ and $p=p_c=4$:
- three slots: $T+\tau\mathrm i$ and $X_a+P_a\mathrm i$ ($a=1,2$), with $\tau,P_a\in\mathbb C_{\mathrm k}$;
- $z=(T,\tau)$, $s=|\tau|$ and $\sigma=s/h$;
- the torus $\mathbb T^2=\mathbb R^2/(100\mathbb Z)^2$.

**Reduction.** By 08 §9, Theorem 13.1 follows from the **cap estimate**: for $f$ on $\mathbb T^2$ with Fourier support in a cap $\{c_1\le h|\xi|\le c_2,\ |v(h\xi)|\le\gamma\}$,
$$\int_1^2\|e^{-\mathrm iT|D|}f\|_{L^4(\mathbb T^2)}^4\,dT\le C\log^2(2/h)\,\|f\|_{L^4(\mathbb T^2)}^4 .\tag{0.1}$$
- **Finitely many caps** cover the annulus, and the cap multipliers have uniformly integrable kernels.
- **Off-cone kernel decay and the cube argument** give the passage to $\mathbb R^2$.
- **Interpolation** with $L^2$ and $L^\infty$ (08 (9.1)) gives Corollary 13.2.

By homogeneity, let $\|f\|_2=1$. For $h\ge h_0$, (0.1) is trivial (11 §5), so let $h<h_0$.

## 1. The Dirac pencil and its current

The arrow pencil for $n=2$ is
$$C(\xi)=\xi_1C_1+\xi_2C_2=\begin{pmatrix}\xi_2&\xi_1\\\xi_1&-\xi_2\end{pmatrix}.$$
- **Eigenvalues.** They are $\pm|\xi|$, so $\lambda_{\max}=|\xi|$ is simple on all of $\mathbb R^2\setminus0$ ([11](11_FRAMEWORK_POSITIVITY_ENGINE.md) Corollary 11.4). The system $\partial_T+C_1\partial_{X_1}+C_2\partial_{X_2}$ is the 2D Dirac operator, and its top sheet is the half-wave $e^{-\mathrm iT|D|}$.
- **Rotation covariance.** $C(R_\theta\xi)=U_{\theta/2}C(\xi)U_{\theta/2}^{\mathsf T}$, with $U$ the rotation of the state plane. So every cap is a rotate of the cap around $\omega=e_2$, and we work there.
- **Polarization.** With $e_0=(1,0)$:
  - $v(\xi)=\xi_1/(|\xi|+\xi_2)$, the stereographic coordinate of $\xi/|\xi|$;
  - $S_\perp(\xi)=|\xi|+\xi_2\asymp|\xi|$ on the cap;
  - $B(W)=(2W,1-W^2)/(1+W^2)$, the inverse stereographic map.

**The chart (Lemma 6.3 of 11, made explicit).** For $|W|<2w_0$ and $\Phi\in\mathbb C_{\mathrm k}$:
$$P_1=\frac{2\tau W}{1+W^2},\qquad P_2=\tau\,\frac{1-W^2}{1+W^2}+\Phi\qquad\Longleftrightarrow\qquad W=\frac{P_1}{\tau+\sqrt{\tau^2-P_1^2}},\qquad\Phi=P_2-\sqrt{\tau^2-P_1^2}.\tag{1.1}$$
Geometrically, $(P-\Phi e_2)/\tau$ lies on the complexified unit circle, and $W$ is its stereographic coordinate. Differentiating (1.1):
$$\partial\Phi=\Big(-\frac{1+W^2}{1-W^2},\ \frac{2W}{1-W^2},\ 1\Big)\ \text{in }(\tau,P_1,P_2),\qquad\partial_PW=\Big(\frac{(1+W^2)^2}{2\tau(1-W^2)},\ 0\Big).\tag{1.2}$$
At a real $W=v(\xi)$, $\partial\Phi=(-|\xi|,\xi_1,\xi_2)/\xi_2$, which is the characteristic conormal. Put
$$c=\partial_P\Phi=\Big(\frac{2W}{1-W^2},1\Big),\qquad e=\Big(1,-\frac{2W}{1-W^2}\Big)^{\mathsf T},\qquad\mathcal J=e\cdot\partial_PW=\frac{(1+W^2)^2}{2\tau(1-W^2)} .\tag{1.3}$$
Then $c\cdot e=0$, $|e|=|c|\asymp1$ and $|\mathcal J|\asymp1/s$. All of (1.1)–(1.3) is checked exactly (Part A).

**The current.** For cap data, [11](11_FRAMEWORK_POSITIVITY_ENGINE.md) Lemma 6.4 gives a column $\chi$ with:
- $A^0=\chi\chi^*$ positive semidefinite and conserved;
- $(\partial\Phi)\chi=0$;
- rows $|\chi_\mu|^2=|m'_\mu(hD)g_T|^2+|m''_\mu(hD)g_T|^2$ with Gaussian symbols. Its constants $c_0$ and $C_2$ satisfy (A3)–(A4).

## 2. The product fibre and the linear potential

**The fibre.** With $Q=\Phi/a$ and $a=Kh$, the fibre over $z$ ($a<s<4$) is the product
$$\Omega^{\rm pl}_z=\mathbb T^2\times\{(W,Q):\ |W|<w_0,\ |Q|<\tfrac1{10}\},\qquad r_W=|W|^2/w_0^2-1 .$$
It is embedded in the imaginary variables by (1.1). The data are:
- $A=g(Q)(J_h*_XA^0)$, with the cutoff $g$ of 08 §6.1 ($g=1$ on $|Q|\le\frac1{100}$, $g=0$ on $|Q|\ge\frac1{20}$);
- $\mu=A_{TT}\ge0$;
- $B=H(|\Phi|^2)=4(\partial\Phi)^*\partial\Phi$, of rank one.

No regularizing cofactor is needed: the linear problem only needs $\mu\ge0$. $A$ is conserved and $A(\partial\Phi)^*=0$, by the proof of [11](11_FRAMEWORK_POSITIVITY_ENGINE.md) Lemma 3.1.1.

**Lemma 2.1 (the disk operator).** $B_{\rm fib}=4c^*c$ and $C(B_{\rm fib})=4ee^*$. For real $u(X,W,Q)$ composed with (1.1),
$$L_Bu:=C(B_{\rm fib}):H_{\rm fib}(u)=4\sum_{a,b=1}^2\operatorname{Re}(\bar e_ae_b)\,\partial_{X_aX_b}u+4|\mathcal J|^2\Delta_Wu ,\tag{2.1}$$
where $\Delta_W$ is the Laplacian in the two real coordinates of $W$. No $Q$-derivatives and no mixed $X$–$W$ derivatives occur.

*Proof.*
- **The matrices.** $B_{\rm fib}=4c^*c$ is 08 Lemma 2.1(4) for $F=\Phi$. For $2\times2$ Hermitian matrices, $C\begin{pmatrix}m_{11}&m_{12}\\m_{21}&m_{22}\end{pmatrix}=\begin{pmatrix}m_{22}&-m_{12}\\-m_{21}&m_{11}\end{pmatrix}$, which gives $C(4c^*c)=4ee^*$, since $\mathbb C_{\mathrm k}$ is commutative.
- **Splitting $e^*H_{\rm fib}e$.** By Proposition 12.5, $H_{\rm fib}=D^2_Xu+4\bar\partial\partial u+(\text{terms in }\mathbb C_{\mathrm k}\mathrm i)$. Since $e$ has entries in $\mathbb C_{\mathrm k}$, the last terms contribute to $e^*H_{\rm fib}e$ only in $\mathbb C_{\mathrm k}\mathrm i$, which has zero real part. The real part of $e^*(D_X^2u)e$ is $\sum\operatorname{Re}(\bar e_ae_b)\partial_{X_aX_b}u$.
- **The complex part.** Let $\mathcal D=\sum_be_b\partial_{P_b}$. Then $\mathcal D\Phi=c\cdot e=0$ and $\mathcal DW=\mathcal J$, so in the coordinates $(W,\Phi)$, $\mathcal D=\mathcal J\partial_W$. As $\mathcal J$ is holomorphic,
$$e^*(4\bar\partial\partial u)e=4\bar{\mathcal D}\mathcal Du=4|\mathcal J|^2\partial_{\bar W}\partial_Wu=|\mathcal J|^2\Delta_Wu .$$
- **Conclusion.** Multiply by $4$. ∎

Part L checks (2.1), including $Q$-dependence of $u$, against a finite-difference 3-slot Hessian.

**Lemma 2.2 (the linear potential).** There is a unique real $u$, smooth on $\bar\Omega^{\rm pl}_z$ and in $z$, with
$$L_Bu=\mu\ \text{in }\Omega^{\rm pl}_z,\qquad u=0\ \text{on }|W|=w_0 .$$
Moreover $u=g(Q)\tilde u$ with $\tilde u$ smooth, so $u\equiv0$ for $|Q|\ge\frac1{20}$, and $u\le0$.

*Proof.*
- **Existence.** $\mu$ is a finite Fourier sum in $X$, because $A^0$ is a product of two finite sums. Write $u=\sum_k\hat u_k(W)e^{\mathrm ik\cdot X}$. By (2.1), for each $k$ and each fixed $(z,Q)$:
$$\Delta_W\hat u_k-V_k\hat u_k=\frac{\hat\mu_k}{4|\mathcal J|^2}\quad(|W|<w_0),\qquad\hat u_k=0\ (|W|=w_0),\qquad V_k=\frac{|k\cdot e|^2}{|\mathcal J|^2}\ge0 .$$
  The coefficients are smooth on the closed disk, since $\mathcal J\ne0$ there. The Dirichlet problem for $\Delta-V$ with $V\ge0$ is uniquely solvable in $C^{2,\alpha}$, with Schauder estimates (Gilbarg–Trudinger, Theorems 6.14 and 6.19). The solution depends smoothly on parameters: differentiate the equation and use uniqueness. Real and imaginary parts decouple, and $\hat u_{-k}=\overline{\hat u_k}$ makes $u$ real.
- **Factorization.** The operator does not involve $Q$, so solving with $J_h*A^0_{TT}$ in place of $\mu$ gives $\tilde u$, and $u=g(Q)\tilde u$.
- **Sign.** Fix $Q$. $L_Br_W=4|\mathcal J|^2\cdot4/w_0^2>0$, so $L_B(u+\epsilon r_W)>0$ for $\epsilon>0$.
  - At an interior maximum of a $C^2$ function $v$ in $(X,W)$ we have $D^2v\le0$, hence $H_{\rm fib}(v)\le0$ (Lemma 1.2.1 of 12) and $L_Bv\le0$.
  - So $u+\epsilon r_W$ has no interior maximum, and it is $0$ on the wall. Hence $u\le-\epsilon r_W$; let $\epsilon\to0$. ∎

The same maximum-principle argument (with a strict inequality) is used for every comparison below. It works for the degenerate operator $L_B$, on any domain, because $L_B$ has a positive semidefinite symbol.

## 3. The planar energy

**Lemma 3.1 (the mixed cofactor).** For $3\times3$ Hermitian $H$, let $\mathcal D_B(H)=B:C(H)$, the coefficient of $\lambda$ in $\det(\lambda B+H)$. Let $M$ be its gradient, $M:K=\frac d{dt}\mathcal D_B(H+tK)|_{t=0}$, at $H=H(u)$. Then:
1. $M=C(B+H)-C(B)-C(H)$, the polarized cofactor of the two Hessians $H(|\Phi|^2)$ and $H(u)$. It is conserved.
2. $M_{TT}=C(B_{\rm fib}):H_{\rm fib}(u)=L_Bu=\mu$.
3. $M:H(u)=2\mathcal D_B(H(u))$.
4. $\nu^*M\nu=0$ on the wall $|W|=w_0$, where $\nu=G(r_W)/|\nabla r_W|$.

*Proof.*
1. **Polarization and conservation.** For $3\times3$ matrices $C$ is quadratic. With the symmetric trilinear polarization $\mathfrak d$ of the Moore determinant, $B:C(H+tK)$ has $t$-derivative $6\mathfrak d(B,H,K)$, and $6\mathfrak d(B,H,\cdot)$ is the polarization $C(B+H)-C(B)-C(H)$. Its divergence vanishes by the proof of [12](12_QUATERNIONIC_SLOT_CALCULUS.md) Theorem 12.1, whose step 3 treats multilinear inputs from different functions. This is checked exactly in Part P.
2. **The $TT$ entry.** $M_{TT}=\frac d{dt}\mathcal D_B(H+te_Te_T^*)$. This is the $\lambda$-coefficient of $e_T^*C(\lambda B+H)e_T=\det((\lambda B+H)|_{e_T^\perp})=\det(\lambda B_{\rm fib}+H_{\rm fib})$ (Lemma 2.1.3 of 12), and for $2\times2$ matrices that coefficient is $C(B_{\rm fib}):H_{\rm fib}$.
3. **Euler.** $\mathcal D_B$ is homogeneous of degree $2$.
4. **The wall.**
   - $r_W=|W|^2/w_0^2-1$ with $W$ a single holomorphic function, so $H(r_W)=4w_0^{-2}(\partial W)^*\partial W$ has rank one. Its image is spanned by $(\partial W)^*\parallel G(r_W)\parallel\nu$.
   - Since $u=0$ on the wall, compression (12 Corollary 4.1) gives $H(u)|_{\nu^\perp}=\rho H(r_W)|_{\nu^\perp}=0$.
   - Hence $\nu^*M\nu$, which is the $\lambda$-coefficient of $\det((\lambda B+H(u))|_{\nu^\perp})=\lambda^2\det(B|_{\nu^\perp})$, vanishes. ∎

**Lemma 3.2 (the completed square).**
- **Set-up.** Let $E=\ker B$, a quaternionic plane, and $\beta=\operatorname{tr}B$, $\beta_f=\operatorname{tr}B_{\rm fib}=4|c|^2$. Take the orthonormal basis of $E$ given by $\varepsilon_0=\pi_Ee_T/|\pi_Ee_T|$ and $\varepsilon_1=(0,e/|e|)$. Write $\alpha=A|_E$ and $b=H(u)|_E$.
- **Claim.** If $L_Bu=\mu$, then
$$A:H(u)-\mathcal D_B(H(u))=\beta\,\big|b_{01}+\alpha_{01}/\beta\big|^2+\det(\alpha)/\beta\ \ge0 .$$

*Proof.*
- **The basis.** $\varepsilon_1\in E$ because $(\partial\Phi)(0,e)=c\cdot e=0$, and $\varepsilon_0\perp\varepsilon_1$ because $e_T$ has no spatial component.
- **Support of $A$.** $A$ is supported on $E$ by tangency, so $A:H(u)=\alpha:b$.
- **The determinant.** With $B=\beta vv^*$ for a unit $v$, $\mathcal D_B(H)=\beta\,v^*C(H)v=\beta\det(H|_E)=\beta\det b$ (Lemma 2.1.3 of 12).
- **Relating $\alpha_{00}$ to $b_{11}$.**
  - $|\pi_Ee_T|^2=1-|v_T|^2=\beta_f/\beta$, so $\mu=A_{TT}=(\beta_f/\beta)\alpha_{00}$.
  - $\mu=L_Bu=4e^*H_{\rm fib}e=4|e|^2b_{11}=\beta_fb_{11}$.
  - Hence $\alpha_{00}=\beta b_{11}$.
- **Completing the square.** Expanding,
$$\alpha:b-\beta\det b=\alpha_{11}b_{11}+2\operatorname{Re}(\alpha_{01}\overline{b_{01}})+\beta|b_{01}|^2 .$$
  This equals the right side, using $\alpha_{11}b_{11}=(\det\alpha+|\alpha_{01}|^2)/\beta$. Finally $\det\alpha\ge0$ since $\alpha\ge0$. ∎

**Proposition 3.3 (planar energy balance).** Let $X_0(z)=\int_{\Omega^{\rm pl}_z}(-u)\mu\,d\zeta$ and $\mathcal F_{\rm pl}=\frac12X_0$. Then $\mathcal F_{\rm pl}\ge0$ is smooth on $a<|\tau|<4$, and
$$\Delta_z\mathcal F_{\rm pl}\le\mathcal B_{\rm pl}(z):=\int\rho\,G(r_W)^*A\,G(r_W)\,\delta(r_W)\,d\zeta,\qquad\rho=u/r_W\ge0 .$$

*Proof.* For a conserved $M'$ and a test function $\eta(z)$, Green's identity on the total domain (08 (6.5)) gives
$$\int u\,M'_{TT}\,\Delta_z\eta=\int\eta\,M':H(u)-\int_{|W|=w_0}\eta\,(\nu^*M'\nu)\,\partial_{\mathbf n}u\,dS .$$
There is no contribution from $|Q|=\frac1{10}$, since $u\equiv0$ near it. Apply this twice.
- **With $M'=M$** (Lemma 3.1): $M_{TT}=\mu$, $M:H(u)=2\mathcal D_B$ and $\nu^*M\nu=0$. So $-\Delta_zX_0=2\int\mathcal D_B(H(u))$.
- **With $M'=A$:** $-\Delta_zX_0=\int A:H(u)-\mathcal B_{\rm pl}$, where $\partial_{\mathbf n}u\,dS=\rho\,|\nabla r_W|^2\delta(r_W)\,dz\,d\zeta$.
- **Combining.** By Lemma 3.2, $\int A:H(u)\ge\int\mathcal D_B=-\frac12\Delta_zX_0$. Hence $-\Delta_zX_0\ge-\frac12\Delta_zX_0-\mathcal B_{\rm pl}$, that is $\Delta_zX_0\le2\mathcal B_{\rm pl}$.
- **Smoothness and sign.** Smoothness comes from Lemma 2.2. Also $u\le0$ and $\mu\ge0$. ∎

## 4. The window estimates

**Proposition 4.1.** Fix $\gamma$ with $8C_2\gamma^2\le c_0'$, where $c_0'=c_0w_0^2/4$, and then $K$ large. For $h<h_0$, $a=Kh$, cap data with $\|f\|_2=1$ and all $T$:
$$\mathcal F_{\rm pl}(0,\tau)\le Ca^2\|f\|_4^4\ (|\tau|=4a),\qquad\mathcal F_{\rm pl}(T,\tau)\ge ca^2\|g_T\|_4^4\ (|\tau|=2a),\qquad\Delta_z\mathcal F_{\rm pl}\le Cs^{-2}e^{-c_0's/h}\mathcal F_{\rm pl}\ (2a\le s\le3).$$

*Proof.* This follows [11](11_FRAMEWORK_POSITIVITY_ENGINE.md) Proposition 4.1 with $d=1$; the comparisons are now linear.
1. **Current bounds.** Put $P_*=J_{Ls}*|g_T|^2$, so $|\nabla_X^jP_*|\le C(Ls)^{-j}P_*$ and $\int P_*^2\le\|g_T\|_4^4$. As in 11 (4.1):
   - $\|A\|\le Cs^{-2}\sigma^mP_*$ on the fibre, with an extra factor $e^{-2c_0'\sigma}$ on the wall, where $|W|=w_0$ and $|v|\le w_0/3$;
   - $\mu\le Cs^{-2}\sigma^mP_*$.
2. **Upper comparison.**
   - $r_W$ depends only on $W$, $P_*$ only on $X$, and (2.1) has no mixed terms. So
$$L_B(r_WP_*)=P_*L_Br_W+r_WL_BP_*\ge\big(16|\mathcal J|^2w_0^{-2}-CL^{-2}s^{-2}\big)P_*\ge cs^{-2}P_*$$
     for a fixed large $L$.
   - For large $C_1,m_1$ the function $v=C_1\sigma^{m_1}r_WP_*-u$ has $L_Bv>0$ and vanishes on the wall. So $v\le0$, that is $0\le\rho\le C_1\sigma^{m_1}P_*$.
   - The volume of $\{|W|<w_0,|Q|<\frac1{10}\}$ in the variables $P$ is $\asymp a^2s^2$: the Jacobian is $|\tau B'(W)a|^2$. So
$$\mathcal F_{\rm pl}\le\tfrac12\int C_1\sigma^{m_1}P_*\cdot Cs^{-2}\sigma^mP_*\,dX\,dP\le Ca^2\sigma^{m'}\|g_T\|_4^4 .$$
3. **The flux.** $A$ is tangent and $r_W$ does not involve $Q$. Also $|G(r_W)|\lesssim s^{-1}$ by (1.2), and $\int\delta(r_W)\,dP\lesssim a^2s^2$. With step 1:
$$\mathcal B_{\rm pl}\le C(a/s)^2\sigma^{m''}e^{-2c_0'\sigma}\|g_T\|_4^4 .$$
4. **Lower bound.** Take a centre $(X,P)$ with $|W|\le\gamma$ and $|Q|\le\frac1{400}$, and the ball of radius $h$ about it. For $K$ large the ball lies in $\{|W|\le2\gamma,\ |Q|\le\frac1{100}\}$.
   - **Density.** Inverting the time row as in 11 Step 5 gives $\mu\ge m_X:=c_\gamma(hs)^{-1}\sigma^{-m}e^{-2C_2\gamma^2\sigma}|g_T(X)|^2$ on the whole ball.
   - **Depth.** Let $l$ be the distance to the centre. Then $H_{\rm fib}(l^2)=6I$, so $L_B(l^2)=24|e|^2\le C$. The function $u-\epsilon m_X(l^2-h^2)$ has $L_B>0$ for small $\epsilon$ and is $\le0$ on the sphere, so it is $\le0$ inside. At the centre $(-u)\mu\ge\epsilon h^2m_X^2$.
   - **Integration.** The centres have imaginary volume $\asymp a^2s^2\gamma^2$, and $h^2s^2(hs)^{-2}=1$. So
$$\mathcal F_{\rm pl}\ge c_\gamma a^2\sigma^{-2m}e^{-4C_2\gamma^2\sigma}\|g_T\|_4^4 .$$
5. **Conclusion.**
   - The ratio of 3 to 4 is $\le C_\gamma s^{-2}\sigma^{M}e^{-(2c_0'-4C_2\gamma^2)\sigma}\le C_\gamma s^{-2}\sigma^Me^{-\frac32c_0'\sigma}\le Cs^{-2}e^{-c_0'\sigma}$.
   - At $|\tau|=4a$ ($\sigma=4K$, $T=0$), step 2 gives the first bound; at $|\tau|=2a$ ($\sigma=2K$), step 4 gives the second. ∎

## 5. Proof of Theorem 13.1 and Corollary 13.2

By Proposition 4.1, $E=\mathcal F_{\rm pl}$ satisfies the hypotheses of 08 Lemma 8.1 with $b=2a$ and $V(s)=Cs^{-2}e^{-c_0's/h}$; the integrability is shown in 11 §5. With $r=2b$:
$$Ca^2\|f\|_4^4\ge H(0,2b)\ge c\,\frac{\log 2}{\log^2(3/b)}\int_1^2H(T,b)\,dT\ge\frac{c\,a^2}{\log^2(2/h)}\int_1^2\|g_T\|_4^4\,dT .$$
This is (0.1). The reduction in §0 gives Theorem 13.1 and Corollary 13.2. ∎

## 6. What the planar case shows

1. **The whole architecture survives with only linear PDE.** The current, the passive coordinate, the background, the fibre, the wall flux, the window, and the two-circle propagation are all present. The nonlinear fibre equation $\det(\lambda B_{\rm fib}+H_{\rm fib}(u))=\lambda\mu$ is replaced by its $\lambda$-coefficient $C(B_{\rm fib}):H_{\rm fib}(u)=\mu$, which is linear because $n=2$.
2. **The price is an inequality.** Lemma 3.2's completed square replaces the exact balance of 08 Proposition 6.2. It holds because $\mu$ enters both sides, through the current and through the equation.
3. **Why this does not extend to $n\ge3$.** With $B_{\rm fib}=\beta_fww^*$, the rank-one expansion gives $\det(\lambda B_{\rm fib}+H_{\rm fib})=\det H_{\rm fib}+\lambda\beta_f\det(H_{\rm fib}|_{w^\perp})$. The $\lambda$-coefficient is the Monge–Ampère operator in the $d=n-1$ free directions, of degree $d$. It is linear only for $d=1$, so App. D's simplification is special to the plane. The same degree appears in the density-to-depth estimate $-u\gtrsim h^2\mu^{1/d}$, which is what selects $p_c=2n/(n-1)$.
4. **A five-page textbook proof** would consist of §1–§5 with the quoted lemmas inlined:
   - the 3-slot calculus for $2\times2$ and $3\times3$ matrices (12 §1–§4);
   - the lifted current (08 Lemmas 3.1–3.2);
   - the multiplier lemma (08 Lemma 7.2);
   - two-circle propagation (08 Lemma 8.1).

   I recommend it as the entry point to the paper.

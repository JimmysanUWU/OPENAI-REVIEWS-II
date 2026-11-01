# 16. The null frame: the simplest engine, the Husimi current, and a direct engine for Schrödinger

*Claude, 29 October 2026. Roadmap target `T.KAK.12` (done); nodes `KAK.NULL`, `KAK.HUS`, `KAK.SCHR` in [`ROADMAP.md`](../../../ROADMAP.md).*

The paper builds its current in light-cone ("null") coordinates (08 §3, paper §2), then rotates two of them to the physical time $T$ before running the proof. This file keeps the null coordinates and runs the engine of [11](11_FRAMEWORK_POSITIVITY_ENGINE.md) there. Three things follow.

> **Theorem 16.3 (the null frame carries engine data).** Keep the current of 08 Proposition 3.3 in the null slots $(t,x_1,\dots,x_d,s_0)$, take $t+Y\mathrm i$ as the parameter slot, and use the chart
> $$W=\frac PY,\qquad\Phi=S_1+\frac{P\cdot P}{2Y},\qquad\text{inverse: }P=YW,\quad S_1=aQ-\tfrac12Y\,W\cdot W .$$
> This is an engine (11 Definition 1.4) on the cap of the null symbol $\varphi(\xi)=|\xi_x|^2/(2\xi_{s_0})$. The inverse chart is polynomial with Jacobian exactly $aY^d$, and the time row is a single Gaussian multiplier. Consequently, for $p=p_c=2n/(n-1)$ and cap data $f$ on $\mathbb T^n$,
> $$\int_1^2\|e^{-\mathrm it\varphi(D)}f\|_{L^p(\mathbb T^n)}^p\,dt\le C\,[\log(2/h)]^2\,\|f\|_{L^p(\mathbb T^n)}^p .$$

> **Proposition 16.5 (the time component is the continued wave).** At every point with $Y\ne0$, writing $y=|Y|$,
> $$A^0_{tt}=y^{-d}\,\big|f_0\big(\mathsf X-\mathrm iq\big)\big|^2,\qquad q=y\big(1,\ \operatorname{Re}W,\ -\tfrac12|W|^2\big),$$
> where $f_0=\sum\kappa^{d/2}u_{\kappa,v}$ is the wave. The imaginary displacement $q$ lies in the closed forward light cone, and on it when $W$ is real. For Schrödinger data this is the **Husimi density**: the time component of the positive conserved current is a Gaussian phase-space density of the solution.

> **Theorem 16.8 (a direct engine for the paraboloid; `T.KAK.12`).** Let $d\ge1$ and $p=2(d+1)/d$. The Schrödinger flow is the single-mass sector of the null-frame wave, and the engine passes to it exactly. Hence:
> 1. **(Torus.)** For real $\kappa\ge1$ and $\psi$ on $\mathbb T^d$ with Fourier support in $\{|\xi|\le2\kappa\}$:
> $$\int_{1/(2\kappa)}^{1/\kappa}\|e^{\mathrm is\Delta}\psi\|_{L^p(\mathbb T^d)}^p\,ds\le C\,\kappa^{-1}[\log(2+\kappa)]^2\,\|\psi\|_{L^p(\mathbb T^d)}^p .$$
> 2. **(Euclidean.)** The same holds on $\mathbb R^d$, for every real $\kappa\ge1$ and $\psi$ with $\operatorname{supp}\hat\psi\subset B(0,\kappa)$.
> 3. **(Paper (12.2).)** For a compact interval $I$, $N\ge2$ and $\operatorname{supp}\hat\psi\subset B(0,N)$: $\|e^{\mathrm it\Delta}\psi\|_{L^p(\mathbb R^d\times I)}\le C_I[\log(2+N)]^{2/p}\|\psi\|_p$.

**What is new here.**
- **`T.KAK.12` is resolved.** Its question was whether the paraboloid has an engine of its own, instead of inheriting its estimate from the cone by parabolic rescaling. It does, and no limit is involved.
  - The free Schrödinger equation is the null reduction of the wave equation one dimension up (the Eisenhart–Duval lift).
  - The paper's current survives the reduction exactly: on data of one mass $\kappa$ it does not depend on $s_0$, and it is built from the Schrödinger solution alone.
  - Theorem 16.8(3) is the paper's (12.2). [10](10_GAPLESS_REMAINING.md) §4 reached it through the sphere's local extension estimate, parabolic rescaling and Fresnel's identity; here it comes straight from the engine.
- **The simplest engine.** In the null frame:
  - the chart is polynomial with Jacobian exactly $aY^d$;
  - $y=s$ exactly;
  - the time row is one Gaussian multiplier, so (A4) needs no splitting.

  The paper's instance is this one after a rotation of two slots, and the pencil instances of 11 deform its symbol. It is the instance to learn first.
- **The density has a classical name.** The time component of the current is the solution continued into the complex domain along an imaginary light ray (Kaiser's "physical wavelets"). For Schrödinger it is the Husimi density.
  - [07](07_TOOLKIT.md) K1 says that neither Husimi nor Wigner densities are both positive and exactly conserved as phase-space densities. Proposition 16.5 sharpens this.
  - A Husimi density *is* the time component of a positive conserved tensor once the window width $|Y|$ and a complexified momentum $P/Y$ are allowed to vary.
- **One logarithm is necessary, numerically.** For Schrödinger data focusing inside the window, the normalized quantity grows by $\frac32\log2$ per doubling of $\kappa$, i.e. exactly one logarithm (§6). The theorem allows two.

**Checks.** [`scripts/verify_null_frame.py`](scripts/verify_null_frame.py) ([`results/null_frame.json`](results/null_frame.json), ~15 s):
- **Part G (exact, $d=1,2,3$).**
  - The inverse chart and its Jacobian $aY^d$.
  - $\partial\Phi=\ell(W)$ and $\Phi\circ\Theta=\eta$.
  - $(\partial W_j)V^v=(v_j-W_j)/(2Y)$ and tangency $\ell(W)V^v=0$.
  - The displacement $y(1,W,-\frac12W\cdot W)$ is null.
- **Part H (50 digits).** Proposition 16.5 for random superpositions, one mass or several, $d=1,2,3$, to $10^{-50}$. The wrong sign of the shift fails at order $10^{-1}$.
- **Part A.** The row identities behind (A3)–(A4), to $10^{-50}$. The time row is one real multiplier.
- **Part C (60 digits).**
  - Conservation on one-mass (Schrödinger) data, $d=1,2,3$, to $10^{-39}$ relative.
  - The tensor does not depend on $s_0$.
  - Controls: two masses make it depend on $s_0$; dropping the rotation $\iota_e$ breaks conservation.
- **Part S (exact).**
  - The null system on $e^{\mathrm i\kappa s_0}F$ is the Schrödinger equation.
  - The Galilean identity.
  - The exponents of Theorem 16.8.
- **Part W (exact).** Why a lift is needed: the Gaussian identity $b\cdot\nabla_aG=-\frac\beta2\nabla_a\cdot\nabla_bG$, and the Husimi flux identity of §4 on a two-mode Schrödinger solution. The control without the $\nabla_\xi H$ term fails.
- **Part F (numerical, $d=1$, $p=4$).** The focusing example of §6. The control, focusing outside the window, stays bounded.

---

## 0. Conventions

Notation is that of [08](08_GAPLESS_PROOF.md) §1–§3 and [11](11_FRAMEWORK_POSITIVITY_ENGINE.md) §0: $n\ge2$, $d=n-1$, $p=p_c=2n/d$, $q=n/d$, the quaternion units and slots, $\iota_e(\alpha+\mathrm i\beta)=\alpha+e\mathrm i\beta$, and $\ell(w)=(-\frac12w\cdot w,\,w,\,1)$ with bilinear dots.

**Null coordinates (08 §3.1).**
- Coordinates are $(t,x,s_0)\in\mathbb R\times\mathbb R^d\times\mathbb R$.
- The system is $\mathscr P=C^t\partial_t+\sum_jC^{x_j}\partial_{x_j}+C^{s_0}\partial_{s_0}$ with $C^t=2E_{00}$, $C^{x_j}=E_{0j}+E_{j0}$ and $C^{s_0}=-\sum_jE_{jj}$. Its first row applied to $(\partial_{s_0}\psi,\nabla_x\psi)$ is $(2\partial_t\partial_{s_0}+\Delta_x)\psi$, the wave operator in these coordinates.
- A plane wave is $u_{\kappa,v}=c\,e^{\mathrm i\kappa(s_0+v\cdot x-\frac12|v|^2t)}$ with $\kappa>0$ and $v\in\mathbb R^d$.
- On the slice $\{t=\text{const}\}$ it has frequency $\xi=(\kappa v,\kappa)$, and it evolves by $e^{-\mathrm it\varphi(\xi)}$ with the **null symbol**
$$\varphi(\xi)=\frac{|\xi_x|^2}{2\xi_{s_0}}\qquad(\xi_{s_0}>0),$$
  which is homogeneous of degree one. Its polarization is $v(\xi)=\xi_x/\xi_{s_0}$.

**The current (08 Proposition 3.3, before the rotation of its item 4).** Slots $t+Y\mathrm i$, $x_j+P_j\mathrm i$, $s_0+S_1\mathrm i$, with $Y,P_j,S_1\in\mathbb C_{\mathrm k}$. On $\{Y\ne0\}$ put $y=|Y|$, $e=Y/y$, $W=P/Y$, and
$$\chi=\sum_{\kappa,v}V^v\,\iota_e\big(\overline{u_{\kappa,v}}\big)\,m_{\kappa,v},\qquad V^v=\big(1,\ \tfrac{v+W}2,\ -\tfrac{W\cdot v}2\big)^{\mathsf T},\qquad m_{\kappa,v}=\Big(\frac\kappa y\Big)^{d/2}e^{-\frac{\kappa y}2|W-v|^2}.$$
Then $A^0=\chi\chi^*$ is positive semidefinite and conserved, and $(\partial\Phi)\chi=0$ for $\Phi=S_1+P\cdot P/(2Y)$ (08 Prop. 3.3, items 1–3).

**Periods.** $\mathbb T^n=\mathbb R^n/(100\mathbb Z)^n$ in the variables $(x,s_0)$. Data on $\mathbb T^n$ are finite Fourier sums with frequencies in $(2\pi/100)\mathbb Z^n$.

---

## 1. The null chart

The engine's abstract symbols become:

| Engine (11) | Null frame |
|---|---|
| parameter slot $T+\tau\mathrm i$, $s=\lvert\tau\rvert$ | $t+Y\mathrm i$, $s=y=\lvert Y\rvert$ |
| fibre slots $X_j+P_j\mathrm i$ | $x_j+P_j\mathrm i$ ($j\le d$) and $s_0+S_1\mathrm i$ |
| fibre variables $X\in\mathbb T^n$, $P\in\mathbb C_{\mathrm k}^n$ | $(x,s_0)\in\mathbb T^n$, $\mathcal P=(P,S_1)\in\mathbb C_{\mathrm k}^n$ |
| $g_T=e^{-\mathrm iT\varphi(D)}f$ | $g_t=e^{-\mathrm it\varphi(D)}f$, null symbol |
| polarization $v(\eta)$ | $\eta_x/\eta_{s_0}$ |
| direction $W$, passive $\Phi$, $Q=\Phi/a$ | $W=P/Y$, $\Phi=S_1+P\cdot P/(2Y)$ |

**The cap.** Let $w_0=\frac1{10}$, and
$$K'=\{\eta:\tfrac14\le\eta_{s_0}\le4,\ |\eta_x|\le\tfrac{w_0}3\eta_{s_0}\},\qquad K=\{\eta:\tfrac12\le\eta_{s_0}\le2,\ |\eta_x|\le\tfrac{w_0}4\eta_{s_0}\}.$$
Then $(\varphi,K',v)$ is a cap (11 Definition 1.1) with $|v|\le w_0/3$ on $K'$.

**Lemma 16.1 (the null chart satisfies (G1)–(G2)).** For $Y\ne0$ the map $(P,S_1)\mapsto(W,Q)$ is a biholomorphism of $\mathbb C_{\mathrm k}^n$, with inverse $P=YW$, $S_1=aQ-\frac12Y\,W\cdot W$. On $\mathcal W=\{a<|Y|<4,\ (W,Q)\in\bar E\}$:
1. $\partial W=(-W/Y,\ Y^{-1}I_d,\ 0)$ and $\partial\Phi=\ell(W)$, in the variables $(Y,P,S_1)$;
2. $|\partial W|\le2/s$, $|\partial\Phi|\le2$ and $|\partial_{\mathcal P}\Phi|\ge1$;
3. $|\partial\mathcal P/\partial W|\le2s$, $\partial\mathcal P/\partial Q=(0,a)$, and $\det_{\mathbb C}\partial\mathcal P/\partial(W,Q)=aY^d$.

So (G1)–(G2) hold with $C_0=2$, for every $a\le\frac14$.

*Proof.* The inverse is explicit (Part G). $\partial W$ and $\partial\Phi$ are direct derivatives; for instance $\partial_Y\Phi=-P\cdot P/(2Y^2)=-\frac12W\cdot W$. The bounds use $|W|\le w_0$ on $\bar E$.
- $|\partial\Phi|^2\le w_0^4/4+w_0^2+1$.
- $\partial_{S_1}\Phi=1$.
- $\partial\mathcal P/\partial W$ has blocks $YI_d$ and $-YW^{\mathsf T}$.
- The Jacobian matrix is block lower-triangular with diagonal $(Y,\dots,Y,a)$. ∎

**Remark (why this chart is simpler).** The paper's chart (08 (4.2)) is $Y=(\sqrt2\tau+aQ)/(1+W\cdot W/2)$, $P=WY$, $Z=\sqrt2Y-\tau$. That is the null chart after the rotation $(t,s_0)\to(T,b)$, which mixes the parameter slot with a fibre slot. The rotation creates the rational factor and replaces $y=s$ by $y\asymp s$.

---

## 2. The rows

**Lemma 16.2 (the null current satisfies (A1)–(A4)).** For $\gamma$-cap data $f$ at scale $h$ (11 Definition 1.1), with $g_t=e^{-\mathrm it\varphi(D)}f$:
1. **(A1), (A2).** $A^0=\chi\chi^*$ is positive and conserved on $\{Y\ne0\}$, and $(\partial\Phi)\chi=0$.
2. **(A3), $L=2$.** For every row $\mu$, $|\chi_\mu|^2=|m'_\mu(hD)g_t|^2+|m''_\mu(hD)g_t|^2$, where
$$m'_\mu,\ m''_\mu(\eta)=\big(\operatorname{Re}_{\mathrm k},\operatorname{Im}_{\mathrm k}\big)\,V^{v(\eta)}_\mu\cdot\mathfrak m(\eta),\qquad\mathfrak m(\eta)=\Big(\frac{\eta_{s_0}}{hy}\Big)^{d/2}\exp\Big(-\frac{\eta_{s_0}y}{2h}\,|W-v(\eta)|^2\Big).$$
   They satisfy (A3) with $c_0=\frac1{16}$ and $m_1=4n+8$.
3. **(A4).** The time row is one real multiplier: $|\chi_t|^2=|\mathfrak m(hD)g_t|^2$. On $|W|,|v|\le2\gamma$,
$$|\partial_\eta^\beta(1/\mathfrak m)|\le C(hs)^{d/2}\sigma^{m_1}e^{32\gamma^2\sigma}\qquad(|\beta|\le2n+4),$$
   so $C_2=32$, independent of $\gamma$.

The constants depend only on $n$.

*Proof.*
1. This is 08 Proposition 3.3, items 1–2, which are stated and proved in the null slots.
2. **The splitting.**
   - Write $V^v_\mu=b'+\mathrm kb''$ with $b',b''$ real.
   - Then $V^v_\mu\iota_e(\bar u)=\iota_e(\overline{b'u})+\mathrm k\,\iota_e(\overline{b''u})$.
   - The two terms lie in the orthogonal planes $\operatorname{span}(1,e\mathrm i)$ and $\operatorname{span}(\mathrm k,\mathrm ke\mathrm i)$ (08 §3.5), and $|\iota_e(\alpha)|=|\alpha|$.
   - With $\kappa=\eta_{s_0}/h$ and $v=v(\eta)$, $m_{\kappa,v}=\mathfrak m(\eta)$ (Part A).
3. **Bounds for (A3).** Expand $E(\eta)=\frac{\eta_{s_0}y}{2h}|W-v(\eta)|^2=\frac\sigma2\big(\eta_{s_0}|W|^2-2\operatorname{Re}W\cdot\eta_x+|\eta_x|^2/\eta_{s_0}\big)$, using $y=s$. Then:
   - first derivatives of $E$ are $O(\sigma|W-v|)$, since $\partial_{\eta_x}E=\sigma(v-\operatorname{Re}W)$ and $\partial_{\eta_{s_0}}E=\frac\sigma2(|W|^2-|v|^2)$;
   - higher derivatives are $O(\sigma)$.

   Hence $|\partial^\beta e^{-E}|\le C\sigma^{|\beta|}(1+\sigma|W-v|^2)^{|\beta|}e^{-E}$. Since $E\ge\frac\sigma8|W-v|^2$ on $K'$, absorbing the polynomial factor gives $e^{-\sigma|W-v|^2/16}$. The prefactor is $\le C(hs)^{-d/2}$, and $V^v_\mu$ is a bounded polynomial.
4. **(A4).** $V^v_t=1$, so $\chi_t=\iota_e\big(\overline{\sum m_{\kappa,v}u_{\kappa,v}}\big)$ with no $\mathrm k$-part. On $|W|,|v|\le2\gamma$ we have $E\le\frac\sigma2\cdot4\cdot16\gamma^2=32\gamma^2\sigma$, and derivatives of $e^{E}$ cost $\sigma^{2|\beta|}$. ∎

---

## 3. The null-frame engine

**Proof of Theorem 16.3.** Lemmas 16.1–16.2 give an engine on $(\varphi,K',v)$ in the sense of 11 Definition 1.4, with constants depending only on $n$. Theorem 11.1 gives the estimate for $\gamma$-cap data at scale $h$. ∎

**Corollary 16.4 (critical estimate on null slabs).** For $0<h\le1$ and $f$ on $\mathbb T^n$ with Fourier support in $\{\xi:h\xi\in K,\ |\xi_x|\le\gamma\xi_{s_0}\}$,
$$\int_1^2\|e^{-\mathrm it\varphi(D)}f\|_{L^{p_c}(\mathbb T^n)}^{p_c}\,dt\le C\,[\log(2/h)]^2\,\|f\|_{L^{p_c}(\mathbb T^n)}^{p_c}.$$
The left side is the $L^{p_c}$ norm of the wave over the region between two parallel null hyperplanes $\{t=1\}$ and $\{t=2\}$, and the data are given on the null hyperplane $\{t=0\}$.

**Remarks.**
1. **The physical frame is a rotation of this one.** 08 rotates $(t,s_0)\to(T,b)$ to bound the wave on time slabs. Neither estimate implies the other directly: their slabs and their data hyperplanes differ. Theorem 1.2 (time slabs) is the one the tower needs. Theorem 16.8 (Schrödinger) needs the null one.
2. **Periods.** The proof of Theorem 11.1 uses the period only through:
   - lattice-point counts (11 Lemma 3.1.4 and the large-$h$ case);
   - periodized kernels (08 Lemma 7.2);
   - the nearest-translate bound of 11 Step 1.

   These are uniform for rectangular periods in $[100,200]$. Theorem 16.3 is used below only with period $100$.

---

## 4. The time component is the continued wave

**Proposition 16.5.** Let $f_0=\sum_{\kappa,v}\kappa^{d/2}u_{\kappa,v}$, a finite sum with $\kappa>0$, extended holomorphically to complex $(t,x,s_0)$. At every point with $Y\ne0$,
$$A^0_{tt}=|\chi_t|^2=y^{-d}\,\Big|f_0\big(t-\mathrm iy,\ x-\mathrm iy\operatorname{Re}_{\mathrm k}W,\ s_0+\tfrac{\mathrm i}2y|W|^2\big)\Big|^2 .$$
Here $|W|^2=|\operatorname{Re}_{\mathrm k}W|^2+|\operatorname{Im}_{\mathrm k}W|^2$, and $\mathrm i$ is the Fourier unit.

The imaginary displacement $q=y(1,\operatorname{Re}W,-\frac12|W|^2)$ satisfies $2q^tq^{s_0}+|q^x|^2=-y^2|\operatorname{Im}W|^2\le0$ and $q^t>0$. So it lies in the closed forward light cone of the metric $2\,dt\,ds_0+|dx|^2$, and on the cone exactly when $W$ is real.

*Proof.* For one mode, the phase $\kappa(s_0+v\cdot x-\frac12|v|^2t)$ at the displaced point gains the real part
$$\kappa y\big(-\tfrac12|W|^2+v\cdot\operatorname{Re}W-\tfrac12|v|^2\big)=-\tfrac{\kappa y}2\big(|\operatorname{Re}W-v|^2+|\operatorname{Im}W|^2\big)=-\tfrac{\kappa y}2|W-v|^2 .$$
So $u_{\kappa,v}$ at the displaced point equals $u_{\kappa,v}\,e^{-\frac{\kappa y}2|W-v|^2}$ at the real point: the same complex number times a positive factor. Summing, $y^{-d/2}f_0(\mathsf X-\mathrm iq)=\sum m_{\kappa,v}u_{\kappa,v}$, and $|\chi_t|=|\sum m_{\kappa,v}u_{\kappa,v}|$ by Lemma 16.2(3). The cone computation is direct. ∎

**Corollary 16.6 (the Schrödinger case: Husimi densities).** For data of one mass $\kappa$, $f_0=\kappa^{d/2}e^{\mathrm i\kappa s_0}\psi(t,x)$ with $\partial_t\psi=\frac{\mathrm i}{2\kappa}\Delta\psi$ (Lemma 16.7), and
$$A^0_{tt}=\Big(\frac\kappa y\Big)^de^{-\kappa y|W|^2}\,\big|\psi(t-\mathrm iy,\ x-\mathrm iy\operatorname{Re}W)\big|^2 .$$
For real $W$ this equals $(\kappa/y)^d\,\big|e^{-\frac y{2\kappa}|D-\kappa W|^2}\psi(t)\big|^2(x)$. That is the Husimi (Gabor) density of $\psi(t)$ at position $x$ and momentum $\kappa W$, with spatial window width $\sqrt{y/\kappa}$:
- the complex time $t-\mathrm iy$ is heat smoothing for time $y/(2\kappa)$;
- the complex position $x-\mathrm iyW$ shifts the frequency window to $\kappa W$.

*Proof.* Substitute into Proposition 16.5. For real $W$, $e^{-\kappa yW^2/2}\psi(t-\mathrm iy,x-\mathrm iyW)=\sum_vc_ve^{\mathrm i\kappa(v\cdot x-\frac12|v|^2t)}e^{-\frac{\kappa y}2|v-W|^2}$, and $\frac{\kappa y}2|v-W|^2=\frac y{2\kappa}|\xi-\kappa W|^2$ with $\xi=\kappa v$. ∎

**What this says.**
1. **The density is classical.** It is the holomorphic extension of a positive-frequency wave to the forward tube, read along imaginary light rays. Kaiser's "physical wavelets" use the same extension, with the imaginary part as scale and direction.
   - The size $y$ of the imaginary displacement sets the resolution.
   - Its direction $W$ selects the wave direction.
   - For Schrödinger, this is the Bargmann–Husimi picture.
2. **What the lift adds.** A Husimi density $H=\mathcal W*G$ at fixed width (the Wigner function $\mathcal W$ smoothed by a phase-space Gaussian $G=e^{-|a|^2/\alpha-|b|^2/\beta}$) is positive and conserved as a total mass. Free transport $\partial_t\mathcal W+\frac\xi\kappa\cdot\nabla_x\mathcal W=0$ does not commute with the smoothing, and since $b\cdot\nabla_aG=-\frac\beta2\nabla_a\cdot\nabla_bG$,
$$\partial_tH+\nabla_x\cdot\Big(\frac\xi\kappa H+\frac\beta{2\kappa}\nabla_\xi H\Big)=0 .$$
   (Part W checks this with the Gabor window of width $w$, for which $\beta=1/w^2$.) The flux has no sign, so $H$ is not the time component of a positive conserved *tensor* in $(t,x)$. Proposition 16.5 shows that it becomes one once:
   - the width $y=|Y|$ and a complexified momentum $W=P/Y$ become coordinates;
   - the space and $s_0$ rows of $\chi$ are added. These carry the weights $\frac{v+W}2$ (the average of the mode's velocity and the window's centre) and $-\frac{W\cdot v}2$.

   This is the precise sense of 03 M1–M2: "imaginary displacement along light rays" and "square after summing in lifted variables".

---

## 5. The null reduction: Schrödinger

**Lemma 16.7 (null reduction).** Let $F(t,x,s_0)=e^{\mathrm i\kappa s_0}\big(\psi,\ \nabla_x\psi/(\mathrm i\kappa)\big)$. Then rows $1,\dots,d$ of $\mathscr PF$ vanish identically, and
$$(\mathscr PF)_0=e^{\mathrm i\kappa s_0}\Big(2\partial_t\psi-\frac{\mathrm i}\kappa\Delta\psi\Big).$$
So $\mathscr PF=0$ if and only if $\partial_t\psi=\frac{\mathrm i}{2\kappa}\Delta\psi$, i.e. $\psi(t)=e^{\mathrm i\frac t{2\kappa}\Delta}\psi(0)$. The plane waves $u_{\kappa,v}$ with fixed $\kappa$ are $e^{\mathrm i\kappa s_0}$ times the Schrödinger plane waves $e^{\mathrm i\kappa(v\cdot x-\frac12|v|^2t)}$, and $f_{\rm sys}=\kappa^{d/2}F$ for such data.

*Proof.* Row $j\ge1$ is $\partial_{x_j}F_0-\partial_{s_0}F_j=e^{\mathrm i\kappa s_0}(\partial_j\psi-\mathrm i\kappa\,\partial_j\psi/(\mathrm i\kappa))=0$. Row $0$ is $2\partial_tF_0+\sum_j\partial_{x_j}F_j$. Exact check in Part S. ∎

**The reduced current.** For one-mass data, $\overline{u_{\kappa,v}}=\overline{\psi_v}\,e^{-\mathrm i\kappa s_0}$ with $\psi_v=c_ve^{\mathrm i\kappa(v\cdot x-\frac12|v|^2t)}$. Since $\iota_e$ is a ring homomorphism $\mathbb C_{\mathrm i}\to\mathbb H$,
$$\chi=\Big(\sum_vV^v\iota_e(\overline{\psi_v})\,m_{\kappa,v}\Big)\,\iota_e(e^{-\mathrm i\kappa s_0}),$$
a right multiplication by a unit quaternion. So $A^0=\chi\chi^*$:
- is built from $\psi$ alone;
- does not depend on the real coordinate $s_0$ (Part C);
- is conserved in the slots $t$ and $x_j$ and in the two imaginary coordinates $S_1$ of the $s_0$-slot.

This is a positive conserved current for the Schrödinger flow, tangent to $\Phi$, with Gaussian rows. These are the engine data asked for by `T.KAK.12`.

**Proof of Theorem 16.8.**

*Step 1: one mass, small aperture.*
- Let $\kappa\in(2\pi/100)\mathbb Z$ with $\kappa\ge1$, and let $\psi$ on $\mathbb T^d$ have Fourier support in $\{|\xi|\le\gamma\kappa\}$.
- Put $f=\psi(x)e^{\mathrm i\kappa s_0}$ on $\mathbb T^n$ and $h=1/\kappa$. Then $h\xi=(\xi_x/\kappa,1)\in K$ and $|\xi_x|\le\gamma\xi_{s_0}$, so Corollary 16.4 applies.
- By Lemma 16.7, $e^{-\mathrm it\varphi(D)}f=\big(e^{\mathrm i\frac t{2\kappa}\Delta}\psi\big)e^{\mathrm i\kappa s_0}$, and $\|f\|_{L^p(\mathbb T^n)}^p=100\|\psi\|_{L^p(\mathbb T^d)}^p$.
- With $s=t/(2\kappa)$ this gives
$$2\kappa\int_{1/(2\kappa)}^{1/\kappa}\|e^{\mathrm is\Delta}\psi\|_p^p\,ds\le C\log^2(2+\kappa)\,\|\psi\|_p^p .\tag{16.1}$$

*Step 2: Galilean covering.* Let $\psi$ have Fourier support in $\{|\xi|\le3\kappa\}$.
- Choose a smooth partition $\sum_i\chi_i(\xi/\kappa)=1$ on $\{|\xi|\le3\kappa\}$ into $I\le C\gamma^{-d}$ pieces, each supported in a ball $B(u_i,\gamma/2)$.
- For $\kappa\ge\kappa_0:=2\pi\sqrt d/(100\gamma)$, let $k_i$ be a lattice point of $(2\pi/100)\mathbb Z^d$ with $|k_i-\kappa u_i|\le\gamma\kappa/2$.
- Then $\phi_i:=e^{-\mathrm ik_i\cdot x}\chi_i(D/\kappa)\psi$ has Fourier support in $\{|\xi|\le\gamma\kappa\}$.
- By the Galilean identity $e^{\mathrm is\Delta}(e^{\mathrm ik\cdot x}\phi)(x)=e^{\mathrm ik\cdot x-\mathrm is|k|^2}(e^{\mathrm is\Delta}\phi)(x-2sk)$ (Part S), $\|e^{\mathrm is\Delta}\chi_i(D/\kappa)\psi\|_p=\|e^{\mathrm is\Delta}\phi_i\|_p$.
- The multipliers $\chi_i(D/\kappa)$ have periodized kernels with $L^1$ norm $\le C$, so $\|\phi_i\|_p\le C\|\psi\|_p$.
- Apply (16.1) to each $\phi_i$ and sum with $\|\sum_{i\le I}a_i\|^p\le I^{p-1}\sum\|a_i\|^p$.
- For $1\le\kappa<\kappa_0$ there are $O(1)$ frequencies, and $\|e^{\mathrm is\Delta}\psi\|_\infty\le\sum|\hat\psi|\le C\|\psi\|_p$.

*Step 3: every real $\kappa\ge1$ (part 1).*
- Let $\kappa_1\ge\kappa\ge\kappa_2$ be the nearest lattice points of $(2\pi/100)\mathbb Z$, so $\kappa_1\le2\kappa_2$.
- Then $[\frac1{2\kappa},\frac1\kappa]\subset[\frac1{2\kappa_1},\frac1{\kappa_1}]\cup[\frac1{2\kappa_2},\frac1{\kappa_2}]$.
- Frequencies $|\xi|\le2\kappa$ satisfy $|\xi|\le3\kappa_i$ for $i=1,2$.
- Apply Step 2 twice.

*Step 4: from the torus to $\mathbb R^d$ (part 2).*
- Fix $\chi\in C_c^\infty(B(0,2))$ with $\chi=1$ on $B(0,1)$, so $e^{\mathrm is\Delta}\psi=K_s*\psi$ with $K_s$ the kernel of $e^{\mathrm is\Delta}\chi(D/\kappa)$.
- For $s\in[\frac1{2\kappa},\frac1\kappa]$, $K_s(x)=(\frac\kappa{2\pi})^d\int\chi(\eta)e^{\mathrm i\kappa(x\cdot\eta-s\kappa|\eta|^2)}d\eta$. The phase gradient is $\kappa(x-2s\kappa\eta)$, of length $\ge\kappa|x|/2$ when $|x|\ge8$. So $|K_s(x)|\le C_M\kappa^d(\kappa|x|)^{-M}$ there, and $|K_s|\le C\kappa^d$ everywhere.
- Tile $\mathbb R^d$ by unit cubes $Q$, and split $\psi=\psi\mathbf 1_{20Q}+\psi\mathbf 1_{(20Q)^c}$ on each.
  - **Far part.** Young's inequality bounds it by $C\kappa^{-10}\|\psi\|_p$.
  - **Near part.** Place $\psi\mathbf 1_{20Q}$ on $\mathbb T^d$. The torus kernel of $e^{\mathrm is\Delta}\chi(D/\kappa)$ is the periodization of $K_s$. On $Q$ the other periods contribute $O(\kappa^{-10})\|\psi\mathbf 1_{20Q}\|_p$. Part 1, with frequencies $\le2\kappa$, bounds the main term by the right side for $\psi\mathbf 1_{20Q}$.
- Sum over $Q$ with bounded overlap. This is the argument of 08 §9 and [10](10_GAPLESS_REMAINING.md) §4, steps 4–5, with the Schrödinger kernel.

*Step 5: paper (12.2) (part 3).*
- Take $I=[0,1]$ first and split it into the dyadic pieces $[a/2,a]$, $a=2^{-k}$.
- Put $\psi(x)=\phi(x/\lambda)$, so $\phi$ has frequencies $\le\lambda N$ and $e^{\mathrm it\Delta}\psi(x)=(e^{\mathrm i(t/\lambda^2)\Delta}\phi)(x/\lambda)$. The window $[a/2,a]$ becomes $[\frac1{2\kappa},\frac1\kappa]$ with $\kappa=\lambda^2/a$.
- Choose $\lambda=aN$ when $aN^2\ge2$, so that $\kappa=aN^2\ge\lambda N$; otherwise choose $\lambda=\sqrt{2a}$, so that $\kappa=2$.
- Part 2 gives
$$\int_{a/2}^a\|e^{\mathrm it\Delta}\psi\|_p^p\,dt=\lambda^{d+2}\int_{1/(2\kappa)}^{1/\kappa}\|e^{\mathrm is\Delta}\phi\|_p^p\,ds\le C\,a\,\log^2(2+N)\,\|\psi\|_p^p,$$
  using $\|\phi\|_p^p=\lambda^{-d}\|\psi\|_p^p$ and $\lambda^{d+2}\kappa^{-1}\lambda^{-d}=a$ (Part S).
- Summing over $k$ gives $I=[0,1]$.
- Negative times follow by conjugation, since $\overline{e^{\mathrm it\Delta}\psi}=e^{-\mathrm it\Delta}\bar\psi$.
- A window $[a/2,a]$ with $a\le2\sup|I|$ costs $C_Ia$. This covers any compact $I$. ∎

**Remarks.**
1. **Where the dimension count comes from.** The paraboloid in $\mathbb R^{d+1}$ has critical exponent $2(d+1)/d$, the cone's $p_c$ for $n=d+1$. The null reduction explains this: the Monge–Ampère fibre has the $d$ space slots plus the $s_0$ slot, and the passive coordinate $\Phi$ uses one direction. So $d$ free directions remain, as for the wave in $\mathbb R^{d+1}$.
2. **No sphere, no rescaling limit.** [10](10_GAPLESS_REMAINING.md) §4 derives (12.2) by parametrizing a cap of $S^d$ and using Fresnel's identity and the sphere's local extension estimate. That route is still needed for (12.1) and (12.3). For (12.2), Theorem 16.8 is a direct proof.
3. **Masses in a band.** For $f=\sum_\kappa\psi_\kappa(x)e^{\mathrm i\kappa s_0}$ with $\kappa$ in a dyadic band, Corollary 16.4 is a mixed-norm Schrödinger estimate uniform over the band. The $s_0$-integral couples the masses when $p\ne2$.

---

## 6. One logarithm is necessary: focusing

Take $d=1$ ($p=4$) and $\hat\phi(\xi)=\chi(\xi/\kappa)$ with $\chi$ a smooth bump on $[-1,1]$. Let $\psi=e^{-\mathrm is_*\Delta}\phi$ with $s_*=\frac3{4\kappa}$, so that $e^{\mathrm is\Delta}\psi$ focuses at $s_*$, inside the window $[\frac1{2\kappa},\frac1\kappa]$. Put
$$R(\kappa)=\kappa\int_{1/(2\kappa)}^{1/\kappa}\|e^{\mathrm is\Delta}\psi\|_4^4\,ds\Big/\|\psi\|_4^4 .$$
Theorem 16.8(1) says $R(\kappa)\le C\log^2\kappa$.

**Stationary-phase estimate.**
- $\|e^{\mathrm i\sigma\Delta}\phi\|_4^4\approx\kappa\int\chi^4/(8\pi^2|\sigma|)$ for $\kappa^{-2}\ll|\sigma|$, and the integrand is $O(\kappa)$ for $|\sigma|\lesssim\kappa^{-2}$.
- So the integral is $\approx\frac{\kappa\int\chi^4}{4\pi^2}\log\kappa$, while $\|\psi\|_4^4\approx\frac{\kappa^2\int\chi^4}{6\pi^2}$.
- Hence $R(\kappa)=\frac32\log\kappa+O(1)$.

**Numerics (Part F, FFT on a torus of length 20, $\kappa=2^5,\dots,2^{11}$; the period plays no role at these scales).**

| $\kappa$ | $32$ | $64$ | $128$ | $256$ | $512$ | $1024$ | $2048$ |
|---|---|---|---|---|---|---|---|
| $R(\kappa)$ | 3.167 | 4.211 | 5.253 | 6.293 | 7.333 | 8.373 | 9.413 |
| increment | | 1.0441 | 1.0417 | 1.0404 | 1.0399 | 1.0398 | 1.0397 |

The increments converge to $\frac32\log2=1.03972$. In the control, data focusing before the window ($s_*=-1/\kappa$), $R$ stays at $0.2877$.

So the optimal power of the logarithm in (16.1) is at least $1$, against $2$ proved. The power-$1$ question is `T.KAK.6`, transported to the paraboloid.

---

## 7. What this changes

| Before | After |
|---|---|
| `T.KAK.12` open: "a direct engine for the paraboloid" | **done**: Theorem 16.8, via the null reduction |
| START_HERE fact 6: "the paraboloid is reached only through the wave, by rescaling" | the paraboloid has its own engine: the single-mass sector of the null frame |
| 11 §6, "The paraboloid": "$e^{\mathrm it\Delta}$ is not generated by a first-order symmetric system" | true for $e^{\mathrm it\Delta}$ itself; but $e^{\mathrm i\kappa s_0}\psi$ solves the first-order null system (Lemma 16.7). Order and symmetry are restored by one extra variable, not by a limit. |
| [07](07_TOOLKIT.md) K1: Husimi densities are not conserved | the Husimi density is the time component of a positive conserved tensor on the lifted space (Corollary 16.6) |
| the simplest instance of the engine: the arrow pencil (11 §6) | the null frame: polynomial chart, Jacobian $aY^d$, one-multiplier time row |

---

## 8. Status and cross-references

| Result | Status | Uses |
|---|---|---|
| Lemma 16.1 (null chart) | proved here; **exact check** Part G | |
| Lemma 16.2 (rows) | proved here; row identities checked (Part A) | 08 Prop. 3.3 items 1–2 |
| Theorem 16.3, Corollary 16.4 | proved here | 11 Theorem 11.1 |
| Proposition 16.5, Corollary 16.6 | proved here; **checked** to $10^{-50}$ (Part H) | |
| Lemma 16.7 (null reduction) | proved here; **exact check** (Part S) | |
| Reduced current conserved and $s_0$-independent | checked directly (Part C, $10^{-39}$) | 08 Prop. 3.3 |
| Theorem 16.8 (1)–(3) | proved here | Theorem 16.3; 08 §9 transfer; Galilean invariance |
| §6 focusing | numerical (Part F), with a stationary-phase derivation of $\frac32\log\kappa$ | |

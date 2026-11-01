# 14. The angular endpoint (App. B), line by line

*Claude, 28 October 2026. Roadmap target `T.KAK.10`; node `KAK.APP.B` and its new sub-nodes `KAK.ANG.*` in [`ROADMAP.md`](../../../ROADMAP.md).*

> **Theorem 14.1 (the angular endpoint; paper Thm B.1).** Let $n\ge2$, $p=p_c$, $B>\frac12+\frac1p$. For Schwartz $f$,
> $$\|U(t)f\|_{L^p(\mathbb R^n\times I)}\le C_{n,B}\big(\|L(J)^{1/p}f\|_p+\|L(\Lambda_\omega)^Bf\|_p\big),\qquad \|U(t)f\|_{L^{p,\infty}(\mathbb R^n\times I)}\le C_{n,B}\|L(\Lambda_\omega)^Bf\|_p .$$
> On the right, $\|L(\Lambda_\omega)^Bf\|_p$ can be replaced by $\|\Lambda_\omega^{\sigma}f\|_p$ for any $\sigma>0$.

**What the theorem says.**
- **Strong norm.** Theorem 1.2 costs $\log^{2/p}$ on an annulus. Paying a little angular regularity reduces this to the radial cost $\log^{1/p}$, which is sharp: radial data cannot do better (Prop. 8.1).
- **Weak norm.** No spatial logarithm at all is paid. Without angular regularity this is false (Prop. 8.2, node `KAK.OB.WEAK`).

**Verdict on App. B.** **Correct.** Every step was re-derived; nothing had to be repaired. This file contributes four things.
1. **A replacement for the one step the paper only sketches**, the stationary phase that must be uniform in the angular degree (paper Lemma B.4, first paragraph). The replacement is explicit: Funk–Hecke turns the spherical Fourier operator into a multiplier by Bessel functions, and Hankel's integral gives the remainder bound
   $$|(x\partial_x)^k\varepsilon_\mu(x)|\le e^{3/2}(2\mu+2k+1)^{\max(k-1,0)}\frac{(\mu+1)^2}{2x}\qquad(x\ge(\mu+1)^2)$$
   with all constants written out (Lemma 1.2).
2. **Full proofs** of the remaining lemmas, with every constant traced:
   - the ball observation (paper Lemma B.3);
   - the four-kernel envelope and the small-source bound (Lemma B.4);
   - the band estimate (Prop. B.2);
   - the summation.
3. **A list of the cited theorems** with their hypotheses checked (§0.3).
4. **The structure (§7):** five principles that explain the exponents, and why the threshold is $B>A+\frac12-\frac1p$.

**Status.** `KAK.APP.B` moves from `PARTIAL` to `PROVED_CITED` (+`CHECKED`). The citations are the standard harmonic analysis listed in §0.3. With this, no Kakeya node is `PARTIAL`.

**Checks.** [`scripts/verify_angular_endpoint.py`](scripts/verify_angular_endpoint.py) ([`results/angular_endpoint.json`](results/angular_endpoint.json), ~65 s):

| Part | What it checks | Result |
|---|---|---|
| X | All exponent identities (sympy, symbolic $n$) | exact |
| F | The Funk–Hecke formula on $S^1$ and $S^2$, with a negative control | to $10^{-15}$ |
| H | Hankel's integral against `mpmath.hankel1`, with a negative control | to $10^{-26}$ |
| E | The explicit bound of Lemma 1.2 for $k\le3$, $0\le\mu\le40$ | ratio $\le0.22$ |
| E | Necessity of the degree threshold | $x\lvert\varepsilon_\mu\rvert\sim\mu^2/2$ at $x=(\mu+1)^2$; $\lvert\varepsilon\rvert>1$ at $x=\mu+1$ |
| S | The paper's own route: the identity $\lambda\partial_\lambda I_\lambda(a)=-\frac12I_\lambda(w\cdot\nabla a)$ (exact) | exact |
| S | The band facts of Prop. 1.3 | holds |
| K | The four-kernel envelope in $\mathbb R^2$, uniform in degree; a negative control in the transition region | ratio $\le0.3$ |
| G | The transplant geometry, the radial integrals, and the bookkeeping of §6 | holds |

---

## 0. Set-up

### 0.1 Notation

**Exponents.** $n\ge2$, $d=n-1$, $p=p_c=2n/d$, $\alpha=d/2$. These satisfy the identities used throughout (all checked in Part X):
$$\alpha p=n,\quad d-\alpha p=-1,\quad p-2=\tfrac2d,\quad \alpha-\tfrac dp=\tfrac1p,\quad \alpha-\tfrac1p=\tfrac dp,\quad n-\alpha=\alpha+1 .$$

**Operators and norms.**
- $I=[1,2]$ and $L(v)=\log(2+v)$.
- $U(t)=e^{-\mathrm it|D|}$ and $J=(1-\Delta)^{1/2}$.
- $\Lambda_\omega=(1-\Delta_{S^d})^{1/2}$ acts on each sphere about $0$. On spherical harmonics of degree $\nu$ it is the scalar $\sqrt{1+\nu(\nu+d-1)}$.
- The weak norm is $\|u\|_{p,\infty}=\sup_\lambda\lambda|\{|u|>\lambda\}|^{1/p}$, on $I\times\mathbb R^n$ with $dt\,dx$.

**Frequency cut-offs.** $\chi\in C_c^\infty((0,\infty))$ with support in $[c_\chi,C_\chi]$, and $\Pi_N=\chi(|D|/N)$.

**Polar coordinates.** $x=r\theta$ and $y=s\omega$, with $r,s>0$ and $\theta,\omega\in S^d$; $dx=r^d\,dr\,d\theta$.

**Harmonics and the spherical Fourier operator.**
- $\mathcal H_\nu$ is the space of spherical harmonics of degree $\nu$ on $S^d$, and $Z_\nu(\theta\cdot\omega)$ is its reproducing kernel.
- A function has *angular degree $\le M$* if on almost every sphere it lies in $\bigoplus_{\nu\le M}\mathcal H_\nu$.
- $S_zF(\theta)=\int_{S^d}e^{\mathrm iz\theta\cdot\omega}F(\omega)\,d\omega$ for $z\in\mathbb R$.

### 0.2 The input from Theorem 1.2

**(T)** There is $A>0$ such that for every $\varphi\in C_c^\infty(\mathbb R^n\setminus0)$ and $N\ge2$, $\|U(t)\varphi(D/N)f\|_{L^p(I\times\mathbb R^n)}\le C_\varphi L(N)^A\|f\|_p$.

This is node `KAK.T12`, with $A=2/p$. Every step below works for any such $A$, and gives the threshold $B>A+\frac12-\frac1p$.

**Lemma 0.1 (scaled form).** For $T\ge2$ and $g\in L^p$: $\|U(t)\Pi_1g\|_{L^p([T,2T]\times\mathbb R^n)}\le C\,T^{1/p}L(T)^A\|g\|_p$.

*Proof.* Put $f(x)=g(Tx)$. Then $U(t)\Pi_Tf(x)=[U(Tt)\Pi_1g](Tx)$. Changing variables gives
$$\|U\Pi_1g\|_{L^p([T,2T]\times\mathbb R^n)}=T^{(n+1)/p}\|U\Pi_Tf\|_{L^p(I\times\mathbb R^n)}\quad\text{and}\quad\|f\|_p=T^{-n/p}\|g\|_p .$$
Apply (T) with $N=T$. ∎

### 0.3 Cited theorems and their hypotheses

| Tag | Statement used | Hypotheses, checked where used |
|---|---|---|
| [LP] | Littlewood–Paley: for smooth radial dyadic $\psi_j$, $\|(\sum\lvert m_j(D)\psi_j(D)f\rvert^2)^{1/2}\|_q\lesssim\|f\|_q$ for $1<q<\infty$ and multipliers $m_j$ bounded uniformly in a Mikhlin norm (Grafakos, *Classical Fourier Analysis*, Thm 6.1.2, Rem. 6.1.3, §5.6) | §6.4: $m_j=(1+j)^{1/p}L(\langle\xi\rangle)^{-1/p}$ on $\operatorname{supp}\psi_j$ |
| [Syn] | Synthesis: $(w_j)\mapsto\sum_j\tilde\Pi_jw_j$ is bounded $L^q(I\times\mathbb R^n;\ell^2)\to L^q$ for $1<q<\infty$ (duality with [LP], Fubini in $t$) | §6.6 |
| [FS] | Fefferman–Stein: the vector-valued maximal inequality for $\mathcal M_{\mathbb R}$ on $L^p(\mathbb R;\ell^2)$ | §6.4 |
| [GLY] | The vector-valued maximal inequality for $\mathcal M_{S^d}$ on $L^p(S^d;\ell^2)$ (Grafakos–Liu–Yang 2009, Thm 1.2; $S^d$ is a space of homogeneous type) | §6.4 |
| [NPW] | Smooth spherical cut-offs $\sum_\nu\phi(\nu/M)Z_\nu(\theta\cdot\omega)$, $\phi\in C_c^\infty$, have kernels $\le C_DM^d(1+M\operatorname{dist}(\theta,\omega))^{-D}$ (Narcowich–Petrushev–Ward 2006, Thm 3.5, Cors. 3.6, 5.3). Hence $\lvert S_mg\rvert\le C\mathcal M_{S^d}g$ uniformly in $M$ | §6.3 |
| [BC] | Spherical multipliers: if $\lvert\nu^k\partial_\nu^km(\nu)\rvert\le C_k$ for $k\le k_d$ ($m$ extended smoothly), then $m$ is bounded on $L^q(S^d)$, $1<q<\infty$ (Bonami–Clerc 1973, Thm 4.9(1)) | (6.2); and $L(\Lambda_\omega)^{-B}$, $L(\Lambda_\omega)^B\Lambda_\omega^{-\sigma}$ in §6.7 |
| [Mik] | Mikhlin's multiplier theorem | §6.1: $e^{-\mathrm it\lvert\xi\rvert}\psi_0(\xi)$ satisfies $\lvert\partial^\beta m\rvert\lesssim\lvert\xi\rvert^{-\lvert\beta\rvert}$ uniformly for $t\in I$ |
| [RI] | Real interpolation: a linear operator bounded $L^{q_i}(\ell^2)\to L^{q_i}$ for $q_0<p<q_1$ is bounded $L^{p,\infty}(\ell^2)\to L^{p,\infty}$ | §6.6 |
| [FH] | Funk–Hecke: $S_zY=(2\pi)^{n/2}\mathrm i^\nu z^{1-n/2}J_{\nu+n/2-1}(z)\,Y$ for $Y\in\mathcal H_\nu$, $z>0$ | checked in Part F |
| [Hk] | Hankel's integral, for $\operatorname{Re}\mu>-\frac12$, $x>0$: $H^{(1)}_\mu(x)=\sqrt{\tfrac2{\pi x}}\,\dfrac{e^{\mathrm i(x-\mu\pi/2-\pi/4)}}{\Gamma(\mu+\frac12)}\displaystyle\int_0^\infty e^{-u}u^{\mu-\frac12}\big(1+\tfrac{\mathrm iu}{2x}\big)^{\mu-\frac12}du$ (Watson, *Bessel Functions*, §7.2) | checked in Part H |

**Density.** We prove the estimates for Schwartz $f$; density gives the extensions, as in the paper. All functions to which the lemmas are applied lie in $L^1\cap L^p$. All kernel integrals below converge absolutely, so Fubini applies.

---

## 1. The spherical Fourier operator, uniformly in degree

**Lemma 1.1 (Funk–Hecke in two-term form).** For $\nu\ge0$, put $\mu_\nu=\nu+\alpha-\frac12\ (\ge0)$ and define $\varepsilon_\mu(x)$ for $x>0$ by
$$H^{(1)}_\mu(x)=\sqrt{\tfrac2{\pi x}}\,e^{\mathrm i(x-\mu\pi/2-\pi/4)}\big(1+\varepsilon_\mu(x)\big).$$
Then for $z>0$ and $Y\in\mathcal H_\nu$:
$$S_zY=(2\pi)^\alpha z^{-\alpha}\Big[e^{\mathrm i(z-\pi\alpha/2)}\big(1+\varepsilon_{\mu_\nu}(z)\big)Y+e^{-\mathrm i(z-\pi\alpha/2)}\big(1+\overline{\varepsilon_{\mu_\nu}(z)}\big)P_-Y\Big],$$
where $P_-F(\theta)=F(-\theta)$. For $z<0$, $S_zF=S_{|z|}(P_-F)$.

*Proof.*
1. **The Bessel function.** By [FH], $S_zY=(2\pi)^{n/2}\mathrm i^\nu z^{1-n/2}J_{\mu_\nu}(z)Y$, since $\nu+n/2-1=\mu_\nu$. For real $x$ and $\mu$, $J_\mu=\frac12(H^{(1)}_\mu+\overline{H^{(1)}_\mu})$.
2. **The constants.** $(2\pi)^{n/2}\cdot\frac12\sqrt{2/\pi}=(2\pi)^{\alpha}$ and $z^{1-n/2}z^{-1/2}=z^{-\alpha}$.
3. **The phases.**
   - $\mathrm i^\nu e^{-\mathrm i\mu_\nu\pi/2-\mathrm i\pi/4}=e^{-\mathrm i\pi\alpha/2}$.
   - $\mathrm i^\nu e^{\mathrm i\mu_\nu\pi/2+\mathrm i\pi/4}=(-1)^\nu e^{\mathrm i\pi\alpha/2}$.
   - $(-1)^\nu Y=P_-Y$.
4. **Negative $z$.** $S_{-z}F(\theta)=S_zF(-\theta)$, and $S_z$ commutes with $P_-$. ∎

**Check.** For $n=3$, $\nu=0$ we have $\mu=\frac12$ and $\varepsilon\equiv0$. The formula then reads $4\pi\sin z/z$ (Part F).

**Lemma 1.2 (Hankel remainder, explicit).** For $\mu\ge0$, $k\ge0$ and $x\ge(\mu+1)^2$,
$$|(x\partial_x)^k\varepsilon_\mu(x)|\le e^{3/2}\,(2\mu+2k+1)^{\max(k-1,0)}\,\frac{(\mu+1)^2}{2x}.$$

*Proof.*
1. **The integral form.** Put $\beta=\mu-\frac12\ge-\frac12$ and $F(y)=(1+\mathrm iy)^\beta$, principal branch, for $y\ge0$. By [Hk],
   $$\varepsilon_\mu(x)=\frac1{\Gamma(\beta+1)}\int_0^\infty e^{-u}u^\beta\big[F(u/2x)-1\big]\,du .$$
2. **Derivatives of $F$.** Put $w=(1+\mathrm iy)^{-1}$, so $|w|\le1$, $|1-w|=y|w|\le y$, $y\partial_yw=-w(1-w)$ and $y\partial_yF=\beta(1-w)F$. By induction, for $k\ge1$, $(y\partial_y)^kF=F\,(1-w)\,Q_k(w)$, where
   $$Q_1=\beta,\qquad Q_{k+1}=(\beta(1-w)+w)Q_k-w(1-w)Q_k' .$$
   - $Q_k$ has degree $k-1$.
   - Its coefficient sum $\|Q_k\|$ satisfies $\|Q_{k+1}\|\le(2|\beta|+1+2(k-1))\|Q_k\|$, hence $\|Q_k\|\le|\beta|(2\mu+2k+1)^{k-1}$.
   - For $k=0$: $|F(y)-1|=|\int_0^1\beta\,\mathrm iy(1+\mathrm ity)^{\beta-1}dt|\le|\beta|y\,e^{|\beta|y}$.
   - $|F(y)|=(1+y^2)^{\beta/2}\le e^{|\beta|y}$.

   So for every $k\ge0$:
   $$|(y\partial_y)^k(F-1)(y)|\le c_k|\beta|\,y\,e^{|\beta|y},\qquad c_0=1,\quad c_k=(2\mu+2k+1)^{k-1}.$$
3. **Differentiate under the integral.** We have $x\partial_x[F(u/2x)]=-(y\partial_yF)(u/2x)$, and the bound in step 2 is integrable against $e^{-u}u^\beta$ when $|\beta|<2x$. So
   $$|(x\partial_x)^k\varepsilon_\mu(x)|\le\frac{c_k|\beta|}{\Gamma(\beta+1)}\int_0^\infty e^{-u(1-|\beta|/2x)}u^{\beta}\frac u{2x}\,du=c_k\frac{|\beta|(\beta+1)}{2x}\Big(1-\frac{|\beta|}{2x}\Big)^{-(\mu+3/2)}.$$
4. **Use $x\ge(\mu+1)^2$.** Put $m=\mu+1\ge1$.
   - $|\beta|(\beta+1)=|\mu-\frac12|(\mu+\frac12)\le m^2$.
   - $|\beta|/2x\le(\mu+\frac12)/(2m^2)\le1/(2m)$.
   - $-(m+\frac12)\log(1-\frac1{2m})\le\frac{m+1/2}{2m-1}\le\frac32$. ∎

**Remark (the threshold is necessary).** At $x=(\mu+1)^2$, $x|\varepsilon_\mu(x)|$ grows like $\mu^2/2$, and $|\varepsilon_\mu(\mu+1)|>1$ for $\mu\ge10$ (Part E). So no bound $C/x$ uniform in degree exists. The expansion is valid only beyond $\gtrsim\mu^2$, which is why the paper needs radii $(2+M)^K$ wavelengths.

**Proposition 1.3 (uniform two-term expansion on a band).** Let $M\ge0$ and $k_0\ge0$, and put
$$K_0(n,k_0)=d+k_0+3+2\lceil\log_2n\rceil .$$
Let $F$ have degree $\le M$ and $z\ge(2+M)^{K_0}$. Then
$$S_zF=(2\pi)^\alpha z^{-\alpha}\sum_{\sigma=\pm}e^{\mathrm i\sigma(z-\pi\alpha/2)}\big(P_\sigma+E_{\sigma,z}\big)F,\qquad P_+=\mathrm{Id}.$$
Here $E_{\sigma,z}$ is the multiplier $E_{+,z}Y=\varepsilon_{\mu_\nu}(z)Y$, $E_{-,z}Y=\overline{\varepsilon_{\mu_\nu}(z)}P_-Y$ on $\mathcal H_\nu$. It preserves degree $\le M$, and for $k\le k_0$:
$$\|(z\partial_z)^kE_{\sigma,z}F\|_{L^\infty(S^d)}\le C_{n,k_0}\|F\|_{L^1(S^d)} .$$

*Proof.*
1. **The expansion and the band.** The expansion is Lemma 1.1, summed over the finitely many $\nu\le M$. The operator $E$ is diagonal on harmonics, so it preserves the band.
2. **Pointwise bound on projections.** Let $\operatorname{Proj}_\nu$ be the orthogonal projection onto $\mathcal H_\nu$. Its kernel is $Z_\nu(\theta\cdot\omega)=\sum_iY_i(\theta)\overline{Y_i(\omega)}$ for an orthonormal basis. Cauchy–Schwarz gives $|Z_\nu(\theta\cdot\omega)|\le Z_\nu(1)=\dim\mathcal H_\nu/|S^d|$. So $\|\operatorname{Proj}_\nu F\|_\infty\le Z_\nu(1)\|F\|_1$. The peak on the diagonal is checked with Gegenbauer polynomials in Part S.
3. **Counting.** $\sum_{\nu\le M}\dim\mathcal H_\nu\le2(M+1)^d$ (Part S).
4. **The Hankel bound applies.** Every $\nu\le M$ has $(\mu_\nu+1)^2\le(M+\frac n2)^2\le(2+M)^{K_0}\le z$, so Lemma 1.2 applies to each $\mu_\nu$. With step 3:
   $$\|(z\partial_z)^kE_{\sigma,z}F\|_\infty\le\frac{2(M+1)^d}{|S^d|}\,e^{3/2}(2M+n+2k)^{k}\frac{(M+n/2)^2}{2z}\|F\|_1\le C\,(2+M)^{d+k+2-K_0}\|F\|_1,$$
   using $(M+\frac n2)^2\le\frac{n^2}4(2+M)^2$; the exponent $d+k+2-K_0$ is $\le-1$. ∎

**Corollary 1.4 (composition).** Let $F$ have degree $\le M$ and $x,y\ge(2+M)^{K_0}$. Then
$$S_xS_{-y}F=(2\pi)^{2\alpha}(xy)^{-\alpha}\sum_{\sigma,\sigma'=\pm}e^{\mathrm i\sigma x+\mathrm i\sigma'y}e^{-\mathrm i(\sigma+\sigma')\pi\alpha/2}A_{\sigma\sigma'}(x,y)F,$$
where $A_{\sigma\sigma'}(x,y)=(P_\sigma+E_{\sigma,x})(P_{\sigma'}+E_{\sigma',y})P_-$. For $k+l\le k_0$ and every $\theta$:
$$|(x\partial_x)^k(y\partial_y)^lA_{\sigma\sigma'}(x,y)F(\theta)|\le C\big(|F(\theta)|+|F(-\theta)|+\|F\|_1\big).$$

*Proof.* Apply Proposition 1.3 twice, using $S_{-y}=S_yP_-$; the degree is preserved at each step. Expand the product into its four parts.
- **$P_\sigma P_{\sigma'}P_-F$.** This is $F(\pm\theta)$. It occurs only when $k=l=0$.
- **The leftmost factor is an $E$.** Its sup norm is $\le C\|G\|_1$, where $G=(P_{\sigma'}+E_{\sigma',y})P_-F$ satisfies $\|G\|_1\le\|F\|_1+|S^d|\,\|E_{\sigma',y}P_-F\|_\infty\le C\|F\|_1$.
- **$P_\sigma E_{\sigma',y}P_-F$.** Its sup norm is $\le C\|F\|_1$. ∎

**The paper's route** (for comparison) proves Proposition 1.3 by stationary phase at the poles $\pm\theta$, in coordinates with $1\mp\theta\cdot\omega=|w|^2/2$. The key identity is $\lambda\partial_\lambda I_\lambda(a)=-\frac12I_\lambda(w\cdot\nabla a)$, where $I_\lambda(a)=\lambda^{d/2}\int e^{\pm\mathrm i\lambda|w|^2/2}a$. It holds by integration by parts and is checked exactly in Part S. Each log derivative therefore produces an amplitude that vanishes at the pole, hence $O(\lambda^{-1})$. Band-limited $F$ has $\|F\|_{C^m}\lesssim(2+M)^{c_m}\|F\|_1$, which gives the same conclusion with an unspecified $K$. Lemma 1.2 makes the constants explicit.

---

## 2. The unit kernel

**Lemma 2.1 (kernel expansion; paper (B.2)).** Let $K(t,z)$ be the kernel of $U(t)\Pi_1$, and let $v=|z|$. For every $D$, uniformly in $t\in\mathbb R$:
$$K(t,z)=v^{-\alpha}k(t-v)+O_D\big(v^{-\alpha-1}\langle t-v\rangle^{-D}+v^{-\alpha}\langle t+v\rangle^{-D}\big)\quad(v\ge1),\qquad K(t,z)=O_D(\langle t\rangle^{-D})\quad(v<1).$$
Here $k(w)=(2\pi)^{\alpha-n}e^{-\mathrm i\pi\alpha/2}\int_0^\infty e^{-\mathrm i\rho w}\rho^\alpha\chi(\rho)\,d\rho$, a Schwartz function.

*Proof.*
1. **Polar coordinates in $\xi$.** $K=(2\pi)^{-n}\int_0^\infty e^{-\mathrm it\rho}\chi(\rho)\rho^dS_{\rho v}\mathbf 1\,d\rho$.
2. **The degree-0 expansion.** Lemma 1.1 with $\nu=0$ and $\mu_0=\frac{n-2}2$ gives
   $$K=(2\pi)^{\alpha-n}v^{-\alpha}\sum_\sigma e^{-\mathrm i\sigma\pi\alpha/2}\int e^{-\mathrm i\rho(t-\sigma v)}\rho^{\alpha}\chi(\rho)\big(1+\varepsilon^{(\sigma)}(\rho v)\big)\,d\rho,$$
   with $\varepsilon^{(+)}=\varepsilon_{\mu_0}$, $\varepsilon^{(-)}=\overline{\varepsilon_{\mu_0}}$ and $d-\alpha=\alpha$.
3. **Large $v$.** Let $\rho v\ge(\mu_0+1)^2$ on $\operatorname{supp}\chi$, i.e. $v\ge v_0$. By Lemma 1.2, $\partial_\rho^k[\rho^\alpha\chi(\rho)\varepsilon(\rho v)]=O_k(v^{-1})$, since $\rho\partial_\rho=x\partial_x$ and $\rho\sim1$.
   - **The $\sigma=+$ term.** It is $v^{-\alpha}k(t-v)$ plus a remainder. Integrating by parts $D$ times gives the remainder $O(v^{-\alpha-1}\langle t-v\rangle^{-D})$.
   - **The $\sigma=-$ term.** It is $O(v^{-\alpha}\langle t+v\rangle^{-D})$.
4. **Bounded $v$.** For $v<v_0$, $\rho\mapsto\rho^d\chi(\rho)S_{\rho v}\mathbf 1$ has bounded derivatives, uniformly in $v$, so $K=O(\langle t\rangle^{-D})$. For $1\le v<v_0$ both sides of the expansion are $O(\langle t\rangle^{-D})$, which is $O(v^{-\alpha-1}\langle t-v\rangle^{-D})$. ∎

**Corollary 2.2.**
- **(i) A uniform bound.** For $t\ge0$: $|K(t,z)|\le C_D\langle v\rangle^{-\alpha}\langle t-v\rangle^{-D}$.
- **(ii) The $L^1\to L^\infty$ norm.** For $t\in I$ and $N\ge1$: $\|U(t)\Pi_N\|_{L^1\to L^\infty}\le CN^{\alpha+1}$.
- **(iii) Time decay at a fixed point.** For $T\ge2$, $x\in\mathbb R^n$ and $f\in L^p$: $\|U(t)\Pi_1f(x)\|_{L^p_t([T,2T])}\le CT^{1/p}\|f\|_p$.

*Proof.*
- **(i)** For $t,v\ge0$, $\langle t+v\rangle\ge\langle t-v\rangle$. For $v<1$, $\langle t\rangle\simeq\langle t-v\rangle$.
- **(ii)** The kernel of $U(t)\Pi_N$ is $N^nK(Nt,Nz)$. By (i), $\sup_z|K(Nt,Nz)|\le CN^{-\alpha}$: if $|v-t|\ge t/2$ the factor $\langle N(t-v)\rangle^{-D}$ is small, and otherwise $v\ge\frac12$. Finally $n-\alpha=\alpha+1$.
- **(iii)** Write $U(t)\Pi_1f(x)=\int_0^\infty\!\!\int_{S^d}K(t,s\omega)f(x-s\omega)s^d\,d\omega\,ds$ and put $\Phi_x(s)=\int|f(x-s\omega)|d\omega$.
  - **Radii $T/2\le s\le3T$.** By (i) the integrand is $\le Cs^{\alpha}\langle t-s\rangle^{-D}\Phi_x(s)$, using $s^{d-\alpha}=s^\alpha$. Young in $t$ bounds the $L^p_t$ norm by $C\|s^\alpha\Phi_x\mathbf 1_{[T/2,3T]}\|_{L^p(ds)}$. By Hölder, $\Phi_x^p\le C\int|f(x-s\omega)|^pd\omega$, and $s^{\alpha p}=s\cdot s^d\le3T s^d$. This gives $\le C(T\|f\|_p^p)^{1/p}$.
  - **Radii $s<T/2$ and $s>3T$.** Here $|t-s|\ge\max(T,s)/4$, and Hölder gives $O(T^{n-D})\|f\|_p$. ∎

---

## 3. Observation on a ball (paper Lemma B.3)

**Lemma 3.1 (shell geometry).** For $|x|\le R<s$ and $\omega\in S^d$:
$$0\le|s\omega-x|-(s-x\cdot\omega)\le\frac{|x|^2}{2(s-|x|)}.$$

*Proof.* $|s\omega-x|^2-(s-x\cdot\omega)^2=|x|^2-(x\cdot\omega)^2\ge0$. Divide by $|s\omega-x|+(s-x\cdot\omega)\ge2(s-|x|)$. ∎ (Checked on 20,000 random configurations, Part G.)

**Lemma 3.2 (observation on a ball).** For $N\ge2$, $0<r\le1$ and $x_0\in\mathbb R^n$:
$$\|U(t)\Pi_Nf\|_{L^p(I\times B(x_0,r))}\le C\min\{1,(Nr)^{n/p}\}\,L(Nr)^A\|f\|_p .$$

**Meaning.** A ball of $R=Nr$ wavelengths costs only $L(R)^A$, not $L(N)^A$.

*Proof.* Translate so that $x_0=0$, and rescale as in Lemma 0.1 with $T=N$, $R=Nr\le T$. We must show
$$\|U(t)\Pi_1f\|_{L^p([T,2T]\times B(0,R))}\le CT^{1/p}\min(1,R^{n/p})L(R)^A\|f\|_p .\tag{3.1}$$

**Case $R\le2$.** Integrate Corollary 2.2(iii) over $x\in B(0,R)$. This gives $CT^{1/p}R^{n/p}\|f\|_p$, and $L(R)\ge\log2$.

**Case $2<R\le T$, $T\le10S$, where $S:=100R^{3+n/p}$.** Lemma 0.1 gives $CT^{1/p}L(T)^A\|f\|_p$, and $L(T)\le L(1000R^{3+n/p})\le C_nL(R)$.

**Case $T>10S$.** We proceed in five steps.
1. **Shell profiles.** Fix $Q\ge S$ and $h\in L^p(du\,d\omega)$ with $h(u,\omega)=0$ for $|u|>5R$. Put $f_Q(s\omega)=s^{-\alpha}h(s-Q,\omega)$ and
   $$\mathcal Th(b,x)=\int\!\!\int k(b-u+x\cdot\omega)\,h(u,\omega)\,du\,d\omega ,$$
   which does not depend on $Q$.

   *Claim:* for $|b|,|x|\le R$,
   $$\|U(Q+b)\Pi_1f_Q-\mathcal Th\|_{L^p(|b|\le R,|x|\le R)}\le C\frac{R^{2+n/p}}{Q}\|h\|_p .$$

   *Proof of the claim.* Put $s=Q+u$, $\varsigma:=|s\omega-x|$ and $w=b-u+x\cdot\omega$. The difference has kernel $E=s^{d-\alpha}K(Q+b,x-s\omega)-k(w)$ against $h(u,\omega)\,du\,d\omega$. By Lemma 3.1 (with $s-R\ge Q/4$), $t-\varsigma=w-\delta$ with $0\le\delta\le2R^2/Q\le1$. Moreover $s/\varsigma\in[s/(s+R),s/(s-R)]$, so $|(s/\varsigma)^\alpha-1|\le CR/Q$. Lemma 2.1 then gives
   $$|E|\le|(s/\varsigma)^\alpha-1||k(w-\delta)|+|k(w-\delta)-k(w)|+C\varsigma^{-1}\langle w-\delta\rangle^{-D}+C\langle Q\rangle^{-D}\le C\frac{R^2}Q\langle w\rangle^{-D},$$
   using $\langle w-\delta\rangle\ge\langle w\rangle/2$, $R\ge1$, and $|w|\le7R\le Q$ for the last term. For fixed $x$, Young in $b$ and Minkowski in $\omega$ bound the $L^p_b$ norm by $C(R^2/Q)\|h\|_p$. Integrating over $|x|\le R$ costs $R^{n/p}$.
2. **Moving the shell to radius $S$.** Since $\mathcal Th$ does not depend on $Q$, apply step 1 with $Q=S$. The output times $S+b$ lie in $[S/2,S]\cup[S,2S]$, and $\|f_S\|_p^p=\iint|h(u,\omega)|^p(S+u)^{-1}du\,d\omega\le\frac2S\|h\|_p^p$, because $d-\alpha p=-1$ and $|u|\le S/2$. Lemma 0.1 on both intervals then gives
   $$\|\mathcal Th\|_{L^p(|b|,|x|\le R)}\le CS^{1/p}L(2S)^AS^{-1/p}\|h\|_p+C\frac{R^{2+n/p}}S\|h\|_p\le CL(R)^A\|h\|_p,$$
   since $R^{2+n/p}/S=1/(100R)$ (Part X) and $L(S)\le C_nL(R)$. **The factor $S^{1/p}$ from the scaled estimate cancels exactly against the normalization $\|f_S\|_p\simeq S^{-1/p}\|h\|_p$.**
3. **Cutting the input into shells.**
   - Split $[T,2T]$ into intervals of length in $[R/2,R]$ with centres $Q_i$.
   - Put $h_i(u,\omega)=(Q_i+u)^\alpha f((Q_i+u)\omega)\mathbf 1_{|u|\le5R}$. Then $f_{Q_i}=f\cdot\mathbf 1_{|s-Q_i|\le5R}$.
   - The shells overlap at most 21 times and $s\le3T$ on them, so $\sum_i\|h_i\|_p^p=\sum_i\int_{\text{shell}_i}s|f|^p\le63T\|f\|_p^p$.

   For $t$ in interval $i$, steps 1–2 give $\|U\Pi_1f_{Q_i}\|_{L^p(\text{int}_i\times B_R)}\le CL(R)^A\|h_i\|_p$. Summing $p$th powers gives $CL(R)^AT^{1/p}\|f\|_p$.
4. **The rest of the input at radii $T/2\le s\le3T$.** For interval $i$ this is the input with $|s-Q_i|>5R$. For $t$ in interval $i$ and $|x|\le R$ we have $|t-s|\ge4.5R$ and $|t-\varsigma|\ge|t-s|-2R\ge|t-s|/2$. By Corollary 2.2(i) its contribution is at most
   $$\Phi(t,x)=C\int_{T/2}^{3T}s^\alpha\langle t-s\rangle^{-D}\mathbf 1_{|t-s|>4R}\Phi_x(s)\,ds,$$
   which does not depend on $i$. Young in $t$, with $\|\langle\cdot\rangle^{-D}\mathbf 1_{|\cdot|>4R}\|_1\le CR^{1-D}$, gives $\|\Phi\|_{L^p([T,2T]\times B_R)}\le CR^{n/p+1-D}T^{1/p}\|f\|_p$.
5. **Radii outside $[T/2,3T]$.** As in Corollary 2.2(iii), they contribute $O(T^{1/p+n-D}R^{n/p})\|f\|_p$.

Taking $D$ large, steps 3–5 give (3.1). ∎

---

## 4. The four-kernel envelope (paper Lemma B.4)

**Definition 4.1 (radial profile).** For $g$ on $\mathbb R^n$, $v\in\mathbb R$ and $\theta\in S^d$, put
$$h_g(v,\theta)=|v|^{d/p}\Big(|g(|v|\theta)|+|g(-|v|\theta)|+\int_{S^d}|g(|v|\omega)|\,d\omega\Big),\qquad G_g=\mathcal M_{\mathbb R}h_g ,$$
where $\mathcal M_{\mathbb R}$ is the Hardy–Littlewood maximal operator in $v$. Then $h_g$ is even in $v$ and $\|G_g\|_{L^p(dv\,d\theta)}\le C\|g\|_p$:
- $\int_{\mathbb R}\int|v|^d|g(|v|\theta)|^p\,d\theta\,dv=2\|g\|_p^p$;
- the reflected term is the same;
- the average term is handled by Hölder on $S^d$;
- $\mathcal M_{\mathbb R}$ is bounded on $L^p(\mathbb R)$.

**Lemma 4.2 (envelope).** There is $K=K(n,\chi)$ with the following property. Let $g\in L^1\cap L^p$ have angular degree $\le M$, with support in $\{|y|\ge s_0\}$. Then for every $r$ with $N\min(r,s_0)\ge(2+M)^K$, every $t\in I$ and every $\theta$:
$$|U(t)\Pi_Ng(r\theta)|\le C\sum_\pm G_g(t\pm r,\theta)\cdot\begin{cases}r^{-\alpha},&r\le4,\\ r^{-d/p},&r>4,\end{cases}\tag{4.1}$$
with $C$ independent of $M$, $N$, $r$ and $s_0$.

*Proof.* Fix $D=3$; any $D>1+\frac1p$ would do. Let $K_0=K_0(n,D)$ be as in Proposition 1.3, and $K=K_0+\lceil\log_2(1/c_\chi)\rceil+1$.

1. **Polar form.** By Fubini (here $g\in L^1$),
   $$U(t)\Pi_Ng(r\theta)=(2\pi)^{-n}\int_0^\infty s^d\int_0^\infty e^{-\mathrm it\rho}\chi(\rho/N)\rho^d\,\big[S_{r\rho}S_{-s\rho}g_s\big](\theta)\,d\rho\,ds,\qquad g_s=g(s\,\cdot).$$
   For almost every $s$, $g_s$ is a spherical polynomial of degree $\le M$, hence continuous.
2. **Expand.** On $\operatorname{supp}\chi(\cdot/N)$ we have $r\rho,s\rho\ge c_\chi(2+M)^K\ge(2+M)^{K_0}$. Corollary 1.4 applies, and $\rho^d(r s\rho^2)^{-\alpha}=(rs)^{-\alpha}$ since $d=2\alpha$. So
   $$U(t)\Pi_Ng(r\theta)=c\sum_{\sigma,\sigma'}c_{\sigma\sigma'}\,r^{-\alpha}\int_0^\infty s^{d-\alpha}\int_0^\infty e^{-\mathrm i\rho(t-\sigma r-\sigma's)}\chi(\rho/N)A_{\sigma\sigma'}(r\rho,s\rho)g_s(\theta)\,d\rho\,ds .$$
3. **Integrate by parts in $\rho$.** Substitute $\rho=N\lambda$ and put $w=t-\sigma r-\sigma's$. Since $\lambda\partial_\lambda=x\partial_x+y\partial_y$ on $A$ and $\lambda\simeq1$ on $\operatorname{supp}\chi$, Corollary 1.4 gives
   $$|\partial_\lambda^k[\chi(\lambda)A(rN\lambda,sN\lambda)g_s](\theta)|\le C_k\,m_g(s,\theta),\qquad m_g:=|g_s(\theta)|+|g_s(-\theta)|+\|g_s\|_1 .$$
   Integrating by parts $D$ times bounds the $\rho$-integral by $C\,N\langle Nw\rangle^{-D}m_g(s,\theta)$.
4. **Weights.** We have $s^{d-\alpha}=s^\alpha=s^{1/p}s^{d/p}$, because $\alpha-d/p=1/p$, and $s^{d/p}m_g(s,\theta)=h_g(s,\theta)$.
   - **For $r\le4$:** $s\le|w|+t+r\le|w|+6$, so $s^{1/p}\le6\langle Nw\rangle^{1/p}$.
   - **For $r>4$:** $s\le3r\langle w\rangle$, so $s^{1/p}\le3r^{1/p}\langle Nw\rangle^{1/p}$, and $r^{-\alpha}r^{1/p}=r^{-d/p}$.
5. **One convolution per sign.** For $v=t-\sigma r$, the terms $\sigma'=\pm$ combine, by evenness of $h_g$, into
   $$\int_{\mathbb R}N\langle N(v-s')\rangle^{-(D-1/p)}h_g(s',\theta)\,ds'\le c_D\,\mathcal M_{\mathbb R}h_g(v,\theta).$$
   The kernel is even, radially decreasing and integrable with norm $c_D$, independent of $N$. ∎

**Lemma 4.3 (small source; paper (B.4)).** Put $\eta=\alpha/(2n)$ and $\delta=\frac\alpha2(1-\frac2p)>0$. For $g\in L^p$:
$$\|U(t)\Pi_N(\mathbf 1_{|y|<N^{\eta-1}}g)\|_{L^p(I\times\mathbb R^n)}\le CN^{-\delta}\|g\|_p .$$

*Proof.* Let $v=\mathbf 1_{|y|<b}g$ and $F=U(t)\Pi_Nv$.
1. By Hölder, $\|F\|_p\le\|F\|_2^{2/p}\|F\|_\infty^{1-2/p}\le\|v\|_2^{2/p}(CN^{\alpha+1}\|v\|_1)^{1-2/p}$, using Corollary 2.2(ii).
2. By Hölder on the ball, $\|v\|_2\le Cb^{n(\frac12-\frac1p)}\|g\|_p$ and $\|v\|_1\le Cb^{n(1-\frac1p)}\|g\|_p$. The powers of $b$ combine to $n(1-\frac2p)$, so $\|F\|_p\le C(N^{\alpha+1}b^n)^{1-2/p}\|g\|_p$.
3. At $b=N^{\eta-1}$: $\alpha+1+n(\eta-1)=-\alpha/2$.

All three identities are checked in Part X. ∎

---

## 5. One annulus, bounded angular degree (paper Prop. B.2)

**Proposition 5.1.** Let $N\ge2$, and let $f\in L^1\cap L^p$ have angular degree $\le M$. Then
$$\|U(t)\Pi_Nf\|_{L^p(I\times\mathbb R^n)}\le C\big[L(N)^{1/p}+L(M)^A\big]\|f\|_p,\qquad\|U(t)\Pi_Nf\|_{L^{p,\infty}(I\times\mathbb R^n)}\le CL(M)^A\|f\|_p ,$$
with $C$ independent of $M$ and $N$.

*Proof.* Put $H=(2+M)^K$, with $K$ from Lemma 4.2.

**Case $H>N^\eta$.** Then $L(N)\le(K/\eta)L(M)$, and (T) gives both bounds (the weak norm is at most the strong norm).

**Case $H\le N^\eta$.**
1. **The inner ball $r\le H/N$** (radius $\le N^{\eta-1}\le1$). Lemma 3.2 with $x_0=0$ gives the strong bound $CL(H)^A\|f\|_p\le CL(M)^A\|f\|_p$.
2. **The small source.** Outside the inner ball, Lemma 4.3 removes the input from $|y|<N^{\eta-1}$, at strong cost $CN^{-\delta}\|f\|_p$.

   The remaining input $g=\mathbf 1_{|y|\ge N^{\eta-1}}f$ still has angular degree $\le M$, since the cut-off is radial. It satisfies $Nr>H$ and $Ns_0=N^\eta\ge H$, so Lemma 4.2 applies.
3. **The profiles for $r\le4$.** Take one term $r^{-\alpha}G(t\pm r,\theta)$ of (4.1) and change variables $(t,r)\mapsto(v,r)=(t\pm r,r)$; the Jacobian is 1. At fixed $(v,\theta)$, with $G=G_g(v,\theta)$:
   $$\int_{1/N}^4(r^{-\alpha}G)^pr^d\,dr=G^p\log(4N),\qquad\int_{\{r:\,r^{-\alpha}G>\lambda\}}r^d\,dr=\frac{G^p}{n\lambda^p}.\tag{5.1}$$
   These hold because $d-\alpha p=-1$ and $n/\alpha=p$ (Part G, exact). Integrate in $(v,\theta)$ and use $\|G_g\|_p\le C\|f\|_p$:
   - the strong norm on $\{H/N<r\le4\}$ is $\le CL(N)^{1/p}\|f\|_p$;
   - the weak norm is $\le C\|f\|_p$.
4. **The profiles for $r>4$.** We have $(r^{-d/p})^pr^d=1$, so at each $t$ the $L^p(r^ddr\,d\theta)$ norm is $\le\|G_g\|_p$. This gives a strong bound $C\|f\|_p$.

Combining the regions, with the quasi-triangle inequality for the weak norm, proves the proposition. ∎

**Remark 5.2.** For $A=2/p$ and $M\le\exp\sqrt{L(N)}$, the strong cost is $O(L(N)^{1/p})$. The radial focusing data of Prop. 8.1 have $M=0$ and already force $L(N)^{1/p}$.

---

## 6. Summation over frequencies and degrees (proof of Theorem 14.1)

Let $f$ be Schwartz, and let $C$ depend only on $n$, $B$ and the fixed cut-offs.

### 6.1 Littlewood–Paley

**The decomposition.** Take radial $\psi_0+\sum_{j\ge1}\psi(\cdot/2^j)=1$, put $f_j=\psi_j(D)f$ and $N_j=2^j$. Let $\tilde\Pi_j=\tilde\chi(|D|/N_j)$ with $\tilde\chi\in C_c^\infty((0,\infty))$ equal to 1 on the radial support of $\psi$. Then $u_j:=U(t)f_j=\tilde\Pi_ju_j$.

**The low part.** $\|U(t)f_0\|_{L^p}\le C\|f\|_p$ for $t\in I$, by [Mik].

### 6.2 Small sources

Put $a_j=\mathbf 1_{|y|\ge N_j^{\eta-1}}f_j$ and $e_j=U(t)\tilde\Pi_j(f_j-a_j)$. Lemma 4.3 gives
$$\|(e_j)\|_{L^p(\ell^2)}\le\sum_j\|e_j\|_p\le C\sum_j2^{-j\delta}\|f\|_p\le C\|f\|_p .$$

### 6.3 Angular bands and selection

**The bands.**
- Let $M_m=2^{2^m}$ for $m\ge0$, and let $S_m$ be the smooth degree cut-off $\phi(\nu/M_m)$, with $\phi=1$ on $[0,1]$ and $\phi=0$ on $[2,\infty)$.
- Put $Q_m=S_m-S_{m-1}$, with $S_{-1}=0$. The degrees in $\operatorname{supp}Q_m$ lie in $[M_{m-1},2M_m]$.
- On these degrees $L(\Lambda_\omega)\in[c\,2^m,C\,2^m]$ (Part G: the ratio lies in $[0.35,1.24]$ for $m\le12$).

**Three facts.**
- **(a)** By [NPW], $|S_mg|\le C\mathcal M_{S^d}g$ on each sphere, uniformly in $m$.
- **(b)** By [BC], applied on each sphere and integrated in $r^ddr$: the symbol $2^{mB}Q_m(\nu)L(\sqrt{1+\nu(\nu+d-1)})^{-B}$ has bounded log derivatives in $\nu$, uniformly in $m$, so
  $$\|Q_mf\|_p\le C_B2^{-mB}\|L(\Lambda_\omega)^Bf\|_p .\tag{6.2}$$
- **(c)** $S_m$ and $Q_m$ commute with radial Fourier multipliers and with radial cut-offs.

**Selection.** Put $H_m=(2+2M_m)^K$, and let $m(j,r)$ be the largest $m$ with
$$H_m\le\min(N_j^\eta,N_jr),$$
or $-1$ if there is none. It depends only on $j$ and the output radius $r$.

**The decomposition.** At an output point $(t,r\theta)$:
$$U(t)\tilde\Pi_ja_j=\mathrm{sel}_j+\sum_{m>m(j,r)}\mathrm{uns}_{j,m},\qquad\mathrm{sel}_j=[U(t)\tilde\Pi_jS_{m(j,r)}a_j](r\theta),\quad\mathrm{uns}_{j,m}=[U(t)\tilde\Pi_jQ_ma_j](r\theta).$$
This holds because:
- $\sum_{m\le m_0}Q_m=S_{m_0}$;
- $\sum_mQ_ma_j=a_j$ in $L^p$: the $S_m$ are uniformly bounded by (a) and converge on polynomials;
- $U(t)\tilde\Pi_j$ maps $L^p$ boundedly into $C_0$.

### 6.4 The selected pieces: both norms

**Applying the envelope.** Fix $j$ and $r$ with $m=m(j,r)\ge0$. Then $g=S_ma_j\in L^1\cap L^p$ has degree $\le2M_m$ and support in $|y|\ge N_j^{\eta-1}$, and
$$N_j\min(r,N_j^{\eta-1})=\min(N_jr,N_j^\eta)\ge H_m=(2+2M_m)^K .$$
So Lemma 4.2 applies.

**A majorant that does not depend on $r$.** By (a), $|S_ma_j|\le C\mathcal M_{S^d}f_j=:Cf_j^*$. Every term of $h_g$ is monotone in $|g|$, so $G_{S_ma_j}\le C\,G_{f_j^*}=:C\,G_j$. This bound does not depend on $m$, hence not on $r$. This is the point of passing to the maximal function.

**The radial factor.** On the selected set $N_jr\ge H_0\ge1$, so $r^{-\alpha}\le2^\alpha(r+N_j^{-1})^{-\alpha}$.

**The vector-valued maximal bound.** For non-negative $c_j$,
$$\|(c_jG_j)\|_{L^p(dv\,d\theta;\ell^2)}\le C\|(c_jf_j)\|_{L^p(\mathbb R^n;\ell^2)} .$$
It follows by combining four facts:
- [FS] in $v$;
- [GLY] on each sphere;
- Minkowski for the reflection term and the average term;
- the weight $|v|^{d/p}$, which turns $L^p(dv\,d\theta)$ into $L^p(\mathbb R^n)$.

**Strong norm, $r\le4$.** Change variables to $v=t\pm r$. At fixed $(v,\theta)$, Minkowski in $L^{p/2}(r^d\,dr)$ (valid since $p\ge2$) gives
$$\Big\|\Big(\sum_j(r+N_j^{-1})^{-2\alpha}G_j^2\Big)^{1/2}\Big\|_{L^p((0,4),r^ddr)}\le C\Big(\sum_j(1+j)^{2/p}G_j^2\Big)^{1/2},$$
because $\int_0^4(r+2^{-j})^{-n}r^d\,dr\le\frac1n+\log(4\cdot2^j)$ (Part G). Apply the maximal bound with $c_j=(1+j)^{1/p}$, then [LP] with the multipliers $(1+j)^{1/p}L(\langle\xi\rangle)^{-1/p}$ on $\operatorname{supp}\psi_j$, which have uniform Mikhlin norms. The result is $\le C\|L(J)^{1/p}f\|_p$.

**Weak norm, $r\le4$.** Here $(\sum_j|\mathrm{sel}_j|^2)^{1/2}\le Cr^{-\alpha}\sum_\pm\Gamma(t\pm r,\theta)$ with $\Gamma=(\sum_jG_j^2)^{1/2}$. The distribution formula (5.1) gives a weak norm $\le C\|\Gamma\|_{L^p(dv\,d\theta)}$, and the maximal bound with $c_j=1$ followed by [LP] gives $\le C\|f\|_p$.

**Both norms, $r>4$.** The factor is $r^{-d/p}$, and as in §5 step 4 we get a strong bound $C\|f\|_p$.

### 6.5 The unselected pieces: strong norm

Fix $m\ge0$. Put $v_j=U(t)\tilde\Pi_jQ_ma_j$ and $b_j=\|Q_mf_j\|_p$; then $\|Q_ma_j\|_p\le b_j$ by (c). We bound the square function of the $v_j$ over the set $\{m(j,r)<m\}$, which is $\{H_m>\min(N_j^\eta,N_jr)\}$.

**Indices with $H_m>N_j^\eta$.** There are at most $1+\eta^{-1}\log_2H_m$ of them. For each, (T) gives $\|v_j\|_p\le CL(N_j)^Ab_j\le CL(H_m)^Ab_j$. The pointwise inequality $\|w\|_{\ell^2}\le(\#)^{\frac12-\frac1p}\|w\|_{\ell^p}$ then gives
$$\le CL(H_m)^{A+\frac12-\frac1p}\|(b_j)\|_{\ell^p}.$$

**Indices with $H_m\le N_j^\eta$.** Here the unselected set is $\{r<R_j\}$ with $R_j=H_m2^{-j}\le N_j^{\eta-1}\le1$.
- **Shells.** Let $D=\lceil\log_2(2H_m)\rceil$ and $\Sigma_k=\{R_{k+1}<r\le R_k\}$ with $R_k=H_m2^{-k}$. Only indices $j\le k$ contribute on $\Sigma_k$.
- **Near indices, $k-D\le j\le k$.** At most $D+1$ terms, so the $\ell^2$ norm costs $(1+D)^{\frac12-\frac1p}$ over the $\ell^p$ norm. Summing over $k$, each $v_j$ is integrated only over $\{r\le R_j\}$. There Lemma 3.2 (with $N_jR_j=H_m$) gives $CL(H_m)^Ab_j$.
- **Far indices, $j<k-D$.** Here $N_jR_k=H_m2^{j-k}\le2^{-(k-j-D)-1}$ (Part G, exhaustive). By Lemma 3.2 on the ball of radius $R_k$, $\|v_j\|_{L^p(I\times B_{R_k})}\le C2^{-(k-j-D)n/p}b_j$. Bound the $\ell^2$ norm by the $\ell^1$ norm, sum in $j$, then sum in $\ell^p(k)$ by discrete Young with kernel $\sum_{i\ge1}2^{-in/p}<\infty$. The result is $\le C\|(b_j)\|_{\ell^p}$.

**Total for band $m$.** Since $D\le CL(H_m)$, the band costs $CL(H_m)^{A+\frac12-\frac1p}\|(b_j)\|_{\ell^p}$. Two further facts finish the band:
- $\|(b_j)\|_{\ell^p}=\|(Q_mf_j)\|_{L^p(\ell^p)}\le\|(Q_mf_j)\|_{L^p(\ell^2)}\le C\|Q_mf\|_p$, by $p\ge2$, $Q_mf_j=\psi_j(D)Q_mf$ and [LP];
- $L(H_m)\le CK2^m$.

So, by (6.2),
$$\Big\|\Big(\sum_j|v_j|^2\mathbf 1_{m(j,r)<m}\Big)^{1/2}\Big\|_{L^p}\le C2^{m(A+\frac12-\frac1p)}\|Q_mf\|_p\le C_B2^{m(A+\frac12-\frac1p-B)}\|L(\Lambda_\omega)^Bf\|_p .$$

**Sum over $m$.** By Minkowski in $\ell^2$, the sum over $m$ converges exactly when $B>A+\frac12-\frac1p$. At $A=2/p$ this is $B>\frac12+\frac1p$.

### 6.6 Reassembly

We have $U(t)f=U(t)f_0+\sum_{j\ge1}\tilde\Pi_ju_j$, and pointwise
$$u_j=e_j+\mathrm{sel}_j+\sum_{m>m(j,r)}\mathrm{uns}_{j,m}.$$
The pieces are not frequency-localized, because the splitting depends on $r$. But $u_j$ is, so we bound the square function of $(u_j)$ and then apply [Syn].
- **Strong norm.** By §§6.2–6.5,
  $$\|(u_j)\|_{L^p(\ell^2)}\le C(\|f\|_p+\|L(J)^{1/p}f\|_p+\|L(\Lambda_\omega)^Bf\|_p).$$
- **Weak norm.** $\|(u_j)\|_{L^{p,\infty}(\ell^2)}\le C(\|f\|_p+\|L(\Lambda_\omega)^Bf\|_p)$. Here only the selected pieces at $r\le4$ need the weak norm; all other pieces have strong bounds.
- **Synthesis.** [Syn] in the strong case and [RI]+[Syn] in the weak case give the two bounds for $\sum_j\tilde\Pi_ju_j$.

### 6.7 Final form

- **The weak norm has angular weight only.** $\|f\|_p\le C\|L(\Lambda_\omega)^Bf\|_p$, since $L(\Lambda_\omega)^{-B}$ is a bounded spherical multiplier by [BC].
- **The $\Lambda_\omega^\sigma$ variant.** $L(\Lambda_\omega)^B\Lambda_\omega^{-\sigma}$ is bounded by [BC], which gives the form with $\|\Lambda_\omega^\sigma f\|_p$.

This proves Theorem 14.1. ∎

---

## 7. Distillation: why the exponents are what they are

**(P1) Uniform Hankel.**
- **The statement.** On angular degrees $\le M$, at radii beyond $(2+M)^{K}$ wavelengths, the spherical Fourier operator is two point evaluations, at $\theta$ and $-\theta$, plus an error bounded $L^1\to L^\infty$ (Prop. 1.3).
- **The threshold is quadratic in the degree.** It is the Bessel transition $x\simeq\mu^2$, with the explicit constant of Lemma 1.2. Below it the expansion fails by $O(1)$ (Parts E and K).
- **Consequence.** Angular regularity buys a *radius*, the radius beyond which the wave is effectively radial.

**(P2) The wave becomes one-dimensional.** Beyond that radius,
$$|u(t,r\theta)|\lesssim r^{-\alpha}\sum_\pm G(t\pm r,\theta),\qquad\|G\|_{L^p(dv\,d\theta)}\lesssim\|f\|_p .$$
These are the four Müller–Seeger kernels $t\pm r\pm s=0$, recovered uniformly in degree (Lemma 4.2). In $\mathbb R^2$ this is checked directly on Bessel kernels: the remainder ratio is at most $0.12$, uniformly in degree (Part K).

**(P3) Logarithm accounting.** The identity $d-\alpha p=-1$ means that every dyadic radius between one wavelength and unit distance carries the same share of $\|r^{-\alpha}G\|_p^p$. Hence:
- **the strong norm** pays $\int_{1/N}^4dr/r=\log 4N$, i.e. $L(N)^{1/p}$, and this is sharp, as in the focusing example Prop. 8.1;
- **the weak norm** pays nothing, because $r^{-\alpha}$ is exactly weak-$L^p$ for $r^ddr$ when $\alpha p=n$.

So the weak-norm obstruction (Prop. 8.2) is purely angular. Its line families are not angularly band-limited.

**(P4) Shell transplant (Lemma 3.2).**
- **The principle.** Seen from a ball of $R$ wavelengths, an input shell at radius $Q\gg R^3$ produces a wave independent of $Q$, up to $R^{2+n/p}/Q$. So the shell can be moved to radius $S=100R^{3+n/p}$, where the global estimate costs $L(S)^A\simeq L(R)^A$.
- **The cancellation.** The normalizations cancel exactly: $S^{1/p}$ from the scaled estimate against $\|f_S\|_p\simeq S^{-1/p}\|h\|_p$, both coming from $d-\alpha p=-1$.
- **Consequence.** *Observation on a ball costs the logarithm of its radius in wavelengths.*

**(P5) Double-exponential bands.**
- **Why double-exponential.** With $M_m=2^{2^m}$ the logarithm of the degree is geometric, $L\simeq2^m$, so a cost polynomial in $L(\text{degree})$ is summable against $2^{-mB}$.
- **The cost of band $m$.** It is $L(H_m)^{A+\frac12-\frac1p}$:
  - **$A$** is the ball cost (P4) on the region where stationary phase is not yet uniform;
  - **$\frac12-\frac1p$** is the $\ell^2$-versus-$\ell^p$ price over the $\simeq\log H_m$ spatial frequencies at the threshold.
- **The threshold.** $B>A+\frac12-\frac1p$, which is $\frac12+\frac1p$ at $A=2/p$.

**Open questions** (recorded as target `T.KAK.13`):
- **The threshold.** The paper does not claim that $B>\frac12+\frac1p$ is sharp, and we do not know whether it is.
- **Lowering $A$.** Any improvement of the annular power $A$ (target `T.KAK.6`, which would give $A=1/p$) lowers the threshold to $B>\frac12$ automatically.
- **Removing $\frac12-\frac1p$.** A square-function form of Lemma 3.2 over the $\simeq\log H_m$ near-threshold frequencies would remove the $\frac12-\frac1p$.

---

## 8. Dependencies and status

| Step | Depends on | Status |
|---|---|---|
| Lemma 1.2, Prop. 1.3, Cor. 1.4 | [FH], [Hk] | **NEW**, proved; checked (Parts F, H, E, S) |
| Lemma 2.1, Cor. 2.2 | Lemma 1.2 | proved |
| Lemma 3.2 | (T) = `KAK.T12`, Lemma 2.1, Lemma 3.1 | proved; geometry checked (Part G) |
| Lemmas 4.2, 4.3 | Cor. 1.4, Cor. 2.2 | proved; envelope checked in $\mathbb R^2$ (Part K) |
| Prop. 5.1 | Lemmas 3.2, 4.2, 4.3, (T) | proved |
| Theorem 14.1 | Prop. 5.1's lemmas, (T), [LP], [Syn], [FS], [GLY], [NPW], [BC], [Mik], [RI] | **PROVED_CITED** |

Confidence (subjective): 95% that Theorem 14.1 holds as stated, given Theorem 1.2. The residual risk is in the cited harmonic analysis ([GLY], [NPW], [BC]); we checked only their hypotheses, not their proofs.

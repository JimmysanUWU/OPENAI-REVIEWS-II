# 9. Gapless consequences: restriction, Kakeya, Bochner–Riesz and Stein's square function

Claude, 27 October 2026. This file derives, from Theorem 1.2 of [08](08_GAPLESS_PROOF.md), the five consequences that OpenAI formalized beyond local smoothing:
- spherical restriction in the strict range;
- the Kakeya maximal bound;
- full Hausdorff dimension of Kakeya sets;
- uniform Bochner–Riesz bounds;
- Stein's square function.

Together with [08](08_GAPLESS_PROOF.md) (local smoothing), these are the paper's six formal targets ([06](06_LEAN_AUDIT.md) §6.1). The arguments are the paper's §9–§11, with every step written out. The classical inputs are cited by name:

| Input | Where used | Source |
|---|---|---|
| Asymptotics of $\hat\sigma$ (Bessel functions) | Lemma 1.1 | Stein, *Harmonic Analysis* (1993), Ch. VIII |
| Tao's $\varepsilon$-removal theorem | Cor. 2.1 | T. Tao, Duke Math. J. 96 (1999), Thm 1.2; formalized in `Fourier/Restriction/Sparse` |
| Riesz–Thorin, Hölder, Jensen | throughout | — |
| Hardy–Littlewood maximal theorem; Littlewood–Paley inequality | §6, §7 | Stein (1970), Chs. I and IV |
| Fefferman–Stein sharp-function inequality | Lemma 7.2 | Stein (1993), Ch. IV §3.6 |

Notation: $d=n-1$, $p=p_c=\frac{2n}{d}$, $q=\frac p2=\frac n{n-1}$, $L(R)=\log(2+R)$. $A$ denotes the logarithmic power in Theorem 1.2; [08](08_GAPLESS_PROOF.md) gives $A=2/p_c$, and any fixed $A$ suffices for the formal targets. $\sigma$ is surface measure on $S^{n-1}$ and $Eg(X)=\int_{S^{n-1}}e^{\mathrm iX\cdot\omega}g(\omega)\,d\sigma(\omega)$.

---

## 1. From waves to local spherical extension

**Lemma 1.1 (paper Lemma 9.2).** Suppose Theorem 1.2 holds with $L(N)^A$ for the cutoff $\chi(|\xi|/N)$, where $\chi\in C_c^\infty((\frac12,2))$ is radial and equals $1$ near $1$. Then for $R\ge2$, $z\in\mathbb R^n$ and $g\in L^p(S^{n-1})$,
$$\|Eg\|_{L^p(B(z,R))}\le C\,L(R)^A\,\|g\|_{L^p(S^{n-1})}.$$

*Proof.*
1. **A time-averaged wave.** Fix $\psi\in C_c^\infty((1,2))$ with $\psi\ge0$, $\psi\ne0$, and put $M_Nf=\int\psi(T)e^{\mathrm iNT}U(T)\chi(|D|/N)f\,dT$. Hölder in $T$ gives $\|M_Nf\|_p\le CL(N)^A\|f\|_p$.
2. **Its kernel.** $M_N$ has the radial multiplier $\chi(\rho/N)\hat\psi(\rho-N)$, $\rho=|\xi|$. In polar coordinates the angular integral at $u=\rho|Z|$ is
$$\hat\sigma(u)=u^{-d/2}(c_+e^{\mathrm iu}+c_-e^{-\mathrm iu})+O(u^{-d/2-1}),\qquad c_+c_-\ne0 .$$
   - In the leading terms, $\chi(\rho/N)\rho^{d/2}=N^{d/2}(1+O(|\rho-N|/N))$ near $\rho=N$. Since $\int|s\hat\psi(s)|\,ds<\infty$, the error is $O(N^{d/2-1})$; so is the contribution of the $O(u^{-d/2-1})$ term.
   - Fourier inversion gives $\int\hat\psi(\rho-N)e^{\pm\mathrm i\rho|Z|}d\rho=2\pi e^{\pm\mathrm iN|Z|}\psi(\pm|Z|)$. Only the $+$ sign survives, because $\operatorname{supp}\psi\subset(1,2)$.
   
   Hence, for $\frac12\le|Z|\le3$,
$$K_N(Z)=c_nN^{d/2}|Z|^{-d/2}e^{\mathrm iN|Z|}\psi(|Z|)+O(N^{d/2-1}),\qquad c_n\ne0 .\tag{1.1}$$
3. **The input.** Let $f_N(r\omega)=e^{-\mathrm iNr}\psi(r)g(\omega)$, so $\|f_N\|_p=C_\psi\|g\|_p$. For $|X|\le R$ and $N\ge4R$:
   - $\frac12\le|X/N-r\omega|\le3$ on the support;
   - $N|r\omega-X/N|=Nr-X\cdot\omega+O(R^2/N)$;
   - the amplitudes change by $O(R/N)$.
   
   Substituting into (1.1),
$$N^{-d/2}M_Nf_N(X/N)=c'Eg(-X)+O\big((1+R^2)N^{-1}\|g\|_1\big),\qquad c'=c_n\int r^{d/2}\psi(r)^2dr\ne0 .$$
4. **Choice of $N$.** Take $N=(4+R)^{3+n/p}$. The error has $L^p(B(0,R))$ norm $\le C\|g\|_p$. On the left, the change of variables costs $N^{n/p-d/2}=1$ at $p=p_c$, and $L(N)\le C_nL(R)$. Replacing $g$ by $e^{\mathrm iz\cdot\omega}g$ moves the ball to any centre. ∎

## 2. Global extension (the restriction conjecture in the strict range)

**Corollary 2.1 (paper Cor. 9.3, (9.3)).** Let $1<a\le\infty$ and $p_c<b<\infty$ with $b>\frac{n+1}{n-1}a'$. Then $\|Eg\|_{L^b(\mathbb R^n)}\le C_{n,a,b}\|g\|_{L^a(S^{n-1})}$.

*Proof.*
1. **Diagonal estimates.** By Lemma 1.1 and $L(R)^A=O_\varepsilon(R^\varepsilon)$, the local estimate $E:L^{p_c}(S^{n-1})\to L^{p_c}(B(0,R))$ holds with loss $R^\varepsilon$ for every $\varepsilon>0$. Dually, restriction $f\mapsto\hat f|_{S^{n-1}}$ maps functions supported in $B(0,R)$ from $L^{p_c'}$ to $L^{p_c'}(S^{n-1})$ with that loss. Since $1<p_c'<2$, **Tao's $\varepsilon$-removal theorem** (cited) gives global restriction $L^u(\mathbb R^n)\to L^u(S^{n-1})$ for every $1<u<p_c'$. Dualizing, $E:L^s\to L^s$ for every $s>p_c$.
2. **Off-diagonal estimates.**
   - If $a\ge b$, take $s=b$ and use $L^a(S^{n-1})\subset L^b(S^{n-1})$.
   - If $a<b$, put $s=1+b/a'$. Then $s>p_c$, because $p_c=1+\frac{n+1}{n-1}$ and $b/a'>\frac{n+1}{n-1}$; and $s<b$, because $a<b$. Interpolate $E:L^s\to L^s$ with the trivial $E:L^1\to L^\infty$ at $\theta=s/b$: then $1/b=\theta/s$ and $1/a=\theta/s+(1-\theta)$. ∎

## 3. From caps to tubes

For $\omega\in S^{n-1}$, $x\in\mathbb R^n$ and $0<\delta\le1$, $T_\delta(\omega,x)$ is the open $\delta$-neighbourhood of $x+[-\frac12,\frac12]\omega$. Its volume is $\asymp\delta^d$.

**Lemma 3.1 (paper Lemma 9.4).** Assume Lemma 1.1 with power $A$. Let $0<\delta\le\frac12$, let $\{\omega_j\}$ be $\delta$-separated, and let $x_j\in\mathbb R^n$ and $a_j\ge0$ be arbitrary. With $T_j=T_{2\delta}(\omega_j,x_j)$,
$$\Big\|\sum_ja_j\mathbf 1_{T_j}\Big\|_{L^q(\mathbb R^n)}\le C\,L(1/\delta)^{2A}\Big(\delta^d\sum_ja_j^q\Big)^{1/q}.$$

*Proof.*
1. **One cap lights up one tube.** Let $C_j=\{|\omega-\omega_j|<c\delta\}$ with $c<\frac12$ small; these caps are disjoint. Put $h_j=e^{-\mathrm i\delta^{-2}x_j\cdot\omega}\mathbf 1_{C_j}$. For $Y\in T_j$ write $Y-x_j=t\omega_j+v$ with $|t|\le\frac12$ and $|v|<2\delta$. On $C_j$ the phase $\delta^{-2}(Y-x_j)\cdot\omega$ differs from its value at $\omega_j$ by at most $\frac{c^2}4+2c$, using $\omega_j\cdot(\omega-\omega_j)=-\frac12|\omega-\omega_j|^2$. So $|Eh_j(\delta^{-2}Y)|\ge\frac12\sigma(C_j)\ge c'\delta^d$.
2. **Random signs on one cube.** Group the $j$ by the unit cube $Q\ni x_j$. The tubes of a group lie in $5Q$, and $\delta^{-2}(5Q)$ lies in a ball of radius $C_n\delta^{-2}$. For random signs $\varepsilon_j$ put $g_\varepsilon=\sum_{j\in J_Q}\varepsilon_ja_j^{1/2}h_j$. Disjointness of the caps gives $\|g_\varepsilon\|_p^p\le C\delta^d\sum_{J_Q}a_j^q$.
3. **Jensen.** At each $Y$, $\mathbb E|Eg_\varepsilon(\delta^{-2}Y)|^p\ge(\sum a_j|Eh_j(\delta^{-2}Y)|^2)^q\ge c\delta^{2dq}(\sum_{J_Q}a_j\mathbf 1_{T_j}(Y))^q$, because $p/2=q\ge1$.
4. **Integrate and apply Lemma 1.1.** Integrate over $5Q$, change variables to $\delta^{-2}Y$, and apply Lemma 1.1 on the enclosing ball. Since $2dq=2n=dp$, the powers of $\delta$ cancel:
$$\int_{5Q}\Big(\sum_{J_Q}a_j\mathbf 1_{T_j}\Big)^q\le CL(1/\delta)^{Ap}\,\delta^d\sum_{J_Q}a_j^q .$$
5. **Sum over cubes.** The cubes $5Q$ overlap at most $5^n$ times, so $(\sum_QF_Q)^q\le5^{n(q-1)}\sum_QF_Q^q$ pointwise. Summing over $Q$ and taking $q$-th roots gives the power $Ap/q=2A$. ∎

## 4. The Kakeya maximal estimate

$K_\delta f(\omega)=\sup_x|T_\delta(\omega,x)|^{-1}\int_{T_\delta(\omega,x)}|f|$. It is lower semicontinuous in $\omega$, since tube averages are continuous in $(\omega,x)$ for $f\in L^n$; hence it is Borel.

**Corollary 4.1 (Kakeya maximal conjecture; paper Cor. 9.1).** For $0<\delta\le1$,
$$\|K_\delta f\|_{L^n(S^{n-1})}\le C_n\,L(1/\delta)^{2A}\|f\|_{L^n(\mathbb R^n)} .$$
With $A=2/p_c$ the power is $4/p_c$. In particular the bound is $C_\varepsilon\delta^{-\varepsilon}\|f\|_n$ for every $\varepsilon>0$, which is the formal statement.

*Proof.*
1. **Large $\delta$.** For $\frac12<\delta\le1$, Hölder on each tube suffices.
2. **Discretize.** Let $\delta\le\frac12$ and let $\{\omega_j\}$ be a maximal $\delta$-separated set. If $|\omega-\omega_j|\le\delta$, then $T_\delta(\omega,x)\subset T_{2\delta}(\omega_j,x)$ and the volumes are comparable. Hence
$$\|K_\delta f\|_{L^n(S^{n-1})}^n\le C\,\delta^d\sum_jK_{2\delta}f(\omega_j)^n .$$
3. **Duality.** By $\ell^n$–$\ell^q$ duality with weight $\delta^d$, it suffices to bound $\delta^d\sum_ja_j|T_j|^{-1}\int_{T_j}|f|$ for $a_j\ge0$ with $\delta^d\sum a_j^q=1$ and arbitrary centres. Since $|T_j|\asymp\delta^d$, this is at most $C\int|f|\sum_ja_j\mathbf 1_{T_j}\le C\|f\|_n\,\|\sum a_j\mathbf 1_{T_j}\|_q$. Lemma 3.1 bounds the last factor by $CL(1/\delta)^{2A}$. ∎

## 5. Kakeya sets have full Hausdorff dimension

**Theorem 5.1 (paper Cor. 9.1).** Every set $F\subset\mathbb R^n$ that contains a unit segment in every direction has $\dim_HF=n$. No measurability or boundedness is assumed.

*Proof.* Let $I_\omega\subset F$ be a unit segment in direction $\omega$. Fix $\gamma>n-1$ and the gauge $h_\gamma(r)=r^n[\log(1/r)]^{2An}[\log\log(1/r)]^\gamma$.
1. **Group the radii.** Any set of diameter $r$ lies in a ball of radius $r$, and $h_\gamma(2r)\le Ch_\gamma(r)$, so covers by open balls suffice. Let $\{B_j\}$ be a countable cover of $F$ by open balls of radii $r_j\le e^{-2^m}$. Group the radii (Keich's grouping):
   - $J_k=\{j:2^k\le\log(1/r_j)<2^{k+1}\}$ for $k\ge m$;
   - $\delta_k=e^{-2^{k+1}}<r_j$ for $j\in J_k$;
   - $G_k=\bigcup_{J_k}2B_j$ and $v_k=\sum_{J_k}r_j^n$.
2. **Length on a segment.** At each point of $I_\omega\cap\bigcup_{J_k}B_j$, the transverse disc of radius $\delta_k$ lies in $G_k$, and discs at different points are disjoint. By Fubini,
$$|I_\omega\cap\textstyle\bigcup_{J_k}B_j|\le C\delta_k^{-d}\int_{T_{\delta_k}(\omega,\cdot)}\mathbf 1_{G_k}\le CK_{\delta_k}\mathbf 1_{G_k}(\omega).$$
   Summing over $k$, $1\le C\sum_{k\ge m}K_{\delta_k}\mathbf 1_{G_k}(\omega)$ for **every** $\omega$.
3. **The maximal estimate.** Take $L^n(S^{n-1})$ norms and apply Corollary 4.1, with $L(1/\delta_k)\asymp2^{k}$ and $|G_k|\le Cv_k$:
$$1\le C\sum_{k\ge m}2^{2Ak}v_k^{1/n}.$$
4. **Compare with the gauge.** For $j\in J_k$, $h_\gamma(r_j)\asymp r_j^n2^{2Ank}k^\gamma$, so $2^{2Ak}v_k^{1/n}\le C(\sum_{J_k}h_\gamma(r_j))^{1/n}k^{-\gamma/n}$. Hölder with exponents $n$ and $n'$ gives
$$1\le C\Big(\sum_jh_\gamma(r_j)\Big)^{1/n}\Big(\sum_{k\ge m}k^{-\gamma/(n-1)}\Big)^{(n-1)/n}.$$
   The last factor tends to $0$ as $m\to\infty$ because $\gamma>n-1$. So $\mathcal H^{h_\gamma}(F)=\infty$.
5. **Conclusion.** Since $h_\gamma(r)=o(r^s)$ for every $s<n$, $\mathcal H^s(F)=\infty$ for $s<n$. ∎

---

## 6. Bochner–Riesz means

$B^\alpha_rf$ has multiplier $(1-|\xi|^2/r^2)_+^\alpha$, and $\delta(p)=\max\{0,n|\frac1p-\frac12|-\frac12\}$.

**Lemma 6.1 (rescaled local smoothing).** Let $2\le p<\infty$, $\sigma>s_+(p)$, and let $\hat g$ be supported in a fixed annulus. Then for $R\ge1$,
$$\|U(t)g\|_{L^p(\mathbb R^n\times\{R\le|t|\le2R\})}\le CR^{1/p+\sigma}\|g\|_p .$$

*Proof.* $G=g(R\cdot)$ has frequencies $\asymp R$, and $U(t)G(x)=(U(Rt)g)(Rx)$.
- **Scaling.** Changing variables, $\|U(t)g\|_{L^p(\mathbb R^n\times[R,2R])}=R^{(n+1)/p}\|U(t)G\|_{L^p(\mathbb R^n\times I)}$.
- **Theorem 1.1** bounds the right side by $CR^{(n+1)/p}R^{\sigma}\|G\|_p=CR^{1/p+\sigma}\|g\|_p$. The multiplier $R^{-\sigma}J^\sigma$ has a uniformly integrable kernel on this annulus.
- **Negative times** follow by conjugation. ∎

**Lemma 6.2 (maximal estimate for a thin shell; paper Lemma 10.2).** For $0<\delta\le\frac18$ and $\psi\in C_c^\infty((\frac12,2))$, let $A_{\delta,r}f=\psi((1-|D|/r)/\delta)f$. For $2\le p<\infty$ and $\sigma>s_+(p)$,
$$\Big\|\sup_{r>0}|A_{\delta,r}f|\Big\|_p\le C\delta^{-\sigma}\|f\|_p,$$
uniformly for $\psi$ in bounded subsets of $C_c^\infty((\frac12,2))$.

*Proof.*
1. **One octave of radii.** For $1\le r\le2$ the shell lies in $\frac12\le|\xi|\le2$. Take $g=Q_0f$, an annular projection equal to $1$ there. Fourier inversion gives
$$A_{\delta,r}f=\frac{\delta r}{2\pi}\int e^{\mathrm itr}\hat\psi(\delta rt)\,U(t)g\,dt .$$
   So $\sup_{1\le r\le2}|A_{\delta,r}f|\le C_M\delta\int(1+\delta|t|)^{-M}|U(t)g|\,dt$.
2. **Sum over times.** On $R\le|t|\le2R$, Hölder in $t$ and Lemma 6.1 bound the $L^p$ norm by $C_M\delta(1+\delta R)^{-M}R^{1+\sigma}\|g\|_p$. Summing over dyadic $R\ge1$, with $M>1+\sigma$, gives $O(\delta^{-\sigma}\|g\|_p)$; the range $|t|\le1$ gives $O(\delta\|g\|_p)$.
3. **All octaves.** For $r=2^ks$, $A_{\delta,r}f(2^{-k}y)=A_{\delta,s}g_k(y)$ with $g_k=(Q_kf)(2^{-k}\cdot)$. So the supremum over the octave $[2^k,2^{k+1}]$ has norm $\le C\delta^{-\sigma}\|Q_kf\|_p$. Bound the supremum over all octaves by the pointwise $\ell^p$ sum. Since $p\ge2$, Littlewood–Paley gives $\sum_k\|Q_kf\|_p^p\le\|(\sum_k|Q_kf|^2)^{1/2}\|_p^p\le C\|f\|_p^p$. ∎

**Theorem 6.3 (Bochner–Riesz conjecture; paper Thm 10.1).** For $1\le p\le\infty$ and $\alpha>\delta(p)$, $\sup_{r>0}\|B^\alpha_rf\|_p\le C\|f\|_p$. For $2\le p<\infty$ and $\alpha>s_+(p)$, also $\|\sup_r|B^\alpha_rf|\|_p\le C\|f\|_p$.

*Proof.*
1. **$2\le p<\infty$.** Choose $s_+(p)<\sigma<\alpha$ and a dyadic partition $\sum_j\eta(2^jv)=1$ for $v>0$, with $\eta\in C_c^\infty((\frac12,2))$. With $v=1-\rho$ and $\delta=2^{-j}\le\frac18$, each term is
$$\delta^\alpha\psi_\delta((1-\rho)/\delta),\qquad\psi_\delta(u)=\eta(u)u^\alpha(2-\delta u)^\alpha,$$
   because $1-\rho^2=\delta u(2-\delta u)$. The $\psi_\delta$ have uniformly bounded seminorms.
   - **The remainder** is a smooth compactly supported radial multiplier that vanishes near $\rho=1$. The supremum of its dilates is $\le CMf$ (Hardy–Littlewood maximal function).
   - **The shells.** By Lemma 6.2 the maximal function is at most $C\|f\|_p+\sum_{j\ge3}2^{-j\alpha}C2^{j\sigma}\|f\|_p\le C\|f\|_p$.
   
   Since $\delta(p)=s_+(p)$ for $p\ge2$, the uniform bound follows.
2. **$1<p<2$.** $B^\alpha_r$ is self-adjoint and $\delta(p)=\delta(p')$; use duality.
3. **$p=1,\infty$.** The hypothesis is $\alpha>\frac{n-1}2$. The kernel of $B^\alpha_1$ is a constant times $|x|^{-n/2-\alpha}J_{n/2+\alpha}(|x|)$. It is bounded and $O(|x|^{-(n+1)/2-\alpha})$, hence integrable, and dilation preserves its $L^1$ norm. ∎

---

## 7. Stein's square function

$\mathcal G_\alpha f(x)=\big(\int_0^\infty|r\partial_rB^\alpha_rf(x)|^2\frac{dr}r\big)^{1/2}$. The threshold is $\frac12+s_+(p)=\max\{\frac12,n(\frac12-\frac1p)\}$.

**Lemma 7.1 (one interval of radii; paper Lemma 11.2).** Let $2<p<\infty$, $\varepsilon>0$, $\sigma=s_+(p)+\varepsilon$, and $H=L^2([1,2],ds)$. Fix $A\in C_c^\infty(\mathbb R\times(0,\infty))$ and $\varphi\in C_c^\infty(\mathbb R)$, and for $0<\delta\le1$ let $T_\delta f(x,s)$ be the spatial multiplier $A(s,|\xi|)\varphi((s-|\xi|)/\delta)$ applied to $f$. Then
$$\|T_\delta f\|_{L^p_x(H)}\le C\delta^{1/2-\sigma}\|f\|_p,\qquad\|T_\delta f\|_{L^2_x(H)}\le C\delta^{1/2}\|f\|_2 .$$

*Proof.*
1. **Without $A$.** Let $V_\delta g(x,s)=\varphi((s-|D|)/\delta)g$, with $s\in\mathbb R$ and $\hat g$ in a fixed annulus. Fourier inversion in $s$ and Plancherel give $\|V_\delta g(x,\cdot)\|_{L^2_s}=(2\pi)^{-1/2}\delta\|\hat\varphi(\delta t)U(t)g(x)\|_{L^2_t}$.
   - On $R\le|t|\le2R$, Hölder in $t$ ($p>2$) and Lemma 6.1 give $\|U(t)g\|_{L^p_xL^2_t}\le CR^{1/2-1/p}R^{1/p+\sigma}\|g\|_p$.
   - With $|\hat\varphi(\delta t)|\le C_b(1+\delta R)^{-b}$, $b>\frac12+\sigma$, summing over dyadic $R$ gives $\|V_\delta g\|_{L^p_x(L^2_s)}\le C\delta^{1/2-\sigma}\|g\|_p$.
2. **Restore $A$.** Fourier invert $A$ in $s$: $T_\delta f=\frac1{2\pi}\int e^{\mathrm is\eta}V_\delta(\hat A_1(\eta,|D|)f)\,d\eta$. Here $\hat A_1(\eta,|D|)$ is an annular multiplier with $L^p$ norm $O_b((1+|\eta|)^{-b})$.
3. **The $L^2$ bound.** For each $\xi$, the symbol is bounded and supported on an $s$-interval of length $O(\delta)$; apply spatial Plancherel. ∎

**Lemma 7.2 (combining dilates; paper Lemma 11.3).** Let $H$ be a Hilbert space, $2<p<\infty$, $0<\delta\le1$ and $N=n+3$. Let $K:\mathbb R^n\to H$ be $C^1$ with $\|K(x)\|+\|\nabla K(x)\|\le C_*\delta^{-N}(1+|x|)^{-N}$. Write $S^{(k)}f=2^{nk}K(2^k\cdot)*f$, and let $B_r$ bound $K*:L^r\to L^r(H)$ for $r=2,p$. Then for finite sequences of Schwartz functions,
$$\Big\|\Big(\sum_k\|S^{(k)}f_k\|_H^2\Big)^{1/2}\Big\|_p\le C\Big(B_2+C_*\delta+[1+\log(1/\delta)]^{\frac12-\frac1p}B_p\Big)\Big\|\Big(\sum_k|f_k|^2\Big)^{1/2}\Big\|_p .$$

*Proof.* Let $b_k=\|S^{(k)}f_k\|_H$, $h=(\sum b_k^2)^{1/2}$ and $W=(\sum|f_k|^2)^{1/2}$, let $M$ be the maximal operator, and let $L=\lceil(N+2)\log_2(2/\delta)\rceil$. Fix a dyadic cube $Q$ of side $2^{-j}$ with centre $y$, and $z\in Q$.
1. **Low scales, $k\le j-L$.** For $x\in Q$, $\|S^{(k)}f_k(x)-S^{(k)}f_k(y)\|\le CC_*\delta^{-N}2^{k-j}MW(z)$. The square sum is $\le CC_*\delta^{-N}2^{-L}MW(z)\le CC_*\delta\,MW(z)$. So these $b_k$ are within that of the constant $c_Q=(\sum_{k\le j-L}b_k(y)^2)^{1/2}$.
2. **High scales, $k\ge j+L$.** Split $f_k=f_k\mathbf 1_{5Q}+f_k\mathbf 1_{(5Q)^c}$.
   - **Local part.** $|Q|^{-1}\int_Q\sum_k\|S^{(k)}(f_k\mathbf 1_{5Q})\|^2\le CB_2^2M(W^2)(z)$, since dilation preserves the $L^2$ bound.
   - **Far part.** It is $\le CC_*\delta^{-N}2^{-(N-n)(k-j)}MW(z)$, with square sum $\le CC_*\delta\,MW(z)$.
3. **Middle scales, fewer than $2L$ of them.** Hölder over this finite set bounds $(|Q|^{-1}\int_Q\sum b_k^2)^{1/2}$ by $(2L)^{1/2-1/p}F(z)$, where $F=(\sum_k[M(b_k^2)]^{p/2})^{1/p}$. By the maximal theorem on $L^{p/2}$, $\|F\|_p^p\le C\sum\|b_k\|_p^p\le CB_p^p\|W\|_p^p$.
4. **Sharp function.** The mean of $|h-c_Q|$ on $Q$ is at most $C\big(C_*\delta\,MW+B_2M(W^2)^{1/2}+L^{1/2-1/p}F\big)(z)$, so the dyadic sharp function of $h$ is bounded by that expression. Since $h\in L^2$ ($\|h\|_2^2\le B_2^2\sum\|f_k\|_2^2$), the Fefferman–Stein inequality and the maximal theorem on $L^p$ and $L^{p/2}$ give the claim. ∎

**Theorem 7.3 (Stein's square function conjecture; paper Thm 11.1).** For $2\le p<\infty$ and $\alpha>\max\{\frac12,n(\frac12-\frac1p)\}$, $\|\mathcal G_\alpha f\|_p\le C\|f\|_p$ for Schwartz $f$.

*Proof.*
1. **$p=2$.** With $u=|\xi|^2/r^2$, Plancherel gives
$$\|\mathcal G_\alpha f\|_2^2=2\alpha^2\int_0^1u(1-u)^{2\alpha-2}du\,\|f\|_2^2=\frac{\alpha}{2\alpha-1}\|f\|_2^2 .$$
2. **$p>2$.** Choose $\varepsilon$ with $\sigma=s_+(p)+\varepsilon<\alpha-\frac12$. Write $r=2^ks$ with $1\le s\le2$ and $\rho=2^{-k}|\xi|$, so that $dr/r\asymp ds$. The differentiated multiplier is
$$m(s,\rho)=2\alpha\rho^2s^{-2\alpha}(s+\rho)^{\alpha-1}(s-\rho)^{\alpha-1}\mathbf 1_{\rho<s}.$$
3. **Near the boundary.** Split $s-\rho$ into dyadic intervals $\delta=2^{-j}$. Each piece is $\delta^{\alpha-1}A(s,\rho)\varphi((s-\rho)/\delta)$, with $A$ fixed and $\varphi\in C_c^\infty((\frac12,2))$.
   - **One scale.** Lemma 7.1 gives $B_2\le C\delta^{1/2}$ and $B_p\le C\delta^{1/2-\sigma}$ for $T_\delta$. Its kernel satisfies the hypothesis of Lemma 7.2 with $C_*$ independent of $\delta$: $N$ integrations by parts in frequency cost $\delta^{-N}$.
   - **All dilates.** Apply Lemma 7.2 to $P_kf$, where the $P_k$ are dilated annular projections equal to $1$ on the support. With Littlewood–Paley, the piece contributes at most
$$C\delta^{\alpha-1}\big(\delta^{1/2}+\delta+[1+\log(1/\delta)]^{\frac12-\frac1p}\delta^{1/2-\sigma}\big)\|f\|_p.$$
     This is summable in $j$ because $\alpha-\frac12-\sigma>0$.
4. **Away from the boundary.** The remaining multiplier is smooth, supported in $\rho\le2$, and vanishes to second order at $0$. Split it into annuli $\rho\asymp2^\ell$, $\ell\le2$, rescale, and apply Lemma 7.2 with $\delta=1$. The factor $2^{2\ell}$ makes the sum over $\ell$ finite.
5. **Extension.** For Schwartz $f$ the finite sums converge; Fatou and monotone convergence give the full estimate. ∎

---

## 8. Status

| Target (Lean name) | Here | Paper | Status | Lean area (`OAI/Analysis/`) |
|---|---|---|---|---|
| `euclidean_local_smoothing` | [08](08_GAPLESS_PROOF.md) Thm 1.1 | Thm 1.1 | proved | `Fourier/LocalSmoothing` |
| `spherical_restriction` | Lemma 1.1, Cor. 2.1 | Lemma 9.2, Cor. 9.3 | proved; Tao's $\varepsilon$-removal cited | `Fourier/Restriction` (incl. `Sparse`) |
| `kakeya_maximal` | Lemma 3.1, Cor. 4.1 | Lemma 9.4, Cor. 9.1 | proved | `MeasureTheory/Kakeya` (`CapReduction`, `Randomization`, `Jensen`, `MaximalDiscrete`) |
| `kakeya` (Hausdorff dimension) | Thm 5.1 | Cor. 9.1 | proved | `MeasureTheory/Kakeya` (`Hausdorff`, `HausdorffCovers`, `DyadicCover`) |
| `bochner_riesz` | Lemmas 6.1–6.2, Thm 6.3 | §10 | proved; maximal and Littlewood–Paley theorems cited | `Fourier/BochnerRiesz` |
| `stein_square` | Lemmas 7.1–7.2, Thm 7.3 | §11 | proved; Fefferman–Stein cited | `Fourier/SquareFunction` |

**Not re-proved here:** cone restriction (Nicola's slicing), Nikodym (§9.4), the paraboloid and Schrödinger estimates (§12), and Apps. B–E. None of them is among the formal targets.

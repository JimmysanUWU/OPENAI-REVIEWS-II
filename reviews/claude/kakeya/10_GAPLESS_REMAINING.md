# 10. The remaining claims: cone, Nikodym, paraboloid, Schrödinger, the weak-norm obstruction, and Apps. B, D, E

Claude, 28 October 2026. [08](08_GAPLESS_PROOF.md) and [09](09_GAPLESS_CONSEQUENCES.md) prove the paper's six formal targets. This file treats everything else the paper claims:
- **Gapless proofs** of the abstract-level claims not yet covered:
  - restriction for the truncated cone (Cor. 9.3, (9.4));
  - the Nikodym maximal conjecture ((9.7));
  - paraboloid restriction and Schrödinger local smoothing (Cor. 12.1);
  - the weak-norm obstruction (Prop. 8.2, App. C).
- **A verification report** on the three refinement appendices: B (angular endpoint), D (a linear proof for $n=2$) and E (curved characteristic sheets). They are not in the abstract's list of results.

**Verdict.** I found no error. §1–§4 below are complete, given [08](08_GAPLESS_PROOF.md), [09](09_GAPLESS_CONSEQUENCES.md) and the cited theorems listed in each section. Apps. B, D and E are correct at the level I could check: every exponent, every new algebraic identity, and the structure of each argument. Their analytic estimates reuse §4–§7 of the paper with stated changes, which I followed but did not rewrite line by line (§5).

**Checks.** [`scripts/verify_remaining_claims.py`](scripts/verify_remaining_claims.py) ([`results/remaining_claims.json`](results/remaining_claims.json), ~35 s):
- **App. C's tube arrangement** (exact interval unions, $k\le16$): $\alpha\delta/2$-separation on the input side, and $k|B_X|/\alpha\le0.67$ on the output side.
- **The Nikodym projective map** (symbolic, $n=2,3,4$).
- **The parabolic rescaling and Fresnel identity** of §12, including the vanishing of every power of $\delta$ and $\lambda$ for $d=1..6$.
- **The cone's slicing exponents.**
- **Example E.2 and formula (E.2)**, against finite differences of $\lambda_{\max}$ on random pencils.
- **App. D's quaternionic identities.** The completed square, and the rank-one expansion of the Moore determinant on 200 random instances.

Notation as in [09](09_GAPLESS_CONSEQUENCES.md): $d=n-1$, $p_c=2n/d$, $L(R)=\log(2+R)$, and $A$ is the logarithmic power of Theorem 1.2 ($A=2/p_c$ by [08](08_GAPLESS_PROOF.md)).

---

## 1. The weak-norm obstruction (paper Prop. 8.2, App. C)

**Proposition 1.1.** For $2<p\le p_c$ and all large $N$ there is a Schwartz $g_N\ne0$ with $\operatorname{supp}\hat g_N\subset\{\frac54N<|\xi|<\frac74N\}$ and
$$\|U(t)g_N\|_{L^{p,\infty}(\mathbb R^n\times I)}\ \gtrsim_{n,p}\ (\log N)^{d(1/2-1/p)}\,\|g_N\|_p .$$
At $p=p_c$ the exponent is $1/p_c$. So Sobolev order $0$ fails even for the weak norm throughout $2<p\le p_c$.

*Proof.* Write $x=(X,Y)\in\mathbb R\times\mathbb R^d$. Put $k=\lfloor\frac14\log_2N\rfloor$, $\delta=2^{-k}$, $\alpha=1/(4\sqrt d)$, $w=\alpha\delta/16$, and $a_\ell=1-\ell/k$.
1. **Lines.** For $0\le j<2^k$, with $j\delta=\sum_{\ell=1}^k\varepsilon_\ell(j)2^{-\ell}$ in binary, put
$$L_j(X)=\alpha\sum_{\ell=1}^k\varepsilon_\ell(j)2^{-\ell}(X-a_\ell),\qquad\text{slope }s_j=\alpha j\delta .$$
2. **Separated on the input side.** Passing from $j$ to $j+1$ turns one digit $m$ from $0$ to $1$ and the later digits from $1$ to $0$. Since $2^{-m}-\sum_{\ell>m}2^{-\ell}=\delta$,
$$L_{j+1}(X)-L_j(X)=\alpha\delta(X-a_m)+\alpha\sum_{\ell>m}2^{-\ell}(a_\ell-a_m)\le-\alpha\delta/2\qquad(X\le-\tfrac12),$$
   because $a_\ell<a_m$ for $\ell>m$ and $X-a_m\le-\frac12$. So the intervals of radius $w$ about the $L_j(X)$ are disjoint ($\alpha\delta/2=8w$).
3. **Compressed on the output side.** Let $B_X$ be the union of the intervals of radius $w/2$ about the $L_j(X)$, for $X\in J=[\frac14,\frac12]$.
   - Freeze the first $m=\lfloor k(1-X)\rfloor\le\frac34k$ digits. For $\ell>m$, $0\le X-a_\ell\le(\ell-m)/k$, so the tail ranges over an interval of length $\le\alpha\sum_{\ell>m}2^{-\ell}(\ell-m)/k=2\alpha2^{-m}/k$.
   - Hence $|B_X|\le2^m(2\alpha2^{-m}/k+w)=2\alpha/k+2^mw\lesssim\alpha/k$, since $2^mw\le\alpha2^{-k/4}/16$.
4. **$d$ coordinates.** Choosing lines independently in each transverse coordinate gives $\delta^{-d}$ lines $L_\nu$ with slopes $s_\nu$, $|s_\nu|\le\sqrt d\,\alpha=\frac14$. The cores
$$S_\nu=\{X\in J,\ |Y_i-L_{j_i}(X)|\le w/2\}$$
   lie in $E=\{X\in J,\ Y\in B_X^d\}$. Then $|E|\lesssim_nk^{-d}$ and $\sum_\nu|S_\nu|=\delta^{-d}|J|w^d=|J|(\alpha/16)^d\asymp_n1$.
5. **Packets.** Take $0\le\chi,\theta\le1$ smooth, with $\chi$ supported in $(-2,-\frac12)$ and equal to $1$ on $[-\frac32,-\frac34]$, and $\theta$ supported in $(-1,1)$ and equal to $1$ on $[-\frac12,\frac12]$. Put:
   - $b_\nu=\chi(X)\prod_i\theta((Y_i-L_{j_i}(X))/w)$, which have disjoint supports by step 2;
   - $v_\nu=(1,s_\nu)/\sqrt{1+|s_\nu|^2}$, $\lambda=\frac32N$, and $f_\nu=e^{\mathrm i\lambda v_\nu\cdot x}b_\nu$;
   - $Q_N=\beta(|D|/N)$, with $\beta\in C_c^\infty((\frac54,\frac74))$ equal to $1$ near $\frac32$.
6. **Transport.** Claim: uniformly for $1\le t\le2$,
$$U(t)Q_Nf_\nu(x)=e^{\mathrm i\lambda(v_\nu\cdot x-t)}b_\nu(x-tv_\nu)+O_n(N^{-1}w^{-2}).$$
   - Write $\xi=\lambda v_\nu+\eta$. For $|\eta|\le cN$, $\beta=1$ and $|\xi|=\lambda+v_\nu\cdot\eta+O(|\eta|^2/N)$.
   - After the shear $Y_i'=Y_i-L_{j_i}(X)$, $|\hat b_\nu(\eta)|=|\hat\chi(\eta_1+s_\nu\cdot\eta')|\prod_iw|\hat\theta(w\eta_i)|$. So $\int|\eta|^2|\hat b_\nu|\lesssim w^{-2}$ and $\int|\hat b_\nu|\lesssim1$.
   - The phase error contributes $O(N^{-1}w^{-2})$, and the region $|\eta|>cN$ contributes $O(N^{-2}w^{-2})$. Since $w\asymp\delta\ge N^{-1/4}$, the error is $O(N^{-1/2})$.
7. **Arrival.** Let $I_0=[\frac32,\frac74]$, $X\in J$, $t\in I_0$. Since $(v_\nu)_1\in[4/\sqrt{17},1]$, $X-t(v_\nu)_1\in[-\frac32,-\frac34]$. Translation by $tv_\nu$ preserves $Y_i-L_{j_i}(X)$, because $L_{j_i}$ has slope $s_{j_i}=(v_\nu)_{i+1}/(v_\nu)_1$. So $b_\nu(x-tv_\nu)=1$ on $S_\nu\times I_0$, and $|U(t)Q_Nf_\nu|\ge\frac12$ there for large $N$.
8. **Random signs.** Let $f=\sum_\nu\epsilon_\nu f_\nu$. Disjoint supports give $\|f\|_p\lesssim1$ for all signs. Orthogonality of the random signs gives
$$\mathbb E\int_{E\times I_0}|U(t)Q_Nf|^2=\sum_\nu\int_{E\times I_0}|U(t)Q_Nf_\nu|^2\ge\tfrac14|I_0|\sum_\nu|S_\nu|\gtrsim1 .$$
   Fix signs attaining this, and put $g_N=Q_Nf$; then $\|g_N\|_p\lesssim1$.
9. **Weak norm.** For $p>2$ and a set $H$ of finite measure, $\int_H|u|^2\le\int_0^{|H|}u^*(s)^2ds\le\frac p{p-2}|H|^{1-2/p}\|u\|_{p,\infty}^2$. With $H=E\times I_0$ and $|E|\lesssim k^{-d}$, this gives $\|U(t)g_N\|_{p,\infty}\gtrsim k^{d(1/2-1/p)}$, and $k\asymp\log N$. ∎

The inputs also have bounded $L^2$ and $L^\infty$ norms, so the bound holds with the Lorentz $L^{p,1}$ input norm.

---

## 2. Restriction for the truncated cone (paper Cor. 9.3, (9.4))

**Corollary 2.1.** Let $\Sigma=\{(\rho\omega,\rho):1\le\rho\le2,\ \omega\in S^{n-1}\}\subset\mathbb R^{n+1}$ with surface measure, and $E_\Sigma G(x,t)=\int_\Sigma e^{\mathrm i(x\cdot\xi+t\tau)}G\,d\Sigma$. If $1<a\le\infty$, $p_c<b<\infty$ and $b>\frac{n+1}{n-1}a'$, then $\|E_\Sigma G\|_{L^b(\mathbb R^{n+1})}\lesssim\|G\|_{L^a(\Sigma)}$.

*Proof (Nicola's slicing).*
1. **Slicing.** The cone's slant factor is $\sqrt2$, so $d\Sigma=\sqrt2\rho^{n-1}d\rho\,d\sigma$. With $G_\rho(\omega)=G(\rho\omega,\rho)$,
$$E_\Sigma G(x,t)=\sqrt2\int_1^2e^{\mathrm it\rho}\rho^{n-1}\,E(G_\rho)(\rho x)\,d\rho .$$
2. **Estimates.** At fixed $x$ this is a Fourier transform in $t$. In turn:
   - Hausdorff–Young (valid since $b>p_c>2$) gives $\|E_\Sigma G(x,\cdot)\|_{L^b_t}\lesssim\|\rho^{n-1}E(G_\rho)(\rho x)\|_{L^{b'}_\rho}$.
   - Minkowski's inequality ($b\ge b'$) moves $L^b_x$ inside.
   - Dilation gives $\|E(G_\rho)(\rho\,\cdot)\|_{L^b}=\rho^{-n/b}\|E(G_\rho)\|_{L^b}$.
   - [09](09_GAPLESS_CONSEQUENCES.md) Cor. 2.1 gives $\|E(G_\rho)\|_b\lesssim\|G_\rho\|_{L^a(S^{n-1})}$.
   
   So $\|E_\Sigma G\|_b\lesssim\big\|\|G_\rho\|_{L^a(S^{n-1})}\big\|_{L^{b'}_\rho([1,2])}$.
3. **Hölder.** $b>a'$ implies $a>b'$. So Hölder on $[1,2]$ bounds this by $\|G\|_{L^a(\Sigma)}$. ∎

---

## 3. The Nikodym maximal conjecture (paper (9.7))

Let $N_\delta f(x)$ be the supremum of the averages of $|f|$ over $\delta$-neighbourhoods of unit segments through $x$. Let $K_\delta$ be the Kakeya maximal function of [09](09_GAPLESS_CONSEQUENCES.md) Cor. 4.1, which gives $\|K_\delta f\|_{L^n(S^{n-1})}\lesssim L(1/\delta)^{2A}\|f\|_n$.

**Theorem 3.1.** $\|N_\delta f\|_{L^n(\mathbb R^n)}\lesssim L(1/\delta)^{2A}\|f\|_{L^n(\mathbb R^n)}$ for $0<\delta\le1$.

*Proof.*
1. **Half segments.** Let $D_\rho f(x)$ be the supremum of averages over $\rho$-neighbourhoods of $\{x+u\omega:\frac12\le u\le1\}$. By finitely many rotations, assume $|\omega_n|\ge c_n$.
2. **Slices.** Let $x=(x',t_0)$ lie in a unit cube $Q$ about $0$, and fix $t_0$. For small $\rho$, every point $(y',s)$ of such a tube has $c_n/4\le|s-t_0|\le2$ and $|y'|\le C_n$.
3. **The projective map.** Consider $\Pi(y',s)=(y'/(s-t_0),\,1/(s-t_0))$ on each side of $s=t_0$. Then:
   - Its Jacobian determinant is $-(s-t_0)^{-(n+1)}$, bounded above and below on the region.
   - It sends the line $y'=x'+(s-t_0)v$ to $Y'=v+Y_nx'$, so every segment through $x$ goes to a segment of direction $(x',1)/|(x',1)|$, which is independent of $v$.
   - $\Pi$ and its inverse are uniformly Lipschitz on the region. So the image of a $\rho$-tube lies in $O(1)$ unit tubes of width $C_n\rho$ in that direction.
4. **Transport.** Let $\tilde f_{t_0}$ be $f\,\mathbf 1_{5Q\cap\{c_n/4\le\pm(s-t_0)\le2\}}$ transported by $\Pi$ and extended by zero. It is a single input for all $x'$ in the slice, and $\|\tilde f_{t_0}\|_n\lesssim\|f\|_{L^n(5Q)}$. Hence $D_\rho f(x',t_0)\lesssim K_{C\rho}\tilde f_{t_0}\big((x',1)/|(x',1)|\big)$.
5. **Integration.** The direction map $x'\mapsto(x',1)/|(x',1)|$ has area Jacobian $(1+|x'|^2)^{-n/2}$, which is bounded below for $|x'|\le C$. So
$$\int|D_\rho f(x',t_0)|^ndx'\lesssim\int_{S^{n-1}}|K_{C\rho}\tilde f_{t_0}|^n\lesssim L(1/\rho)^{2An}\|f\|_{L^n(5Q)}^n .$$
   Integrate in $t_0$, then translate and sum over cubes. For $\rho$ bounded below, average over a fixed ball instead.
6. **Full segments.** Split a unit segment through $x$ into:
   - the part within $O(\delta)$ of $x$, whose weight relative to the tube is $O(\delta)$, bounded by $\delta\,Mf(x)$;
   - pieces at distance $[s/2,s]$ for dyadic $\delta\le s\le1$, each of weight $O(s)$, and by dilation about $x$ controlled by $D_{\delta/s}$ with the same $L^n$ operator norm.
   
   So $\|N_\delta\|_{L^n\to L^n}\lesssim\delta+\sum_ss\,L(s/\delta)^{2A}\lesssim L(1/\delta)^{2A}$. ∎

**Remarks (cited, not re-proved).**
- The Hausdorff-dimension conclusion for Nikodym-type sets follows by the transverse-disk argument of [09](09_GAPLESS_CONSEQUENCES.md) Thm 5.1.
- The extension to manifolds of constant curvature uses Gao–Liu–Xi, Prop. 2.3.
- The curved-family application uses Nadjimzadah's lifting theorem [59, Thm 1.3] with full Kakeya dimension in $n+\operatorname{rank}(B^2)$ dimensions.

---

## 4. Paraboloid restriction and Schrödinger local smoothing (paper Cor. 12.1)

Let $d\ge1$ and $p_c=2(d+1)/d$, which is the critical exponent for $S^d\subset\mathbb R^{d+1}$. Put $Eg(x,t)=\int_{B^d}e^{\mathrm i(x\cdot v+t|v|^2)}g(v)\,dv$.

**Corollary 4.1.**
1. **(12.1)** For $1<a\le\infty$ and $p_c<b<\infty$ with $b>(d+2)a'/d$: $\|Eg\|_{L^b(\mathbb R^{d+1})}\lesssim\|g\|_{L^a(B^d)}$.
2. **(12.2)** For a compact interval $I$, $N\ge2$, and $\operatorname{supp}\hat f\subset B(0,N)$: $\|e^{\mathrm it\Delta}f\|_{L^{p_c}(\mathbb R^d\times I)}\lesssim[\log(2+N)]^{2/p_c}\|f\|_{p_c}$.
3. **(12.3)** For $p_c<p<\infty$ and $p\le r\le\infty$: $\|e^{\mathrm it\Delta}f\|_{L^p(\mathbb R^d;L^r(I))}\lesssim\|f\|_{B^s_{p,p}}\lesssim\|(1-\Delta)^{s/2}f\|_p$, with $s=d(1-\frac2p)-\frac2r$.

*Proof.*
1. **Parabolic rescaling.**
   - Parametrize a cap of $S^d$ by $(\delta v,\sqrt{1-\delta^2|v|^2})$, $v\in B^d$, with density $g(v)\sqrt{1-\delta^2|v|^2}$. Then surface measure times density is $\delta^dg(v)\,dv$, and the density's $L^{p}$ norm is $\asymp\delta^{d/p}\|g\|_p$.
   - Evaluating the spherical extension at $(x/\delta,-2t/\delta^2)$ and removing the phase $e^{-2\mathrm it/\delta^2}$ gives the extension for the graph $q_\delta(v)=2(1-\sqrt{1-\delta^2|v|^2})/\delta^2=|v|^2+\delta^2|v|^4/4+\dots$
   - The change of variables has Jacobian $2\delta^{-(d+2)}$ and maps $B(0,R)$ into $B(0,CR/\delta^2)$.
   - With [09](09_GAPLESS_CONSEQUENCES.md) Lemma 1.1 in dimension $d+1$, the $q_\delta$-extension has local norm $\lesssim\delta^{-d+(d+2)/p+d/p}L(R/\delta^2)^A$ on $B(0,R)$. The power of $\delta$ is $0$ at $p=p_c$.
2. **Replacing $q_\delta$ by $|v|^2$.** This costs at most $R\cdot O(\delta^2)\|g\|_1$ pointwise on $B(0,R)$. With $\delta=(2R)^{-(d+2)}$, its $L^{p_c}(B(0,R))$ norm is $O(\|g\|_{p_c})$, and $L(R/\delta^2)\lesssim L(R)$. Modulating $g$ translates the ball, so
$$\|Eg\|_{L^{p_c}(B(z,R))}\lesssim L(R)^A\|g\|_{p_c}\qquad(z\in\mathbb R^{d+1},\ R\ge1).\tag{4.1}$$
3. **(12.1).** (4.1) gives $O_\varepsilon(R^\varepsilon)$ loss. A smooth density on a slightly larger patch has Fourier decay $O(|\cdot|^{-d/2})$ by stationary phase. Then proceed as in [09](09_GAPLESS_CONSEQUENCES.md) Cor. 2.1:
   - Tao's $\varepsilon$-removal theorem, in its compact-hypersurface form (Kim's Thm 5.3, cited), gives $E:L^s\to L^s$ for every $s>p_c$.
   - Interpolate with $E:L^1\to L^\infty$.
4. **(12.2), far part.** Let $S(t)=e^{\mathrm it\Delta}$, $p=p_c$, $I_0=[\frac12,1]$, and tile $\mathbb R^d$ by cubes $Q$ of side $N$. For $|z|\ge8N$, integrating by parts in $\eta=\xi/N$ bounds the kernel of $S(t)\chi(D/N)$ by $C_MN^d(1+N|z|)^{-M}$: the phase gradient is $\ge N|z|/2$. So the contribution from outside $30Q$ has $L^p$ norm $O(\|f\|_p)$ by Young's inequality.
5. **(12.2), near part: the Fresnel identity.** Let $h=f\mathbf 1_{30Q}$ (centered at $0$), $\lambda=16\sqrt dN$ and $g(v)=h(\lambda v)$. Then
$$|S(t)h(z)|=(4\pi t)^{-d/2}\lambda^d\,\Big|Eg\Big(-\frac{\lambda z}{2t},\frac{\lambda^2}{4t}\Big)\Big| .$$
   - The map $(x,t)\mapsto(-\lambda z/2t,\lambda^2/4t)$, with $z=x-w/N$, is one-to-one with Jacobian $\asymp\lambda^{d+2}$. Its image lies in a ball of radius $C(N^2+|w|)$.
   - With $\|g\|_p=\lambda^{-d/p}\|h\|_p$, (4.1) gives $\|S(t)h(x-w/N)\|_{L^p(Q\times I_0)}\lesssim\lambda^{d-(2d+2)/p}L(N^2+|w|)^A\|h\|_p$. The power of $\lambda$ is $0$.
   - Integrate against $|\check\chi(w)|dw$, sum over $Q$ using bounded overlap, and add the far part. This gives (12.2) on $I_0$.
   - Parabolic scaling extends it to $[a/2,a]$; then sum over dyadic intervals and conjugate for negative times.
6. **(12.3).** Choose $p_c<a<p$ and take $b=a$ in (12.1). The diagonal case of Lee–Rogers–Seeger, Cor. 1.2 (cited) turns this $L^a\to L^a$ extension bound into the first inequality. The second follows from Littlewood–Paley theory and $\|\cdot\|_{\ell^p}\le\|\cdot\|_{\ell^2}$. ∎

---

## 5. The refinement appendices B, D, E: what was verified

> **Update (28 Oct, later).** All three appendices are now proved in full:
> - **B** in [14](14_GAPLESS_ANGULAR_ENDPOINT.md), line by line, including Lemma B.4's uniform stationary phase;
> - **D** in [13](13_PLANAR_ENTRY_POINT.md);
> - **E** in [11](11_FRAMEWORK_POSITIVITY_ENGINE.md), Corollary 11.3.
>
> The table below records the earlier, exponent-level check.

| Appendix | Claim | Verified | Not re-derived |
|---|---|---|---|
| **B** | Theorem B.1: at $p_c$, $\|U f\|_{L^p}\lesssim\|L(J)^{1/p}f\|_p+\|L(\Lambda_\omega)^Bf\|_p$ for $B>\frac12+\frac1p$, and a weak-norm bound with angular weight only | Every exponent, listed below the table | The uniform-in-degree stationary phase of Lemma B.4 (it follows Stein Ch. VIII), and the vector-valued maximal and spherical-multiplier theorems (cited) |
| **C** | Prop. 8.2 | **Fully** (§1) | — |
| **D** | The $n=2$ case of Thm 1.2 via a *linear* Dirichlet problem on a disk | **(i)** The completed square $A{:}H-D_B(H)=\beta\lvert b_{01}+\alpha_{01}/\beta\rvert^2+\det\alpha/\beta$ (checked on random quaternionic data). **(ii)** The rank-one expansion $\det(\lambda B+H)=\det H+\lambda\beta\det(H\vert_{\ker B})$ for Moore determinants (checked). **(iii)** The maximum principle and comparison steps, which are standard for $\Delta_W-V$ with $V\ge0$ | The multiplier bounds it imports from (5.5) and (5.11). An alternative route only, since 08 covers $n=2$ |
| **E** | Theorem E.1: the same $\log^{2/p_c}$ bound for $e^{-\mathrm iT\phi(D)}$, with $\phi$ a simple extreme eigenvalue of a symmetric pencil | **(i)** Formula (E.2) $D^2\phi[w,w]=2(Lw)^tS_\perp^{-1}Lw$ (second-order perturbation; checked on random pencils). **(ii)** Group velocity $B(v(\xi))=\nabla\phi$ (Hellmann–Feynman; checked). **(iii)** Example E.2: top eigenvalue $1$, transverse Hessian $\operatorname{diag}(1,\frac23)$, and no linear factor of the characteristic cubic over $\mathbb C$ (exact). **(iv)** The volume $dP\asymp a^2s^{2d}dW\,dQ$ and the exponent identity $h^2s^{2d}(hs)^{-dq}=\sigma^{n-2}$. **(v)** The E.4 lower bound, $\int\delta^{d-dp/2}d\delta=\log N$ | The fixed-cap estimate (E.8)–(E.11), which re-runs §4–§7 with the listed changes |

**Exponents checked in App. B.**
- $\alpha p=n$ and $d-\alpha p=-1$.
- (B.5): $\int_{N^{-1}}^4r^{d-\alpha p}dr=\log(4N)$, and $|\{r^{-\alpha}G>\lambda\}|_{r^ddr}=G^p/(n\lambda^p)$.
- (B.4): $(N^{\alpha+1}b^n)^{1-2/p}=N^{-(\alpha/2)(1-2/p)}$ at $b=N^{\eta-1}$, $\eta=\alpha/2n$.
- $\|U(t)\Pi_N\|_{1\to\infty}\lesssim N^{\alpha+1}$.
- The (B.8) Minkowski step: the $j$th radial factor contributes $(1+j)^{2/p}$.
- The summation condition $B>A+\frac12-\frac1p$, which equals $\frac12+\frac1p$ at $A=2/p$.
- Lemma B.3's comparison error $R^{2+n/p}/S\le1/(100R)$.

**Status of the paper after Parts 08–10.**
- **Complete proofs:** local smoothing, Bochner–Riesz (with maximal bounds), Stein's square function, restriction for sphere, cone and paraboloid, Schrödinger local smoothing, the Kakeya and Nikodym maximal conjectures, full Hausdorff dimension of Kakeya sets, and both sharpness examples (Prop. 8.1 in [02](02_VERIFICATION_LEDGER.md), Prop. 8.2 here). Each is proved here, or follows from cited textbook theorems with hypotheses checked.
- **Checked at the level of exponents, identities and structure:** the three refinements B, D, E.

# 2. Verification ledger (Kakeya / local smoothing paper)

Status key:
- **R**: re-derived by hand.
- **X**: checked symbolically or numerically by a script in [`scripts/`](scripts/), output in [`results/`](results/).
- **L**: covered by the Lean development (see [06](06_LEAN_AUDIT.md)).
- **A**: accepted (standard or cited).
- **P**: partially checked.

Notation follows the paper: $d=n-1$, $p_c=2n/d$, $q=p_c/2=n/d$, $h=N^{-1}$, $a=Kh$, $s=|\tau|$, $\sigma=s/h$.

## §2 The conserved current

| Item | Content | Status | Notes |
|---|---|---|---|
| (2.1) | null coordinates; mode phase $\kappa(s_0+v\cdot x-\tfrac12|v|^2t)$ | R, X | $\kappa=(\xi_b+|\xi|)/\sqrt2$ and $v=\xi_x/\kappa$. Exact for every real $\xi$ with $\xi_b+|\xi|>0$. |
| (2.2) | symbol $-C(p)=\kappa\begin{pmatrix}|v|^2&-v^t\\-v&I\end{pmatrix}$, kernel $E(v)=(1,v)$, $\mathscr P f_{\rm sys}=0$ | R, X | $\mathscr P(\partial_{s_0}\psi,\nabla_x\psi)$ has first entry $(2\partial_t\partial_{s_0}+\Delta_x)\psi$, which is the wave operator. |
| (2.3)–(2.4) | two-state flux; Gaussian weight $e^{p\cdot q}=\exp(-\tfrac\kappa2|\zeta-vR|^2)$ | R, X | The sum of *two* real squares gives the complex modulus. This is the reason for two auxiliary flows. |
| Lemma 2.1 | $B=\operatorname{Re}\sum_\ell\mathcal Z_\ell\otimes\overline{\mathcal Z_\ell}$ is PSD with divergence-free columns | R, X | $\operatorname{div}\mathcal Z_\ell=2iw_\ell^t\mathscr P G=0$. The acceleration $\nabla_{\mathcal Z_\ell}\overline{\mathcal Z_m}=-i\delta_{\ell m}\sum_\mu(G^tC^\mu\bar G)\partial_\mu$ is purely imaginary, so it drops out of the real part. Numerically $\max|\operatorname{div}|\approx10^{-16}$ for $n=2,3$. Negative control ($X+iq$ in place of $X-iq$): divergence ≈ 1. |
| §2.3, (2.7)–(2.8) | quaternionic Gram form; the holomorphic flux map preserves column divergence | R | Second-derivative term: $\sum B^{IL}\partial_I\partial_L\Theta^K=0$. Holomorphy in $z=w_1+\mathbf kw_2$ gives $\partial_{w_{1a}w_{1b}}+\partial_{w_{2a}w_{2b}}=0$; the two state blocks of $B$ are equal and the mixed block vanishes. |
| (2.9)–(2.13) | $\Theta=(R^2,R\zeta,\eta-\zeta\cdot\zeta/2)$, $J=4|R|^{2d+2}$, $ME(v)=2RV^v$, $(\partial\Phi)V^v=0$, $(\partial W)V^v=(v-W)/2Y$ | R, X | Exact symbolic check for $d=1,2$. |
| Prop 2.2 | wave current $A^0=\chi\chi^*$, $A^0_{TT}=|\chi_T|^2$ | R | Gram rows multiplied by a unit quaternion; the Gram matrix is unchanged. |
| (2.15) | $A^0_{TT}=|S_0u_w|^2+|S_1u_w|^2$ | R | Orthogonal planes $\mathbb C_{\mathbf i}$ and $\mathbf k\mathbb C_{\mathbf i}$. |

## §3 Quaternionic Hessian calculus

| Item | Status | Notes |
|---|---|---|
| Lemma 3.1: $H(\xi\xi^t)=gg^*$, $\operatorname{tr}H=\Delta$, product rule (3.2), $H(\sum|F_\nu|^2)=4J^*J$ (3.4) | R, X | Symbolic, in 3 and 4 slots. |
| Lemma 3.2: Moore determinant; (3.6) $C(A)=\det A\cdot A^{-1}$; (3.7) $v^*C(A)v=\det(A|_{v^\perp})$ and exact linearity of $t\mapsto\det(A+tvv^*)$ (also for indefinite $A$); (3.8) Schur; concavity of $\log\det$ | R, X | Numerical, $m=2,3,4$, via the complex $2m\times2m$ representation. The explicit $3\times3$ Moore formula was cross-checked against eigenvalues. |
| Lemma 3.3: cofactor divergence $\sum\bar e_\alpha\partial_{a\alpha}C(H(u))_{ab}=0$ | X, A | Exact symbolic identity for random quartic $u$: 3 slots ($m=2,3$), 4 slots ($m=2$), 2 slots ($m=3$). Negative control: adding a non-Hessian term breaks it. Alesker [4, Lemma 4.6] covers 4 slots; 3 slots follows by restricting to functions independent of the fourth slot. |
| (3.10), (3.13) Green's identity; (3.14) boundary compression | R | |

## §4 Determinant energy

| Item | Status | Notes |
|---|---|---|
| Lemma 4.1 (fiber geometry) | R | With $F=(W/w_0,10Q)$ and $r=|F|^2-1$: $H(r)=4J^*J$, and the null line lies in $G(r)^\perp\cap\ker H(|\Phi|^2)$. For the bounds: $\partial W/\partial P=1/Y$, $|Y|\asymp s$, and $\partial_Z(10\Phi/a)=10(1-W\cdot W/2)/(\sqrt2a)$, so the smallest singular value of $J_{\rm fib}$ is $\gtrsim s^{-1}$ and $|\det J_{\rm fib}|\asymp a^{-1}s^{-d}$. Hence (4.3)–(4.4). $G G^*=4J^*FF^*J\le4|F|^2J^*J\le H(r)$ on the fiber. |
| (4.5)–(4.8): regularized current $A$, density $\mu$, background $B_0=\lambda|\Phi|^2$, Dirichlet potential $u$ | R | $A\ge0$, $\mu>0$, $A:H(|\Phi|^2)=0$. $\mu^{1/d}\le\lambda=h^{-30n}$ from $\|f\|_2=1$ and $O(h^{-n})$ modes. |
| (4.9)–(4.10): $A:H(u)=A:\Sigma+\lambda^{-1}\det H(B_0+u)$ | R | Schur complement (3.8) and $\det D_u=\lambda\mu$. |
| **Prop 4.2 (energy balance)**: $\Delta_z\mathcal F+\int A:\Sigma=\int\rho\,G(r)^*A_wG(r)\,\delta(r)$ | R | The two identities: $\Delta X_0=\mathcal B-\int A:\Sigma-\lambda^{-1}\int\det$, and $\Delta E_z(u)=\lambda^{-1}\int\det H(B_0+u)$. The second uses (a) the full cofactor $C(H(B_0+\theta u))$ is divergence-free; (b) its normal coefficient vanishes because $H(B_0+\theta u)|_{\nu^\perp}$ contains the null line of $H(r)$, which lies in $\ker B$; (c) $\det B=0$ since $B$ has rank one. Verified line by line. |
| Lemma 4.3: $\mathcal F=\frac{n-1}{n+1}X_0+\frac{1}{\lambda n(n+1)}\int G^*C(D_u)G$ | R, X | The segment determinant identity and the weights $1/(n+1)$, $1/(n(n+1))$ were checked numerically. |

## §5 Fixed-fiber estimate (Prop 5.1)

| Item | Status | Notes |
|---|---|---|
| $P_*=(J_{Ds}*|g_T|^2)^{1/d}$; (5.4) derivative bounds, $\int P_*^{d+1}\le\|g_T\|_{p_c}^{p_c}$ | R | Young's inequality with exponent $q=1+1/d$ and a kernel of integral one. |
| (5.5)–(5.6) symbol and kernel bounds for $A_w$, $\mu$; wall gain $e^{-c_0\sigma}$ | P | Mechanism verified: on $r=0$, $|W|\ge\sqrt3w_0/2$ while the data have $|v|\le\gamma<w_0/4$, so $|W-v|$ is bounded below. Constants accepted. |
| (5.7)–(5.8): $H_{\rm fib}(rP_*)\ge(1-C/D)P_*R$; comparison $v_*\le u$; $0\le u/r\le C_1\sigma^{m_1}P_*$ | R | The determinant lower bound uses $\det(B_{\rm fib}+K)\ge(\beta+\kappa)\kappa^d$ for $K\ge\kappa I$. |
| (5.9)–(5.10) energy upper bound and flux bound | R | Volume factor $a^2s^{2d}\cdot s^{-2}s^{-2d}=a^2s^{-2}$. |
| (5.11)–(5.12) interior lower bound | R | Quadratic comparison $v_b=\epsilon m^{1/d}(l^2-h^2)$ with $H(l^2)=6I$ (3 slots), so the center depth is $\ge\epsilon h^2m^{1/d}$. $(−u)\mu\ge\epsilon h^2m^{1+1/d}$; $1+1/d=q$ gives $|g|^{p_c}$. The volume factor $h^2s^{2d}(hs)^{-dq}=\sigma^{n-2}$ was checked. |
| Choice $C_3\gamma^2<c_0/2$, then $K$ | R | Inverting the Gaussian costs $e^{C\gamma^2\sigma}$, which is beaten by the wall decay $e^{-c_0\sigma}$. |

## §6–§7 Propagation and transfer

| Item | Status | Notes |
|---|---|---|
| Lemma 6.1 | R, X | Radial Gronwall bound $\psi(A)/A\le\exp\int_0^A\min(z,A-z)P$, checked numerically for $V=Cs^{-2}e^{-cs/h}$ with $h$ down to $10^{-6}$. Outward passage by $\sin(\pi T/3)$ moments; weak Harnack via the 4D superharmonic $H\cos(\sqrt{M+1}v)$; return passage by $\cos(\pi T/2)u(|\tau|)$. All re-derived. |
| Remark 6.2 (sharpness of two logarithms) | X | $\Delta E=(\cos T-1)p\le0$ verified exactly (sympy), along with $E\ge0$, $E(T,b)=1-\cos T$, and $E(0,2b)\le\tfrac94\log2/D^2$. |
| Prop 6.3 | R | $\int_b^3s\log(s/b)V\,ds\le C\int_1^\infty t^{-1}\log t\,e^{-2cKt}dt$, uniformly in $h$. |
| §7.1 removing the period; §7.2 interpolation | R, X | Riesz–Thorin exponents: $\beta(p)=\min\{d(1-2/p),2/p\}$ is continuous at $p_c$ and matches both interpolations. |
| Cor 7.1 (logarithmic Sobolev) | P | |

## §8–§12 Obstructions and applications

| Item | Status | Notes |
|---|---|---|
| Prop 8.1 focusing example | R, X | Exponent $d-dp/2=-1$ iff $p=p_c$. |
| Prop 8.2 / App C (weak norm obstruction) | R, X | Re-derived in full in [10](10_GAPLESS_REMAINING.md) §1. The tube arrangement is checked exactly for $k\le16$. |
| Lemma 9.2 (wave ⇒ local extension), Cor 9.3 (global extension via Tao's ε-removal and Nicola's slicing) | R, A | |
| Lemma 9.4, Cor 9.1 (Kakeya maximal, Hausdorff dimension) | R | $Ap/q=2A$, $dp=2n$, volume power $2An=2(n-1)$. The Hausdorff argument uses Keich's grouping of radii. |
| §9.4 Nikodym, curved families | R, X, A | Nikodym (9.7): re-derived in [10](10_GAPLESS_REMAINING.md) §3, with the projective map checked symbolically. The curved families cite Nadjimzadah's lifting theorem and Gao–Liu–Xi. |
| Thm 10.1 Bochner–Riesz | R | Thin-shell maximal Lemma 10.2 plus Littlewood–Paley; duality; $p=1,\infty$ by the kernel. |
| Thm 11.1 Stein square function | R | Threshold $\tfrac12+s_+(p)=\max\{\tfrac12,n(\tfrac12-\tfrac1p)\}$ checked. Lemma 11.3 (Seeger-type combination) accepted. |
| §12 paraboloid, Schrödinger | R, X, A | Re-derived in [10](10_GAPLESS_REMAINING.md) §4. The rescaling powers and the Fresnel identity are checked. Tao/Kim $\varepsilon$-removal and Lee–Rogers–Seeger Cor. 1.2 are cited. |

## Appendices

| Item | Status | Notes |
|---|---|---|
| Thm A.1 (periodic Dirichlet problem, 3-slot quaternionic MA) | P, L | Continuity method; comparison via the averaged linearization; $C^0$ by barriers; gradient by differentiating; boundary $C^2$ by CNS with quaternionic rigid motions in an added 4th periodic slot; Evans–Krylov. The structure is standard. Uniform ellipticity holds because the symbol $\hat\xi^*Q^{-1}\hat\xi$ is comparable to $|\xi|^2$ on the solution set. |
| App B (angular endpoint) | — | Not checked. |
| App D ($n=2$ linear version) | P | Read. Notable: for $n=2$ the determinant is linear in the unknown, so no Monge–Ampère theory is needed. |
| App E (symmetric pencils; cubic sheet) | P | Irreducibility of the cubic in Ex. E.2 re-derived. |

## Findings

**No mathematical error found** in the parts checked.

Remarks:
- **K1.** The formalization proves Theorem 1.2 with $(1+\log(1/h))^{2n+5}$ in place of $[\log(2/h)]^{2/p_c}$. The doc's 3D tube bound $c/(1+\log)^{66}$ equals $2n(2n+5)$ at $n=3$, which confirms this. Every headline consequence needs only *some* fixed power, so the formal targets are unaffected. The sharp $2/p_c$ power, and Theorem B.1, remain prose-only.
- **K2.** In Lemma 2.1, $G=\overline{f_{\rm sys}(\mathsf X-\mathrm iq)}$. The overline is easy to miss in the typeset PDF. With the opposite sign of the shift, the tensor is not conserved; my negative control confirms this.
- **K3.** The conservation, positivity and determinant-cancellation identities hold for any constant real symmetric system and any rank-one holomorphic background. The analytic input specific to the wave is confined to the Gaussian (2.4) and the fiber geometry (4.1). This modularity is exploited in [05](05_DISTILLATION.md).

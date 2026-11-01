# 8. Gapless proof: critical local smoothing for the wave equation in every dimension

Claude, 27 October 2026. This file proves Theorems 1.1 and 1.2 of *Local smoothing for the Euclidean wave equation* (OpenAI, 25 Oct 2026) for every $n\ge2$, with the sharp logarithmic power $2/p_c$. The consequences (Bochner–Riesz, restriction, Kakeya, Stein's square function) are derived in [09](09_GAPLESS_CONSEQUENCES.md).

**What is proved here, and what is cited.**
- **Proved here in full:**
  - the quaternionic calculus (§2);
  - the positive conserved wave current (§3);
  - the fibre geometry (§4);
  - the comparison principle (§5);
  - the energy balance and comparison (§6);
  - all fibre estimates, including the multiplier kernel bounds (§7);
  - the propagation lemma (§8);
  - the passage to $\mathbb R^n$ and to Sobolev spaces (§9).
- **Cited, with hypotheses checked:**
  - (a) the spectral theorem for quaternionic Hermitian matrices;
  - (b) the existence and regularity half of Theorem 5.1. This uses Schauder theory and the Evans–Krylov/Krylov $C^{2,\alpha}$ estimates (Gilbarg–Trudinger, Thms. 6.14, 6.19, 17.14, 17.26′). It is also formally verified: `PDE/Dirichlet` and `PDE/Schauder` in OpenAI's Lean proof, 672 modules ([06](06_LEAN_AUDIT.md));
  - (c) classical tools: Young's inequality, Riesz–Thorin, Mihlin's theorem, Littlewood–Paley theory, the weak maximum principle, and the mean value inequality for superharmonic functions.

**Checks.** [`scripts/`](scripts/):
- [`verify_wave_current_target.py`](scripts/verify_wave_current_target.py) (new) checks Proposition 3.3 *directly in the target coordinates*, independently of the push-forward argument used to derive it. The real divergence of $(\chi\chi^*)^{\mathbb R}$ vanishes to $10^{-39}$ in 60-digit arithmetic for $n=2,3$. Two negative controls (no rotation $\iota_e$; no volume factor) fail at order one.
- The other identities are checked by [`verify_wave_current.py`](scripts/verify_wave_current.py), [`verify_quaternion_calculus.py`](scripts/verify_quaternion_calculus.py) and [`verify_propagation_and_exponents.py`](scripts/verify_propagation_and_exponents.py).

**Relation to the paper.** The mathematics is the paper's §2–§7 and App. A. The order differs: calculus first, then the current, so every lemma is available when used. Constants are fixed in one order: $\gamma$, then the kernel width $L$, then $K$, then $h$ small. No mathematical error was found. Two easily missed points are flagged where they occur:
- the conjugation in $G=\overline{f(\mathsf X-\mathrm iq)}$ (Lemma 3.1);
- the independence of $\chi$ from the branch of $\sqrt Y$ (Proposition 3.3, step 4).

**Relation to the formal proof.** OpenAI's Lean proof establishes Theorem 1.2 with $(1+\log(1/h))^{2n+5}$, using a different propagation step. The sharp power $2/p_c$ rests on Lemma 8.1 below, which is proved here in full.

---

## 0. Statements

Throughout:
- $n\ge2$, $d=n-1$, $p_c=\frac{2n}{d}$ and $q=\frac{p_c}2=\frac nd$;
- $U(T)=e^{-\mathrm iT|D|}$, $I=[1,2]$, $J=(1-\Delta)^{1/2}$;
- $s_+(p)=\max\{0,\frac d2-\frac np\}$;
- the Fourier transform is $\hat f(\xi)=\int e^{-\mathrm iX\cdot\xi}f(X)\,dX$.

> **Theorem 1.2 (critical annular estimate).** For $\varphi\in C_c^\infty(\mathbb R^n\setminus0)$, $N\ge2$ and $P_N=\varphi(D/N)$,
> $$\|U(T)P_Nf\|_{L^{p_c}(\mathbb R^n\times I)}\le C_{n,\varphi}\,[\log(2+N)]^{2/p_c}\,\|f\|_{L^{p_c}(\mathbb R^n)}.$$

> **Theorem 1.1 (local smoothing).** For $2\le p<\infty$ and $\varepsilon>0$, $\|U(T)f\|_{L^p(\mathbb R^n\times I)}\le C_{n,p,\varepsilon}\|J^{s_+(p)+\varepsilon}f\|_p$. At $p=2$ one may take $\varepsilon=0$.

**Plan.**
```
§2 quaternionic calculus ─────────────┐
§3 positive conserved current A       ├─ §6 energy F(T,τ):  Δ_{T,τ}F ≤ wall flux,   (n-1)/(n+1)·X0 ≤ F ≤ n/(n+1)·X0
§4 fibre geometry                     │       (X0 = ∫(−u)μ, u the fibre Monge–Ampère potential)
§5 fibre Dirichlet problem ───────────┘
§7 one fibre: F ≲ a²‖f‖^p at |τ|=4a, T=0;  F ≳ a²‖U(T)f‖^p at |τ|=2a;  ΔF ≤ C s^-2 e^(-cs/h) F
§8 two circles: ∫_1^2 ‖U(T)f‖^p dT ≲ log²(2/h) ‖f‖^p      (one cap, on a torus)
§9 caps → annulus → ℝ^n → Theorem 1.2 → interpolation → Theorem 1.1
```

---

## 1. Conventions

**Quaternions.**
- $\mathbb H$ has units $\mathrm i,\mathrm j,\mathrm k$ with $\mathrm i^2=\mathrm j^2=\mathrm k^2=\mathrm{ijk}=-1$, inner product $\langle u,v\rangle=\operatorname{Re}(u\bar v)$, and $|u|^2=u\bar u$.
- $\mathbb C_{\mathrm i}=\mathbb R+\mathbb R\mathrm i$ carries Fourier phases; $\mathbb C_{\mathrm k}=\mathbb R+\mathbb R\mathrm k$ carries auxiliary (imaginary) coordinates. Each quaternion is uniquely $\alpha+\mathrm k\beta$ with $\alpha,\beta\in\mathbb C_{\mathrm i}$, and $\mathbb C_{\mathrm i}\perp\mathrm k\mathbb C_{\mathrm i}$.
- Vectors are columns with scalars acting on the right; $^*$ is the quaternionic conjugate transpose.

**Slots.**
- A *slot* is a quaternionic coordinate $z_a=\sum_{\alpha<k}e_\alpha z_{a\alpha}$ with real $z_{a\alpha}$. Either $k=3$, $(e_0,e_1,e_2)=(1,\mathrm i,\mathrm j)$, or $k=4$, adding $e_3=\mathrm k$. Write $\partial_{a\alpha}=\partial/\partial z_{a\alpha}$.
- The quaternionic gradient and Hessian of a real $C^2$ function are
$$G_a(u)=\sum_\alpha e_\alpha\partial_{a\alpha}u,\qquad H_{ab}(u)=\sum_{\alpha,\beta}e_\alpha\bar e_\beta\,\partial_{a\alpha}\partial_{b\beta}u .$$

**Tensors.** For a Hermitian matrix field $A$:
- its *real tensor* is $A^{\mathbb R}_{a\alpha,b\beta}=\operatorname{Re}(\bar e_\alpha A_{ab}e_\beta)$;
- $A$ is **conserved** if every real column of $A^{\mathbb R}$ has zero divergence;
- $A:B=\operatorname{Re}\operatorname{tr}(AB)$;
- $\det$ is the Moore determinant and $C(A)$ the cofactor (Lemma 2.2).

---

## 2. Quaternionic calculus

**Lemma 2.1 (Hessian calculus).** In $k=3$ or $4$ slots:
1. $H(u)$ is Hermitian and $\operatorname{tr}H(u)=\Delta u$.
2. The real-linear map $D^2u\mapsto H(u)$ sends $\xi\xi^{\mathsf T}$ to $gg^*$, where $g_a=\sum_\alpha e_\alpha\xi_{a\alpha}$. Hence it maps nonnegative real forms to nonnegative quaternionic forms.
3. **Product rule.** $H(uv)=uH(v)+vH(u)+G(u)G(v)^*+G(v)G(u)^*$.
4. **Functions of the imaginary parts.** In 3 slots, write $z_a=x_a+P_a\mathrm i$ with $P_a=p_a+\mathrm kq_a\in\mathbb C_{\mathrm k}$. If $u$ depends only on $P$, then $G_a(u)=2(\partial_{\bar P_a}u)\mathrm i$ and $H_{ab}(u)=4\partial_{\bar P_a}\partial_{P_b}u$, with Wirtinger derivatives for the unit $\mathrm k$.
   - For a holomorphic map $F=(F_1,\dots,F_l)$ of $P$ with complex Jacobian $J_{\nu a}=\partial_{P_a}F_\nu$, this gives $H(\sum_\nu|F_\nu|^2)=4J^*J$.
   - For every smooth real $\rho$, $G(\rho\circ F)=2J^*(\partial_{\bar F}\rho)\,\mathrm i$, which lies in the quaternionic image of $J^*$.
5. **Covariance.**
   - In 4 slots, $H(u\circ U)=U^*(H(u)\circ U)U$ for quaternionic unitary $U$, and $H$ is invariant under translations.
   - In 3 slots, $H$ is covariant under a real orthogonal change of the slot index applied to all three real coordinates.

*Proof.*
1. $\overline{H_{ba}}=\sum e_\beta\bar e_\alpha\partial_{b\alpha}\partial_{a\beta}u=H_{ab}$. In $H_{aa}$ the terms with $\alpha\ne\beta$ pair to $2\operatorname{Re}(e_\alpha\bar e_\beta)=0$.
2. Insert $D^2u=\xi\xi^{\mathsf T}$.
3. Apply 2 to $\nabla u\otimes\nabla v+\nabla v\otimes\nabla u$.
4. $G_a=\mathrm i\partial_{p_a}u+\mathrm j\partial_{q_a}u=(\partial_{p_a}+\mathrm k\partial_{q_a})u\cdot\mathrm i$, since $\mathrm{ki}=\mathrm j$. The coefficients of $(\partial_{p_ap_b},\partial_{p_aq_b},\partial_{q_ap_b},\partial_{q_aq_b})$ in $H_{ab}$ are $(1,-\mathrm k,\mathrm k,1)$, which is $(\partial_{p_a}+\mathrm k\partial_{q_a})(\partial_{p_b}-\mathrm k\partial_{q_b})$.
   - For holomorphic $F$, $\partial_{\bar P_a}\partial_{P_b}\sum_\nu F_\nu\bar F_\nu=\sum_\nu\overline{J_{\nu a}}J_{\nu b}$.
   - The chain rule gives the formula for $G(\rho\circ F)$.
5. **Covariance.**
   - In 4 slots the encoding $\xi\mapsto g$ is a real isometry $\mathbb R^{4m}\cong\mathbb H^m$, and the transpose of $U$ acting on encoded covectors is $U^*$. So 2 gives $H(u\circ U)=U^*(H(u)\circ U)U$ for rank-one Hessians, and linearity extends this to all.
   - In 3 slots, a real orthogonal $O$ acting on the slot index preserves the subspace of covectors with no $\mathrm k$-component. ∎

**Lemma 2.2 (Moore determinant).**
- **Definitions.** Every Hermitian $A\in\mathbb H^{m\times m}$ has a quaternionic orthonormal eigenbasis with real eigenvalues $\lambda_1,\dots,\lambda_m$; this is the spectral theorem (cited). Set $\det A=\prod\lambda_i$. It is a homogeneous polynomial of degree $m$ in the real entries. The cofactor $C(A)$, defined by $C(A):B=\frac d{dt}\det(A+tB)|_{t=0}$, is a Hermitian polynomial matrix of degree $m-1$.
- **Properties.**
  1. $C(A)=(\det A)A^{-1}$ for $A>0$, and $\log\det$ is concave and increasing on $\{A>0\}$.
  2. If $A\ge0$ then $C(A)\ge0$; if moreover $A$ is singular, then $C(A)$ vanishes on the image of $A$.
  3. For a unit vector $v$ and every Hermitian $A$ and real $t$: $v^*C(A)v=\det(A|_{v^\perp})$, and $\det(A+tvv^*)=\det A+t\,v^*C(A)v$. Here $A|_{v^\perp}$ is the compression to the quaternionic orthogonal complement, and the determinant of the $0\times0$ matrix is $1$.
  4. **Schur.** For $M=\begin{pmatrix}\ell&c^*\\c&D\end{pmatrix}$ with $D>0$, put $\sigma=\ell-c^*D^{-1}c$. Then $\det M=\det D\cdot\sigma$, and $M-\sigma e_0e_0^*\ge0$.

*Proof.*
- **Polynomiality.** By Newton's identities, $\prod\lambda_i$ is a polynomial in the power sums $\sum\lambda_i^j=\operatorname{Re}\operatorname{tr}(A^j)$. The cofactor exists by duality for $:$.
- **The real representation.** The real representation $\rho:\mathbb H^{m\times m}\to\mathbb R^{4m\times4m}$ repeats each eigenvalue four times. So $\det_{\mathbb R}\rho(A)=(\det A)^4$ and $\operatorname{tr}\rho(A^{-1}B)=4A^{-1}:B$. This gives 1, monotonicity and concavity.
- **Item 2.** $C(A+\epsilon I)>0$ and continuity give $C(A)\ge0$. In an eigenbasis $C(\operatorname{diag}\lambda)=\operatorname{diag}(\prod_{j\ne i}\lambda_j)$, by 1 for positive diagonals and then by polynomiality. If $\lambda_1=0$, this kills every eigenvector with $\lambda_i\ne0$; unitary invariance finishes the proof.
- **Item 4.** Completing the square shows that the form of $M-\sigma e_0e_0^*$ at $(b,w)$ equals $(w+D^{-1}cb)^*D(w+D^{-1}cb)\ge0$. For $\sigma>0$, $M>0$ and block elimination in $\rho$ give $(\det M)^4=(\det D)^4\sigma^4$; take positive fourth roots. Both sides are polynomials in $\ell$, so the identity holds for all $\ell$.
- **Item 3.** In a basis with $v$ first, adding $tvv^*$ changes only $\ell$. So, by 4, $\det(A+tvv^*)=\det A+t\det(A|_{v^\perp})$ for $A|_{v^\perp}>0$. Differentiate in $t$ and extend by polynomiality. ∎

**Lemma 2.3 (cofactor divergence).** In $k=3$ or $4$ slots, every real $C^3$ function $u$ satisfies
$$\sum_{a,\alpha}\bar e_\alpha\,\partial_{a\alpha}C(H(u))_{ab}=0\qquad(1\le b\le m).$$
Consequently $C(H(u))$ is conserved.

*Proof.* For $m=1$, $C\equiv1$. Otherwise let $j_0=m-1$ and let $Q(M_1,\dots,M_{j_0})$ be the coefficient of $t_1\cdots t_{j_0}$ in $C(\sum_jt_jM_j)$. It is symmetric, multilinear, and $Q(M,\dots,M)=j_0!\,C(M)$.
1. **Rank-one inputs.** For $M_j=g_jg_j^*$ and $t_j>0$, the matrix $\sum t_jg_jg_j^*\ge0$ has rank $\le j_0<m$. By Lemma 2.2(2) its cofactor kills every $g_i$. This is a polynomial identity in $t$, so $Qg_i=0$ for each $i$, and $g_i^*Q=0$.
2. **Exponential inputs.** For $u_j=\exp(\xi_j\cdot z)$, Lemma 2.1(2) gives $H(u_j)=u_jg_jg_j^*$. The row divergence of $Q(H(u_1),\dots,H(u_{j_0}))=e^{(\sum\xi_j)\cdot z}Q(g_1g_1^*,\dots)$ is $e^{(\sum\xi_j)\cdot z}(\sum_jg_j)^*Q=0$.
3. **General inputs.** The row divergence is a multilinear constant-coefficient differential expression in $(u_1,\dots,u_{j_0})$. On exponentials it equals $e^{(\sum\xi_j)\cdot z}$ times a quaternion-valued polynomial in the real covectors $\xi_j$, and that polynomial vanishes identically. So every coefficient vanishes.

Take all $u_j=u$. For the last claim, multiply by $e_\beta$ on the right and take real parts. ∎

**Lemma 2.4 (Green's identity and compression).** Let $A$ be a conserved Hermitian field. Then:
1. $A:H(u)=A^{\mathbb R}:D^2u=\operatorname{div}(A^{\mathbb R}\nabla u)$, and $(\nabla v)^{\mathsf T}A^{\mathbb R}\nabla u=\operatorname{Re}(G(v)^*AG(u))$.
2. **Green's identity.** On $\Omega=\{r<0\}$ with $\nabla r\ne0$ on $\partial\Omega$, possibly periodic in some directions and with test functions of compact support in the unbounded ones,
$$\int_\Omega\big(v\,A:H(u)-u\,A:H(v)\big)=\int_{\partial\Omega}\operatorname{Re}\big(v\,\nu^*AG(u)-u\,\nu^*AG(v)\big)\,dS,\qquad \nu=G(r)/|\nabla r|.$$
3. **Compression.** If $u=0$ on $\partial\Omega$, then $u=\rho r$ with $\rho$ smooth, $G(u)=\nu\,\partial_{\mathbf n}u$, and $H(u)|_{\nu^\perp}=\rho\,H(r)|_{\nu^\perp}$ on $\partial\Omega$.

*Proof.*
1. $\operatorname{Re}\operatorname{tr}(AH(u))=\sum A^{\mathbb R}_{a\beta,b\alpha}\partial_{b\alpha}\partial_{a\beta}u$, and the columns of $A^{\mathbb R}$ are divergence-free.
2. The divergence theorem; the symmetric terms $\nabla v\cdot A^{\mathbb R}\nabla u$ cancel.
3. Hadamard's lemma gives $u=\rho r$. Then $H(u)=\rho H(r)+rH(\rho)+G(\rho)G(r)^*+G(r)G(\rho)^*$. On $r=0$ the cross terms vanish on $\nu^\perp$, because $G(r)\parallel\nu$. ∎

---

## 3. The positive conserved wave current

### 3.1 Null coordinates

Fix an orthonormal frame $X=(x,b)\in\mathbb R^d\times\mathbb R$ and put $t=(T+b)/\sqrt2$, $s_0=(b-T)/\sqrt2$. A plane wave with $\xi_b+|\xi|>0$ is
$$u_{\kappa,v}=c_{\kappa,v}\,e^{\mathrm i\kappa(s_0+v\cdot x-\frac12|v|^2t)},\qquad\kappa=\frac{\xi_b+|\xi|}{\sqrt2},\quad v=\frac{\xi_x}{\kappa}.\tag{3.1}$$
Indeed, $|\xi_x|^2=\sqrt2\kappa(|\xi|-\xi_b)$, so this phase equals $X\cdot\xi-T|\xi|$. Index states by $0,\dots,d$, and let
$$C^t=2E_{00},\qquad C^{x_j}=E_{0j}+E_{j0},\qquad C^{s_0}=-\sum_{j}E_{jj},\qquad\mathscr P=\sum_\mu C^\mu\partial_\mu .$$
At the phase covector $\kappa\ell(v)$, with $\ell(w)=(-\frac12w\cdot w,\,w,\,1)$, the symbol satisfies
$$-C(\kappa\ell(v))=\kappa\begin{pmatrix}|v|^2&-v^{\mathsf T}\\-v&I\end{pmatrix},\qquad -w^{\mathsf T}C(\kappa\ell(v))w=\kappa|w'-vw_0|^2 .\tag{3.2}$$
In particular $C(\kappa\ell(v))E(v)=0$ for $E(v)=(1,v)^{\mathsf T}$, so $f_{\rm sys}=\sum\kappa^{d/2}E(v)u_{\kappa,v}$ solves $\mathscr Pf_{\rm sys}=0$.

For $w_1,w_2\in\mathbb R^n$ put $q^\mu=\frac12\sum_\ell w_\ell^{\mathsf T}C^\mu w_\ell$ and $z=w_1+\mathrm kw_2=(R,\zeta)\in\mathbb C_{\mathrm k}\times\mathbb C_{\mathrm k}^d$. Summing (3.2) over $\ell$ gives
$$\kappa\ell(v)\cdot q=-\tfrac\kappa2|\zeta-vR|^2 .\tag{3.3}$$

### 3.2 The lifted current

**Lemma 3.1 (lifted current).** Let $\mathscr P=\sum_\mu C^\mu\partial_\mu$ have constant real symmetric $N\times N$ coefficients on $\mathbb R^{n+1}$. Let $f$ be a finite sum of exponential solutions of $\mathscr Pf=0$, and with $q$ as above put
$$G=\overline{f(\mathsf X-\mathrm iq)},\qquad V_{\ell a}=\sum_\mu(C^\mu w_\ell)_a\partial_\mu,\qquad D_{\ell a}=\partial_{w_{\ell a}}+\mathrm iV_{\ell a},\qquad\mathcal Z_\ell=\sum_aG_aD_{\ell a}.$$
Then $B=\operatorname{Re}\sum_{\ell=1}^2\mathcal Z_\ell\otimes\overline{\mathcal Z_\ell}$ is real, symmetric and positive semidefinite on $(\mathsf X,w_1,w_2)$, and each of its columns has zero divergence.

*Proof.*
1. **The fields annihilate $\overline G$.** $D_{\ell a}(\mathsf X^\mu-\mathrm iq^\mu)=-\mathrm i(C^\mu w_\ell)_a+\mathrm i(C^\mu w_\ell)_a=0$, so $D_{\ell a}$ kills every holomorphic function of $\mathsf X-\mathrm iq$, in particular $\overline G$. Conjugating gives $\partial_{w_{\ell a}}G=\mathrm iV_{\ell a}G$.
2. **Divergence of $\mathcal Z_\ell$.** $\operatorname{div}\mathcal Z_\ell=\sum_a(\partial_{w_{\ell a}}G_a+\mathrm iV_{\ell a}G_a)=2\mathrm iw_\ell^{\mathsf T}\mathscr PG=0$, since $\mathscr PG=\overline{(\mathscr Pf)(\mathsf X-\mathrm iq)}=0$.
3. **Acceleration.** $\operatorname{div}(\mathcal Z_\ell\otimes\overline{\mathcal Z_\ell})=(\operatorname{div}\mathcal Z_\ell)\overline{\mathcal Z_\ell}+\nabla_{\mathcal Z_\ell}\overline{\mathcal Z_\ell}$. With the flat connection, $\mathcal Z_\ell\overline{G_b}=0$, and $\mathcal Z_\ell$ applied to the coefficient $(C^\mu w_m)_b$ gives $\delta_{\ell m}(C^\mu G)_b$. So
$$\nabla_{\mathcal Z_\ell}\overline{\mathcal Z_m}=-\mathrm i\,\delta_{\ell m}\sum_\mu(\overline G^{\mathsf T}C^\mu G)\,\partial_\mu ,$$
which is $\mathrm i$ times a real field and vanishes after taking real parts.
4. **Positivity.** If $\mathcal Z_\ell=A_\ell+\mathrm iB_\ell$ with real fields, its contribution to $B$ is $A_\ell\otimes A_\ell+B_\ell\otimes B_\ell$. ∎

**Remark.** The conjugation in $G$ matters. With $f(\mathsf X+\mathrm iq)$ the divergence does not vanish; [`verify_wave_current.py`](scripts/verify_wave_current.py) runs this negative control.

### 3.3 Push-forward

**Lemma 3.2 (push-forward of a conserved tensor).** Let $B$ be a smooth symmetric tensor with divergence-free columns on an open $U\subset\mathbb R^D$, and let $\Theta:U\to\Theta(U)$ be a diffeomorphism with $J=|\det D\Theta|$. Define $\tilde B$ by $\tilde B\circ\Theta=J^{-1}D\Theta\,B\,D\Theta^{\mathsf T}$. If $\sum_{I,L}B^{IL}\partial_I\partial_L\Theta^K=0$ for every $K$, then $\tilde B$ has divergence-free columns.

*Proof.* For $\varphi\in C_c^\infty(\Theta(U))$, change variables:
$$\int\tilde B^{KL}\partial_L\varphi=\int(\partial_I\Theta^K)B^{IL}\partial_L(\varphi\circ\Theta)=-\int\varphi\circ\Theta\,\big[(\operatorname{div}B)^I\partial_I\Theta^K+B^{IL}\partial_I\partial_L\Theta^K\big]=0.\qquad∎$$

### 3.4 The wave current

**Target coordinates.** Each spacetime index $\mu\in\{t,x_1,\dots,x_d,s_0\}$ gets a slot $\mathsf X^\mu+\omega_\mu\mathrm i$ with $\omega_\mu\in\mathbb C_{\mathrm k}$. Since $\mathrm{ki}=\mathrm j$, its real coordinates $(\mathsf X^\mu,\operatorname{Re}\omega_\mu,\operatorname{Im}_{\mathrm k}\omega_\mu)$ carry the units $(1,\mathrm i,\mathrm j)$. Write $\omega=(Y,P,S_1)\in\mathbb C_{\mathrm k}\times\mathbb C_{\mathrm k}^d\times\mathbb C_{\mathrm k}$. On $\{Y\ne0\}$ put:
- $y=|Y|$, $e=Y/y$, $W=P/Y$, and $\Phi=S_1+\dfrac{P\cdot P}{2Y}$ (all dots bilinear);
- $V^v=\big(1,\tfrac{v+W}2,-\tfrac{W\cdot v}2\big)^{\mathsf T}$;
- $m_{\kappa,v}=(\kappa/y)^{d/2}\exp\!\big(-\tfrac{\kappa y}2|W-v|^2\big)$;
- $\iota_e(\alpha+\mathrm i\beta)=\alpha+e\mathrm i\beta$.

**Proposition 3.3 (the wave current).** For a finite superposition $u_{\rm w}=\sum u_{\kappa,v}$ with $\kappa>0$, set
$$\chi=\sum_{\kappa,v}V^v\,\iota_e\big(\overline{u_{\kappa,v}}\big)\,m_{\kappa,v},\qquad A^0=\chi\chi^* .\tag{3.4}$$
On $\{Y\ne0\}$:
1. $A^0$ is smooth, positive semidefinite and conserved.
2. $(\partial\Phi)\chi=0$, hence $(\partial\Phi)A^0=0$, where $\partial\Phi=\ell(W)$ is the complex differential in $(Y,P,S_1)$.
3. $(\partial W_j)V^v=(v_j-W_j)/(2Y)$.
4. The same holds after the orthogonal change $(t,s_0)\to(T,b)$ applied to the two slots. There $V^v_T=(1+W\cdot v/2)/\sqrt2$ and $A^0_{TT}=|\chi_T|^2$.

*Proof.*
1. **The flux map.** Apply Lemma 3.1 to the wave system, extended by zero to a passive $\eta\in\mathbb C_{\mathrm k}$. Let
$$\Theta(\mathsf X,z,\eta)=(\mathsf X,\omega),\qquad\omega_\mu=\tfrac12z^{\mathsf T}C^\mu z+\delta_{\mu s_0}\eta,\qquad\text{i.e.}\qquad(Y,P,S_1)=\big(R^2,R\zeta,\eta-\tfrac12\zeta\cdot\zeta\big).$$
   Its complex Jacobian is $2R^{d+1}$, so $J=4|R|^{2d+2}$, and on $Y\ne0$ there is a local inverse $R=\sqrt Y$, $\zeta=P/R$, $\eta=\Phi$.
2. **Second derivatives.** Since $\Theta$ is holomorphic in $z=w_1+\mathrm kw_2$, $(\partial_{w_{1a}w_{1b}}+\partial_{w_{2a}w_{2b}})\Theta=0$. The two state blocks of $B$ equal $\operatorname{Re}(G_a\overline{G_b})$, the mixed state block vanishes, and $\Theta$ is affine in $\mathsf X$ and $\eta$. So Lemma 3.2 applies.
3. **Gram rows.** $B^{IJ}=\operatorname{Re}(Q_I\overline{Q_J})$ with $Q_I=\mathcal Z_1^I+\mathrm k\mathcal Z_2^I$; the cross terms are orthogonal to $\mathbb R$. Explicitly:
   - $Q_{\mathsf X^\mu}=\sum_aM_{\mu a}\mathrm iG_a$, where $M_{\mu a}=(C^\mu z)_a=\partial_{z_a}\omega_\mu$;
   - $Q_{w_{1a}}=G_a$ and $Q_{w_{2a}}=\mathrm kG_a$.
   
   The pushed rows $\sum_I(\partial_I\Theta^K)Q_I$ over slot $\mu$ are $\big(Q,\ \sum_a\overline{M_{\mu a}}G_a,\ \mathrm k\sum_a\overline{M_{\mu a}}G_a\big)$, with $Q=Q_{\mathsf X^\mu}$. Since $-\mathrm ic\,\mathrm i=\bar c$ for $c\in\mathbb C_{\mathrm k}$, these rows are $(Q,-\mathrm iQ,-\mathrm jQ)$. These are the rows $\bar e_\alpha\tilde\chi_\mu$ of $(\tilde\chi\tilde\chi^*)^{\mathbb R}$, with $\tilde\chi_\mu=J^{-1/2}\sum_aM_{\mu a}\mathrm iG_a$. So $\tilde B=(\tilde\chi\tilde\chi^*)^{\mathbb R}$, which is conserved by Lemma 3.2.
4. **Evaluation.**
   - $\overline{u(\mathsf X-\mathrm iq)}=\overline{u}\,e^{\kappa\ell(v)\cdot q}=\overline u\,e^{-\frac\kappa2|\zeta-vR|^2}$ by (3.3), and $|\zeta-vR|^2=y|W-v|^2$.
   - $ME(v)=(2R,\ \zeta+Rv,\ -\zeta\cdot v)^{\mathsf T}=2RV^v$.
   - $J^{-1/2}=\frac12|R|^{-d-1}$.
   
   Hence $\tilde\chi_\mu=\sum V^v_\mu\,g\,\mathrm i\,\overline u\,m_{\kappa,v}$ with $g=R/|R|\in\mathbb C_{\mathrm k}$. Multiply all rows on the right by the unit quaternion $(g\mathrm i)^{-1}=-\mathrm i\bar g$; this leaves every Gram matrix unchanged. Since $\mathrm i\bar g=g\mathrm i$ for $g\in\mathbb C_{\mathrm k}$, we get $g\mathrm i\,\overline u\,(-\mathrm i\bar g)=g\overline u\bar g=\iota_{g^2}(\overline u)$, and $g^2=Y/y=e$. Thus $\tilde\chi(g\mathrm i)^{-1}=\chi$, independently of the branch of $\sqrt Y$. Conservation is local and $\chi$ is smooth on $Y\ne0$, so item 1 holds.
5. **Items 2–4.**
   - Directly, $\ell(w)V^v=-\frac12(w-W)\cdot(w-v)$, so $\ell(W)V^v=0$. Also $\partial\Phi=(-\tfrac12W\cdot W,\ W,\ 1)=\ell(W)$, and $W=P/Y$ gives 3.
   - An orthogonal change of the slot index preserves everything.
   - Finally $T=(t-s_0)/\sqrt2$. ∎

([`verify_wave_current_target.py`](scripts/verify_wave_current_target.py) checks items 1, 2 and 4 directly in the target coordinates, independently of steps 1–4.)

### 3.5 The time component on a cap

At a fixed auxiliary point write $V^v_T=b_0+\mathrm kb_1$, with $b_0=(1+\operatorname{Re}W\cdot v/2)/\sqrt2$ and $b_1=\operatorname{Im}_{\mathrm k}W\cdot v/(2\sqrt2)$. Let $\mathsf S_\ell$ be the scalar multiplier $u_{\kappa,v}\mapsto b_\ell m_{\kappa,v}u_{\kappa,v}$. Then
$$\chi_T=\iota_e\big(\overline{\mathsf S_0u_{\rm w}}\big)+\mathrm k\,\iota_e\big(\overline{\mathsf S_1u_{\rm w}}\big),\qquad A^0_{TT}=|\mathsf S_0u_{\rm w}|^2+|\mathsf S_1u_{\rm w}|^2 ,\tag{3.5}$$
because $\operatorname{span}(1,e\mathrm i)\perp\operatorname{span}(\mathrm k,\mathrm ke\mathrm i)$.

**Cap data.** Work on the torus $\mathbb T^n=\mathbb R^n/(100\mathbb Z)^n$. For $0<h\le1$ let $\hat f$ be supported in $c_1/h\le|\xi|\le c_2/h$ and in the cap $|v|\le\gamma$ in some frame. Here $u_{\rm w}=e^{-\mathrm iT|D|}f$, and $\|f\|_2=1$. Then:
- $\kappa\asymp h^{-1}$;
- the cap contains $O(h^{-n})$ lattice frequencies, so $\sum|c_{\kappa,v}|\lesssim h^{-n/2}$ by Cauchy–Schwarz and Parseval.

The construction is periodic in the original spatial variables.

---

## 4. The fibre geometry

Fix $w_0=\frac1{10}$, an aperture $0<\gamma<w_0/8$ (chosen in §7), $K\ge1$, $a=Kh$ with $4a\le1$, and $s=|\tau|\in(a,4)$.
- **Parameters.** Slot 0 is $T+\tau\mathrm i$ with $\tau\in\mathbb C_{\mathrm k}$, and $z=(T,\tau)\in\mathbb R^3$.
- **Spatial slots.** These are $x_j+P_j\mathrm i$ and $b+Z\mathrm i$. In physical auxiliary coordinates, $Y=(\tau+Z)/\sqrt2$ and $S_1=(Z-\tau)/\sqrt2$.
- **The fibre at $z$.** $X$ runs over $\mathbb T^n$ and $(P,Z)\in\mathbb C_{\mathrm k}^d\times\mathbb C_{\mathrm k}$ over
$$\Omega_z:\quad\mathrm r:=\frac{|W|^2}{w_0^2}+100|Q|^2-1<0,\qquad Q=\Phi/a .\tag{4.1}$$
- **The inverse chart.**
$$Y=\frac{\sqrt2\,\tau+aQ}{1+W\cdot W/2},\qquad P=WY,\qquad Z=\sqrt2\,Y-\tau .\tag{4.2}$$
  It is holomorphic and depends smoothly on $\tau$. On $\bar\Omega_z$, $|1+W\cdot W/2|\ge\frac12$ and $|aQ|\le s/10$, so $y=|Y|\asymp s$.
- **Notation.** $H,G$ are the full $(n+1)$-slot operators; $H_{\rm fib},G_{\rm fib}$ are the spatial parts; $d\zeta=dX\,dP\,dZ$.

**Lemma 4.1 (fibre geometry).**
1. $H_{\rm fib}(\mathrm r)>0$. The full $H(\mathrm r)\ge0$ has a one-dimensional quaternionic kernel contained in $G(\mathrm r)^\perp\cap\ker H(|\Phi|^2)$.
2. On the closed fibres, uniformly in $a<s<4$:
$$H_{\rm fib}(\mathrm r)\ge cs^{-2}I,\qquad G_{\rm fib}(\mathrm r)G_{\rm fib}(\mathrm r)^*\le H_{\rm fib}(\mathrm r),\qquad dP\,dZ\asymp a^2s^{2d}\,dW\,dQ,\qquad\det H_{\rm fib}(\mathrm r)\asymp a^{-2}s^{-2d}.$$
3. $H_{\rm fib}(|\Phi|^2)$ has rank one, with nonzero eigenvalue between two absolute constants. Every first derivative of $W$ is $O(s^{-1})$; in particular $|G(|W|^2)|\lesssim s^{-1}$.

*Proof.* Let $F=(W/w_0,10Q)$, holomorphic in $(\tau,P,Z)$, with complex Jacobian $J$ (spatial block $J_{\rm fib}$). Since $\mathrm r=|F|^2-1$, Lemma 2.1(4) gives $H(\mathrm r)=4J^*J$ and $G(\mathrm r)=2J^*F\mathrm i$.
1. **Positivity and kernel.** $J_{\rm fib}$ is invertible (the chart (4.2)), so $H_{\rm fib}(\mathrm r)>0$. $J$ has rank $n$ on $n+1$ variables; its kernel is the direction in which $\tau$ varies with $(W,Q)$ fixed. The quaternionic span $N_z$ of that kernel is $\ker H(\mathrm r)$, and it is $\perp G(\mathrm r)\in J^*(\cdot)$. Moreover $\partial\Phi$ is proportional to the last row of $J$, so $H(|\Phi|^2)=4(\partial\Phi)^*\partial\Phi$ vanishes on $N_z$.
2. **Bounds.**
   - From (4.2), derivatives of $(P,Z)$ in $W$ are $O(s)$ and in $Q$ are $O(a)\le O(s)$. So $\|J_{\rm fib}^{-1}\|\lesssim s$ and $H_{\rm fib}(\mathrm r)\ge cs^{-2}I$.
   - $G_{\rm fib}G_{\rm fib}^*=4J_{\rm fib}^*FF^*J_{\rm fib}\le4|F|^2J_{\rm fib}^*J_{\rm fib}\le H_{\rm fib}(\mathrm r)$, since $|F|\le1$ on $\bar\Omega_z$.
   - Wedging $dP=Y\,dW+W\,dY$ with $dZ=\sqrt2\,dY$ gives the complex Jacobian $\sqrt2aY^d/(1+W\cdot W/2)$ of $(W,Q)\mapsto(P,Z)$. Its squared modulus gives the volume, and $\det(4J_{\rm fib}^*J_{\rm fib})=4^n|\det_{\mathbb C}J_{\rm fib}|^2$.
3. **The background direction.**
   - $\partial_P\Phi=W$ and $\partial_Z\Phi=(1-W\cdot W/2)/\sqrt2$ is bounded away from $0$, so $H_{\rm fib}(|\Phi|^2)=4(\partial_{\rm fib}\Phi)^*\partial_{\rm fib}\Phi$ has the stated eigenvalue.
   - $W=P/Y$ with $|Y|\asymp s$. ∎

On these fibres $\int\delta(\mathrm r)\,dP\,dZ\le Ca^2s^{2d}$, by coarea on the fixed ellipsoid in $(W,Q)$.

---

## 5. The Dirichlet problem on a periodic fibre

Let $\Lambda\subset\mathbb R^m$ be a lattice, $V\subset\mathbb R^{2m}$ a bounded smooth domain, and $\Omega=(\mathbb R^m/\Lambda)\times V$. The periodic variables occupy the real parts of $m$ three-real-coordinate slots, and $V$ occupies the $\mathrm i,\mathrm j$ parts. $H$ is the 3-slot Hessian.

**Lemma 5.1 (comparison).** Let $\Omega'\subset\Omega$ be open, and let $\psi_1,\psi_2\in C^2(\Omega')\cap C(\bar\Omega')$ have $Q_i=H(\psi_i)>0$, with continuous extensions to $\bar\Omega'$. If $\det Q_1\ge\det Q_2$ in $\Omega'$ and $\psi_1\le\psi_2$ on $\partial\Omega'$, then $\psi_1\le\psi_2$ in $\Omega'$.

*Proof.* By Lemma 2.2(1),
$$\log\det Q_1-\log\det Q_2=\Big(\int_0^1((1-\theta)Q_2+\theta Q_1)^{-1}d\theta\Big):H(\psi_1-\psi_2)=:L(\psi_1-\psi_2)\ge0.$$
The coefficient matrix of $L$ is positive definite on $\bar\Omega'$. As a real operator, $Lw=\mathcal M^{\mathbb R}:D^2w$ with $\xi^{\mathsf T}\mathcal M^{\mathbb R}\xi=\hat\xi^*\mathcal M\hat\xi>0$, where $\hat\xi_a=\sum_\alpha e_\alpha\xi_{a\alpha}$. So $L$ is elliptic without zeroth-order term, and the weak maximum principle gives $\psi_1-\psi_2\le\max_{\partial\Omega'}(\psi_1-\psi_2)\le0$. ∎

**Theorem 5.2 (periodic Dirichlet problem; paper Thm A.1).** Assume:
- $\mathrm r\in C^\infty(\bar V)$ is negative in $V$, with $d\mathrm r\ne0$ on $\partial V$ and $H(\mathrm r)>0$ on $\bar\Omega$;
- $B_0,g\in C^\infty(\bar\Omega)$ are real with $H(B_0)\ge0$ and $g>\det H(B_0)$.

Then there is a unique $u\in C^\infty(\bar\Omega)$ with
$$\det H(B_0+u)=g,\qquad H(B_0+u)>0,\qquad u|_{\partial\Omega}=0,$$
and there are constants $b\ge b'>0$ with $b\,\mathrm r\le u\le b'\mathrm r$. For smooth finite-dimensional families of data and domains, $b,b'$ can be chosen locally uniformly, and $u$ depends smoothly on the parameters up to the boundary, after a smooth identification of the domains with a fixed one.

*Proof.* **Status:** steps 1–3 are proved here; steps 4–7 use cited elliptic theory, whose hypotheses are checked. The whole theorem is also in the Lean proof.
1. **Uniqueness and barriers.**
   - Uniqueness is Lemma 5.1 applied both ways.
   - **Barriers.** Choose $0<b'<1<b$ with $\det H(B_0+b'\mathrm r)<g_t<\det H(B_0+b\mathrm r)$ on $\bar\Omega$, for every $g_t$ in the continuity path below. The left inequality holds for small $b'$ by the uniform gap $g_t-\det H(B_0)>0$. The right holds for large $b$ since $H(B_0+b\mathrm r)\ge b\,cI$. Lemma 5.1 then gives $b\mathrm r\le u\le b'\mathrm r$.
   - Writing $u=\alpha\mathrm r$ near $\partial\Omega$ gives $b'\le\alpha\le b$ and $\nabla u=\alpha\nabla\mathrm r$ on $\partial\Omega$.
2. **The continuity path.** Let $g_t=(1-t)\det H(B_0+\mathrm r)+tg$, $0\le t\le1$. Then $u=\mathrm r$ solves at $t=0$. Also $g_t>\det H(B_0)$, since $\det(H(B_0)+H(\mathrm r))>\det H(B_0)$ by strict monotonicity.
3. **Openness.** At a solution $\psi=B_0+u$, with $Q=H(\psi)$ positive definite on $\bar\Omega$, the linearization of $\log\det H(B_0+u)-\log g_t$ is $Lw=Q^{-1}:H(w)$. It is uniformly elliptic with smooth coefficients and no zeroth-order term, with symbol $\hat\xi^*Q^{-1}\hat\xi$.
   - The Dirichlet problem $Lw=f$, $w|_{\partial\Omega}=0$, is uniquely solvable in $C^{k+2,\alpha}$ for $f\in C^{k,\alpha}$. **Cited:** Schauder theory, GT Thms. 6.14 and 6.19, applied on a finite cover of the flat quotient; the method of continuity from the flat Laplacian, solved by Fourier series in the periodic variables.
   - The implicit function theorem gives openness.
4. **$C^1$ bound.** Differentiating along a constant real unit vector $e$ gives $L\psi_e=(\log g_t)_e$ and, by concavity, $L\psi_{ee}\ge(\log g_t)_{ee}$. Since $L\mathrm r\ge c\operatorname{tr}Q^{-1}\ge c\,m\,g_t^{-1/m}\ge c_0>0$, the maximum principle applied to $\pm\psi_e+C\mathrm r$ bounds $|\nabla\psi|$ by its boundary values, which are bounded by step 1. Applied to $\psi_{ee}+C\mathrm r$, it reduces the upper bound for $D^2\psi$ to the boundary.
5. **Boundary $C^2$.** Lift to $\tilde\Omega=\Omega\times(\mathbb R/\mathbb Z)^m$, adding a constant fourth coordinate in each slot. On functions independent of the new coordinates the 4-slot Hessian equals the 3-slot one, and 4-slot $H$ is covariant under quaternionic rigid motions (Lemma 2.1(5)). At a boundary point:
   - **Tangential–tangential:** $D^2u(v,w)=\alpha D^2\mathrm r(v,w)$.
   - **Tangential–normal** (Caffarelli–Nirenberg–Spruck). Let $v=P_{\tan}D^2\mathrm r\,e/|\nabla\mathrm r|$ and $A=v\nu^*-\nu v^*-\nu(\nu^*v)\nu^*$, a quaternionic skew-Hermitian matrix. The field $X=(e+Ax)\cdot\nabla$ generates rigid motions, and $X\mathrm r$ vanishes to second order at the point.
     - Covariance gives $L(X\psi)=X\log g_t$, so $|L(Xu)|\le C\operatorname{tr}Q^{-1}$.
     - In a small chart ball the barrier $w=C_1|x|^2-C_2\mathrm r$ satisfies $|Xu|\le w$ on the boundary and $L(w\pm Xu)\le0$. Here $H(|x|^2)=8I$ and $L\mathrm r\ge c\operatorname{tr}Q^{-1}$.
     - So $|Xu|\le w$, which bounds $|\partial_\nu Xu|=|D^2u(\nu,e)|$ at the point.
   - **Normal–normal.** In a quaternionic orthonormal basis starting with $\nu$, the compression of $Q$ to $\nu^\perp$ is $\ge b'H(\mathrm r)|_{\nu^\perp}\ge cI$, by Lemma 2.4(3). The Schur formula gives $\psi_{\nu\nu}=g_t/\det(\cdot)+z^*(\cdot)^{-1}z-c_1$, with $z,c_1$ already bounded.
   - **Global.** With the interior reduction and $\Delta\psi=\operatorname{tr}Q>0$, this bounds $|D^2\psi|$ globally.
6. **$C^{2,\alpha}$ bound.** Now all eigenvalues of $Q$ lie in a fixed $[c,C]\subset(0,\infty)$, and $M\mapsto\log\det\mathcal H(M)$ is concave and uniformly elliptic on the convex set of such Hessians. **Cited:** the interior Evans–Krylov estimate and Krylov's boundary estimate (GT Thms. 17.14, 17.26′, hypotheses (i)′–(ii)′) give a uniform $C^{2,\alpha}$ bound. Schauder bootstrapping gives smoothness up to the boundary, so the path is closed and reaches $t=1$.
7. **Descent and parameters.**
   - Uniqueness on $\tilde\Omega$ makes the solution independent of the added coordinates, so it descends to $\Omega$.
   - **Parameters.** Pull back to a fixed domain. The linearization is conjugate to the Dirichlet isomorphism of step 3, so the implicit function theorem in $C^{k+2,\alpha}_0$ gives smooth dependence. The strict barrier inequalities persist locally. ∎

---

## 6. The determinant energy

Keep the cap data of §3.5 and the fibres of §4. Put $\lambda=h^{-30n}$.

**6.1 The data of the Dirichlet problem.**
- **Smoothing kernel.** $J_\ell$ is the periodization of $c_n\ell^{-n}(1+|X/\ell|^2)^{-M/2}$, with $M=2n+4$ and integral one.
- **Cutoff.** $0\le g\le1$ is smooth on $\mathbb C_{\mathrm k}$, equal to $1$ on $|Q|\le\frac1{100}$ and to $0$ on $|Q|\ge\frac1{20}$.
- **Current, density and background.**
$$A_{\rm w}=g(Q)\,(J_h*_XA^0),\qquad A=A_{\rm w}+h^{50n}C(H(\mathrm r)),\qquad\mu=A_{TT},\qquad B_0=\lambda|\Phi|^2,\qquad B=H(B_0).\tag{6.1}$$

**Lemma 6.1 (properties of the data).**
1. $A\ge0$ is conserved on the total domain $\Omega=\{(z,\zeta):a<|\tau|<4,\ \zeta\in\Omega_z\}$, and $\mu\ge h^{50n}\det H_{\rm fib}(\mathrm r)>0$.
2. $A:B=0$ and $\det B=0$. $B_{\rm fib}$ has rank one, with nonzero eigenvalue $\beta\asymp\lambda$.
3. On the regularizer, $\nu^*C(H(\mathrm r))\nu=0$ on the wall.
4. For $h\le h_0(n,K)$, $\mu^{1/d}\le\lambda$.

*Proof.*
1. **Conservation.**
   - $X$-convolution preserves conservation.
   - $g(Q)$ is a function of the holomorphic $\Phi$, so $G(g(Q))$ lies in the image of $(\partial\Phi)^*$ (Lemma 2.1(4)), which $J_h*A^0$ annihilates (Prop. 3.3(2)). The product rule then preserves conservation.
   - $C(H(\mathrm r))$ is conserved and $\ge0$ by Lemmas 2.2–2.3.
   
   Its image is $\ker H(\mathrm r)$, and its $TT$ entry is $\det H_{\rm fib}(\mathrm r)$, by Lemma 2.2(3) with $v=e_T$.
2. **The background.**
   - $H(|\Phi|^2)=4(\partial\Phi)^*\partial\Phi$ has rank one, so $\det B=0$ (as $n+1\ge2$).
   - $A_{\rm w}H(|\Phi|^2)=0$ by tangency, and $H(|\Phi|^2)C(H(\mathrm r))=0$ because $\ker H(\mathrm r)\subset\ker H(|\Phi|^2)$ (Lemma 4.1). Hence $A:B=0$.
   - Lemma 4.1(3) gives $\beta$.
3. The image of $C(H(\mathrm r))$, which is $\ker H(\mathrm r)$, is $\perp G(\mathrm r)\parallel\nu$.
4. Each $m_{\kappa,v}\le(\kappa/y)^{d/2}\lesssim h^{-d}$ (as $y\gtrsim s\ge h$), $V^v$ is bounded, and $\sum|c|\lesssim h^{-n/2}$. So $\|A^0\|\lesssim h^{-n-2d}$. The regularizer adds $O(h^{50n}a^{-2}s^{-2d})$. Hence $\mu\le h^{-3n}$ and $\mu^{1/d}\le h^{-3n}\le\lambda$. ∎

**The potential.** By Lemma 4.1 and Lemma 6.1, Theorem 5.2 applies on each fibre: $m=n$, the lattice of the torus in the rotated frame, the ellipsoid image $V_z$, the defining function $\mathrm r$, $B_0$, and $g=\lambda\mu>0=\det B_{\rm fib}$. It gives
$$D_u:=H_{\rm fib}(B_0+u)>0,\qquad\det D_u=\lambda\mu,\qquad u|_{\partial\Omega_z}=0,\qquad u\le0,\tag{6.2}$$
with $u$ smooth on the closure of the total domain. Write $\rho=u/\mathrm r\ge0$, which is smooth up to the wall.

**The energy.** Put $D_v=H_{\rm fib}(B_0+v)$ and
$$X_0(z)=\int_{\Omega_z}(-u)\mu\,d\zeta,\qquad E_z(u)=\frac1\lambda\int_{\Omega_z}u\int_0^1\det D_{\theta u}\,d\theta\,d\zeta,\qquad\mathcal F(z)=X_0(z)+E_z(u).\tag{6.3}$$

**Proposition 6.2 (energy balance).** $\mathcal F$ is smooth on $a<|\tau|<4$, and
$$\Delta_z\mathcal F+\int_{\Omega_z}A:\Sigma\,d\zeta=\mathcal B(z):=\int\rho\,G(\mathrm r)^*A_{\rm w}G(\mathrm r)\,\delta(\mathrm r)\,d\zeta .$$
Here $\Sigma\ge0$ is defined below, and $\delta(\mathrm r)\,d\zeta=dS_{\rm fib}/|\nabla_{\rm fib}\mathrm r|$ on one fibre. In particular $\Delta_z\mathcal F\le\mathcal B$.

*Proof.*
1. **The Schur split.** In the slot order $(0,{\rm fib})$, write
$$H(B_0+u)=\begin{pmatrix}\ell&c^*\\c&D_u\end{pmatrix}=\sigma_0e_Te_T^*+\Sigma,\qquad\sigma_0=\ell-c^*D_u^{-1}c,\quad\Sigma\ge0 .$$
   By Lemma 2.2(4) and (6.2), $\det H(B_0+u)=\lambda\mu\sigma_0$. With $A:B=0$ this gives
$$A:H(u)=A:H(B_0+u)=\mu\sigma_0+A:\Sigma=\lambda^{-1}\det H(B_0+u)+A:\Sigma .\tag{6.4}$$
2. **Green's formula in the parameters.** Let $\eta\in C_c^\infty$ depend only on $z$, so $H(\eta)$ has the single entry $H_{TT}(\eta)=\Delta_z\eta$. For a conserved $M$, Lemma 2.4(2) on the total domain (no end terms) and $G(u)=\nu\partial_{\mathbf n}u$ give
$$\int_\Omega uM_{TT}\Delta_z\eta=\int_\Omega\eta\,M:H(u)-\int_{\partial_{\rm lat}\Omega}\eta\,(\nu^*M\nu)\,\partial_{\mathbf n}u\,dS .\tag{6.5}$$
   On the wall, $\partial_{\mathbf n}u=\rho|\nabla\mathrm r|$ and $dS=|\nabla\mathrm r|\delta(\mathrm r)\,dz\,d\zeta$.
3. **$M=A$.** By Lemma 6.1(3) only $A_{\rm w}$ contributes on the wall. With (6.4):
$$\Delta_zX_0=\mathcal B-\int A:\Sigma-\lambda^{-1}\int\det H(B_0+u).$$
4. **$M=C(H(B_0+\theta u))$.** It is conserved by Lemma 2.3, $M_{TT}=\det D_{\theta u}$, and $M:H(u)=\partial_\theta\det H(B_0+\theta u)$. On the wall, Lemma 2.4(3) gives $H(B_0+\theta u)|_{\nu^\perp}=(B+\theta\rho H(\mathrm r))|_{\nu^\perp}$. This contains the common null line of $B$ and $H(\mathrm r)$ (Lemma 4.1), which lies in $\nu^\perp$. So $\nu^*M\nu=\det(\cdot|_{\nu^\perp})=0$ by Lemma 2.2(3), whatever the sign of the full Hessian. Integrating in $\theta\in[0,1]$, dividing by $\lambda$ and using $\det B=0$:
$$\Delta_zE_z(u)=\lambda^{-1}\int\det H(B_0+u).$$
5. **Conclusion.** Add the two identities. They hold in distributions, and the chart (4.2) together with the smooth dependence in Theorem 5.2 makes them pointwise. ∎

**Lemma 6.3 (energy comparison).** On every fibre,
$$\mathcal F=\frac{n-1}{n+1}X_0+\frac1{\lambda n(n+1)}\int_{\Omega_z}G_{\rm fib}(u)^*C(D_u)G_{\rm fib}(u)\,d\zeta,\qquad\text{so}\qquad\frac{n-1}{n+1}X_0\le\mathcal F\le\frac n{n+1}X_0 .$$

*Proof.*
1. **The segment.** $D_{\theta u}=(1-\theta)B_{\rm fib}+\theta D_u$ and $B_{\rm fib}$ has rank one. By Lemma 2.2(3),
$$\det D_{\theta u}=\theta^n\det D_u+(1-\theta)\theta^{n-1}C(D_u):B_{\rm fib}.$$
   Integrating in $\theta$ gives $\mathcal F=\frac n{n+1}X_0-\frac{I}{\lambda n(n+1)}$, where $I=\int(-u)\,C(D_u):B_{\rm fib}\ge0$.
2. **Euler's identity.** $C(D_u):D_u=n\det D_u=n\lambda\mu$. Integrating by parts on the fibre ($C(D_u)$ is conserved and $u=0$ on the wall) gives
$$n\lambda X_0-I=\int(-u)\,C(D_u):H_{\rm fib}(u)=\int G_{\rm fib}(u)^*C(D_u)G_{\rm fib}(u)\ge0 .\qquad∎$$

---

## 7. Estimates on a fixed fibre

**Proposition 7.1 (paper Prop. 5.1).** Fix an annulus $c_1\le h|\xi|\le c_2$. There are $\gamma>0$, $K\ge1$ and $c,C>0$, depending only on $n$ and the annulus, with the following property. Let $a=Kh$ with $0<h<h_0$, let $f$ be as in §3.5 with $\|f\|_2=1$, and put $g_T=e^{-\mathrm iT|D|}f$. Then, for all $T\in\mathbb R$,
$$\mathcal F(0,\tau)\le Ca^2\|f\|_{p_c}^{p_c}\ \ (|\tau|=4a),\qquad\mathcal F(T,\tau)\ge ca^2\|g_T\|_{p_c}^{p_c}\ \ (|\tau|=2a),\qquad\Delta_z\mathcal F\le V(s)\mathcal F\ \ (2a\le s\le3),\tag{7.1}$$
where $V(s)=Cs^{-2}e^{-cs/h}$.

Throughout, $\sigma=s/h\ge2K$ and $2a\le s\le3$. Two regions of the fibre matter.
- **The wall region** $\mathrm r=0$, $|Q|\le\frac1{20}$, where $g(Q)\ne0$. There $|W|\ge\frac{\sqrt3}2w_0$.
- **The centre region** $|W|\le\gamma$, $|Q|\le\frac1{400}$.

**Lemma 7.2 (multiplier kernels).** Let $\mathfrak a(\eta)$ be smooth, supported in a fixed compact subset of $\mathbb R^n\setminus0$, with $|\partial^\beta\mathfrak a|\le\mathsf B$ for $|\beta|\le2n+4$. On $\mathbb T^n$ the multiplier $\mathfrak a(hD)$ satisfies
$$|\mathfrak a(hD)g|^2\le C\mathsf B^2\,J_h*|g|^2 .$$
Also:
- $J_h*J_h\le CJ_h$;
- $J_h\le(Ls/h)^nJ_{Ls}$ for $Ls\ge h$;
- $\inf J_{Ls}\ge c_Lh^{M-n}$.

*Proof.*
1. **The Euclidean kernel.** It is $K(X)=(2\pi h)^{-n}\int\mathfrak a(\eta)e^{\mathrm iX\cdot\eta/h}d\eta$. Applying $(1-\Delta_\eta)^{n+2}$ gives $|K|\le C\mathsf Bh^{-n}(1+|X/h|^2)^{-n-2}$.
2. **The torus kernel.** By Poisson summation the torus kernel is the periodization of $K$, hence $\le C\mathsf BJ_h$. Cauchy–Schwarz against the probability density $J_h$ gives the main bound.
3. **The kernel facts.**
   - On $\mathbb R^n$, $(j_h*j_h)(X)\le2j_h(X/2)\le2^{M+1}j_h(X)$, because $|Y|$ or $|X-Y|$ is $\ge|X|/2$; periodization preserves this.
   - Comparing $(1+|X/h|^2)^{-M/2}\le(1+|X/(Ls)|^2)^{-M/2}$ gives the second.
   - The nearest lattice term on a torus of fixed diameter gives the third, using $s\ge h$. ∎

**Step 1: the spatial average.** For a constant $L\ge1$ fixed in Step 3, let $P_*=(J_{Ls}*|g_T|^2)^{1/d}$. Then
$$|\nabla^j_XP_*|\le C(Ls)^{-j}P_*\ (j=1,2),\qquad\int P_*^{d+1}dX\le\|g_T\|_{p_c}^{p_c}.\tag{7.2}$$
The first holds because $|\nabla^jj_\ell|\le C\ell^{-j}j_\ell$ termwise. The second is Young's inequality with exponent $q$ and a kernel of integral one. Also $P_*>0$.

**Step 2: bounds on the current.** By (3.4), each real scalar part of $\chi$ is a multiplier applied to $g_T$, with symbol $b(\xi)(\kappa/y)^{d/2}e^{-\Psi}$, where $\Psi=\frac{\kappa y}2|W-v|^2$ and $b$ is a component of $V^v$. Insert a fixed annular cutoff in $\eta=h\xi$ that equals $1$ on the data and is supported in $|v|<w_0/3$.
- **Derivatives.** On its support $h\kappa$ and $v$ are smooth in $\eta$, each $\eta$-derivative of $\Psi$ is $O(\sigma)$, and $(\kappa/y)^{d/2}\asymp(hs)^{-d/2}$. So the derivatives through order $2n+4$ are $\le C(hs)^{-d/2}\sigma^me^{-\Psi}$.
- **Gain on the wall.** In the wall region, $|W-v|\ge(\frac{\sqrt3}2-\frac13)w_0$, so $e^{-\Psi}\le e^{-c_0\sigma}$, with $c_0$ independent of $\gamma$.

$\|A^0\|=|\chi|^2$ is the sum of the squared scalar parts. Lemma 7.2 and $(hs)^{-d}(s/h)^n=s^{-2d}\sigma^{n+d}$ then give
$$\|A_{\rm w}\|\le Cs^{-2d}\sigma^mP_*^d\ \ (\mathrm r\le0),\qquad\|A_{\rm w}\|\le Cs^{-2d}\sigma^me^{-c_0\sigma}P_*^d\ \ (\text{wall}),\qquad\mu\le Cs^{-2d}\sigma^mP_*^d\ \ (\mathrm r\le0).\tag{7.3}$$
The last inequality also uses $h^{50n}\det H_{\rm fib}(\mathrm r)\le Ch^{50n}a^{-2}s^{-2d}\le s^{-2d}c_Lh^{M-n}\le s^{-2d}P_*^d$ for small $h$.

**Step 3: an upper comparison and the energy bound.** Let $R=H_{\rm fib}(\mathrm r)$. By Lemma 2.1(3), (7.2), $R\ge cs^{-2}I$ and $G_{\rm fib}G_{\rm fib}^*\le R$, for every quaternionic vector $e$:
$$e^*H_{\rm fib}(\mathrm rP_*)e\ \ge\ (1-C/L)\,P_*\,e^*Re .\tag{7.4}$$
- The cross terms are $\le2C(e^*Re)^{1/2}\cdot(Ls)^{-1}P_*|e|$, and $|e|\le Cs(e^*Re)^{1/2}$.
- The term $\mathrm r\,H(P_*)$ costs $CL^{-2}P_*e^*Re$, since $|\mathrm r|\le1$.

Fix $L$ with $1-C/L\ge\frac12$, and let $v_*=C_1\sigma^{m_1}\mathrm rP_*$. Then $H_{\rm fib}(v_*)\ge\kappa I$ with $\kappa=c\,C_1\sigma^{m_1}P_*s^{-2}$. Since $\det(B_{\rm fib}+K)\ge(\beta+\kappa)\kappa^d\ge c\lambda\kappa^d$ for $K\ge\kappa I$, (7.3) gives
$$\det\big(B_{\rm fib}+H_{\rm fib}(v_*)\big)\ \ge\ \lambda\mu\qquad\text{for }C_1,m_1\text{ large}.$$
Both potentials vanish on the wall, so Lemma 5.1 gives $v_*\le u\le0$, that is,
$$0\le u/\mathrm r\le C_1\sigma^{m_1}P_*\quad(\mathrm r\le0).\tag{7.5}$$
By Lemma 6.3, $\mathcal F\le X_0$. With (7.3), (7.5), $dP\,dZ\asymp a^2s^{2d}dW\,dQ$ and (7.2):
$$\mathcal F(T,\tau)\le Ca^2\sigma^m\|g_T\|_{p_c}^{p_c}.\tag{7.6}$$

**Step 4: the flux.** $A_{\rm w}G(|\Phi|^2)=0$ by tangency, so $G(\mathrm r)^*A_{\rm w}G(\mathrm r)=w_0^{-4}G(|W|^2)^*A_{\rm w}G(|W|^2)\le Cs^{-2}\|A_{\rm w}\|$. With the wall bound in (7.3), (7.5) and $\int\delta(\mathrm r)\,dP\,dZ\le Ca^2s^{2d}$, Proposition 6.2 gives
$$\Delta_z\mathcal F\le\mathcal B\le C(a/s)^2\sigma^me^{-c_0\sigma}\|g_T\|_{p_c}^{p_c}.\tag{7.7}$$

**Step 5: recovering the wave in the interior.** Consider a centre $(X,P,Z)$ in the centre region, and the ball of radius $h$ about it in the $3n$ real spatial coordinates.
1. **The ball stays inside.** On the ball $W=P/Y$ moves by $\le Ch/s\le C/(2K)$ and $Q$ by $\le C/K$. So for $K\ge C/\gamma$ the ball lies in $|W|\le2\gamma$, $|Q|\le\frac1{100}$, where $g(Q)=1$ and $\mathrm r<0$.
2. **Inverting the window.** At every auxiliary point of the ball, invert $\mathsf S_0$ of (3.5) on the frequencies of $g_T$. Its symbol is $b_0(\kappa/y)^{d/2}e^{-\Psi}$, with $b_0\ge c>0$. Multiply the reciprocal by a cutoff equal to $1$ on the data and supported in $|v|\le2\gamma$; there $|W-v|\le4\gamma$, so $\Psi\le C'\gamma^2\sigma$. Its derivatives through order $2n+4$ are
$$\le B_{\rm inv}=C_\gamma(hs)^{d/2}\sigma^me^{C'\gamma^2\sigma},$$
   where $C'$ does not depend on $\gamma$.
3. **Density bound at a point.** By Lemma 7.2 and (3.5), $|g_T(X)|^2\le CB_{\rm inv}^2\,J_h*A^0_{TT}(\cdot,P',Z')(X)$ for each auxiliary point $(P',Z')$ of the ball. Since $J_h(X'-Y)\asymp J_h(X-Y)$ for $|X'-X|\le h$, every point of the ball has
$$\mu\ \ge\ m_X:=c_\gamma(hs)^{-d}\sigma^{-m}e^{-C_2\gamma^2\sigma}|g_T(X)|^2 .\tag{7.8}$$
4. **Depth at the centre.** Let $l$ be the Euclidean distance to the centre and $v_b=\epsilon m_X^{1/d}(l^2-h^2)$, with fibre Hessian $6\epsilon m_X^{1/d}I$. Since $m_X^{1/d}\le\mu^{1/d}\le\lambda$ (Lemma 6.1(4)),
$$\det(B_{\rm fib}+6\epsilon m_X^{1/d}I)\le(C\lambda+6\epsilon\lambda)(6\epsilon)^dm_X\le\lambda m_X\le\lambda\mu$$
   for small $\epsilon=\epsilon(n)$. Also $u\le0=v_b$ on the sphere. Lemma 5.1 on the ball gives $u\le v_b$, hence $(-u)\mu\ge\epsilon h^2m_X^{1+1/d}$ at the centre. Here $1+\frac1d=q$ and $m_X^q\propto|g_T(X)|^{p_c}$.
5. **Integrate over the centres.** The imaginary volume of the centre region is $\asymp a^2s^{2d}\gamma^{2d}$, and $h^2s^{2d}(hs)^{-dq}=\sigma^{n-2}\ge1$. With Lemma 6.3,
$$\mathcal F(T,\tau)\ \ge\ \tfrac{n-1}{n+1}X_0\ \ge\ c_\gamma a^2\sigma^{-m}e^{-C_3\gamma^2\sigma}\|g_T\|_{p_c}^{p_c}.\tag{7.9}$$

**Step 6: choice of constants and conclusion.**
1. **Constants.** Choose $\gamma<w_0/8$ with $C_3\gamma^2<c_0/2$; the constants $c_0,C_3$ do not depend on $\gamma$. Then choose $K$ large enough for the ball inclusions, then $h_0$ for Lemma 6.1(4) and Step 2.
2. **The circles.** At $|\tau|=4a$ ($\sigma=4K$), (7.6) with $T=0$ gives the first bound in (7.1). At $|\tau|=2a$ ($\sigma=2K$), (7.9) gives the second.
3. **The potential.** Dividing (7.7) by (7.9):
$$\frac{\Delta_z\mathcal F}{\mathcal F}\le C_\gamma s^{-2}\sigma^{2m}e^{-(c_0-C_3\gamma^2)\sigma}\le Cs^{-2}e^{-(c_0/4)s/h}.\qquad∎$$

---

## 8. Propagation across the thin cylinder

**Lemma 8.1 (two-circle propagation; paper Lemma 6.1).** Assume:
- $0<b<\frac14$ and $D=\log(3/b)$;
- $V:[b,3]\to[0,\infty)$ is continuous with $\int_b^3s\log(s/b)V(s)\,ds\le M$ and $\sup_{[1/2,3/2]}V\le M$;
- $E\ge0$ is smooth on $\mathbb R_T\times\{b\le|\tau|\le3\}\subset\mathbb R^3$, with $\Delta_{T,\tau}E\le V(|\tau|)E$ in the interior.

Then the circular mean $H(T,s)$ satisfies, for $b\le r\le1$,
$$H(0,r)\ \ge\ c(M)\,\frac{\log(r/b)}{D^2}\int_1^2H(T,b)\,dT .$$

*Proof.*
1. **Averaging.** Averaging over rotations preserves the inequality, since $V$ is radial. So $\Delta H\le V(|\tau|)H$.
2. **Radial solutions.** For $S\in\{1,3\}$ and fixed $\kappa>0$, let $u$ solve $u''+s^{-1}u'=(V+\kappa^2)u$ on $(b,S)$ with endpoint values $\{0,1\}$ in either order. We claim $u\ge c\cdot(\text{the function affine in }\log s\text{ with the same values})$.
   - With $x=\log(s/b)$ and $A=\log(S/b)$ the equation is $u_{xx}=P(x)u$, where $P=s^2(V+\kappa^2)$.
   - Let $\psi''=P\psi$ with $\psi(0)=0$, $\psi'(0)=1$. Then $\psi\ge x$ and $\psi(x)/x=1+\int_0^x\frac{z(x-z)}xP(z)\frac{\psi(z)}z\,dz$.
   - Since $z(x-z)/x\le\min(z,A-z)$, Gronwall's inequality gives $\psi(A)/A\le\exp\int_0^A\min(z,A-z)P(z)\,dz$.
   - In the variable $s$, this integral is at most $\int_b^Ss\log(s/b)V\,ds+\kappa^2S^2\int_0^\infty te^{-2t}dt\le M+\kappa^2S^2/4$.
   
   So $u=\psi/\psi(A)\ge c\,x/A$. The integral is symmetric under $x\mapsto A-x$, which gives the other order.
3. **Outward passage.** Let $\mathsf m(s)=\int_0^3\sin(\pi T/3)H(T,s)\,dT$. Integrating by parts twice,
$$\int_0^3\sin(\tfrac{\pi T}3)\partial_T^2H\,dT=-(\pi/3)^2\mathsf m+\tfrac\pi3(H(0,s)+H(3,s))\ge-(\pi/3)^2\mathsf m .$$
   So $\mathsf m''+s^{-1}\mathsf m'\le(V+\pi^2/9)\mathsf m$, and $\mathsf m(3)\ge0$. The maximum principle against the radial solution with $S=3$, $u(b)=1$, $u(3)=0$ gives $\mathsf m(s)\ge c\,\mathsf m(b)\log(3/s)/D$. Hence
$$\int_0^3\!\!\int_{3/4\le|\tau|\le5/4}H\,d\tau\,dT\ge c\,\mathsf m(b)/D .$$
4. **Weak Harnack on a fixed region.** On $(-2,4)\times\{\frac12<|\tau|<\frac32\}$ we have $\Delta H\le MH$. Then $H(T,|\tau|)\cos(\sqrt{M+1}\,v)$ is nonnegative and superharmonic in $\mathbb R^4$ on a slab $|v|<\pi/(2\sqrt{M+1})$.
   - **The mean value inequality** on the 4-ball of radius $2\rho$, restricted to $B_\rho(x)\times(-\rho,\rho)$, gives $H(x)\ge c\rho^{-3}\int_{B_\rho(x)}H$ for one small $\rho$ and all $x$ near the connected band $[-1,3]\times\{\frac34\le|\tau|\le\frac54\}$.
   - **Chaining.** Cover the band by finitely many balls of radius $\rho/4$. If two of them, $B$ and $B'$, overlap, then $B'\subset B_\rho(x)$ for all $x\in B$, so $\inf_BH\ge c\rho^{-3}\int_{B'}H$. A finite chain and step 3 give $H(T,1)\ge c(M)\mathsf m(b)/D$ for $|T|\le1$.
5. **Return passage.** Let $u$ be the radial solution with $S=1$, $\kappa=\pi/2$, $u(b)=0$, $u(1)=1$. The function $w=c(M)\mathsf m(b)D^{-1}\cos(\pi T/2)u(|\tau|)$ solves $(-\Delta+V)w=0$ on $|T|<1$, $b<|\tau|<1$. It vanishes on $|T|=1$ and on $|\tau|=b$, and it is $\le H$ on $|\tau|=1$ by step 4.
   - The maximum principle for $-\Delta+V$ with $V\ge0$ gives $H\ge w$.
   - At $T=0$, step 2 gives $H(0,r)\ge c(M)\mathsf m(b)D^{-1}\log(r/b)/\log(1/b)$.
   
   Finally $\mathsf m(b)\ge\sin(\pi/3)\int_1^2H(T,b)\,dT$ and $\log(1/b)\le D$. ∎

**Remark (sharpness).** Remark 6.2 of the paper gives a supersolution with $V=0$ attaining $\log(r/b)/D^2$ up to the constant $\frac94\log2$ at $r=2b$. This is verified exactly in [`verify_propagation_and_exponents.py`](scripts/verify_propagation_and_exponents.py); [07](07_TOOLKIT.md) K4–K5 identifies the constant and explains why only a two-dimensional auxiliary time gives logarithms.

**Proposition 8.2 (one cap on the torus; paper Prop. 6.3).** Let $f$ on $\mathbb T^n$ have Fourier support in $c_1/h\le|\xi|\le c_2/h$ and in one cap of the aperture $\gamma$ of Proposition 7.1, with $0<h\le1$. Then
$$\Big(\int_1^2\|U(T)f\|_{L^{p_c}(\mathbb T^n)}^{p_c}dT\Big)^{1/p_c}\le C\,[\log(2/h)]^{2/p_c}\|f\|_{L^{p_c}(\mathbb T^n)} ,$$
uniformly in $h$ and in the cap.

*Proof.* By homogeneity let $\|f\|_2=1$.
1. **Small $h$.** For $h<h_0$, put $a=Kh$ and $b=2a<\frac14$, and apply Lemma 8.1 to $E=\mathcal F\ge0$. This is smooth on $a<|\tau|<4$, which contains $b\le|\tau|\le3$.
   - With $V(s)=Cs^{-2}e^{-cs/h}$ and $s=bt$, $\int_b^3s\log(s/b)V\le C\int_1^\infty\frac{\log t}te^{-2cKt}dt$, uniformly in $h$. Also $V$ is bounded on $[\frac12,\frac32]$.
   - By (7.1), $H(0,2b)\le Ca^2\|f\|_{p_c}^{p_c}$ and $H(T,b)\ge ca^2\|U(T)f\|_{p_c}^{p_c}$.
   - Take $r=2b$, cancel $a^2$, and use $D=\log(3/(2Kh))\le C\log(2/h)$.
2. **Large $h$.** For $h\ge h_0$ the frequencies lie in a fixed finite set $\mathcal S$, and $\|U(T)f\|_p\le\sum_{\xi\in\mathcal S}|\hat f(\xi)|\,\|e^{\mathrm iX\cdot\xi}\|_p\le C\|f\|_p$. ∎

---

## 9. Euclidean space and the Sobolev estimate

**Proof of Theorem 1.2.** Write $h=N^{-1}$ and $P_h=\varphi(hD)$.
1. **Caps.** Choose a finite smooth partition $\varphi=\sum_\nu\varphi_\nu$ subordinate to caps of aperture $\gamma$. Each $\varphi_\nu(hD)$ has $L^1$-bounded kernel $h^{-n}\check\varphi_\nu(\cdot/h)$, and periodization does not increase $L^1$ norms. So Proposition 8.2 gives, on $\mathbb T^n$,
$$\|U(T)P_hf\|_{L^{p_c}(I\times\mathbb T^n)}\lesssim[\log(2/h)]^{2/p_c}\|f\|_{L^{p_c}(\mathbb T^n)}.$$
2. **Kernel decay.** The kernel $K_h(T,X)=(2\pi h)^{-n}\int e^{\mathrm i(X\cdot\eta-T|\eta|)/h}\varphi(\eta)\,d\eta$ satisfies $|K_h|\le C_Mh^M(1+|X|)^{-M}$ for $T\in I$ and $|X|\ge3$. Integrate by parts along the phase gradient, which has length $\ge|X|/3$.
3. **Localization.** Tile $\mathbb R^n$ by unit cubes $Q$ and let $10Q$ be the concentric cube of side 10.
   - **Far part.** On $Q$, the contribution of $f\mathbf 1_{(10Q)^c}$ is $\le C_Mh^M((1+|\cdot|)^{-M}*|f|)$. Its $L^{p_c}$ norm over all cubes is $O(h^M\|f\|_{p_c})$ by Young.
   - **Near part.** Periodize $f_Q=f\mathbf 1_{10Q}$ on the torus of side 100. On $Q$, the Euclidean and periodic outputs differ by at most $C_Mh^M\|f_Q\|_1$, since the nonzero translates have arguments of length $\ge94.5|k|$. Apply step 1.
   - **Sum.** Take $p_c$-th powers and sum over $Q$; the cubes $10Q$ overlap boundedly.
   
   Finally $\log(2/h)\asymp\log(2+N)$. ∎

**Interpolation.**
- **$L^2$.** By Plancherel, $\|U(T)P_h\|_{L^2\to L^2_{T,X}}\le\|\varphi\|_\infty$.
- **$L^\infty$.** Its kernel has $L^1$ norm $O(h^{-d/2})$. Partition the annulus into $O(h^{-d/2})$ caps of aperture $c\sqrt h$. After the linear phase is removed, the remaining phase $T|z|^2/(\sqrt{\rho^2+h|z|^2}+\rho)$ has bounded derivatives in the cap coordinates $\eta=(\sqrt h\,z,\rho)$, so each cap kernel has bounded $L^1$ norm.
- **Riesz–Thorin** with Theorem 1.2 gives
$$\|U(T)P_hf\|_{L^p(I\times\mathbb R^n)}\lesssim_ph^{-s_+(p)}[\log(2/h)]^{\beta(p)}\|f\|_p,\qquad\beta(p)=\min\{d(1-\tfrac2p),\tfrac2p\}.\tag{9.1}$$
  The log exponent is $\frac2{p_c}\theta$ with $\theta=2n(\frac12-\frac1p)$, which equals $d(1-\frac2p)$, below $p_c$. Above $p_c$ it is $2/p$, and the power of $h^{-1}$ is $\frac d2(1-\frac{p_c}p)=s_+(p)$.

**Proof of Theorem 1.1.** Let $p\ge2$ and $g=J^{s_+(p)+\varepsilon}f$.
- **Dyadic frequencies $2^j$, $j\ge1$.** By (9.1) the contribution is $\le C2^{-j\varepsilon}(1+j)^{\beta(p)}\|g\|_p$. Here $2^{j(s_++\varepsilon)}J^{-s_+-\varepsilon}\tilde P_j$ has a uniformly integrable kernel. Sum over $j$.
- **Low frequencies.** $e^{-\mathrm iT|\xi|}\chi_0(\xi)\langle\xi\rangle^{-s_+-\varepsilon}$ satisfies Mihlin bounds uniformly in $T\in I$.
- **$p=2$.** Use conservation of the $L^2$ norm. ∎

---

## 10. Status summary and cross-references

| Step | Paper | Here | Status | Closest Lean area (`OAI/Analysis/`) | Check |
|---|---|---|---|---|---|
| Hessian calculus, Moore determinant, cofactor divergence | Lemmas 3.1–3.3 | §2 | proved (spectral theorem cited) | `Calculus/Quaternion` (`CofactorDivergence`, `Polarization`, `Gram`) | `verify_quaternion_calculus.py` |
| Lifted current, push-forward, wave current, time component | Lemma 2.1, (2.7)–(2.15), Prop. 2.2 | §3 | proved; **Prop. 3.3 also verified directly** | `Calculus/Current` (`Wave`, `WaveDivergence`, `Holomorphic`) | `verify_wave_current.py`, `verify_wave_current_target.py` |
| Fibre geometry | Lemma 4.1 | §4 | proved | `Calculus/Current/Physical`, `…/Fiber` | — |
| Comparison principle | App. A | Lemma 5.1 | proved | `PDE/Dirichlet` | — |
| Dirichlet existence and regularity | Thm. A.1 | Thm. 5.2 | standard theory cited; hypotheses checked; formally verified | `PDE/Dirichlet` (398), `PDE/Schauder` (274) | — |
| Energy balance and comparison | Prop. 4.2, Lemma 4.3 | §6 | proved | `Calculus/Current/Cancellation`, `…/Regularized*` | `verify_quaternion_calculus.py` (Lemma 4.3) |
| Fibre estimates | Prop. 5.1 | §7 | proved, including kernel bounds | `Calculus/Current/Fiber/*`, `PDE/Concentration` | exponent bookkeeping in `verify_propagation_and_exponents.py` |
| Propagation (sharp $\log^2$) | Lemma 6.1, Prop. 6.3 | §8 | proved | **not used by Lean**, which uses `Concentration/CorrectedEnergy` and `Cylinder/Estimate` with power $2n+5$ | `verify_propagation_and_exponents.py` |
| Euclidean transfer, interpolation, Theorem 1.1 | §7 | §9 | proved (Riesz–Thorin, Mihlin, Littlewood–Paley cited) | `Fourier/LocalSmoothing` | — |

**Net result.** Every step of Theorems 1.1 and 1.2 is now either proved in this file or a cited textbook theorem whose hypotheses are checked here. The one substantial cited block, the existence and regularity theory of Theorem 5.2, is also certified by OpenAI's Lean proof, which I built ([06](06_LEAN_AUDIT.md)). The sharp power $2/p_c$ is prose-only, but its proof (Lemma 8.1, Proposition 8.2) is elementary and complete above.

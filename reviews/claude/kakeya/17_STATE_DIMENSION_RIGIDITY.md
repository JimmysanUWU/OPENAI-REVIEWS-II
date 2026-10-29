# 17. Which systems carry the engine: state-dimension rigidity

*Claude, 29 October 2026. Roadmap target `T.KAK.7` (sharpened and corrected); node `KAK.RIGID` in [`ROADMAP.md`](../../../ROADMAP.md).*

[11](11_FRAMEWORK_POSITIVITY_ENGINE.md) Theorem 11.2 builds engine data for a symmetric pencil whose extremal eigenvalue is simple and whose state dimension $N$ equals the space dimension $n$. 11 §7.3 and target `T.KAK.7` then said that for $N>n$ the pushed-forward current is "still positive, conserved and tangent", so that only the Gaussian bounds (A3)–(A4) would need re-deriving. **The tangency claim is wrong.** This file proves what is true.

> **Theorem 17.2 (tangency is a rank condition).** A lifted current admits a passive coordinate if and only if, for every complex state $z$, the $n+1$ vectors
> $$z,\ C_1z,\ \dots,\ C_nz$$
> are linearly dependent (after compressing to the span $V$ of the extremal eigenvectors). Equivalently, the complexified flux map $z\mapsto(\frac12z^{\mathsf T}C^\mu z)_\mu$ has its image in a hypersurface.

> **Theorem 17.4 (state-dimension rigidity).** Let the extremal sheet have constant multiplicity $m$ and be transversally curved, and let $N'=\dim V$. Then $N'\ge n+m-1$, and a passive coordinate exists **if and only if $N'=n+m-1$**, that is, iff the extremal eigenvectors fill an open subset of $V$.
>
> For a simple sheet ($m=1$) this says $N'=n$. So App. E's hypothesis "state dimension $=n$" is necessary.

> **Proposition 17.6 (the three division algebras).** The bound $N'=n+m-1$ holds for:
> - real symmetric pencils with $N'=n$ (App. E; the wave is the arrow pencil);
> - complex Hermitian $m'\times m'$ pencils in $n=2m'-1$ variables ($m=2$);
> - quaternionic Hermitian $m''\times m''$ pencils in $n=4m''-3$ variables ($m=4$).
>
> The isotropic members, $C(\xi)^2=|\xi|^2I$, are the **Hopf pencils** over $\mathbb R,\mathbb C,\mathbb H,\mathbb O$. Their flux map is the Hopf map, and their polarization is the Hopf fibration. By Hurwitz–Radon these four are the only Clifford pencils meeting the bound (Corollary 17.7).

**What this changes.**
- **A correction.** 11 §7.3 and `T.KAK.7` are wrong about tangency for $N>n$; both are corrected (§6). A simple sheet cannot carry engine data from a single lift unless its eigenvectors span exactly $n$ dimensions.
- **A test anyone can run.** Whether a constant-coefficient system can carry a tangent lifted current is decided by linear algebra: the rank of $[z\,|\,C_1z\,|\cdots|\,C_nz]$ at one random complex $z$.
- **A sharper target.** `T.KAK.7` becomes: prove (A3)–(A4) for the fibre-integrated currents of complex and quaternionic Hermitian pencils that meet the bound. These pass the tangency test (§4).
- **One picture for the known instances.**

  | System | Polarization |
  |---|---|
  | wave (arrow pencil) | stereographic projection (11 Corollary 11.4) |
  | 2D Dirac | complex squaring: the $\mathbb R$-Hopf map |
  | 3D massless Dirac (Weyl), real form | the Hopf fibration $S^3\to S^2$ |
  | the quaternionic and octonionic Hopf pencils | $S^7\to S^4$ and $S^{15}\to S^8$ |
  | Maxwell's equations | fails the test: $N'=6>n+m-1=4$ |

**Checks.** [`scripts/verify_state_dimension.py`](scripts/verify_state_dimension.py) ([`results/state_dimension.json`](results/state_dimension.json), ~12 s):
- **Part R.**
  - The arrow pencils ($n=2,3,4$) land on the light cone (exact) and have rank $n$.
  - Random pencils with $N=n+1$ have $\det[z\,|\,C_1z\,|\cdots]\not\equiv0$ (exact) and rank $n+1$.
  - The arrow pencil coupled to one extra state ($N'=4$) has rank $4$.
  - The uncoupled sum has $N'=n$, and its compression has rank $n$.
  - The fibre identity of Lemma 17.3(4) (exact).
- **Part H**, the Hopf pencils over $\mathbb R,\mathbb C,\mathbb H,\mathbb O$:
  - anticommutation and the Hopf identity (exact);
  - multiplicity $m=k$;
  - rank $n$;
  - $\Phi\circ\Theta=\eta$ and tangency at complex points (40 digits);
  - polarization equal to the Hopf map.
  - Controls that fail: four of the five quaternionic matrices; a perturbed $\mathbb C$-pencil; arrow $\oplus$ arrow with one lift.
- **Part Q.**
  - Random complex Hermitian pencils ($m'=2,3,4$, $n=2m'-1$) and quaternionic ones ($m''=2,3$, $n=4m''-3$) meet $N'=n+m-1$, are transversally curved at a sample direction, and pass with rank $n$.
  - With one direction fewer they fail.
  - Maxwell fails.
- **Part W.** Hurwitz–Radon: $\rho(2(n-1))\ge n$ only for $n=2,3,5,9$ (checked for $n\le40$; the bound $\rho(N)\le2\log_2N+2$ covers the rest).

---

## 0. Setting

**Pencils.**
- $C^0=I_N$ and $C^1,\dots,C^n$ are real symmetric $N\times N$ matrices.
- $C(\xi)=\sum_j\xi_jC^j$, and for $c=(c_0,c')\in\mathbb C^{1+n}$, $\mathcal C(c)=c_0I+C(c')$.
- The system is $\mathscr P=\partial_T+\sum_jC^j\partial_{X_j}$, as in 11 §6.

**The extremal sheet.** On an open cone $\Gamma\subset\mathbb R^n\setminus0$:
- $\varphi=\lambda_{\max}\circ C$ has constant multiplicity $m$;
- $D^2\varphi(\xi)$ is positive definite on $\xi^\perp$ (transversally curved).

$E(\xi)$ is the $\varphi$-eigenspace and $P(\xi)$ the projection onto it; both are real-analytic (Kato).

$V=\operatorname{span}\bigcup_{\xi\in\Gamma}E(\xi)$, $N'=\dim V$, and $C_V^\mu=\Pi_VC^\mu|_V$, the compression. By analyticity, the span over any open subcone is $V$.

**Cap data.** Finite sums $f=\sum_\xi e_\xi u_\xi$ with $\xi$ in a subcone, $e_\xi\in E(\xi)$, and $u_\xi=c_\xi e^{\mathrm i(X\cdot\xi-T\varphi(\xi))}$. They solve $\mathscr Pf=0$ and take values in $V$.

**Definition 17.1 (single lift; passive coordinate).**
1. **The single lift.** Fix an auxiliary state space $U$ with $V\subseteq U\subseteq\mathbb R^N$.
   - Apply 08 Lemma 3.1 to $\mathscr P$, with auxiliary states $w_1,w_2\in U$ and $z=w_1+\mathrm kw_2\in U_{\mathbb C}$.
   - Push forward along $\Theta(X,z,\eta)=(X,\ Q(z)+c\eta)$, where $Q(z)=(\frac12z^{\mathsf T}C^\mu z)_{\mu=0}^n$, $c\in\mathbb C^{1+n}$ is a fixed passive direction, and $\eta\in\mathbb C_{\mathrm k}$.
   - Where $\Theta$ is a submersion, integrate over its fibres.

   The pointwise pushed tensor is $(\tilde\chi\tilde\chi^*)^{\mathbb R}$ with $\tilde\chi_\mu=\sum_a(C^\mu z)_a\,\mathrm iG_a$ and $G=\overline{f(\mathsf X-\mathrm iq)}$ (08 Prop. 3.3, step 3, which uses no property of the $C^\mu$). So $A^0=\int\tilde\chi\tilde\chi^*$ over each fibre, with the coarea weight. When $N=n$ there is no integral, and this is 08 Prop. 3.3 and 11 Lemma 6.4.
2. **A passive coordinate.** This is a holomorphic $\Phi$ with $\partial\Phi\ne0$, on an open set $\mathcal U$ of the target's auxiliary space $\mathbb C_{\mathrm k}^{1+n}$, such that $(\partial\Phi)A^0=0$ for every cap data. We require $\mathcal U\ni Q(z_0)$ for some real $z_0\in E(\xi_0)\setminus0$, $\xi_0\in\Gamma$: the centres of the engine's fibres are such points.

   Axiom (A2) together with $|\partial_P\Phi|\ge1/C_0$ in (G2) gives a passive coordinate.

Tangency is a pointwise property on each fibre. Nothing below depends on whether the fibre integrals converge.

---

## 1. Tangency at a point

**Lemma 17.1.** $\Phi$ is a passive coordinate if and only if, at every lifted point $(z,\eta)$ with $\Theta(z,\eta)\in\mathcal U$,
$$\Pi_V\,\nabla_z(\Phi\circ\Theta)=\Pi_V\,\mathcal C\big(\partial\Phi(\Theta(z,\eta))\big)\,z=0 .$$
In particular, $\Phi\circ\Theta$ must not depend on the $V$-components of $z$.

*Proof.*
- **The gradient.** $\partial_{z_a}(\Phi\circ\Theta)=\sum_\mu\partial_\mu\Phi\,(C^\mu z)_a=(\mathcal C(\partial\Phi)z)_a$, by holomorphy and the symmetry of $C^\mu$. So
$$(\partial\Phi)\tilde\chi=\sum_a(\mathcal C(\partial\Phi)z)_a\,\mathrm iG_a .$$
- **Splitting.** Write $\mathcal C(\partial\Phi)z=\alpha+\mathrm k\beta$ with $\alpha,\beta$ real. The terms $\alpha_a\mathrm iG_a\in\mathbb C_{\mathrm i}$ and $\mathrm k\beta_a\mathrm iG_a\in\mathrm k\mathbb C_{\mathrm i}$ are orthogonal, so $(\partial\Phi)\tilde\chi=0$ iff $\alpha\cdot G=\beta\cdot G=0$.
- **The range of $G$.** At a fixed lifted point, $G=\sum_\xi e_\xi\overline{u_\xi(\mathsf X-\mathrm iq)}$. The amplitudes are arbitrary and every $e_\xi\in E(\xi)$ is allowed, so $G$ ranges over $V\otimes\mathbb C_{\mathrm i}$. Hence the condition is $\Pi_V\alpha=\Pi_V\beta=0$.
- **Integrating.** $(\partial\Phi)A^0=0$ iff $\int|(\partial\Phi)\tilde\chi|^2=0$ over each fibre, iff the integrand vanishes there. ∎

---

## 2. The rank criterion

Put $M_V(z)=[\,z\,|\,C_V^1z\,|\cdots|\,C_V^nz\,]$, an $N'\times(n+1)$ matrix, so that $M_V(z)c=\mathcal C_V(c)z$. The complex Jacobian of $Q_V=Q|_{V_{\mathbb C}}$ at $z$ is $M_V(z)^{\mathsf T}$.

**Theorem 17.2.** Take $U=V$. A passive coordinate exists near $Q(z_0)$ if and only if $\operatorname{rank}M_V(z)\le n$ for every $z\in V_{\mathbb C}$. Then, for a passive direction $c$ transverse to the image, $\Phi$ can be chosen with $\Phi\circ\Theta=\eta$. For $U=\mathbb R^N$ the rank condition is still necessary.

*Proof.*
- **(⇒)** By Lemma 17.1, at every $(z,\eta)$ near $(z_0,0)$ with $z\in V_{\mathbb C}$: $M_V(z)\,\partial\Phi(\Theta(z,\eta))=\Pi_V\mathcal C(\partial\Phi)z=0$ with $\partial\Phi\ne0$. So $\operatorname{rank}M_V\le n$ on an open subset of $V_{\mathbb C}$. Its $(n+1)$-minors are polynomials, so they vanish identically. This argument only uses lifted points with $z\in V_{\mathbb C}$, so it also covers $U=\mathbb R^N$.
- **(⇐)**
  - $\operatorname{rank}DQ_V\le n$ everywhere, and $\operatorname{rank}DQ_V(z_0)\ge n$ by Lemma 17.3(3) below. By lower semicontinuity the rank is $n$ near $z_0$.
  - By the constant rank theorem, $Q_V$ maps a neighbourhood of $z_0$ onto a complex hypersurface $H$ near $Q(z_0)$.
  - Take a defining function $\Psi$ of $H$ and a direction $c\notin T_{Q(z_0)}H$. Define $\Phi$ near $Q(z_0)$ by $\Psi(\omega-c\,\Phi(\omega))=0$ and $\Phi(Q(z_0))=0$; this is possible by the implicit function theorem, since $d\Psi\cdot c\ne0$.
  - Then $\Phi(Q_V(z)+c\eta)=\eta$, so $\nabla_z(\Phi\circ\Theta)=0$, and Lemma 17.1 applies. ∎

**Remark (auxiliary states).** Lift with $U=V$. A larger auxiliary space only adds conditions: $\Phi\circ\Theta$ must also be compatible with the $V^\perp$-components of $z$. For instance, for the arrow pencil with an uncoupled extra scalar state, the lift with $U=V=\mathbb R^n$ works, and that is the natural choice.

---

## 3. Rigidity

**Lemma 17.3 (four facts about an extremal sheet of constant multiplicity).** Let $\xi_0\in\Gamma$ and $z_0\in E(\xi_0)\setminus0$ real.
1. **(Hellmann–Feynman on eigenspaces.)** $P(\xi)C^jP(\xi)=\partial_j\varphi(\xi)P(\xi)$. Hence $Q(e)=\frac12|e|^2(1,\nabla\varphi(\xi))$ for $e\in E(\xi)$.
2. **(The eigen-cone is a submanifold.)** Near $(\hat\xi_0,z_0)$ the map $\iota:(\hat\xi,e)\mapsto e$, from the eigenbundle $\{(\hat\xi,e):e\in E(\hat\xi)\}$ to $V$, is an injective immersion. Its image $\mathcal E$ has dimension $n-1+m$. In particular $N'\ge n+m-1$.
3. $\operatorname{rank}_{\mathbb R}DQ_V(z_0)\ge n$.
4. **(The fibre identity.)** For real $z\in V$,
$$z^{\mathsf T}S_V(\xi_0)z=2\big(\varphi(\xi_0)Q_0(z)-\xi_0\cdot Q'(z)\big),\qquad S_V=\varphi(\xi_0)I-C_V(\xi_0),$$
   where $S_V\ge0$ on $V$ and $\ker S_V=E(\xi_0)$.

*Proof.*
1. Differentiate $C(\xi)P=\varphi P$ in $\xi_j$ and multiply by $P$ on both sides. The terms $PC\,\partial_jP\,P$ and $\varphi P\,\partial_jP\,P$ cancel because $PC=\varphi P$.
2. **Injectivity.** If $\iota(\hat\xi_1,e)=\iota(\hat\xi_2,e)$ with $e\ne0$, then item 1 gives $\nabla\varphi(\hat\xi_1)=\nabla\varphi(\hat\xi_2)$. The Gauss map $\hat\xi\mapsto\nabla\varphi(\hat\xi)$ has injective differential $D^2\varphi|_{\hat\xi^\perp}$, so it is locally injective.

   **Immersion.** Take a curve $(\hat\xi(t),e(t))$ with $e'(0)=0$. Differentiating $Q'(e(t))=\frac12|e|^2\nabla\varphi(\hat\xi(t))$ gives $0=\frac12|e|^2D^2\varphi\,\hat\xi'$, so $\hat\xi'=0$. The remaining tangent vector is vertical, with image $e'(0)=0$.
3. By item 1, $Q_V(re(\hat\xi))=\frac12r^2(1,\nabla\varphi(\hat\xi))$. Its differential has rank $1+(n-1)$: the $r$-direction, plus $D^2\varphi$ on $\hat\xi^\perp$.
4. The identity is the definition of $S_V$. For $z\in V$, $z^{\mathsf T}C(\xi_0)z\le\varphi|z|^2$, since $\varphi$ is the top eigenvalue. Equality forces $z\in E(\xi_0)$, which lies in $V$. ∎

**Theorem 17.4.** Under the standing assumptions, $N'\ge n+m-1$. If a passive coordinate exists (for any auxiliary space $U$), then $N'=n+m-1$. Conversely, if $N'=n+m-1$, i.e. if the extremal eigenvectors fill an open subset of $V$, then a passive coordinate exists for $U=V$.

*Proof.*
- **The bound** is Lemma 17.3(2).
- **If.** Suppose $N'=n+m-1$.
  - Then $\mathcal E$ is open in $V$ near $z_0$ (invariance of domain).
  - Every $z\in\mathcal E$ lies in some $E(\hat\xi)$, so $\mathcal C_V(-\varphi(\hat\xi),\hat\xi)z=0$ and $\operatorname{rank}M_V(z)\le n$ on a real open set.
  - The minors are polynomials, so the bound holds on all of $V_{\mathbb C}$, and Theorem 17.2 (⇐) applies.
- **Only if.** Suppose a passive coordinate exists.
  - By Theorem 17.2 (⇒), $\operatorname{rank}M_V\le n$ everywhere. With Lemma 17.3(3) the real rank of $DQ_V$ is exactly $n$ near $z_0$.
  - So the level set $F=\{z\in V:Q_V(z)=Q_V(z_0)\}$ is a submanifold of dimension $N'-n$ near $z_0$.
  - By Lemma 17.3(4), $z^{\mathsf T}S_Vz$ is constant on $F$, equal to $z_0^{\mathsf T}S_Vz_0=0$. So $F\subset E(\xi_0)$.
  - By Lemma 17.3(1), $F$ is exactly the sphere $\{e\in E(\xi_0):|e|=|z_0|\}$ near $z_0$, of dimension $m-1$.
  - Hence $N'-n=m-1$. ∎

**Corollary 17.5 (simple sheets).** For a simple, extremal, transversally curved eigenvalue, a single lift carries a passive coordinate iff $N'=n$. In that case the compressed pencil $(C_V^j)$ is an App. E pencil of size $n$, and Theorem 11.2 gives full engine data.

So the hypothesis "state dimension $=n$" of 11 Theorem 11.2 is necessary for simple sheets. The uncoupled case $N>n=N'$ reduces to it by compression.

---

## 4. The division-algebra families

**Proposition 17.6.**
1. **(Real.)** $N'=n$, $m=1$: App. E. The arrow pencil is the wave (11 Corollary 11.4), and its flux map lands on the light cone: $\sum_j(z^{\mathsf T}A_jz)^2=(z^{\mathsf T}z)^2$.
2. **(Complex Hermitian.)** Let $H_1,\dots,H_n$ be complex Hermitian $m'\times m'$ matrices, realified to real symmetric $2m'\times2m'$ matrices. A simple (complex) top eigenvalue has real multiplicity $m=2$. Suppose it is transversally curved and its eigenvectors span $\mathbb C^{m'}$. Then $N'=2m'$, and a passive coordinate exists iff $n=2m'-1$.
3. **(Quaternionic Hermitian.)** In the same way, $m=4$ and $N'=4m''$, and a passive coordinate exists iff $n=4m''-3$.
4. **(Hopf pencils.)** Let $\mathbb A\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$ with $k=\dim\mathbb A$, $L_a$ left multiplication, and
$$C(r,a)=\begin{pmatrix}r&L_a^{\mathsf T}\\L_a&-r\end{pmatrix}\quad\text{on }\mathbb A^2=\mathbb R^{2k},\qquad\xi=(r,a)\in\mathbb R^{1+k}.$$
   - $C(\xi)^2=|\xi|^2I$, so $n=k+1$, $N=N'=2k$ and $m=k$: the bound holds with equality.
   - The flux map is the Hopf map $(x,y)\mapsto\frac12\big(|x|^2+|y|^2,\ |x|^2-|y|^2,\ 2(\langle y,L_ex\rangle)_e\big)$, whose image is the light cone.
   - The passive coordinate is $\Phi(\tau,P)=\omega\cdot P-\sqrt{\tau^2-P\cdot P+(\omega\cdot P)^2}$, with $\Phi\circ\Theta=\eta$.
   - A real state $z$ is a top eigenvector of $C(h(z))$, where $h(z)=Q'(z)/Q_0(z)\in S^k$: **the polarization is the Hopf fibration.**

   The case $k=1$ is the arrow pencil of $\mathbb R^2$ (the 2D Dirac operator, whose flux map is complex squaring). The case $k=2$ is the massless 3D Dirac (Weyl) operator in real form. The cases $k=2,4,8$ are the isotropic members of items 2 and 3 (and of the octonionic $2\times2$ analogue).

*Proof.*
- Items 2–3 are Theorem 17.4 with the stated $m$ and $N'$.
- **Item 4.**
  - Block multiplication and $L_a^{\mathsf T}L_a=L_aL_a^{\mathsf T}=|a|^2I$ (a normed division algebra) give $C(\xi)^2=|\xi|^2I$. The two eigenvalues $\pm|\xi|$ have equal multiplicity since $\operatorname{tr}C_j=0$.
  - $z^{\mathsf T}C(\xi)z=r(|x|^2-|y|^2)+2\langle y,L_ax\rangle$. The vectors $L_ex$, for $e$ in an orthonormal basis, are orthogonal of length $|x|$. Hence $\sum_e\langle y,L_ex\rangle^2=|x|^2|y|^2$ and $|Q'|=Q_0$.
  - So every real $z$ satisfies $z^{\mathsf T}C(h)z=|z|^2$ with $|h|=1$, i.e. it is a top eigenvector of $C(h)$, and $\operatorname{rank}M\le n$.
  - $\Phi\circ\Theta=\eta$ follows from the cone identity: at $P=Q'+\omega\eta$ the square root equals $\omega\cdot Q'$. ∎

**Corollary 17.7 (Hurwitz–Radon).** Among Clifford pencils ($C^jC^l+C^lC^j=2\delta_{jl}I$), the bound $N'=n+m-1$ holds exactly for $n\in\{2,3,5,9\}$, $N=2(n-1)$: the four Hopf pencils.

*Proof.*
- **The dimension.** The eigenvalues are $\pm|\xi|$, each with multiplicity $N/2$. Also $V=\mathbb R^N$, since bottom eigenvectors at $\xi$ are top eigenvectors at $-\xi$. So the bound reads $N=2(n-1)$.
- **The obstruction.** $C^1C^2,\dots,C^1C^n$ are $n-1$ anticommuting orthogonal complex structures on $\mathbb R^N$. So $n\le\rho(N)$, the Hurwitz–Radon number ($\rho(2^{4a+b}\cdot\text{odd})=8a+2^b$).
- **The count.** $\rho(2(n-1))\ge n$ holds only for $n=2,3,5,9$. This is checked for $n\le40$ (Part W); beyond that, $\rho(N)\le2\log_2N+2<N/2+1$. ∎

**Numerical evidence (Part Q).**

| Pencil | $n$ | $N'$ | $m$ | $n+m-1$ | rank $M$ | passes |
|---|---|---|---|---|---|---|
| complex Hermitian $2\times2$, $3\times3$, $4\times4$ | 3, 5, 7 | 4, 6, 8 | 2 | 4, 6, 8 | $n$ | yes |
| the same with one direction fewer | 2, 4, 6 | 4, 6, 8 | 2 | 3, 5, 7 | $n+1$ | no |
| quaternionic Hermitian $2\times2$, $3\times3$ | 5, 9 | 8, 12 | 4 | 8, 12 | $n$ | yes |
| the same with one direction fewer | 4, 8 | 8, 12 | 4 | 7, 11 | $n+1$ | no |
| Maxwell (isotropic medium) | 3 | 6 | 2 | 4 | 4 | no |

The passing pencils are random, and transversally curved at a sample direction. Their top sheets are those of Hermitian determinantal hyperbolic polynomials, which are in general not quadrics.

---

## 5. The known instances, sorted

| System | $N$ | $n$ | $m$ | $N'$ | single lift tangent? |
|---|---|---|---|---|---|
| wave = arrow pencil; App. E pencils | $n$ | $n$ | 1 | $n$ | yes (11 Theorem 11.2) |
| wave in null coordinates; Schrödinger ([16](16_NULL_FRAME_AND_SCHRODINGER.md)) | $n$ | $n$ | 1 | $n$ | yes: $\Theta$ is a local biholomorphism (08 Prop. 3.3; here $C^t=2E_{00}$, outside §0's normalization) |
| arrow $\oplus$ uncoupled scalar | $n+1$ | $n$ | 1 | $n$ | yes, after compressing to $V$ |
| arrow coupled to one more state; generic $N=n+1$ | $n+1$ | $n$ | 1 | $n+1$ | **no** |
| arrow $\oplus$ arrow, one lift | $2n$ | $n$ | 2 | $2n$ | **no** (two lifts, added, do work) |
| 3D Weyl/Dirac, real form; quaternionic and octonionic Hopf | $2k$ | $k+1$ | $k$ | $2k$ | yes; (A3)–(A4) open |
| complex Hermitian, $n=2m'-1$; quaternionic Hermitian, $n=4m''-3$ | $2m'$, $4m''$ | | 2, 4 | $=N$ | yes; (A3)–(A4) open |
| Maxwell; four of the five quaternionic Clifford matrices | 6, 8 | 3, 4 | 2, 4 | 6, 8 | **no** |

**Two notes.**
1. **Sums of lifts.** The theorem concerns one lift. Currents of decoupled subsystems can be added. Each needs its own lift, and each summand must satisfy the theorem.
2. **No new estimates from the isotropic cases.** The Hopf pencils' extremal sheet is the round cone, so even with (A3)–(A4) they would reprove the wave estimate. The non-isotropic Hermitian pencils are the interesting ones.

---

## 6. Corrections and the refined target

- **11 §7.3**, row "State dimension $=n$". This row said that for $N>n$ the push-forward is "still positive, conserved and tangent". Corrected: it is still positive and conserved. For a simple sheet it is tangent to some passive coordinate only if $N'=n$ (Corollary 17.5). With multiplicity $m$, the condition is $N'=n+m-1$ (Theorem 17.4).
- **`T.KAK.7`** is restated: *Hermitian pencils.*
  - Complex and quaternionic Hermitian pencils that meet the bound pass the tangency test.
  - Prove (A3)–(A4) for their fibre-integrated currents. Whether the fibre integrals converge with uniform symbol bounds is part of the problem: the fibres of $\Theta$ are complexified spheres $S^{m-1}$, non-compact in their imaginary directions, and the Gaussian windows narrow along them.
  - Then decide whether their sheets reach beyond App. E. This connects to determinantal representations of hyperbolic polynomials, real symmetric versus Hermitian.

---

## 7. Status and cross-references

| Result | Status | Uses |
|---|---|---|
| Lemma 17.1, Theorem 17.2 | proved here; checked (Part R) | 08 Lemma 3.1, Prop. 3.3 step 3 |
| Lemma 17.3, Theorem 17.4, Corollary 17.5 | proved here; checked on instances (Parts R, Q) | Kato (analytic eigenprojection); constant rank theorem |
| Proposition 17.6 | proved here; **exact checks** (Part H), numerical (Part Q) | normed division algebras |
| Corollary 17.7 | proved here (Part W) | Hurwitz–Radon theorem |

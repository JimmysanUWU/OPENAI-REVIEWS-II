# 21. Why $N'=n+m-1$: the flux map is a projection of a rank-one cone

*Claude, 1 November 2026. Distillation of [17](17_STATE_DIMENSION_RIGIDITY.md) (tangency) and [20](20_GAUGE_ORBIT_OBSTRUCTION.md) (divergence) by classical invariant theory. Roadmap node `KAK.QUOTIENT`; target `T.KAK.7`.*

17 proved that a single lift of a symmetric pencil is tangent to a passive coordinate iff $N'=n+m-1$. Here $N'$ is the dimension spanned by the extremal eigenvectors and $m$ their multiplicity. The equality cases were real symmetric, complex Hermitian and quaternionic Hermitian pencils, with
$$n=N'\ (\mathbb R),\qquad n=2m'-1\ (\mathbb C),\qquad n=4m''-3\ (\mathbb H).$$
17 obtained these counts by a rank computation. This note explains them, as one theorem:

> **Theorem 21.1.** Let $\mathbb F\in\{\mathbb R,\mathbb C,\mathbb H\}$, $K=\{\pm1\},U(1),Sp(1)$ its unit group, $m=\dim_{\mathbb R}\mathbb F$. Let the pencil consist of $\mathbb F$-Hermitian matrices, realified on $V=\mathbb R^{N'}$. Then the flux map factors as
> $$Q\;=\;L\circ\pi\colon\ V_{\mathbb C}\xrightarrow{\ \pi\ }X_{\mathbb F}\xrightarrow{\ L\ }\mathbb C^{1+n},$$
> where $\pi$ is the quotient by the complexified gauge group $K_{\mathbb C}$ and $L$ is **linear**. The quotient $X_{\mathbb F}$ is a classical rank-one cone:
>
> | $\mathbb F$ | $K_{\mathbb C}$ | $\pi(z)$ | $X_{\mathbb F}$ | $\dim X_{\mathbb F}$ |
> |---|---|---|---|---|
> | $\mathbb R$ | $\{\pm1\}$ | $zz^{\mathsf T}$ | rank-one symmetric matrices (Veronese cone) | $N'$ |
> | $\mathbb C$ | $\mathbb C^\times$ | $uv^{\mathsf T}$ | rank-one $m'\times m'$ matrices (Segre cone) | $2m'-1=N'-1$ |
> | $\mathbb H$ | $SL(2,\mathbb C)$ | $(\det(z_a,z_b))_{a<b}$ | rank-two skew $2m''\times2m''$ matrices (Plücker cone over $\mathrm{Gr}(2,2m'')$) | $4m''-3=N'-3$ |
>
> In all three cases $\dim X_{\mathbb F}=N'-(m-1)$.
>
> Consequently $\operatorname{rank}DQ\le N'-(m-1)$. The tangency condition $N'=n+m-1$ says exactly that **the quotient has dimension $n$**. Then $L$ maps the cone, generically finitely, onto the image hypersurface of $Q$.

**What this changes.**
1. **The numerology is explained.** The numbers $2m'-1$ and $4m''-3$ are the dimensions of the Segre and Plücker cones. The tangency theorem says: *the flux map is a linear projection of a rank-one cone of the same dimension as the hypersurface it covers.*
2. **The divergence of 20 is explained.** The fibres of $\pi$ are the $K_{\mathbb C}$-orbits: finite for $\mathbb R$, noncompact for $\mathbb C,\mathbb H$. The single lift pushes forward from $V_{\mathbb C}$, so it integrates over these orbits, and that integral diverges (20 Theorem 20.2). For $\mathbb R$ there is nothing to integrate, which is why App. E's engine exists.
3. **`T.KAK.7` gets a precise form.** Build a positive conserved current **on the cone $X_{\mathbb F}$**, where the flux map is linear, instead of pushing one forward from $V_{\mathbb C}$. For $\mathbb R$ the cone is $V_{\mathbb C}/\pm1$, and the paper's current is such a current.

**Checks.** [`scripts/verify_gauge_quotient.py`](scripts/verify_gauge_quotient.py) writes [`results/gauge_quotient.json`](results/gauge_quotient.json) and takes under 1 s.
- **$\mathbb R$.** The factorization through $zz^{\mathsf T}$.
- **$\mathbb C$** ($m'=2,3,4$):
  - each $C^\mu$ pairs the $\pm\mathrm k$ eigenspaces of $J$ and kills each of them;
  - $Q=u^{\mathsf T}Lv$ on random points;
  - the generic rank of $\pi$ is $2m'-1$.
- **$\mathbb H$** ($m''=2,3$):
  - the $K$-invariant quadratic forms span exactly $\binom{2m''}2$ dimensions;
  - every $C^\mu$ lies in their span;
  - the invariant map has generic rank $4m''-3$.
- **Flux rank, both Hermitian families.** The generic flux rank is exactly $n=\dim X_{\mathbb F}$. With one spatial direction removed it stays $n>n-1$, which is 17's failure of tangency.
- **Controls.** A generic real $N=n+1$ pencil has flux rank $n+1$. A random antisymmetric matrix does not commute with the pencil.

---

## 1. The factorization

**Setting.** The three structures are:
- $\mathbb C$: the pencil commutes with a complex structure $J$, which is orthogonal and antisymmetric with $J^2=-1$;
- $\mathbb H$: it commutes with three, $I_1,I_2,I_3$, satisfying the quaternion rules;
- $\mathbb R$: there is no extra structure.

$K$ is the group they generate: $\{\pm1\}$, $\{e^{\alpha J}\}$, or the unit quaternions acting by right multiplication. As in 20 §2, $K_{\mathbb C}=\exp(\mathfrak k\otimes\mathbb C_{\mathrm k})$ preserves $Q$ and acts freely on the nonzero vectors of each extremal eigenspace.

*Proof of Theorem 21.1.*
- **$\mathbb R$.** $Q^\mu(z)=\frac12\langle C^\mu,zz^{\mathsf T}\rangle$. The map $z\mapsto zz^{\mathsf T}$ is the quotient by $\pm1$, with image the rank-$\le1$ symmetric matrices, of dimension $N'$.
- **$\mathbb C$.**
  - Let $U_\pm\subset V_{\mathbb C}$ be the $\pm\mathrm k$-eigenspaces of $J$, each of dimension $m'=N'/2$.
  - **Isotropy.** For $u,u'\in U_+$, $\mathrm ku^{\mathsf T}Cu'=(Ju)^{\mathsf T}Cu'=-u^{\mathsf T}CJu'=-\mathrm ku^{\mathsf T}Cu'$. Hence $u^{\mathsf T}Cu'=0$, and likewise on $U_-$.
  - So $Q^\mu(u+v)=u^{\mathsf T}C^\mu v=\operatorname{tr}(L^\mu\,uv^{\mathsf T})$: linear in $uv^{\mathsf T}$.
  - $K_{\mathbb C}=\{e^{\theta J}\}$ acts by $(u,v)\mapsto(e^{\mathrm k\theta}u,e^{-\mathrm k\theta}v)$. Its invariant polynomials are generated by the $u_av_b$ (weight zero monomials).
  - The quotient is the cone of rank-$\le1$ matrices, of dimension $2m'-1$.
- **$\mathbb H$.**
  - $\mathbb H\otimes\mathbb C\cong M_2(\mathbb C)$, so as a right $SL(2,\mathbb C)$-module $V_{\mathbb C}\cong(\mathbb C^2)^{2m''}$.
  - The $K$-invariant quadratic forms on $V$ are the matrices commuting with $I_1,I_2,I_3$. These are the quaternionic Hermitian forms, of real dimension $m''+4\binom{m''}2=\binom{2m''}2$.
  - By the first fundamental theorem for $SL(2)$ they are spanned by the brackets $[z_a,z_b]=\det(z_a,z_b)$.
  - So $Q^\mu=\sum_{a<b}c^\mu_{ab}[z_a,z_b]$ is linear in the skew matrix $Z^{\mathsf T}JZ$, which has rank $\le2$.
  - This is the affine cone over $\mathrm{Gr}(2,2m'')$, of dimension $2(2m''-2)+1=4m''-3$.

In each case the generic fibre of $\pi$ is a free $K_{\mathbb C}$-orbit, so $\dim X_{\mathbb F}=N'-\dim_{\mathbb C}K_{\mathbb C}=N'-(m-1)$. ∎

**Corollary 21.2 (the "if" half of 17 Theorem 17.4, by invariant theory).** For these pencils, if $N'=n+m-1$ then $\operatorname{rank}DQ_V\le n$ everywhere, so a passive coordinate exists (17 Theorem 17.2).

*Proof.* $DQ=L\circ D\pi$, and $\operatorname{rank}D\pi\le\dim X_{\mathbb F}=n$. ∎

The "only if" half, $N'\ge n+m-1$ always, is 17 Lemma 17.3(2) (the eigen-cone immersion). In the tangency class, $\operatorname{rank}DQ_V(z_0)\ge n$ (17 Lemma 17.3(3)). So $L|_{X_{\mathbb F}}$ has generic rank $n=\dim X_{\mathbb F}$: it is generically finite onto the image hypersurface. For $\mathbb R$ this finite map is the double cover $\pm z$ of the original construction.

---

## 2. Consequences and limits

1. **One picture for 17 and 20.**
   - **Tangency (17):** the quotient has the right dimension.
   - **Divergence (20):** the quotient map has noncompact fibres, and the plane wave sees none of their imaginary directions.

   The two theorems are the two halves of the factorization $Q=L\circ\pi$.
2. **The Hopf pencils.** They are the isotropic members of each row (17 Proposition 17.6), and their sheets are round cones. The octonionic Hopf pencil ($m=8$, $n=9$) passes the tangency test but has no gauge *group*: $S^7$ is not a group. So Theorem 21.1 does not cover it. Its fibres are still complexified 7-spheres (17 Lemma 17.3(4)), and whether its fibre integral diverges is not decided here.
3. **What `T.KAK.7` now asks.**
   - Find a positive conserved current on $X_{\mathbb C}$ (rank-one matrices $M=uv^{\mathsf T}$) or $X_{\mathbb H}$ (rank-two skew matrices), whose pushed time density through the linear map $L$ satisfies (A3)–(A4).
   - On these cones the natural holomorphic coordinates are the entries of $M$. The real displacement $q$, which supplies the Gaussian damping, is **not** a function on the quotient: it is $K$- but not $K_{\mathbb C}$-invariant (20 Theorem 20.2(1)). This is the precise point where a new idea is needed.
4. **What is not claimed.** No new current is constructed, and no new estimate is proved. The theorem reorganizes what is known and isolates the obstruction.

---

## 3. Status

| Result | Status | Uses |
|---|---|---|
| Theorem 21.1 | proved here (classical invariant theory: invariants of tori and the first fundamental theorem for $SL(2)$); checked numerically | 17, 20 |
| Corollary 21.2 | proved here | 17 Theorem 17.2 |
| §2 | discussion | — |

**Subjective confidence:** ≈97% for Theorem 21.1 and Corollary 21.2.

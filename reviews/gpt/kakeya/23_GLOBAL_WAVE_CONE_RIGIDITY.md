# 23. Global curvature at critical state dimension forces a wave cone

GPT, 1 November 2026. **Prose derivation using explicitly cited classical topology and Clifford theory; exact algebraic checks. Not Lean-formalized.** New to this repository; worldwide novelty is not claimed.

The route is a change of language: local eigenvalue geometry becomes a sphere bundle; a short homotopy obstruction forces two spectral levels; their first two moments determine the dispersion and force the matrix multiplication rules. The resulting algebra is a spin factor.

## 0. Global hypotheses

Let $C(\xi)=\sum_{j=1}^n\xi_j C^j$ be a real symmetric $N\times N$ pencil, $n\ge2$. Assume on **every** $\xi\ne0$:

1. $\varphi(\xi)=\lambda_{\max}(C(\xi))$ has the same multiplicity $m$;
2. $D^2\varphi(\xi)$ is positive definite on $\xi^\perp$;
3. $N=n+m-1$.

Apply these conditions after compression if the original system has unused states. They are stronger than the local cap hypotheses of **KAK.RIGID** and [22](22_CURVATURE_AND_QUADRIC_FIBRES.md). In particular, globality must not be inferred from a successful local calculation.

## 1. A sphere bundle whose total space is a sphere

Put $k=n-1$, and form the unit top-eigenvector bundle
$$
\mathcal M=\{(\xi,a)\in S^k\times S^{N-1}:C(\xi)a=\varphi(\xi)a\}.
$$
Constant multiplicity gives a smooth rank-$m$ spectral projection locally, so $\mathcal M\to S^k$ is a smooth $S^{m-1}$-bundle.

**Lemma 23.1 (all unit states are used, exactly once).** Evaluation
$$
\mathcal M\longrightarrow S^{N-1},\qquad(\xi,a)\longmapsto a
$$
is a diffeomorphism. Thus there is a sphere bundle
$$
S^{m-1}\longrightarrow S^{k+m-1}\longrightarrow S^k.
\tag{23.1}
$$
Also $k\ge m$.

*Proof.* Hellmann–Feynman on the eigenspace gives
$$
(a^{\mathsf T}C^j a)_{j=1}^n=\nabla\varphi(\xi).
$$
If a tangent vector to $\mathcal M$ has zero evaluation derivative, differentiating this identity gives $D^2\varphi(\xi)\dot\xi=0$. Since $\dot\xi\perp\xi$, curvature gives $\dot\xi=0$. Evaluation is therefore an immersion; dimensions agree, since
$$
\dim\mathcal M=k+m-1=N-1.
$$
It is a local diffeomorphism. Its image is open and compact, hence all of the connected sphere.

To prove injectivity, suppose the same unit state $a$ is top at two distinct unit directions $\xi,\eta$. The linear functional $\ell_a(\zeta)=a^{\mathsf T}C(\zeta)a$ is at most the convex function $\varphi(\zeta)$, with equality at both endpoints. Convexity forces equality on their segment. If $\eta\ne-\xi$, this contradicts strict transverse curvature along that segment. If $\eta=-\xi$, $a$ is both a top and bottom state of $C(\xi)$; this would make $C(\xi)$ scalar and its top multiplicity $N>m$. Thus evaluation is injective, and a bijective local diffeomorphism is a diffeomorphism.

Finally, the top spaces of $C(\xi)$ and $C(-\xi)$ are the orthogonal top and bottom spaces of $C(\xi)$, each of dimension $m$. They are disjoint, so $N\ge2m$, equivalently $k\ge m$. ∎

## 2. Topology forces two levels; two moments force the pencil

**Theorem 23.2 (global wave-cone rigidity).** Under §0,
$$
\boxed{n=m+1,\qquad N=2m,\qquad m\in\{1,2,4,8\}.}
\tag{23.2}
$$
There are a linear form $\alpha$, a positive definite quadratic form $q(\xi)=\xi^{\mathsf T}Q\xi$, and matrices $D_j=C^j-\alpha_jI$ such that
$$
\boxed{\varphi(\xi)=\alpha(\xi)+\sqrt{q(\xi)},\qquad
D_iD_j+D_jD_i=2Q_{ij}I.}
\tag{23.3}
$$
Both $\alpha$ and $Q$ are recovered from traces:
$$
\alpha_j=\frac{\operatorname{tr}C^j}{N},\qquad
Q_{ij}=\frac{\operatorname{tr}(D_iD_j)}{N}.
\tag{23.4}
$$

*Proof.* First force $k=m$. If $m\ge2$ and $k>m$, the homotopy exact sequence of (23.1) contains
$$
0=\pi_m(S^k)\longrightarrow
\pi_{m-1}(S^{m-1})=\mathbb Z
\longrightarrow\pi_{m-1}(S^{k+m-1})=0,
$$
which is impossible. Here we use only the standard exact sequence for a fibre bundle and the basic homotopy groups below and on the sphere dimension [H]. Lemma 23.1 gives $k\ge m$, so $k=m$. If $m=1$, (23.1) is a connected two-sheeted covering $S^k\to S^k$. For $k\ge2$ the base is simply connected, so such a covering is impossible; hence $k=m=1$. This proves the first two identities in (23.2).

At each $\xi\ne0$, $C(\xi)$ now has exactly two eigenvalues, each with multiplicity $m$: its top eigenvalue and the negative of the top eigenvalue at $-\xi$. Write them as $\alpha(\xi)\pm\beta(\xi)$, $\beta>0$. Taking the trace makes $\alpha$ linear. For $D(\xi)=C(\xi)-\alpha(\xi)I$,
$$
D(\xi)^2=\beta(\xi)^2I,\qquad
\beta(\xi)^2=\frac{\operatorname{tr}D(\xi)^2}{N}=q(\xi).
$$
This is a quadratic form, positive for every $\xi\ne0$. Polarization proves (23.3)–(23.4). Changing spatial variables by $Q^{-1/2}$ makes the $D_j$ a symmetric Clifford family.

For the remaining dimension assertion, [22, Corollary 22.2](22_CURVATURE_AND_QUADRIC_FIBRES.md) supplies a representation of $\mathrm{Cl}(m-1)$ on $\mathbb R^m$, since $k=m$. The classical real Clifford table and its period-eight continuation allow this precisely for $m=1,2,4,8$ [B, §2.3]. Hopf pencils supply examples in each of these four cases. ∎

The converse algebra is immediate: a symmetric Clifford family of size $2m$ in $m+1$ variables, shifted by $\alpha I$ and subjected to an invertible spatial change, has the stated global multiplicity and positive transverse curvature. Thus the possible triples $(n,N,m)$ are
$$
(2,2,1),\quad(3,4,2),\quad(5,8,4),\quad(9,16,8).
$$
No classification up to orthogonal conjugacy is needed or asserted.

## 3. The bare algebra: a spin factor

**Corollary 23.3 (Jordan closure).** The space $\mathbb RI+\operatorname{span}\{D_j\}$ is closed under $X\circ Y=(XY+YX)/2$. In coordinates $sI+D(u)$, its product is
$$
\boxed{(s,u)\circ(t,v)=
(st+u^{\mathsf T}Qv,\;sv+tu).}
\tag{23.5}
$$
Its positive semidefinite cone is $s\ge\sqrt{u^{\mathsf T}Qu}$, a linear image of the Lorentz cone.

*Proof.* Expand the product using (23.3). The eigenvalues of $sI+D(u)$ are $s\pm\sqrt{u^{\mathsf T}Qu}$. ∎

This is the classical **spin factor** [B, §3.3]. Its algebra and its connection to Lorentz geometry are established mathematics. The contribution here is the chain of hypotheses that forces this algebra from a globally curved critical-size pencil. A possible nonlinear-looking dispersion has become a drifted, anisotropic wave cone.

## 4. What fails locally

- The non-Clifford $3\times3$ complex-Hermitian example in [22's checker](scripts/verify_curvature_fibre.py) has $n=5,N=6,m=2$ and positive curvature at the north direction. At the opposite direction its transverse Hessian is
  $$
  \operatorname{diag}(0,0,6,32/3).
  $$
  The global hypothesis fails exactly where required. This is an exact calculation, not a numerical scan.
- The three-state arrow pencil has a simple top at the north pole and a double top at the south pole. It remains a valid local wave representation.
- Repeating a Clifford representation preserves global curvature and doubles multiplicity, but loses $N=n+m-1$.

Consequently this result does not exclude new **local** sheets or solve **T.KAK.7**'s current-construction problem. It says that a globally regular extension at critical state dimension cannot introduce a new dispersion geometry.

## 5. Sources and checks

**[H]** Allen Hatcher, [Algebraic Topology, Chapter 4](https://pi.math.cornell.edu/~hatcher/AT/ATch4.pdf): Theorem 4.41 (fibration exact sequence), Corollaries 4.9 and 4.25 (sphere homotopy groups), Example 4.42 (discrete fibres give coverings). The bundle is smooth and locally trivial here, so its homotopy-lifting hypotheses hold.

**[B]** John C. Baez, *The Octonions*, Bull. Amer. Math. Soc. 39 (2002), 145–205: [§2.3](https://math.ucr.edu/home/baez/octonions/node6.html) (real Clifford representations and period eight), [§3.3](https://math.ucr.edu/home/baez/octonions/node11.html) (spin factors and Lorentz geometry). These are accepted classical ingredients, not new results of this project.

~~~text
python reviews/gpt/kakeya/scripts/verify_global_rigidity.py
~~~

The standard-library script checks exact trace reconstruction, Clifford relations, spin-factor products, and sixteen explicit state-to-frequency selectors for drifted anisotropic Hopf models. It also checks the three hypothesis-removal examples above and the classical dimension-table arithmetic. Results: [global_rigidity.json](results/global_rigidity.json). The topological proof is verified through its cited theorem and hypotheses, not by finite computation.

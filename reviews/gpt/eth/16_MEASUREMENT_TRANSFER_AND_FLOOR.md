# 16. Complete copies, measurement transfer, and an ungraded cost floor

GPT / Astra, 1 November 2026. **Complete prose derivation under stated assumptions; exact finite checks.**

The abstract ungraded design bound is proved in [15](15_UNGRADED_SPLINE_BOUND.md). This note supplies its algorithmic interface and repairs the specific gap identified in [the earlier assessment](../ASSESSMENT.md): a sparsely filled copy need not preserve a measurement row's support. The repair uses sufficiently many **complete copies**, not a claim about every copy.

The resulting floor concerns the published measurement-and-Taylor **enumeration bound**. It is not a lower bound on SAT algorithms, on all implementations of deformation counting, or on the number of distinct attainable measurement vectors after further compression.

## 1. Labels, copies, and the support injection

Use the notation of [Claude 08, §6](../../claude/eth/08_GAPLESS_PROOF.md) (`ETH.I.TRANSFER`). Labels are
$\mathcal B=([r]\times[r])^Q$, of cardinality $L_{\mathcal B}=r^{2Q}$. Incidence $(f,v)$ carries label
$u_{f,v}=((a_v(i),d_f(i)))_{i=1}^Q$. Incidences with a given label are numbered consecutively from zero. Thus copy $\kappa$ contains label $u$ exactly when $\kappa<c_u$, where $c_u$ is that label's occupancy.

The filling map reads the variable at an occupied position and zero at an empty one. A measurement row is a tensor product of local rows
$$
e_{aa},\quad e_{ab}+e_{ba},\quad e_{ab}-e_{ba}\qquad(a<b).
$$
The last type is antisymmetric. Write $j(\psi)$ for its number of antisymmetric factors. The actual group row is
$(\psi^{\mathsf T}\pi_\kappa\iota)|_{S_{i,g}}$, and zero rows are discarded.

**Lemma 1.1 (support injection).** The map sending a label in $\operatorname{supp}\psi$ to its first-coordinate schedule $(a(1),\ldots,a(Q))$ is injective. Consequently no two occupied labels in this support, in one copy, can read the same variable.

*Proof.* In a diagonal local factor the first coordinate determines the only supported pair. In an off-diagonal symmetric or antisymmetric factor, first coordinate $a$ forces pair $(a,b)$, while first coordinate $b$ forces $(b,a)$. Therefore the complete first-coordinate schedule determines every pair in the supported label. All incidences of one variable share that schedule. Distinct supported labels consequently read distinct variables. ∎

In particular, summing entries into variable columns causes **no sign cancellation for one such row**. Empty positions remain a real issue; the lemma does not fill them.

**Lemma 1.2 (one complete copy).** Suppose copy $\kappa$ contains every label. At node $i$, every tensor eigenvector contributes at least one nonzero group row. Each nonzero group row has
$$
\ell^1\text{-norm}\ \ge\ 2^{\max(j(\psi)-1,0)}.
\tag{16.1}
$$

*Proof.* Let $o$ be the total number of off-diagonal local factors, symmetric or antisymmetric. Then $o\ge j$ and the support has $2^o$ labels, each with coefficient $\pm1$. Fixing the first coordinate at node $i$ either leaves all labels (a diagonal factor) or half (an off-diagonal factor). Every nonempty restriction therefore has $2^o$ or $2^{o-1}$ entries. Completeness keeps them all, and Lemma 1.1 prevents their cancellation under filling. This proves (16.1), including $j=0$. ∎

## 2. Aggregate cost and a sufficient filling condition

**Definition 2.1 (the enumeration model).** Fix uniform Taylor budget
$q=\lfloor MD/Q\rfloor$, where $D$ is the deformation-degree bound and $M$ is the incidence/port degree budget. Put
$$
C_q(n)=\sum_{j=0}^{\min(q,n)}2^j\binom nj,
\qquad
\mathcal C=C_q(n)(2\beta+1)^{c_*}2^s,
\qquad \delta=n^{-1}\log_2\mathcal C.
\tag{16.2}
$$
Here $s$ is the largest group size, $c_*$ the largest total number of nonzero group rows at any node, and $\beta$ the largest row $\ell^1$-norm, taken to be at least one. Rows are counted **with the multiplicity used by the transfer construction**. Polynomial factors and the sum over nodes do not reduce this expression.

This is the exponential part of the bound in Claude 08, Theorem 4.6. Selecting a smaller basis of rows, exploiting correlations between their values, using a smaller actual degree budget, or avoiding the Cartesian enumeration can change this model and is not ruled out below.

**Proposition 2.2 (aggregate measurement bound).** Suppose there are $k$ complete copies. If, at some node, at least $L_{\mathcal B}\varepsilon$ measured eigenvectors have $j(\psi)\ge j_0$, where $j_0$ is a nonnegative integer, then
$$
\delta\ge\frac{kL_{\mathcal B}}n\max(j_0,1)\varepsilon.
\tag{16.3}
$$
In particular, if $M\ge n$ and $k\ge M/(2L_{\mathcal B})$, then
$\delta\ge\max(j_0,1)\varepsilon/2$, the earlier cost condition C3.

*Proof.* Lemma 1.2 gives $c_*\ge kL_{\mathcal B}\varepsilon$. If this number is nonzero, it also gives
$\log_2(2\beta+1)\ge\max(j_0,1)$. Multiplying these bounds in (16.2) proves the claim. If there are no such eigenvectors, the right side is zero. ∎

**Proposition 2.3 (variance slack supplies complete copies).** Suppose the schedules satisfy
$$
\sum_{u\in\mathcal B}(c_u-M/L_{\mathcal B})^2<2M(b+w),
\qquad M\ge8(b+w)L_{\mathcal B}^2,
\tag{16.4}
$$
where $b$ and $w$ bound occurrences and widths. Then at least $M/(2L_{\mathcal B})$ copies are complete.

*Proof.* Every label satisfies
$$
c_u>M/L_{\mathcal B}-\sqrt{2M(b+w)}\ge M/(2L_{\mathcal B}).
$$
The first $\min_uc_u$ copies are complete because numbering starts at zero for each label. ∎

The conditional-expectation construction in Claude 08, Lemma 6.2 already supplies the variance inequality in (16.4). Its stated size assumption $M\ge2(b+w)L_{\mathcal B}^2$ alone is weaker than the sufficient slack used here; it should not be silently replaced by (16.4). The paper's prescribed Part I parameters have $L_{\mathcal B}^2=O(\sqrt n)$ and $M\ge n$, with fixed $b,w$, so they satisfy the stronger assumption for all sufficiently large $n$.

This resolves the earlier witness correctly: the last sparse copy can lose or annihilate rows, while the many complete copies still supply the aggregate cost. It does not retract the witness or prove C3 for every filling allowed by the weaker assumptions.

## 3. The floor for arbitrary group-algebra designs

**Theorem 3.1 (ungraded enumeration floor).** Assume:

1. $r\ge2$, $Q\ge1$, $0\le D<Q$, and $r^{2Q}\le n$;
2. $P(T)$ is any group-algebra swap design of degree at most $D$, with $P(0)=I$, at nodes $1,\ldots,Q$;
3. at each node, every eigenvector on which its defect is nonzero is measured;
4. $M\ge n$, and at least $M/(2L_{\mathcal B})$ copies are complete;
5. the cost is the enumeration model (16.2), with $r$ groups partitioning the $n$ variables at each node.

Then, for all sufficiently large $n$,
$$
\boxed{\quad\delta\ge\frac1{1000}
\frac{(\log_2\log_2n)^2}{\sqrt{\log_2n}}.\quad}
\tag{16.5}
$$
No grading assumption is needed. The threshold is not optimized or made explicit here.

*Proof.* Write $L=\log_2n$, $\ell=\log_2L$, and $x=D/Q$. Three lower bounds apply.

**Groups:** $s\ge n/r$, so $\delta\ge1/r$, and hence
$Q\le L/(2\log_2r)$.

**Taylor terms:** for $0<x\le1/e$,
$$
\delta\ge x\log_2(1/x)-O(\log n/n).
\tag{16.6}
$$
Indeed $q\ge nx-1$, and $\binom nk\ge(n/k)^k$ for integer $k\le n$. Rounding changes the logarithm by at most $O(\log n)$. If $x>1/e$, selecting $k=\lfloor n/e\rfloor\le q$ instead gives a positive constant lower bound on $\delta$.

**Measurements:** if $Q\ge512$ and $D\le Q/16$, Theorem 3.4 of [15](15_UNGRADED_SPLINE_BOUND.md) and Proposition 2.2, with $j_0=\lceil Q/32\rceil$, give
$$
\delta\ge c_0e^{-4D^2/Q},\qquad c_0=3e^{-25}/512.
\tag{16.7}
$$
Here the eigenvalue weights are $p_r\in[1/4,1/2)$, and the factor $Q$ from a heavy row's norm cancels the directional averaging loss in (15.8).

Suppose, for contradiction along arbitrarily large $n$, that
$\delta<\ell^2/(1000\sqrt L)$; in particular $\delta\to0$.
If $D=0$, then $P(T)=I$, every defect eigenvector is measured, and Proposition 2.2 with $j_0=0$, $\varepsilon=1$ gives $\delta\ge1/2$, a contradiction. If $D\ge1$ but $Q<512$, or if $D/Q>1/16$, the Taylor bound gives a positive constant, again a contradiction. Thus (16.7) applies.

For sufficiently large $n$, $\delta\le c_0^2$, so (16.7) implies
$$
\frac{D^2}Q\ge\frac18\ln(1/\delta).
$$
Since $\log_2r\ge\log_2(1/\delta)$,
$$
x^2\ge\frac{\ln(1/\delta)}{8Q}
\ge\frac{\ln2}{4L}(\log_2(1/\delta))^2.
$$
Our assumed upper bound on $\delta$ gives $\log_2(1/\delta)\ge\ell/4$ eventually. Hence
$$
x\ge x_0:=\frac{\sqrt{\ln2}}8\frac\ell{\sqrt L}.
$$
The function $x\log_2(1/x)$ increases on $(0,1/e)$, and eventually $\log_2(1/x_0)\ge\ell/4$. Formula (16.6) therefore gives
$$
\delta\ge\frac{\sqrt{\ln2}}{32}\frac{\ell^2}{\sqrt L}
-O(\log n/n),
$$
which contradicts the assumed coefficient $1/1000$. ∎

**Corollary 3.2 (a checkable repair).** Conditions (16.4), together with the other hypotheses of Theorem 3.1, imply the same floor. Thus complete-copy slack repairs the transfer step and extends the enumeration floor from graded to arbitrary group-algebra designs.

## 4. What remains open

- Can complete-copy slack be weakened to an aggregate occupancy hypothesis that still proves the required cost? The variance argument is sufficient, not necessary.
- Can correlated measurements be encoded substantially more cheaply than the Cartesian product in (16.2)? A theorem about that expression is not a barrier to such an improvement.
- Can deformations outside the swap group algebra avoid the parity-trace obstruction? The general `T.ETH.5` question is still open and needs a precisely restricted worst-case model; trivial already-local readout systems cannot satisfy a universal nontrivial floor.
- Can the weighted-jet obstruction in [15, Theorem 5.2](15_UNGRADED_SPLINE_BOUND.md) identify a useful tradeoff for a hybrid interpolation scheme? No improved counting algorithm is claimed here.

[`scripts/verify_complete_copies.py`](scripts/verify_complete_copies.py) checks support injection and transferred norms with shared variables, exact spectral multiplicities, complete-copy row counts, and sparse-copy counterexamples. It uses the standard library and writes [`results/complete_copies.json`](results/complete_copies.json). These are finite checks; Theorem 3.1 is an asymptotic prose proof.

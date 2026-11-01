# Research ledger: mechanisms, failed routes, and precise next questions

GPT / Astra, 1 November 2026. This records the continuation after the first pushed snapshot, as well as the routes that led to it. “Proved” means a complete prose argument with the stated hypotheses. The new work is not Lean-formalized, and worldwide priority is unclaimed. A finite exploration cannot exhaust all mathematical routes; the purpose here is to preserve useful successes, obstructions, and unfinished interfaces rather than present an unlimited search as complete.

## Round 8: the phone continuation

- **[ETH 21](eth/21_PRODUCT_BIAS_ENERGY.md), T.ETH.8 DONE.** Diagonalize P_i^*P_i rather than the non-self-adjoint biased swaps. The exact optimum is a Bernoulli spectral thinning of scalar capacities, with a completely specified optimizer. A full-retention event alone gives a weak Q-dependent comparison; positive Cauchy-Binet sums improve the loss to theta^D. For the actual theta=1-1/r^2, relative gain is at most D/r^2 and vanishes in Part I. Exact original vertex QPs check the reduction independently. This is energy, not support or arithmetic locality.
- **[Kakeya 22](kakeya/22_CURVATURE_AND_QUADRIC_FIBRES.md).** A gauge group is unnecessary for central divergence. The second eigenvalue variation in the critical state dimension is a square orthogonal design. Its determinant computes the complexified eigensphere coarea measure, and the time entry cancels the residue denominator. All m>1 rigid full fibres diverge, including the octonionic case. The same identity forces a local Clifford module. This is a necessary local structure, not a classification or an integrated current on a quotient.
- **Boundaries preserved.** General infinite central density is weaker than ambient non-integrability. Scalar energy does not count exact defects. Neither result supplies a faster counting algorithm or a new dispersive estimate. No new Lean build is claimed. Claude 18 supersedes the old complete-copy qualification, and Claude 20 repairs the old real-slice-to-ambient inference.

## 1. Results worth retaining

| Route | Concrete result | Remaining interface |
|---|---|---|
| Replace a graded Fourier argument by a spline algebra | [ETH 15](eth/15_UNGRADED_SPLINE_BOUND.md): arbitrary swap designs have a binomial-tail defect bound; random restriction gives the biased version | Does a sharper bound hold for the maximum direction? |
| Treat failed local conditions as a dimension budget | [ETH 15, §5](eth/15_UNGRADED_SPLINE_BOUND.md#5-unequal-jet-orders): unequal retained jet orders obey a weighted subset bound | An actual hybrid counting scheme using those jets |
| Identify the missing measurement-transfer hypothesis | [ETH 16](eth/16_MEASUREMENT_TRANSFER_AND_FLOOR.md): complete copies preserve support and give an aggregate enumeration floor | Compression of the actual measurement image is outside that floor |
| Abstract parity completion to multiplication of traces | [ETH 17, §§1–2](eth/17_TRACE_CAPACITY_AND_APPROXIMATE_LOCALITY.md): multiplication-rank certificate and a joint character-completion theorem | An algorithmic source of several traces sharing the same good edges |
| Replace exact zero by divisibility | [ETH 17, §3](eth/17_TRACE_CAPACITY_AND_APPROXIMATE_LOCALITY.md#3-what-survives-reduction-modulo-a-prime): odd-residue unit defects satisfy the same support bound | Colliding nodes, ramified arithmetic, denominators, and characteristic two need different arguments |
| Replace defect support by defect size | [ETH 17, §§4–6](eth/17_TRACE_CAPACITY_AND_APPROXIMATE_LOCALITY.md): sharp uniform-measure energy, a discrete endpoint kernel, and its square-root transition | Biased energy and an arithmetic error/rounding transfer theorem |
| Inspect the noncompact symmetry before estimating a fibre integral | [Kakeya 18](kakeya/18_NONCOMPACT_FIBRE_OBSTRUCTION.md): exact central divergence on the Weyl gauge fibre | A different reduction must preserve positivity and conservation |
| Test whether the central divergence is harmless almost everywhere | [Kakeya 19](kakeya/19_TRANSVERSE_SINGULARITY.md): inverse-square transverse singularity, failure of local integrability, and an explicit logarithmic coefficient | A quotient or altered current, rather than removal of an exceptional axis |

## 2. The most reusable new algebraic mechanism

The basic datum is a filtered algebra of compatible local objects, together with a trace into a finite-dimensional function algebra. In the cube case, the full system realizes exactly the low Fourier characters, while deletion of one scalar compatibility equation enlarges the trace space by at most one dimension.

A new low-degree trace acts by multiplication. Its image outside the old trace space certifies how many compatibility equations must have been removed. With several traces, products can produce additional independent directions. In ETH 17's four-coordinate example, two degree-one traces separately certify two missing constraints each; their joint linear span certifies four, and their product forces a fifth.

This is more informative than saying that the proof “uses dimensions.” It identifies the operator to compute: multiplication followed by quotient by the baseline trace filtration. For character traces it becomes a union of Hamming balls with costs assigned to generators. For arbitrary traces it is an explicit finite matrix. The same pattern can be sought in other spline or interpolation algebras, but a new application must prove its own trace filtration and its own constraint-cost estimate.

No theorem here identifies this algebraic mechanism with the Hermitian-current mechanism. There is a useful shared habit—inspect what the intermediate representation can carry—but a rank quotient and a divergent positive integral are different mathematical objects. A forced unification would obscure rather than extend the results.

## 3. Sharpness experiments and what they did not establish

The exact exhaustive search [explore_directional_extremizers.py](eth/scripts/explore_directional_extremizers.py), with nodes $1,\ldots,Q$, found:

| $Q$ | Degree $D$ | Minimum possible maximum uniform directional defect mass |
|---|---|---|
| 2 | 1 | $1/2$ |
| 3 | 1 | $1/2$ |
| 3 | 2 | $1/4$ |

For these small instances the search enumerates every candidate deleted-edge set below the optimum and solves the remaining equations exactly. The [saved witnesses](eth/results/directional_extremizers.json) include every vertex polynomial. These are finite computational conclusions, not an asymptotic theorem.

The search suggested the explicit one-edge family in [ETH 15, Proposition 2.4](eth/15_UNGRADED_SPLINE_BOUND.md#2-parity-completes-the-trace-space). At degree $Q-1$, the total missing-edge bound is exactly sharp. This does **not** settle whether the factor $Q$ in passing from total defects to the largest directional defect can be removed in the low-degree range. The one-edge construction has its failures concentrated in one direction, and naive averaging of different constructions can turn zeros into nonzeros rather than distribute a fixed support budget.

## 4. Routes that fail without additional structure

### 4.1 “Small modulo a prime” instead of zero

At an odd prime with separated nonzero node residues, reduction to the residue field preserves the full character basis and the spline dimension argument. Some defects must remain units. Thus merely replacing equality by high divisibility does not evade the obstruction within this arithmetic setting.

For nodes $1,\ldots,Q$, a fixed small prime does not satisfy the separation hypothesis when $Q$ grows. The explicit family $n_x=1-T$ at nodes $1+ip$ shows why this is a real boundary. Such a collision can defeat the residue support statement, but it does not establish efficient interpolation or exact counting at the required modulus.

### 4.2 Characteristic-two nilpotence

The swap group algebra is not semisimple in characteristic two. Its elementary augmentation operators can be nilpotent, which initially looks promising for reducing degree. However, a nilpotent operator can still have nonzero scalar readouts. The $2\times2$ example in ETH 17 separates these facts explicitly. A counting scheme would have to control the readouts and recovery loss, not merely show that a matrix power vanishes.

### 4.3 Common rational denominators

Suppose a rational vertex family is $n_x(T)=p_x(T)/q(T)$ with a **common** denominator nonzero at zero and at every node. Normalize $q(0)=1$. Then $p_x(0)=1$, and every bad-edge event is exactly the numerator event $p_x(t_i)+p_{x\oplus e_i}(t_i)\ne0$. ETH 15 applies with the numerator degree. A common denominator therefore does not remove this support bound.

Vertex-dependent denominators change the edge equations to cross-products and destroy that immediate fixed-degree spline identification. This is a possible mathematical escape from the proof, not a constructed algorithm. No recovery or bit-complexity bound was obtained for it.

### 4.4 More than two local character values

The binary complement identity is special: multiplication by full parity changes every supported coordinate to an unsupported one and vice versa. For a larger character group, adding a fixed nonzero frequency need not erase an arbitrary nonzero frequency. Therefore the half-degree completion argument does not transfer verbatim to larger alphabets.

The multiplication-rank formulation is the correct starting point: compute translated low-weight character sets or their spans and account for products of generators. This suggests a coding-theoretic direction, but this round does not claim a sharp general-alphabet theorem.

### 4.5 Energy as a shortcut to support or running time

The exact uniform energy minimum is classical polynomial approximation after Fourier reduction. It gives an amplitude obstruction but does not count nonzero entries. Conversely, the rank obstruction is insensitive to how small nonzero real defects are. Under a biased product law the Fourier orthogonality used in the energy proof changes, so the scalar optimum is not automatically the optimum there.

Neither statement makes small real errors harmless in an exact arithmetic counting procedure. A meaningful approximate-locality algorithm needs a separate error/recovery theorem. This interface remains open.

### 4.6 Compressing correlated measurements

The product of individual row ranges can greatly overcount a linear map's actual image. For example, repeating one row $m$ times gives only the values of that one row, not the $m$-fold Cartesian product. This elementary observation is why ETH 16 names its enumeration model explicitly.

For a real linear map of rank $k$, its image on a Boolean cube has at least $2^k$ points: select $k$ independent columns and vary only those coordinates. This lower bound is often much smaller than a product-of-ranges bound and does not supply a comparable upper bound. In the transfer construction, preservation of each individual row's support does not prove independence of different rows after variables are identified. No rank-preservation theorem or improved image-enumeration algorithm was established here.

### 4.7 Dividing by gauge volume or ignoring the singular axis

On the Weyl fibre the time integrand is not gauge-invariant, so “divide by the infinite gauge volume” is not a definition of the same positive push-forward. A radial cutoff also produces a divergence source. The further calculation in Kakeya 19 shows that the singularity is not removable by redefining the density on a null set: its nearby positive mass diverges logarithmically.

A gauge section or a compact reduction may still work, but the measure, tensor, and conservation law must be constructed together. A proof of tangency alone leaves those requirements untouched.

### 4.8 The broad square-root law and the remaining logarithm

Already-local readout systems defeat an unrestricted universal floor. A meaningful general theorem must specify a nontrivial worst-case class and a cost model. The new swap results provide one such obstruction inside a specified algebra, not the general theorem.

The one-logarithm local-smoothing target also remains open. The previously distilled scalar propagation estimate is sharp for its stated class. A stronger result would need extra structure, an altered quantity, or a different propagation theorem. Neither the spline calculation nor the Weyl obstruction supplies that missing analytic estimate.

## 5. A concrete next programme

1. **Shared traces:** find a counting or interpolation construction that forces two or more trace characters to use the same exceptional edges. Optimize the weighted Hamming-ball union and compare it with separate parity bounds.
2. **Biased energy:** solved exactly in ETH 21, including the actual-bias relative asymptotic. Efficient evaluation of the subset formula or sharper fixed-bias asymptotics would be further questions.
3. **Measurement image:** use actual filling maps from complete schedules to study rank and integer image size. A useful theorem must address correlations across rows, not just the support of each row.
4. **Conserved reduction:** start with the explicit Weyl coordinates and seek a finite positive current on a quotient, tracking the divergence identity through the reduction. The transverse singularity provides an exact test that any candidate must pass.

These are unfinished research questions, intentionally left precise. All completed arguments, checks, and the small search are included in the repository rather than left in chat.


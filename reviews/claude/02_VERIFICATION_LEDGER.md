# 2. Verification ledger

Status key:
- **R**: re-derived independently, line by line.
- **X**: checked in exact arithmetic by a script in [`scripts/`](scripts/), with output in [`results/`](results/).
- **A**: accepted from a standard or cited source without re-proof.
- **P**: partially checked.

"No issue" means I found no mathematical error after the stated check.

## A. The ETH-refuting chain

| Item | Content | Status | Notes |
|---|---|---|---|
| §2.2, (2.2) | $\#F=\sum_x\prod_e\Phi_e$, $\Phi_e=1-\prod(1-L_{e,v})$ | R | Dummy factors $\Phi=1$ give every variable a position. |
| Def. 3.1 | design: partitions, assignments $\sigma_f$, integer measurements, kernel condition | R | Kernel form ⇔ (2.5): $P_{f,v}(t)\vert_{S_{t,h}}\in\operatorname{row}(C_{t,h})$ for $h\neq\sigma_f(t)$. |
| Lemma 3.2 | measurement separation | R, X | Representatives need not satisfy the $E$ constraints, because invariance is over $\mathbb Q^n$. Cost $2(2D+1)^c\sum_g2^{\lvert S_g\rvert}$. |
| (3.7)–(3.9) | shifted inputs separate; $\deg_X[z^j]\mathcal G_x\le j$ | R, X | Every entry of $\Delta P(z)$ is divisible by $z$. |
| Lemma 3.3 | jets from Boolean substitutions (Möbius inversion) | R, X | Applied over $\mathcal R=\mathcal A[z]/(z^J)$, where the degree in $X$ is $\le q$. |
| §3.4 | Hermite recovery; $QJ>d_H$ | R, X | The confluent Vandermonde matrix is invertible over $\mathbb Q$, coordinatewise in a $\mathbb Q$-algebra. |
| Prop. 3.4, (3.14)–(3.16) | three costs | R | $\eta^qC_q\le(1+2\eta)^n$ with $\eta=R/Q$. |
| (4.6)–(4.7) | $W_0$, $V_0$, $f_{\mathcal I}$, $M_0=V_0f_{\mathcal I}(W_0)$ | R, X | $f_{\mathcal I}(A_0+q)=(-1)^{A_0+q}$ for $0\le q<R$, by the binomial theorem. |
| Lemma 4.1 | $c_{\mathcal B}\le D\Pr[W\notin\mathcal I]$; $M_0$ fixes $\ker C^{\mathcal B}$ | R, X | Exact: rank$(I-M_0)=c_{\mathcal B}$ in all 7 cases. Rows have $\ell_1\le2^Q$ (see S1). |
| Lemma 4.2 | support: $M_0[u,y]\neq0\Rightarrow\lvert E_{u,y}\rvert\le R-1$ | R, X | |
| Lemma 4.3 | $P_u(0)=\delta_u$, $\deg\le R-1$, kernel/locality | R, X | Checked per group and jointly (single coefficient vector). |
| Remark 4.4 | worked example $Q=R=r=2$ | X | The displayed $M_0$ row and $P_u(T)$ match exactly. |
| Lemma 4.5 | balanced schedules via conditional expectations | R | $\mathbb E\Phi\le M(b+w)$ and $\mathbb E\Xi=Qn(1-1/r)$. The final potential is $\le2$, which gives both bounds under (4.13). |
| §4.5, (4.14)–(4.16) | copies, $\iota$, $C_{t,g}$, kernel transfer | R, X | $\sum_gc_{t,g}=rhc_{\mathcal B}\le2r\varepsilon M\le2\tau n$ with $\varepsilon=\tau/(rb)$. Repeated variables in one copy are harmless. |
| Lemma 4.6 | common denominator $(R-1)!Q!$, $O(Q\log Q)$-bit coefficients | R | $\lVert\varphi_{u,y}\rVert_1\le Q+1$ (sharper than $2^Q$). |
| Lemma 4.7 | assembly | R | |
| (2.19), Lemma 4.8 | parameters | R, X | All inequalities hold eventually. Exact thresholds are in [`results/parameters.json`](results/parameters.json). |
| §2.8 | 3-CNF via OR gates; Tovey cycles | R | Gate values are forced, so the count is preserved. Each copy has exactly 3 occurrences. |
| Lemma 6.3 / App. A | disjoint sparsification with $B=O_k((\sigma^{-1}\log(2/\sigma))^{k-1})$ | R | Lemma A.1 (partition, survival), Lemma A.2 (the two priorities), Lemma A.3 ($T_j\le4(3M)^{j-1}n$), tree size $(Cn+1)^22^{\sigma n}$. |
| App. B.4 | bit complexity and space, overview algorithm | R | Everything stays polynomial: there are polynomially many denominators, and exponentially many summands add only logarithmically many bits. |

**End-to-end** ([`scripts/verify_counting_pipeline.py`](scripts/verify_counting_pipeline.py)). I ran 17 random bounded-occurrence CNFs, in 8 parameter regimes, with $n\le8$ and $r\in\{2,3\}$, $Q\in\{2,3,4\}$, $R\in\{2,3\}$. The runs cover:
- many copies, with $h$ up to 3;
- nontrivial measurement buckets, up to 32 assignments per bucket;
- truncated inclusion–exclusion, with $q<n$.

Each run builds the design, checks Definition 3.1 exactly, and runs Lemma 3.2 per pair $(E,U)$, then (3.13) and Hermite interpolation. It returned $H(0)=\#F$ in all 17 runs, 11 of them with a positive count. Symbolic $H(T)$ matched every Taylor coefficient, and $\deg H\le d(R-1)$.

As a negative control, I deleted one measurement row per group. This broke the kernel condition in all 3 trials and gave a wrong count in 2 of 3.

## B. The $\rho(n)$ bound for arbitrary density (Prop. 6.4)

| Item | Status | Notes |
|---|---|---|
| Lemma 5.1 (complemented Beigel–Tarui $p_\ell$, division by $h=Y(Y-1)$) | R, X | Congruence checked for $\ell\le10$, $\lvert y\rvert\le200$. Decomposition exact. |
| Lemma 5.2 (cross-equation truncation) | R, X | Congruence checked on random tuples. Retained degree is exactly $m+2(\ell-1)$, versus $m(2\ell-1)$ untruncated. |
| (5.4)–(5.5), Prop. 5.3 | R | $\mathcal S_\ell$ is linear and applied after interpolation. Each $\mathcal S_\ell(\cdot)$ term is an integer at Boolean points. |
| §5.5, proof of Thm 2.2 | R | Branching on $\lvert S\rvert\le n/\sqrt\lambda$ high-incidence variables; $n'\ge n/2$; $b,w\le\lambda^{O(1)}$. |
| Lemma 6.1 (Tanner code, $\varepsilon$-biased Cayley graph) | R | Mixing (6.3), local code by union bound, distance $\ge m/P$, row count $\le3n$. All constants check: $\varepsilon^2D_0/2=4\mu$ and $\delta(\delta-\varepsilon)=1/(20000B^2\mu^2)$. |
| Lemma 6.2 (greedy rows via exact counts) | R | $\mathbb E_\gamma D(A;h_\gamma)\le\frac34D(A)$, and $(27/32)^n<1$. |
| §6.3–6.4, Prop. 6.4 | R | $2^{(\log n)^{O(1)}}$ oracle calls; $\sigma=1/s$, $s=\lfloor\sqrt{\lambda}/\nu^2\rfloor$. |
| App. B.1–B.3 | R | Bareiss elimination; the minors bound the intermediate integers. |
| §2.9 randomized illustration | R | |

## C. Appendix D (sharper exponent)

| Item | Status | Notes |
|---|---|---|
| D.4.1 coefficient ring $\mathbb Z_{(2)}[U]/(f)$ | R | Irreducibility mod 2 gives the valuation and units. |
| D.4.3, Lemma D.4 | R, X | $W(0)=I$; $W(t)$ is shift-by-$t$ supported; $W(t+\alpha z)$ is integral; off-shift entries are even ($a\le3$). |
| (D.10) | R, X | The identity holds exactly and $\ge-6Dh$ ($a\le4$, all $u\le a$). |
| Lemma D.5 (tapered Hermite recovery) | R, X | Recovered $F(0)\bmod2^p$ from random lifts of the data mod $2^{p_t}$ in 46 trials. Negative control: with $p_t=p$, recovery failed in 46/46. |
| (D.12) geometric mode expansion | X | Error valuation $\ge28$ after 12 terms, in all 32 crossing cases tested. |
| (D.13)–(D.15) mode counting $\exp(O(M/\sqrt u))$ | R | Generating-function bound with $\theta=\sqrt{M/(Dp_t)}$. |
| Lemma D.3 packing; D.4.6 children; D.4.7 totals; (D.11) precision | **P** | Followed the logic; did not re-derive every inequality in (D.4)–(D.6) and (D.11). |
| Lemma D.2 (sparse systems as tensor networks) | P | The degree-array encoding with $M=O(n\nu)$ wires is plausible; not fully checked. |

## D. Corollaries and appendices C

| Item | Status | Notes |
|---|---|---|
| Cor. 1.2 (Max-$k$-SAT) | R | Hoeffding plus a union bound over $2^n$; linear-size counting circuit; threshold search. |
| Cor. 1.3 | R | Parsimonious encodings. |
| Cor. C.1 | A | Uses [Wil13, Thm 6.1]: "ETH false ⇒ E^NP ⊄ SIZE(O(n))". Not checked against the source. |
| Cor. 1.4 | P | Uses App. D's fixed-density bound plus the [JMV13] local reduction. The arithmetic $C_W(dn)>\eta n\sqrt{\log_2 n}$ looks consistent; not fully checked. |
| Cor. C.2 | R | |

## Findings

**No mathematical errors found** in parts A or B. The proofs are careful. All constants I checked are consistent.

Sharpenings and remarks. None of these changes a result.

- **S1.** Each row of $C^{\mathcal B}$ is a tensor product of $Q$ local vectors with support 1 or 2. Its absolute row sum is therefore $\le2^Q$, not merely $\le D=r^{2Q}$. Measurements then lie in $[-2^Q,2^Q]$, and the measurement exponent improves from $2\tau n\log(2D+1)$ to $2\tau n(Q+2)$. The improvement is a $\log r$ factor in a term that is not the bottleneck.
- **S2.** $\operatorname{rank}(I-M_0)=c_{\mathcal B}$ exactly, so Lemma 4.1 is tight for this $M_0$. Section 4 of [04](04_MIRACLES_EXPLAINED.md) shows that the interval choice is optimal among $\mathrm{Sym}(Q)$-invariant constructions up to a $\mathrm{poly}(Q)$ factor.
- **S3.** Rows have degree $\le R-1$, but Def. 3.1 records $\le R$. Hence $d_H=d(R-1)$ suffices, and the scripts use both values.
- **S4.** The matrix formula of [05, Theorem B](05_DISTILLATION.md) gives an alternative family of valid rows. It is equivariant and has a closed form, and it was verified exactly.
- **Notation.** $D$ denotes both $\lvert\mathcal B\rvert$ and the measurement $\ell_1$ bound; this is harmless but confusing (see S1). §2 indexes groups in $[r]$, while Remark 4.4 uses $\{0,1\}$; the remark says so.
- **Related work** (not correctness). Earlier "beating brute force for polynomial systems over $\mathbb F_2$" algorithms are not discussed: Lokshtanov–Paturi–Tamaki–Williams–Yu (SODA 2017), Björklund–Kaski–Williams (ICALP 2019), Dinur (SODA 2021). They use parity-counting and polynomial-method ideas related to §5. A survey of this area belongs in a final version.

# Exact-arithmetic checks

The scripts use Python 3.10+ and the standard library only (`fractions`, `itertools`). Each is deterministic and writes JSON to `../results/`.

| Script | Checks | Runtime |
|---|---|---|
| `eth_core.py` | Library: auxiliary array $\mathcal B$, $W_0$, $V_0$, $f_{\mathcal I}$, $M_0$, $C^{\mathcal B}$, rows $P_u(T)$ (Sections 4.1–4.3), random bounded-occurrence CNFs, brute-force counts. | — |
| `verify_replacement_forms.py` | Lemmas 4.1–4.3 and Remark 4.4, for $(r,Q,R)\in\{(2,2,2),(2,3,2),(2,3,3),(2,4,2),(2,4,3),(3,2,2),(3,3,2)\}$. Exact rank of $I-M_0$. Negative controls ($f\equiv1$; wrong interval). | ~6 s |
| `verify_counting_pipeline.py` | The full algorithm of §3–§4 against brute force on 17 CNFs, plus a symbolic $H(T)$ and a deleted-measurement negative control. | ~1 min |
| `verify_parity_lift.py` | Lemmas 5.1 and 5.2 (congruences, decomposition, the degree $m+2(\ell-1)$) and the coefficient sizes of Lemma B.1. | <1 s |
| `verify_parameters.py` | (2.19) and Lemma 4.8: least $\lambda=\lceil\log_2n\rceil$ at which all hypotheses hold and at which the bound beats $2^n$. | <1 s |
| `verify_appendix_d.py` | Lemma D.4, identity (D.10), Lemma D.5 (with a precision negative control), expansion (D.12). | ~2 s |
| `verify_distilled_swap_formula.py` | Theorem B of `05_DISTILLATION.md` (closed form), and Proposition B″ (symmetric optimality). | ~1 s |
| `verify_padic_truncation.py` | TOOLKIT E2: the $p$-adic cross-equation truncation for $p\in\{2,3,5,7\}$, $\ell\le6$ (congruences, decomposition, truncation, degree bound). | <1 s |
| `verify_cost_floor.py` | TOOLKIT E4: minimizes the method's three costs over integer parameters with a rigorous window bound; compares with $(\log\log n)^2/\sqrt{\log n}$. | ~2 min |

Run them from this directory, for example `python3 verify_counting_pipeline.py`.

# Checks for the Kakeya / local smoothing paper

Requirements: Python 3.10+, `numpy`, `sympy`, `scipy`. Each script is deterministic and writes JSON to `../results/`.

| Script | Checks | Runtime |
|---|---|---|
| `verify_wave_current.py` | The plane-wave phase in null coordinates (2.1); the symbol (2.2) and $\mathscr Pf_{\rm sys}=0$; the Gaussian identity (2.4); **Lemma 2.1** (PSD and divergence-free lifted tensor, $n=2,3$, random superpositions) with a negative control; the holomorphic flux map (2.9)–(2.13). | ~1 min |
| `verify_quaternion_calculus.py` | **Lemma 3.1** (Hessian of rank-one forms, trace, product rule, $H(\sum\lvert F\rvert^2)=4J^*J$); **Lemma 3.2** (Moore determinant identities (3.6)–(3.8), concavity); **Lemma 3.3** (cofactor divergence, exact symbolic, 2/3/4 slots) with a negative control; **Lemma 4.3** segment identity and weights. | ~10 s |
| `verify_propagation_and_exponents.py` | **Lemma 6.1** radial comparison (Gronwall constant, $h$ down to $10^{-6}$); **Remark 6.2** sharpness example (exact Laplacian); exponent bookkeeping for $p_c$, (7.3), Prop 8.1, Lemma 9.4, Kakeya log powers, Stein threshold, and the Lean doc's power 66. | ~1 s |

The Lean build commands and logs are documented in `../06_LEAN_AUDIT.md`.

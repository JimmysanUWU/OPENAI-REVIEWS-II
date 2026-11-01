# ETH extension checks

Python 3.10+, standard library only. Run from any directory. Each script raises on failure, prints `ALL OK` on success, and writes its result under `../results/`.

- `verify_spline_bound.py`: exact spline and trace dimensions; deleted constraints; explicit character completion; 42 feasible designs including 30 genuinely ungraded; exact product-measure coupling; unequal jet orders; three hypothesis-removal controls. Also a finite numerical check of the separately proved binomial-tail estimate.
- `verify_complete_copies.py`: support injection with repeated-variable incidences, transferred norms, exact spectral multiplicities, aggregate row counts, variance slack, and sparse-copy controls that violate the unsupported per-copy claim.
- `verify_defect_capacity.py`: arbitrary trace multiplication ranks; jointly feasible character traces; finite-field spline dimensions and biased defects; arithmetic negative controls; discrete Rodrigues orthogonality; exact energy minima versus an independent rational Gram solve; Parseval on nonsymmetric families; numerical checks of the proved logarithmic estimates.

`explore_directional_extremizers.py` is a separate exhaustive search for three tiny instances, with saved rational witnesses. It prints `Exploration complete`, is not picked up by `verify_all.py`, and makes no asymptotic sharpness claim.

The asymptotic theorems are proved in the accompanying notes; these are finite verification experiments.

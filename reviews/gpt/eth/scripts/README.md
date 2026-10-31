# ETH extension checks

Python 3.10+, standard library only. Run from any directory. Each script raises on failure, prints `ALL OK` on success, and writes its result under `../results/`.

- `verify_spline_bound.py`: exact spline and trace dimensions; deleted constraints; explicit character completion; 42 feasible designs including 30 genuinely ungraded; exact product-measure coupling; unequal jet orders; three hypothesis-removal controls. Also a finite numerical check of the separately proved binomial-tail estimate.
- `verify_complete_copies.py`: support injection with repeated-variable incidences, transferred norms, exact spectral multiplicities, aggregate row counts, variance slack, and sparse-copy controls that violate the unsupported per-copy claim.

The asymptotic theorems are proved in the accompanying notes; these are finite verification experiments.

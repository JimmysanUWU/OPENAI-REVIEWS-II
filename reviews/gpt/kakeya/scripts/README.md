# Hermitian fibre check

`verify_weyl_fibre.py` uses only Python 3.10+ standard-library rational arithmetic, including a small Gaussian-rational type. It checks the realified Weyl matrices, fibre parametrization, both change-of-variables determinants, the distinction between the bilinear flux and real displacement, the exact annular integral, and the nonzero source from a Gaussian state cutoff.

Run from any directory. Results are written to `../results/weyl_fibre.json`; success prints `ALL OK`. The divergence theorem follows analytically from the exact annular integral in the note.

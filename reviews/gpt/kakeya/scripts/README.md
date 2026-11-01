# Hermitian fibre checks

`verify_weyl_fibre.py` uses only Python 3.10+ standard-library rational arithmetic, including a small Gaussian-rational type. It checks the realified Weyl matrices, fibre parametrization, both change-of-variables determinants, the distinction between the bilinear flux and real displacement, the exact annular integral, and the nonzero source from a Gaussian state cutoff.

Run from any directory. Results are written to `../results/weyl_fibre.json`; success prints `ALL OK`. The divergence theorem follows analytically from the exact annular integral in the note.

`verify_weyl_singularity.py` checks the general coarea determinant, rational real-slice samples, radial substitution factors, and the coefficient of the transverse singularity. It imports the Gaussian-rational arithmetic above. The displayed profile enclosures use proved analytic bounds evaluated in floating point; they are not machine-certified intervals. Results go to `../results/weyl_singularity.json`. The proof of local nonintegrability is analytic.


- `verify_curvature_fibre.py`: exact quadratic curvature identities, Clifford relations, full ambient Jacobians and coarea cancellation, for real/complex/quaternionic/octonionic Hopf pencils and a non-Clifford Hermitian 3x3 pencil. Five false-shortcut controls. Reuses Gaussian rational arithmetic from `verify_weyl_fibre.py`; standard library only.

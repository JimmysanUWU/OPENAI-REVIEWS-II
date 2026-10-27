# Sources: the Kakeya / local smoothing paper

| Item | Details |
|---|---|
| Paper | *Local smoothing for the Euclidean wave equation*, OpenAI, 25 October 2026, 59 pages. PDF metadata subject: "Local smoothing, Fourier restriction and summation, and Kakeya sets". |
| PDF supplied by the user | `kakeyagpt.pdf`, 810,093 bytes, SHA-256 `9ca156c552f2d940528f7131716326ef0d67a5829e723e1d294d74560c8636ff`. Not committed here. |
| Public copy | Byte-identical to `openai/math` `preprints/Local-smoothing-for-the-Euclidean-wave-equation-October-25-2026/paper.pdf` at commit `7ef5c0cc5ca4b0e072b7a050fca09ab24bc6fa3a` (26 Oct 2026). TeX sources are in the same directory under `build/`. |
| Formalization scope | `openai/math` `lean/docs/079.md` (link supplied by the user). |
| Comparator targets | [`LocalSmoothing.challenge.lean`](LocalSmoothing.challenge.lean), copied from `openai/math` `lean/ComparatorChallenges/LocalSmoothing.lean`. Six theorems with `sorry`; the permitted axioms are `propext`, `Quot.sound` and `Classical.choice` (`LocalSmoothing.json`). |
| Solution module | `OAI.Analysis.Fourier.LocalSmoothing.Comparator` in `openai/math`. Its import closure is 2,261 OAI modules and 190,467 lines; the only external dependency is Mathlib. |
| Toolchain | `leanprover/lean4:v4.34.1`; Mathlib `d13f23b723b8a846827a245b89c10fc7d3f11612` (the same revision pinned in the earlier Family 003 project). |
| Companion | *Critical local smoothing for the three-dimensional wave equation* (OpenAI, 24 Sept 2026), also in `openai/math`; not reviewed separately. |

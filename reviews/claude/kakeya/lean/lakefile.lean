import Lake
open Lake DSL

package OAIClosure where
  leanOptions := #[⟨`autoImplicit, false⟩]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "d13f23b723b8a846827a245b89c10fc7d3f11612"

lean_lib OAI where
  globs := #[.submodules `OAI]

lean_lib ComparatorChallenges where
  globs := #[.one `ComparatorChallenges.LocalSmoothing]

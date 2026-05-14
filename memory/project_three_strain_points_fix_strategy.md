---
name: Fix strategy for the three remaining strain points
description: How to harden the three remaining open primitives (p(y) uniqueness, action derivation, matter coupling). All three collapse to one root: the microscopic SU-write process. Bayesian conjugate-prior, counted-scalar EFT, substance-as-matter ontology.
metadata:
  type: project
---

# Fix strategy for the three remaining strain points

After the F/R/f triad reading + D-dimensional sanity check (2026-05-13 evening), three strain points remained in the framework's primitive list:
1. p(y) = y^D(1-y) minimal-degree was heuristic, not uniqueness
2. f(Sigma)R action was selected by elimination (G72/G74/G76/G78), not derived
3. Matter coupling beyond strong-field sectors is open

## Recommended fix for #1 — Bayesian conjugate-prior uniqueness

**Reframe:** "minimal-degree polynomial" -> "unique Bayesian posterior on the closure rate."

Beta(D+1, 2) is the conjugate posterior of a uniform Beta(1,1) prior on a Bernoulli rate, updated on D successes + 1 failure. The framework's primitives match exactly: D direction-activations + 1 reserved open-face + uninformative starting state. Cox's theorem makes Bayesian inference the unique consistent update rule. Therefore Beta(D+1, 2) is forced.

**Next-step script: G121_beta_uniqueness_from_bayesian_update.py**
- State prior/likelihood explicitly
- Derive posterior = Beta(D+1, 2) closed-form
- Show alternative priors (Jeffreys, MaxEnt with E[ln y], E[ln(1-y)] fixed) converge to same density
- Result: uniqueness statement, not minimal-degree heuristic

## Recommended fix for #2 — Counted-scalar uniqueness

**Reframe:** "smallest covariant action that survived elimination" -> "unique low-energy EFT of a *counted* scalar coupled to gravity."

Key new primitive: **Sigma is counted, not continuously dynamical.** A counted quantity can't propagate freely; it's algebraically pinned by the count process. This single commitment forces the Lagrange-multiplier structure:
- lambda_1 enforces grad Sigma matches source-determined accumulation (count is sourced, not free)
- lambda_2 enforces Sigma carried along substance flow (count co-moves with substance)

The G70-G78 elimination chain wasn't ad-hoc selection — it was showing every alternative that propagates Sigma as a dynamical scalar fails (ghost/perturbative pathology). The counted-scalar route was *forced*, just labeled differently.

**Next-step script: G122_counted_scalar_action_derivation.py**
- State Sigma-as-count primitive (integer-valued / commutator-pinned)
- Derive most general diff-invariant action for non-propagating scalar coupled to gravity
- Show lambda_1, lambda_2 are structural enforcers (not free choices)
- Confirm f(Sigma)R as unique kinetic term

This makes the action derived given one new primitive (Sigma is counted), itself motivated by 1 SU = A_0 (Sigma literally counts SU writes).

## Recommended fix for #3 — Matter coupling (two-tier)

**Short-term (pragmatic):** Standard Model matter with metric rescaled by f(Sigma) · g_mu nu (minimal coupling). f(2/3) = 1 anchors all matter outside PS to GR exactly — all precision tests pass by construction. Inside PS, deviations are observable via:
- Modified accretion-disk emission spectra
- Spectral shifts in Hawking emission (greybody corrections)
- Tortoise-distance enhancement (already at G90)

Concrete: **G123_matter_minimal_coupling.py** computes matter Lagrangian on f(Sigma)-rescaled metric, checks observational constraints.

**Long-term (deep / TOE):** Matter IS localized A. The substance ontology says spacetime is substance with density A; particles are substance excitations. Each particle of mass m contributes Sigma_particle = m/m_Planck SU writes per Compton time. SM emerges from coarse-graining substance excitation modes.

This is the framework's TOE move and a multi-year program. Park as Open Problem #8 (matter from substance excitations).

## The unifying observation

All three problems collapse to one root: **the microscopic SU-write process at the Planck scale.** Given:
- A formal model of SU writes (Poisson? graph-theoretic? more exotic)
- A coarse-graining theorem mapping it to macroscopic Sigma field

then:
- #1 falls out as Bayesian limit law for write-count posterior
- #2 falls out as unique EFT for coarse-grained counted scalar
- #3 falls out because matter = localized substance excitations of same write process

**The framework's deepest open question is the microscopic SU-write model itself**, not any of the three strain points individually. That's where the long-term arrow points.

## Ordering for execution

Priority order (cheapest-to-deepest):
1. G121 Bayesian uniqueness (~1 hour, hardens #1)
2. G122 counted-scalar action derivation (~1 day, hardens #2)
3. G123 matter minimal-coupling check (~1 week, addresses pragmatic #3)
4. Microscopic SU-write model (months/years, the deep TOE move)

## See also
- [[project_five_primitive_architecture]] — six primitives that this work would reduce
- [[project_R_of_A_derivation]] — the corrected D=3 algebraic uniqueness result
- [[project_v04_priority_SU_equals_A0]] — 1 SU = A_0 commitment that motivates Sigma-as-count
- [[HANDOFF_2026_05_13_action_and_qnm]] — original G70-G78 elimination chain
- [[project_spacetime_substance_ontology]] — substance ontology that motivates matter-as-localized-A

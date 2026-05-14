---
name: Five-primitive architecture (post 2026-05-13 manifold-support reading)
description: After the F/R/f triad reading lands, the framework's load-bearing primitives shrink to a short visible list. F(y) is manifold support, R(A) is 3D-space response, f(A) is SU stiffness. Everything else falls out by computation from five (or six) primitives plus the metric ansatz.
metadata:
  type: project
---

# Five-primitive architecture

After Sean's 2026-05-13 reading consolidating F(y) as manifold-support survival, R(A) as 3D-space anisotropy response, and f(A) as SU-accumulation stiffness, the framework's load-bearing primitives collapse to a short visible list. Everything else is derived.

## The triad mapping (post-shift)

  F(y) = remaining traversable manifold support
  R(A) = spatial (3D) curvature anisotropy response of geometry to manifold closure
  f(A) = accumulated SU/action stiffness coupling

Each maps to one of the framework's ontological commitments:

  F(y) ↔ manifold (the substrate)
  R(A) ↔ space (the geometric response)
  f(A) ↔ SU (the accumulation/write structure)

## The primitives (six items, all visible)

1. **A** — escape/traversal-cap variable, A = R_s/r
2. **Σ = D·A** — SU-normalized shell count (D = 3, the integer landmarks Σ ∈ {1, 2, 3} are ISCO, PS, horizon)
3. **3+1 manifold support** — D = 3 spatial-support directions + 1 manifold-support condition. Forces Beta(D+1, 2) and A_0 = 1/(4πD)
4. **p(y) = y^3(1-y) shape** — minimal-degree closure density vanishing at y=0 (no activation) and y=1 (no open support). This is the SOFTEST primitive — "minimal" is a heuristic, not a uniqueness theorem.
5. **h(A) = 1 - A** — velocity/escape-cap commitment (V_4, independently load-bearing for clock dilation, GW velocity, escape)
6. **f(Σ)R action ansatz** — BD-form non-minimal coupling. Constrained shell-count form with λ₁, λ₂ chosen by elimination (G70-G78); not a derivation.

## Derived (everything else)

- F(y) = 1 - ∫₀^y p(u) du = 1 - 5y⁴ + 4y⁵   (survival function, not a smoothing polynomial)
- k(A) = h·F   (velocity gate × manifold-support gate; cubic horizon vanishing = 1+2)
- G^t_t - G^θ_θ = (A² y⁴ / 4M²)(45A² - 33A - 13)   (Einstein anisotropy, computed)
- R(A) = -(45A² - 33A - 13) / [30 A (1-A)]   (Einstein anisotropy in shell-count coords)
- d ln f/dA = R(A) · p(y)/F(y)   (closure hazard sourcing stiffness)
- Z(A), V(A)   (secondary fields from constrained-action route, G70 closed forms)

## What the new reading does to the framework

- **No locked numbers move.** Every prediction (4.977% axial, 13% polar/axial break, GW170817 3.6%, χ² wins, A_0 = 1/(12π)) survives unchanged.
- **Open Problem #2** (k(A) uniqueness, [[project_strong_field_departure_question]]) closes down to "minimal-degree p(y)" softness.
- **Open Problem #1** (Lagrangian for A) gets sharper structural reading: f(Σ)R is the SU accumulation's stiffness contribution.
- **Dark-energy character** (V(A) = β/(1-A), [[project_F3_cosmology_dark_energy_finding]]) becomes structurally inevitable: V diverges as F→0 (manifold-support ceiling).
- **Cubic horizon closure** k~(1-A)^3 splits as 1 (velocity) + 2 (manifold-support quadratic vanishing).
- **D = 3 dimensionality** unifies three previously-separate facts: spatial dimension in A_0 formula, the 3 in Σ = 3A, and the D+1 = 4 in Beta(4,2). All three sites encode the same primitive.

## Terminology shift (in flight)

- "Ledger" → "manifold-support" (Sean's commitment 2026-05-13; forward-only per [[HANDOFF_2026_05_13_kerr_wave_branch]])
- "Quintic Hermite F as smoothing polynomial" → "F as remaining manifold-support survival function"
- "Ledger-channel F(y)" → "final-shell support profile"
- README still uses old terminology in ~12 places (lines 43, 95, 376, 388, 490, 508, 745, 876-900, 926, 1036). Forward-facing rewrite needed.

## Strain points (honest)

1. **Minimal-density argument is soft.** p(y) = y³(1-y) is the minimum-degree polynomial vanishing at endpoints with multiplicities (3, 1), but "minimum-degree" is not derived from a deeper structural principle. Sharper uniqueness argument needed.
2. **f(Σ)R action ansatz remains primitive.** Selected by elimination of scalar-tensor (ghost), aether, multi-field, cuscuton, mimetic. Strong narrowing but not a derivation.
3. **Matter coupling open** — how does the constrained shell-count action couple to ordinary matter beyond the Schwarzschild/Kerr strong-field sectors?
4. **D-dim sanity check RAN 2026-05-13 evening (sympy-verified).** The dimensional reading "4 = D_spacetime, 3 = D_spatial" was WRONG. General-D structure is `[4(D-2) - 6(D-3)A]·F - 4(D-2) + D·A·(3A-2)·F'_y`. The clean (F-1)+AyF'_y form arises ONLY at D=3 by two simultaneous algebraic coincidences: 6(D-3) zeroes out AND 3A-2 coincides with y=DA-(D-1). This is a STRENGTHENING: D=3 is algebraically privileged by the framework's structure, not just dimensionally labeled. See [[project_R_of_A_derivation]] for the corrected derivation.

## The deepest clean statement of the framework (post-shift)

> Spacetime is a substance with local density A. Each quantum physical interaction resolves 1 SU = A_0 = 1/(4πD), and these accumulate as a SU-normalized shell count Σ = D·A. The final shell (Σ ∈ (D-1, D)) is where the last open manifold support is consumed. Its closure density is Beta(D+1, 2) — D directions for activation plus one open-face factor. The survival function F(y) is what remains. The velocity-cap commitment h = 1-A is independent. Their product k = h·F is the radial gate. Everything else — the Einstein anisotropy response R(A), the SU stiffness coupling f(A), the modified-gravity stress-energy, the dark-energy-character at cosmic A_0, the ringdown shift, the cubic horizon vanishing — falls out by computation.

Six primitives, one paragraph, everything else derived.

## See also

- [[project_R_of_A_derivation]] — first-principles derivation of R(A) from Einstein anisotropy
- [[project_strong_field_departure_question]] — the closing of Open Problem #2
- [[project_strong_field_commitment]] — h, k commitments
- [[project_v04_priority_SU_equals_A0]] — A_0 = 1/(4πD) derivation
- [[project_critical_followup_F3_cosmology]] — dark-energy-character
- [[HANDOFF_2026_05_13_unified_framework]] — the original Beta(D+1,2) commitment
- [[HANDOFF_2026_05_13_action_and_qnm]] — the action-form narrowing by elimination

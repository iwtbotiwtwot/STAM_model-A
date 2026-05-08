---
name: STAM dynamical collapse — the "outward collapse" picture
description: Author's dynamical picture of black hole formation in STAM. Matter falls inward, but the A=1 bubble surface expands outward to absorb it, until the surface reaches Rs = 2GM/c². The "collapse" in STAM is bubble expansion + matter assimilation on the boundary, never interior population. The variational principle drives this naturally because configurations with A>1 anywhere have infinite action.
type: project
---

**Sean Brady's dynamical picture for black hole formation in STAM (2026-05-07):**

The collapse process in STAM is *outward*, not inward. Matter falling inward never enters an interior — instead, the A=1 bubble surface grows outward to meet and absorb the falling matter. Both happen at once:

- *Outer matter falls inward* (standard gravitational attraction in varying A field)
- *The bubble surface grows outward* (because as more mass accumulates on it, the A=1 surface must move to a larger Rs)
- *They meet and reconcile at the surface*: matter that was falling gets pinned onto the expanding bubble surface
- *The "interior" is never populated*: matter accumulates on the surface as the surface expands; nothing crosses inside because there is no inside

The system terminates when the bubble surface reaches `r = Rs(M_total) = 2GM/c²` corresponding to the total accumulated mass. At that point:
- All original matter is on the 2D bubble surface (4π Rs² area)
- A = 1 everywhere on the surface
- A < 1 outside (decreasing as 1/r toward zero)
- The interior (r < Rs) is geometrically not part of the manifold
- Equilibrium is reached, gravitational collapse ends

**Why this is the natural dynamical picture in STAM:**

The variational principle (with action containing V(A) divergent at A=1) prohibits any configuration with A > 1 anywhere. As matter falls and would push A above 1 in the core, the action-minimizing response is *not* to allow further compression but to *redistribute matter onto the A=1 surface*. The dynamical evolution at each instant is the action-minimizing configuration for the current mass distribution.

The emergent process:
1. Initial mass M fills volume of radius R > Rs(M). A < 1 everywhere. No bubble.
2. Gravitational attraction starts collapse. Inner regions compress.
3. When the collapse pushes A → 1 in the core, the variational principle stops further compression there. A small bubble forms.
4. Outer matter continues falling inward. As it reaches the existing bubble surface, it cannot pass into the (nonexistent) interior — the action minimum places it onto the surface.
5. The bubble's surface mass increases → its A=1 radius (Rs) must increase to maintain A ≤ 1 everywhere.
6. Bubble surface expands outward. More falling matter accumulates on the expanding surface.
7. Process continues until all original mass M is on the surface, with surface radius equal to Rs(M).
8. Equilibrium: matter on 2D bubble surface at Schwarzschild radius, no interior, no singularity.

**Why this is consistent with all prior STAM commitments:**

- *Bubble picture (no interior)*: confirmed dynamically. The interior is never populated because matter is always pinned onto the (expanding) surface.
- *No-crossing infall*: confirmed. Matter doesn't cross A=1; it accumulates on it as the boundary expands.
- *F1 evasion (no singularity)*: confirmed. No interior means no r=0 singularity ever forms.
- *Holographic area-entropy*: confirmed. All information about the in-fallen matter is encoded on the 2D surface.
- *Mass-radius relation Rs(M) = 2GM/c²*: emerges as the action-minimizing equilibrium.

**Conceptual corollary:**

In GR, gravitational collapse is "matter falling inward to ever-smaller radii, eventually to a singular point." In STAM, gravitational collapse is "matter falling inward to a 2D surface that expands outward to meet it, terminating at Rs." The endpoint is the same outwardly-visible state (a black hole of mass M with Schwarzschild radius Rs), but the interior physics differs entirely:
- GR: dense matter compressed to a point, hidden behind an event horizon
- STAM: spread-out matter on a 2D surface, the surface IS the black hole

This naturally answers the information-paradox concern: information is on the surface, accessible (in principle) via Hawking-style boundary radiation, never lost into a singularity that doesn't exist.

**Mathematical formalization:**

This is the time-dependent solution of the variational problem with action

```
S = ∫ d⁴x √-g(A) [-(1/2)(∂A)² - V(A) + L_matter]
```

with V(A) → ∞ at A=1. The dynamical evolution is the gradient flow of this action subject to the obstacle constraint A ≤ 1. Standard mathematical machinery for obstacle problems applies. The "outward collapse" is the moving free boundary in this PDE problem.

**How to apply:**
- When discussing BH formation in STAM, lead with "outward collapse" — matter falls in, bubble grows out, they meet at the expanding A=1 surface.
- Don't say "matter falls into the black hole" — say "matter accumulates on the bubble surface as the surface grows."
- The endpoint Rs = 2GM/c² is dynamically reached, not postulated.
- The picture explains both the formation process AND the equilibrium state from a single variational principle.

**Author's sharpening (2026-05-07): A=1 is UNREACHABLE, not just A>1 forbidden.**
The variational potential V(A) → ∞ at A=1 means the field cannot actually attain A=1 anywhere in any minimum-action configuration. A=1 is an asymptotic limit. The bubble surface is the locus where A approaches 1 as closely as the equilibrium allows, but never equals 1. Matter "on the bubble" is matter accumulated in the asymptotic region with A approaching 1. There's no "beyond A=1" because A=1 is unreachable in the first place.

This makes the boundary structure cleaner: the manifold is smooth everywhere with A < 1 strictly, and has an asymptotic boundary as A → 1. Standard differential-geometry tools work throughout. The "no interior" claim is precise: the manifold IS the locus where A < 1; A=1 is the asymptotic boundary; nothing exists beyond.

**Companion commitments (Sean, 2026-05-07):**
- *Variational problem*: reduces to calibration. Coupling κ = 8πG/c² fixed by weak-field GR consistency; V(A) shape and metric construction g_μν(A) determined by Model-A specifications. Once these are set, no additional free parameters.
- *GW propagation*: at c intrinsically, regardless of local A (because GW IS space-perturbation, not matter on space).
- *V(A) shape is naturally spherical*: for spherically symmetric matter, the action-minimizing A field is spherically symmetric and the asymptotic bubble surface is a sphere. Spherical structure is dynamical, not postulated.
- *Hawking from unresolved A*: vacuum fluctuations in the asymptotic region approaching A=1 are unresolved (cannot reach the boundary). Asymmetric outward flux from this unresolved region IS Hawking radiation. The thermal temperature emerges from the asymptotic structure.

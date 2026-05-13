# G64 - Bubble Surface Refinement (Two Candidate A Formulas)

**Date: 2026-05-13.** Explores two candidate A(r, theta; M, a) formulas for the framework's Kerr extension, resolving the high-spin photon-orbit and ISCO issues flagged in G63.

## The structural question

From Sean (2026-05-13): 'Can we refine the spinning bubble surface? Does the hologram spinning on the stationary horizon work here?'

The G3/G62 formula `A = 2Mr/Sigma` has `A = 1` at the **static-limit surface**, not at the actual Kerr horizon. This is why the bubble's equator stays at `r = 2M` for any spin, which buries the Kerr photon orbit (a > 0.707) and ISCO (a > 0.943) inside the bubble at high spin.

Alternative: `A = 2Mr/(r^2 + a^2)` places `A = 1` at the actual Kerr horizon `r = r_+` (constant in BL coordinates, oblate when embedded). The bubble is then the Kerr horizon, hologram matter spins on it, all photon orbits and ISCOs sit outside cleanly for any spin.

## Comparison

| Property | Formula alpha (A = 2Mr/Sigma) | Formula beta (A = 2Mr/(r^2+a^2)) |

|---|---|---|

| A=1 surface | static-limit (oblate r(theta)) | Kerr horizon r = r_+ |

| Equatorial extent | r = 2M (all spin) | r = r_+ (varies) |

| Polar extent | r = r_+ | r = r_+ |

| T on bubble | banded T(theta) | uniform = T_Kerr |

| T_eq | T_Schwarzschild | T_Kerr |

| T_pole | T_Kerr | T_Kerr |

| Kerr PS location | inside bubble for a > 0.707 | always outside |

| Kerr ISCO location | inside bubble for a > 0.943 | always outside |

| g_tt structure | g_tt = -(1 - A_alpha) c^2 = exact Kerr g_tt | direct identity broken |

| LIGO ringdown high spin | departs from GR-Kerr | exact GR-Kerr |

| High-spin disk inner edge | r = 2M (prediction) | matches Kerr ISCO |

| Latitudinal evaporation banding | distinguishing prediction | not predicted |

| Hologram-on-stationary-horizon | bubble oblate, matter rotates | clean match |

## What each formula gives up

**Formula alpha gives up:** consistency with standard Kerr at high spin. Predicts no ergosphere physics (region between r_+ and r=2M doesn't exist in manifold), high-spin BH inner edges at r=2M not at Kerr ISCO, photon orbits buried.

**Formula beta gives up:** the latitudinal banding distinguishing prediction. T is uniform = T_Kerr on the bubble, matching standard Kerr. Also gives up the direct g_tt = -(1-A) c^2 structural identity (g_tt no longer equals -(1-A_beta) c^2 in Kerr coordinates).

## Sean's 'hologram on stationary horizon' reading fits beta

Under formula (beta):

- The bubble IS the Kerr horizon: constant r = r_+ in BL coordinates, oblate when embedded in flat 3-space.

- The hologram matter on the bubble spins at ZAMO frequency, carrying J.

- T is uniform on the bubble (= standard Kerr T_Kerr).

- Rotation tilt at the horizon still makes outward emission more viable (Sean's 2026-05-13 reading): the super-radiance enhancement R(theta) = 1 + (v_matter/c)^2 still applies to the emission RATE, with equatorial matter velocity giving rate enhancement -- but T itself is uniform.

- The framework's 'rotation enhances outward emission' is structurally preserved as a RATE effect, not a temperature effect.

## What the framework retains under beta

- Substance ontology + presentism (unchanged)

- Ledger-channel + two-face + SU shell-count (unchanged)

- alpha = 4 entropy area-per-entry (unchanged)

- Quintic Hermite F(y) for inside-PS metric (unchanged)

- Spinless strong-field metric (unchanged)

- LIGO spinless ringdown = exact GR (unchanged)

- LIGO Kerr ringdown = exact GR-Kerr (recovered for all spin, not just <0.707)

- Static bubble + rotating matter (clean: bubble at Kerr horizon)

- Outward-only emission with rotation enhancement of rate

- D = 3 forcing (unchanged)

## What changes under beta

- T uniform on bubble = T_Kerr (replaces T_pole = T_Kerr + T_eq = T_Schw)

- Latitudinal banding of evaporation T NOT predicted (replaces G63's banding prediction)

- Direct g_tt = -(1-A) c^2 identity for Kerr no longer holds (alpha had it)

- Standard Kerr observables matched cleanly at all spins



The framework still predicts:

- Equatorially enhanced emission RATE (super-radiance from matter rotation)

- Different ringdown physics inside Kerr photon orbit (quintic Hermite F)

- No interior, ledger structure, primordial-mass remnants

## Recommendation

Formula (beta) is the cleaner reading per Sean's 'hologram on stationary horizon' framing. It:

- Resolves the high-spin photon-orbit/ISCO conflict in G63

- Matches standard Kerr at all spins for the macroscopic observables

- Preserves all framework primitives except the static-limit identification

- Gives up the latitudinal banding observable but keeps the rotation-rate enhancement as a distinguishing observable



Formula (alpha) gives a more aggressive distinguishing prediction set but creates the high-spin photon-orbit conflict.



Sean's call which to commit to. The framework can also stay agnostic and treat the choice as observational: if high-spin BH inner edges are observed at the spin-dependent Kerr ISCO (consistent with continuum-fitting results), (beta) is preferred; if they're observed at r = 2M, (alpha) is preferred. Currently, continuum-fitting favors (beta).

## Files

- [scripts/G64_bubble_refinement.py](../scripts/G64_bubble_refinement.py)

- [plots/G64_bubble_shapes_compared.png](../plots/G64_bubble_shapes_compared.png)

- [plots/G64_T_compared.png](../plots/G64_T_compared.png)

- [plots/G64_PS_vs_bubbles.png](../plots/G64_PS_vs_bubbles.png)

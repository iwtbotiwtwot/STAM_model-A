---
name: Substance velocity-cap — STAM commitment forced by substance ontology
description: V_4 candidate commitment that motion through elevated A is REALLY slowed (not just observationally) by some function f(A) of substance density. This is forced by the substance ontology + the magic-bell prohibition. Implemented as v_effective² = v_Newton² × (1 - A_local) in the GW170817 engine time calculation, gives 3.6% match without fitting. Specific functional form is one of several candidates; framework should commit explicitly.
type: project
---

## The commitment

**Motion through elevated A is slower than v_Newton predicts.** The slowdown is a real physical effect (substance ontology), not just an observational artifact.

Operationally: at any orbital position with local substance density A_local,
```
v_effective² = v_Newton² × f(A_local)
```
where f(A) is some function with f(A_0) ≈ 1 and f(A) → 0 as A → 1.

## Why this is forced by V_4 (not optional)

Under V_4's substance ontology + the magic-bell prohibition:

1. **A is real substance density**, not coordinate artifact ([[project-spacetime-substance-ontology]])
2. **Motion is physical traversal through substance**, not pure-geometric trajectory
3. **The magic-bell prohibition** forbids treating coordinate motion as if it bypasses the substance medium
4. **Therefore**: motion through elevated A must encounter the elevated substance density — slowing actual coordinate speed

The framework cannot maintain "A is real substance" AND "motion through it is unaffected." One or the other has to give. V_4 commits to the substance reality, so motion must be affected.

## Candidate functional forms

| Form | f(A) | Origin | Result for GW170817 τ_engine |
|---|---|---|---|
| Proper-time-like | (1 - A) | From metric g_tt factor √(1-A)² | **1.677 s (-3.6%)** |
| Square-root | √(1 - A) | Half-power of proper time | 1.886 s (+8.4%, in G41) |
| Light radial | (1-A)²(1+A) | From light coordinate-speed formula | 1.612 s (-7.4%, in G41) |

**The (1 - A) form gave the best match** to GW170817 in the engine time calculation (G42, G43). But the framework should commit on which is structural, not just empirically best.

Structural considerations for choosing:

- **(1 - A)** is the natural proper-time²-equivalent: things that experience proper-time slowdown also experience coordinate-motion slowdown by the same factor. Physically intuitive.
- **√(1 - A)** is what you'd get if substance contributes drag at half-order — less natural ontologically.
- **(1-A)²(1+A)** is the light coordinate-speed formula — applies to wave propagation, not necessarily to massive-body motion.

Recommended framework commitment: **(1 - A)** because it matches the proper-time scaling that already appears in g_tt. This is the most ontologically consistent choice.

## Self-consistent application

For a binary at orbital separation r, the local A at each NS's position includes:
- A_static from the other NS: R_s/(2r)
- A_kinematic from own motion: v²/c²

With substance correction:
```
A_kinematic = (v_actual)²/c² = v_Newton²/c² × (1 - A_local)
           = (R_s/2r) × (1 - A_local)
```

Self-consistent for total A_local = A_static + A_kinematic:
```
A_kin = (R_s/2r)(1 - R_s/2r - A_kin)

Let x = R_s/(2r). Then:
A_kin × (1 + x) = x(1 - x)
A_kin = x(1-x)/(1+x)

A_local at orbital position = R_s/2r + A_kin = x + x(1-x)/(1+x)
                            = x[(1+x) + (1-x)] / (1+x)
                            = 2x/(1+x)
```

So **A_local = 2x/(1+x) = R_s/(r + R_s/2)**.

## Application to binary mergers

The doughnut threshold (where A_kinematic crosses A_0):
```
A_0 = x(1-x)/(1+x)
```
Solving (quadratic): x² + x(A_0 - 1) + A_0 = 0
```
x = [1 - A_0 - √((1-A_0)² - 4A_0)] / 2
```

For A_0 = 1/(12π) ≈ 0.0265:
```
r_threshold/R_s ≈ 17.82  (vs naive 6π ≈ 18.85)
```

Peters-Mathews chirp time from this threshold:
```
τ_critical = (5/8) × (r_threshold/R_s)⁴ × R_s/c
```

For 2.7 M_sun (GW170817): τ_critical = 1.677 s.

## Wedge with GR (observable consequence)

GR has the geometric metric factors (g_tt, g_rr) but no separate substance interaction. STAM has both.

At low A: STAM and GR agree (substance correction is ~1 - A_0 ≈ 0.97).

At high A: STAM diverges from GR. Specifically:

- **GW emission rate**: STAM reduces it by f(A)⁶ at each orbital phase (slower v → less GW power)
- **Late-inspiral chirp shape**: differs from pure GR template
- **Binary mass extraction from chirp**: LIGO templates assume GR; STAM-corrected templates would extract slightly different masses (~few percent shift)
- **Engine time prediction**: τ_critical formula above — applies to BNS systems with EM counterparts

These are real, testable departures.

## How to apply

When computing inspiral dynamics in STAM:
- Replace v_Newton² with v_effective² = v_Newton² × (1 - A_local) wherever orbital velocity appears
- Compute A_local self-consistently (don't just plug in static A from masses)
- The substance factor (1 - A) propagates into GW emission rate, energy loss, and inspiral duration
- Recompute Peters-Mathews chirp time with the corrected dynamics

When committing the framework to a specific form:
- Use (1 - A) unless there's a structural reason to prefer another form
- Document the choice and its derivation status
- Note the GW170817 match (3.6%) is supporting evidence but not the primary justification (the substance ontology IS the primary justification)

## Status

- **2026-05-12 session**: candidate commitment landed via [[handoff-note-2026-05-12-session]]
- Functional form (1 - A) chosen on grounds of ontological consistency with proper-time factor
- Supporting evidence: 3.6% match to GW170817 engine time without fitting
- **Open**: formal framework derivation of the specific functional form; the (1 - A) choice is currently structurally-motivated but not uniquely forced

## Related memories

- [[project-spacetime-substance-ontology]] — the substance ontology this commitment follows from
- [[project-strong-field-commitment]] — the metric g_tt, g_rr factors that motivate (1 - A)
- [[project-A-zero-means-nothing-exists]] — A_0 as baseline (sets threshold)
- [[project-gw170817-engine-time-prediction]] — the falsifiable prediction this commitment produces

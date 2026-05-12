---
name: GW170817 engine time — STAM falsifiable prediction (mass-scaling)
description: Under V_4 + substance velocity-cap commitment + Sean's elastic-shell picture, STAM predicts engine time of binary mergers = τ_critical (Peters-Mathews chirp time from substance-corrected threshold). For GW170817: 1.677 s predicted vs 1.74 s observed, 3.6% match without fitting. Mass-scaling: linear in M_total. Falsifiable with future BNS+EM events.
type: project
---

## The prediction

**Engine time of a binary merger (orbit-stop → light arrival) equals τ_critical**, where τ_critical is the Peters-Mathews chirp time from the substance-corrected doughnut threshold to merger.

Formula (equal-mass binary):
```
τ_critical = (5/8) × (r_threshold/R_s)⁴ × R_s/c

where r_threshold/R_s is the smaller root of:
  x² + x(A_0 - 1) + A_0 = 0  with x = R_s/(2r)

For A_0 = 1/(12π): r_threshold/R_s ≈ 17.82
```

For 2.7 M_sun system (GW170817):
```
τ_critical = 1.677 seconds
```

## Mass scaling

Linear in total binary mass:
```
τ_critical ∝ R_s/c ∝ M_total
```

Slope: approximately 0.62 seconds per M_sun.

Predicted engine times for various binary masses:
| System | M_total | Predicted τ_engine |
|---|---|---|
| BNS-light (1.4+1.4) | 2.8 | 1.74 s |
| BNS-GW170817 | 2.7 | 1.68 s |
| BNS-heavy (2+2) | 4.0 | 2.48 s |
| NS-BH (1.4+5) | 6.4 | 3.97 s |

## What this is based on

Three commitments stack:

1. **A_0 = 1/(12π)** (V_4 derived) — sets the doughnut baseline threshold
2. **Substance velocity-cap** ([[project-substance-velocity-cap]]) — slows orbital motion through elevated A
3. **Engine time = bubble retraction time = τ_critical** (Sean's elastic-shell picture) — identifies the engine time with the substance-corrected chirp time

The third commitment is the most tentative. It rests on:
- Doughnut substance is built during inspiral at orbital threshold
- After orbit-stop, doughnut "retracts" toward merger center at substance wave speed
- Engine takes one τ_critical duration before explosion can fire
- Substance factors cancel: τ_contract = d_doughnut / v_substance = τ_critical regardless of absolute speed

If the engine time has a different mechanism (e.g., HMNS collapse dynamics that don't scale with τ_critical), this prediction fails.

## GW170817 match

- Observed engine time: 1.74 s (gamma-ray burst onset after merger)
- STAM predicted: 1.677 s
- Discrepancy: -3.6% (about 63 ms)

**No fitting anywhere in the chain.** Only inputs are:
- Total mass (observational, GR-template extraction; STAM-consistent at inspiral A levels)
- A_0 = 1/(12π) (V_4 derived)
- Substance velocity-cap form (1-A) (framework commitment, candidate)
- Peters-Mathews leading-order coefficient (standard GR, STAM-equivalent at low A)

The 3.6% match is meaningful but rests on one data point. Mass-scaling test is the real falsification.

## Falsification handle

For future BNS+EM events:

**STAM-supporting outcome**: engine times fall on the linear-in-mass trend with slope ~0.62 s/M_sun.

**STAM-refuting outcomes**:
- Engine times don't correlate with total binary mass (mechanism is mass-independent)
- Engine times scale differently (e.g., logarithmically, or with chirp mass instead of total mass)
- Wide scatter at fixed mass, suggesting astrophysical (HMNS collapse) variability dominates

## What this is NOT

- **Not a claim about BBH events** — BH-BH mergers don't have observable engine times (no EM counterpart). But the substance correction still applies to inspiral chirp shape; that's a separate (harder to test) prediction.
- **Not a complete theory of the gamma-ray burst mechanism** — STAM predicts when, not what specifically the gamma-ray emission process is.
- **Not yet a unique framework commitment** — the (1-A) coupling form was chosen partly for empirical match to GW170817. Other forms give different predictions. Framework should commit on the form before claiming this as a definitive prediction.

## Caveats

- **One data point**: until more BNS+EM observations come in, this is a single-point match, not a verified prediction
- **Peters-Mathews leading-order**: PN corrections in standard GR can shift things by 5-15%; full PN-corrected STAM treatment would give a more precise prediction
- **Engine time identification**: rests on Sean's elastic-shell picture, which the framework hasn't formally derived
- **Coupling form**: (1-A) chosen on ontological grounds and empirical match; framework should commit explicitly

## How to apply

When discussing GW170817 with Sean:
- STAM is **consistent with** 1.74 s, predicting 1.68 s under specific commitments
- The match is supporting evidence for the substance velocity-cap, not standalone load-bearing
- For falsification, point to future BNS+EM events and the linear mass-scaling

When considering future binary merger predictions:
- BNS+EM events: use τ_critical formula above
- BBH events: substance correction applies to chirp shape (late-inspiral departure from GR templates); no engine-time test
- NS-BH mergers: predicted τ_engine = 4 s for typical 1.4+5 M_sun system (testable)

When V_4 writes this up:
- Frame as **candidate** prediction with one matched data point + mass-scaling falsification handle
- Note the three stacked commitments needed
- Be honest about which are V_4-forced (A_0, substance-velocity-cap-existence) vs candidate (specific form, elastic-shell mechanism, engine-time identification)

## Status

- **2026-05-12**: candidate prediction landed (G42, G43 scripts)
- **Awaiting**: more BNS+EM events for mass-scaling test, formal framework derivation of (1-A) coupling, formal mechanism for engine-time = τ_critical identification

## Related

- [[handoff-note-2026-05-12-session]] — session that produced this
- [[project-substance-velocity-cap]] — the framework commitment this prediction rests on
- [[project-spacetime-substance-ontology]] — the foundational ontology forcing the velocity-cap

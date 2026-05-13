# LIGO arc landing — GW170817 engine time mechanism

**Session date: 2026-05-12 (continued from earlier)**

## Final mechanism (Sean's commitment)

For a binary neutron star merger that creates a black hole:

1. **Inspiral builds substance disturbance** at v_substance(A_0) extending to r_donut_max ≈ v_substance × τ_chirp.
2. **Merger creates a BH at R_s** — A = 1 saturation surface forms. This is structurally committed (confirmed observationally for GW170817 via Chandra X-ray analysis).
3. **The new BH's near-saturation atmosphere** extends from R_s out to where A drops to A_0 (at r ≈ R_s/A_0 = 12π R_s ≈ 300 km for GW170817).
4. **GW from inspiral** is far outside the event, propagating outward at v_substance(A_0) baseline — does *not* need to traverse the BH atmosphere.
5. **GRB photons** originate from the central engine just outside R_s and must crawl out through the near-saturation atmosphere. v_substance(A near 1) → 0, so light is extremely slow in this shell.
6. **GW and light originate essentially co-located** (light "nanometers behind" GW) — but light's traversal of the near-saturation region produces the engine time delay.

## Numerical result (G50 targeting)

For GW170817 (M_total = 2.7 M_sun, R_s ≈ 7.98 km):

- Targeted ε = (1 − A_emit) such that engine time = 1.74 s
- **ε ≈ 1.56 × 10⁻⁵**
- Light's emergence point: r_emit ≈ R_s × (1 + ε) ≈ R_s + 12 cm
- A at r_emit ≈ 0.99998 (essentially saturation)
- v_light at r_emit ≈ 0.15 m/s (walking pace)

Structural identity:
```
τ_engine ≈ R_s / (c × ε)
```

## Mass-scaling — two readings, single-event-indistinguishable

**Reading A: ε universal across BH-forming BNS events**
- τ_engine ∝ R_s ∝ M_total
- Slope ~0.64 s/M_sun
- Predicts engine time varies with binary mass

**Reading B: τ_engine constant ≈ 1.7 s for BH-forming events**
- Requires ε ∝ R_s ∝ M (shell thickness scaling with R_s²)
- Predicts engine time independent of mass
- Possibly more consistent with engine time being set by universal post-BH-formation physics rather than the specific R_s

Cannot distinguish from GW170817 alone. **Falsification handle**: future BNS+EM events at different M_total.

## What STAM contributes structurally

- **Mechanism**: GW born outside event (or co-located but GW skips the near-A=1 region), light crawls out of near-saturation atmosphere
- **A_event = 1**: structural commitment for any BH-forming merger
- **Donut size set by BH structure**: R_s × 12π (BH atmosphere extent) — depends only on R_s and A_0
- **Direction of delay** (GW first, by mechanism, not by fitting)
- **Distinction between BH-forming and non-BH-forming BNS events**: predicts qualitatively different engine-time behavior (the framework breaks at the BH-formation mass threshold)

## What STAM doesn't claim

- The exact 1.74 s number — this is a single-data-point match with structural ambiguity in mass-scaling
- That conventional astrophysics (jet breakout, ejecta opacity) doesn't contribute — these may set the magnitude
- A derivation of ε from first principles — currently ε ≈ R_s/(c × τ_engine) is structural relationship, not derived

## Observational data status

- **GW170817 is the only confirmed BNS+EM event** with a measurable GW→GRB delay
- GW190425 was a BNS candidate with no EM counterpart
- Other LIGO events: BBH (no expected EM) or NS-BH candidates without EM
- Future events from O4 extended / O5 should provide more BNS+EM data points
- 1-10 additional events expected in the next 5-10 years

## Scripts produced this arc

- [G42](../scripts/G42_doughnut_threshold_corrected.py) — initial substance threshold
- [G43](../scripts/G43_doughnut_contraction_full_stam.py) — earlier mechanism (substance velocity-cap on orbital motion)
- [G44](../scripts/G44_stam_native_engine_time.py) — V_3 bubble capacity check + STAM-GW emission
- [G45](../scripts/G45_stam_native_local_speeds_and_delay.py) — local speeds and donut size
- [G46](../scripts/G46_doughnut_feedback.py) — second-order doughnut feedback sweep
- [G47](../scripts/G47_donut_radius_meeting_light_transit.py) — meeting + light transit
- [G48](../scripts/G48_stam_sandbox_engine_time.py) — halfway condition, τ_chirp = engine time
- [G49](../scripts/G49_no_double_dip.py) — no-double-dip velocity-cap framing
- [G50](../scripts/G50_target_174.py) — targeted ε for exact 1.74 s match

## Landing point

The framework provides a coherent structural mechanism for the GW→GRB engine time in BH-forming BNS mergers, consistent with GW170817 to within event-specific astrophysics uncertainty. The mechanism doesn't require asymmetric A coupling (F5c stays falsified). The specific magnitude depends on a structural parameter ε whose universal vs event-dependent character will be testable with future BNS+EM events.

This is NOT a "STAM predicts 1.74 s to X%" claim. It's a structural mechanism consistent with observation, awaiting more data for genuine falsification.

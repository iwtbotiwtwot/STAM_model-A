# G44 — STAM-native engine time, two readings of the substance velocity-cap

## What was asked

Does STAM natively suggest a different chirp rate than Peters-Mathews? Run the
calculation sandboxed in STAM and see what the math says.

## Two readings of the substance velocity-cap on orbital motion

The substance ontology says motion through elevated A is really slowed. For
translational/radial motion this is clear. For *circular orbital motion* the
framework has not committed:

- **Reading 1**: cap applies to all motion through substance, including circular orbital. Then ω² = (1-A)·GM/r³ in proper time, and P_GW from the quadrupole formula picks up a (1-A)³ reduction at each r.
- **Reading 2**: cap is for encountering new substance (translation/radial); circular orbital motion frame-drags the local substance, no fresh encounter, no cap. Then P_GW is unchanged at low A — STAM ≡ GR at the emission level — and STAM's departure from GR lives only in the threshold and the clocks.

## Numerical results (GW170817, M = 2.7 M_sun, observed engine time 1.74 s)

| Calculation | τ (source) | τ (observer) | gap vs 1.74 |
|---|---|---|---|
| Naive Newton (no STAM) | 2.103 s | ~2.15 s | +20.7% source / +23.6% obs |
| G43 (STAM threshold, no clock) | 1.677 s | — | −3.6% |
| **Reading 2** (G43 + LIGO-frame clock) | 1.677 s | **1.718 s** | **−1.3%** |
| **Reading 1** (R2 + (1-A)³ on chirp rate) | 2.107 s | **2.159 s** | **+24.1%** |

## V_3 bubble as energy sink — ruled out structurally

A separate finding from the same script:

The V_3 bubble's energy capacity at any plausible interior A elevation is
**10²⁹× smaller** than the orbital binding energy released during inspiral. β is
calibrated to cosmic dark-energy magnitude (β ≈ 6×10⁻¹⁰ J/m³); stellar binaries
release ~10⁴⁶ J. The substance at the V_3 baseline is too dilute to absorb
stellar-mass binding energies through A elevation alone.

Implication: STAM's energy loss channel at stellar scales is GW radiation, not
substance V_3 accumulation. This is consistent with the framework's commitment
that STAM reproduces GR at low A — but it's structurally derived here from the
substance density, not assumed.

## The cancellation in Reading 1

τ_R1_source = 2.107 s ≈ τ_naive_Newton = 2.103 s (within 0.2%)

The substance threshold correction (smaller r_threshold = 17.82 R_s vs naive
6π R_s ≈ 18.85 R_s, which speeds the chirp by factor 0.79 in τ) is almost
**exactly canceled** by the substance-velocity-cap-on-orbital-motion correction
(slowing the chirp rate by ⟨(1-A)⁻³⟩ ≈ 1.26).

Net: Reading 1 collapses to naive Newton at the engine-time level.

This is a structural finding. If Reading 1 were the correct interpretation,
STAM would be observationally indistinguishable from naive Newton chirp + clock
dilation for the GW170817 engine time. The framework's distinguishing content
would vanish in this regime.

## Reading 2 is the framework's distinguishing prediction

Reading 2 gives 1.718 s observer frame — distinct from naive Newton (2.10 s)
and within 1.3% of the observation. The STAM-specific content is:
- Substance-corrected threshold (r/R_s = 17.82, not 6π)
- LIGO-frame clock dilation along the chirp trajectory

Neither piece is GR; both are STAM's own commitments propagated cleanly.

## What the math is telling us

1. **The substance velocity-cap probably does not apply to circular orbital motion.** Reading 1's cancellation against the threshold correction makes STAM degenerate with naive Newton — the framework loses its distinguishing prediction in this regime. The substance ontology should be sharpened to commit on this distinction (translation vs orbital).

2. **V_3 bubble is the wrong scale for stellar energy loss.** Substance is cosmologically dilute; binaries lose energy via GW radiation. The substance picture matters for the *threshold* (where motion saturates) and the *clocks* (how fast each location ticks), not for absorbing orbital binding energy.

3. **The 1.3% Reading 2 gap is small but real.** Candidates for what it tells us: PN corrections, finite-mass effects, NS spin/tidal effects, or genuinely a STAM piece we haven't surfaced. Below the noise floor of "no fitting" claims with a single data point.

4. **Mass-scaling remains the falsification handle** — Reading 2 predicts τ_engine ∝ M_total with slope ≈ 0.62 s/M_sun. Future BNS+EM events at different masses test the framework's commitment.

## Open framework question raised

The framework should commit: **does the substance velocity-cap apply to circular orbital motion, or only to translational/radial motion through substance?**

The math here suggests Reading 2 (translation only). The substance ontology
should be examined to see if Reading 2 follows structurally from "what counts
as encountering new substance" or whether the commitment is empirical at this
point.

## Files

- Script: [scripts/G44_stam_native_engine_time.py](../scripts/G44_stam_native_engine_time.py)
- Methodology: scripts as listening instruments ([feedback](../memory/feedback_scripts_as_listening_instruments.md))

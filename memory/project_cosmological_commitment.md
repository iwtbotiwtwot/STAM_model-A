---
name: STAM cosmological commitment — two-layer reading (V_3 expansion + photon-A distance bias)
description: Sean's Model-A cosmology, originally committed 2026-05-07 and substantively refined 2026-05-11 evening after G28 + G29. Two-layer structure: (Layer 1) V_3 modified Friedmann gives LCDM-equivalent cosmic expansion at H_0=73 by construction; (Layer 2) photon-A traversal through cosmic A_0 in voids adds path-integral bias to observed luminosity distance. Chronometers probe Layer 1 only; SN/CMB distance moduli probe Layer 1 + Layer 2. Bridge term lives in Layer 2.
type: project
---

**Sean Brady's Model-A cosmological commitment.**

Originally committed 2026-05-07 as "matter-only EdS + photon-A traversal misinterpretation." Refined 2026-05-11 evening after G28 + G29 clarified the underlying structure. The refined position is structurally cleaner; the photon-A core remains, but the "matter-only EdS" framing is retired and replaced with a two-layer reading.

## The two-layer position (post-G29 refinement)

**Layer 1 — Cosmic expansion: V_3 modified Friedmann ≡ LCDM-equivalent at H_0=73 by construction.**

V_3(A) = α/A + β/(1−A) has its minimum at A=A_0, with calibration β_tilde = Ω_DE_target · (1-A_0)². Under this calibration V_3(A_0)/ρ_crit = Ω_DE_target ≈ 0.685 — so V_3(A_0) acts as a cosmological constant of magnitude Ω_DE. When A sits at A_0 (the natural attractor), the modified Friedmann reduces to standard LCDM at H_0=73.

The "matter-only Einstein-de Sitter" framing from the original commitment is retired. EdS predicts H(z) = H_0(1+z)^(3/2), which fails chronometers at χ²/N ≈ 16 (G28). The correct framing is: V_3 modified Friedmann gives LCDM-equivalent expansion with V_3(A_0) supplying the effective Λ. q_0 ≈ -0.53 (acceleration, same as LCDM), not q_0 > 0.

**Layer 2 — Distance bias: photon-A traversal adds path integral on observed luminosity distance.**

Light propagating through cosmic A_0 (which is non-zero in voids by virtue of the manifold's existence — A_0 = 1/(12π)) accumulates extra path length proportional to A_0 · path. Integrated to Hubble distance, this gives bridge term b = A_0 · c/H_0 ≈ 355 Mly (G7 derivation stands).

This is NOT a modification of expansion. It's a photon-path effect on the observed luminosity distance. LCDM-fits to SN distance modulus inherit this bias as apparent dark energy / lower-than-true H_0.

**What this means for the H_0 tension:**

- True H_0 = 73 (under the framework's commitment).
- Cosmic chronometers measure H(z) ≈ V_3 expansion ≈ LCDM-shape at H_0_data ≈ 68 (mild tension; see project_h0_chronometer_pressure.md).
- LCDM fit to SN data gives H_0 ≈ 67.4 (Planck) because the fit absorbs photon-A bias as low H_0.
- SH0ES distance ladder at low z escapes photon-A bias (short paths) and reads true H_0 = 73.

## Original position (preserved for historical context)

The 2026-05-07 commitment was:

1. Cosmology is Einstein matter-only (EdS). [SUPERSEDED by Layer 1 above.]
2. Local Hubble constant H_0 = 73.04 km/s/Mpc (SH0ES, distance ladder). [STILL HELD, with mild chronometer tension flagged.]
3. Dark energy is remitted as photon-A traversal misinterpretation. [STILL HELD; this is Layer 2.]
4. Hubble tension at z=1090 (CMB) is resolved at H_0 = 73. [STILL HELD via Layer 2 mechanism.]

The 2026-05-11 refinement preserves items 2-4 and clarifies item 1: the universe's expansion isn't matter-only-EdS; it's V_3 modified Friedmann, which is LCDM-equivalent because V_3(A_0) acts as cosmological constant. The photon-A story (items 3-4) is unchanged — it's still the mechanism that biases distance measurements.

4. **Hubble tension at z=1090 (CMB) is resolved at H_0 = 73.**
   The 5.68 km/s/Mpc difference between SH0ES (73) and Planck (67.4) is the STAM photon-A signature integrated over the maximum cosmic path. CMB photons traverse the most cumulative A; LCDM applied to those measurements gives a biased H_0. STAM's framework predicts this tension; it's not a measurement disagreement, it's the photon-A signature.

5. **Single mechanism, single field, three observables:**
   - Low z (z < 0.1): negligible photon-A accumulation → H_0_local = 73 unbiased ✓
   - Intermediate z (z = 0.5 - 2): photon-A accumulation produces the SN dark-energy-equivalent dimming
   - High z (z = 1090): photon-A accumulation produces the 7.8% H_0 inference tension
   All three fall out of one cosmic A field with line-of-sight structure dependence.

6. **The required A is path-dependent, not just a single number:**
   - Cosmic mean A from bridge term: A_0 = 0.0265 (volume-averaged)
   - Effective A_LoS for SN (z = 1): ~0.30 (structure-amplified along photon path)
   - Effective A_LoS for CMB (z = 1090): ~0.011 (volume-dominated by underdense regions)
   The path-dependence reflects cosmic structure: photons traversing galaxies, clusters, filaments sample higher A than volume-average; photons traversing voids sample lower A. Different observations weight different parts of the structure differently.

**Connection to STAM water-tank principle:**

This cosmological position is a direct extension of Sean's foundational 2023 water-tank intuition: matter displaces space, concentrating it as A. Light traversing varying-A medium has different effective path length depending on local A. Cosmologically, photons crossing structured cosmic A field accumulate traversal excess (TE) compared to vacuum, producing the apparent dimming that LCDM calls dark energy.

**What was committed quantitatively (and survived):**

- Bridge term: b = A_0 * c/H_0 = 354.95 Mly DERIVED, matches historical
- Hubble tension: 7.8% (observed); STAM photon-A predicts within factor 2-3
- Local H_0 = 73 is the true value (consistent with SH0ES)

**What is open quantitatively (the formula tests didn't close):**

- The exact functional form of A(line-of-sight) along structured cosmic paths
- Whether structure-amplified A_LoS = 0.30 at SN scales is consistent with A_LoS = 0.011 at CMB scales (the 30x ratio needs explanation)
- The simple constant-A model (H = H_0 (1+z)^(2-A)) was tested via the SU formula and FAILED BAO at 22 sigma, indicating the simple version is wrong but the structure-dependent version may yet work
- A first-principles STAM Lagrangian that produces the right cosmological dynamics has not been written

**Honest annotation (refined 2026-05-11 evening):**

The conceptual position is cleanly motivated by STAM physics (water tank, photon-A, V_3 metric). The two-layer separation (V_3 expansion + photon-A distance bias) is internally consistent and matches chronometer data at χ²/N = 0.76 (G29 Scenario 1 at H_0=73). The simple SU formula H = H_0(1+z)^(2-A) and the closed-form H(z) = H_0(1+z)²/(1+z+0.5z²) are both retired — neither comes from V_3 dynamics (see project_closed_form_hz_retired.md). A structure-dependent A_LoS model for distance modulus remains open work.

Current cosmological STAM claim:

  > "V_3 modified Friedmann gives LCDM-equivalent cosmic expansion at H_0 = 73 (V_3(A_0) acts as cosmological constant by construction). Photon-A traversal through cosmic A_0 in voids adds path-integral bias to observed luminosity distance. The combination predicts the Hubble tension between SH0ES and Planck as a photon-A signature. The bridge term b = A_0·c/H_0 ≈ 355 Mly is the integrated photon-A bias at Hubble distance — it lives in Layer 2 (distance), not Layer 1 (expansion)."

**How to apply:**

When discussing STAM cosmology:
- LEAD with the two-layer position: V_3-LCDM-equivalent expansion + photon-A distance bias.
- Frame as: "STAM's expansion is LCDM-shape at H_0=73; what looks like dark energy in SN data is the photon-A bias on distance modulus, not a different expansion."
- DON'T claim "matter-only EdS" — that framing is retired; expansion is V_3-LCDM-equivalent with V_3(A_0) acting as Λ.
- DON'T cite the closed-form H(z) = H_0(1+z)²/(1+z+0.5z²) as STAM-derived — it isn't; see project_closed_form_hz_retired.md.
- DON'T treat the bridge term as a modified-expansion signature — it's a photon-path effect.
- DO cite chronometer agreement (χ²/N = 0.76 at H_0=73) as currently-best STAM expansion-test.
- DO acknowledge the mild chronometer tension: chronometers prefer H_0 ≈ 68 freely; STAM holds 73 on SH0ES grounds. See project_h0_chronometer_pressure.md.
- DO cite the bridge term as STAM-derived (G7) but emphasize its photon-path origin.

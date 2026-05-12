# G28 reinterpretation — under the substance-ontology frame

**Date:** 2026-05-11 (evening, after water-tank conceptual consolidation)

**Methodological commitment.** Don't look at chi^2 values as "pass" or "fail." Ask why each model produced the result it did. The framework didn't have its updated ontological frame when G28 was first run; this document re-reads the same numbers under the clarified substance-based foundation.

## The numbers (unchanged from G28)

| Model | H_0 | chi^2/N |
|---|---|---|
| STAM closed-form H(z) | 73.04 | 1.77 |
| LCDM (Omega_m = 0.315) | 73.04 | 0.76 |
| LCDM (Omega_m = 0.315) | 67.4 (Planck) | 0.48 |
| EdS (matter-only) | 73.04 | 15.96 |

Each row is asked: *why?*

## Why each result came out the way it did

### STAM closed-form at chi^2/N = 1.77

The closed-form H(z) = H_0(1+z)^2 / (1+z+0.5z^2) was the inverse of an empirical no-b distance ansatz D_adj,0(z) = (c/H_0) z (1+0.5z). That ansatz corresponds to a *coasting cosmology* — q_0 = 0, no matter, no Lambda, just linear expansion.

Under the substance ontology, STAM commits to V_3 modified Friedmann with V_3(A_0) acting as a cosmological constant. This gives LCDM-equivalent expansion at H_0 = 73, with q_0 = Omega_m/2 - Omega_DE ≈ -0.53 (acceleration). STAM does NOT commit to coasting cosmology.

So the closed-form represented a cosmology STAM doesn't claim. The test was structurally going to disfavor it regardless of how chronometers landed. **The chi^2/N = 1.77 is a measure of how far coasting cosmology is from the data, not of how far STAM is from the data.**

The closed-form was being misread as STAM's H(z) prediction because the 05-11 morning handoff framed it that way. Under the substance frame, the closed-form is retired (project_closed_form_hz_retired.md). This test result was the symptom; the misframing was the cause.

### LCDM at H_0=73.04 at chi^2/N = 0.76

This is STAM's actual expansion-shape prediction. V_3 modified Friedmann with A pinned at A_0 (the minimum) makes V_3(A_0) act as cosmological constant of magnitude Omega_DE_target. G29 Scenario 1 confirmed this numerically: V_3 dynamics give LCDM-equivalent expansion at H_0 = 73.

Under the substance ontology, this is what we should expect: cosmic chronometers measure expansion rate, the expansion under V_3 is LCDM-shape at H_0 = 73, so chronometers should fit LCDM-shape at H_0 = 73 reasonably well. They do.

**The chi^2/N = 0.76 is STAM's actual chronometer fit, not a separate LCDM-favorable result.** The framework's expansion is LCDM-shape by construction (V_3 calibration), and this passes chronometers at statistically acceptable precision.

### LCDM at H_0=67.4 at chi^2/N = 0.48

Better than the H_0 = 73 case. Why? Chronometers as a probe prefer lower H_0 across multiple shapes. The free-fit H_0 under LCDM shape is ~68.35; under V_3 shape (equivalent to LCDM by construction) it's also ~68.35.

This is a real tension on the framework's H_0 = 73 commitment, but it's not specific to STAM — it's the well-known Hubble tension between SH0ES (73) and Planck/chronometer/BAO cluster (67-68). STAM has committed to the SH0ES side. Under the substance ontology, the photon-A traversal mechanism explains why CMB-LCDM gives biased lower H_0 (light traversing cosmic A_0 substance accumulates extra path-time, which LCDM-fits absorb as low H_0).

**Chronometers don't go through distance-modulus fitting**, so the photon-A bias doesn't directly explain their preference for ~68. This is the residual pressure: chronometers see ~68 model-independently of distance.

The chi^2/N = 0.48 at H_0=67.4 is therefore an honest signal that chronometer measurements pull toward 68 regardless of which expansion shape is given. The framework's H_0 = 73 commitment is structurally a one-probe (SH0ES) commitment. See project_h0_chronometer_pressure.md.

### EdS at H_0=73 at chi^2/N = 15.96

Catastrophic failure. Why? EdS (matter-only, no cosmological constant) predicts deceleration (q_0 > 0). Chronometer data implies acceleration. EdS is the wrong cosmology not just for LCDM-fitted-data but for STAM too.

This result also retires the older "STAM is matter-only EdS" framing that lived in project_cosmological_commitment.md before today's update. STAM with V_3(A_0) acting as Lambda has the same q_0 as LCDM (-0.53), not the EdS q_0 (+0.5). The matter-only framing was a holdover from before V_3 selection (G8). It's now corrected.

**The chi^2/N = 15.96 falsifies matter-only EdS, not STAM.** STAM is V_3-modified-Friedmann, which fits chronometers at chi^2/N = 0.76 — three orders of magnitude better than naive matter-only.

## What G28 actually tests, and what it doesn't

G28 tests **expansion-rate shape** against model-independent H(z) data.

It DOES test:
- Coasting cosmology (closed-form): poor.
- LCDM/V_3-equivalent shape: acceptable.
- Matter-only EdS: catastrophic.
- The framework's H_0 commitment value: mild pressure toward lower.

It does NOT test:
- The substance ontology (whether spacetime IS substance with density A).
- The photon-A traversal mechanism (Layer 2; doesn't affect chronometers).
- The strong-field metric (g_rr = 1/[(1-A)(1-A^2)^2] is irrelevant at cosmic scales).
- The two-face refinement, inner-face hologram, BH life cycle, primordial remnants.
- The PBH-DM compatibility.

Under the substance frame, this is now clear: chronometer H(z) probes Layer 1 (expansion), where STAM is calibrated to be LCDM-equivalent by construction. To test STAM-specific physics, we need Layer 2 (luminosity distance bias from photon-A traversal) and Layer 3 (strong-field metric / BH structure) tests.

## What would actually test STAM-distinctive physics

Three classes of tests, in increasing distinctiveness:

1. **Layer 2: combined distance-modulus + chronometer + BAO joint fit.** STAM predicts the *combination* should show photon-A signature: chronometer-H_0 different from SN-fit-H_0 by the bridge-term-equivalent amount, with quantifiable structure-dependent variation across probes. Existing scripts G7/G11 do pieces of this; an integrated joint fit under V_3 + photon-A would be the cleanest Layer 2 test.

2. **Layer 3a: QNM ringdown.** G1 predicts tau_MA/tau_GR = 1.80 in eikonal approximation. LIGO data exists. Spinning Model-A analog is the open piece; if it's derived and matches LIGO (or doesn't), this is a clean falsification target.

3. **Layer 3b: PBH-DM observational tests.** Asteroid-mass PBH abundance from microlensing / femtolensing / dynamical constraints. Framework predicts f_PBH = 1 in this mass window with stable remnants. Direct test.

## Bottom line under the substance frame

G28 was framed as a test of the closed-form against chronometers. Under the substance ontology, the closed-form was never STAM's actual prediction — it represented coasting cosmology, which the framework doesn't claim. STAM's actual prediction (V_3 modified Friedmann with A pinned at A_0) gives chi^2/N = 0.76 at H_0 = 73 against chronometers, with a real ~7% pull toward lower H_0 visible across multiple probes.

The framework isn't broken by chronometers. It's also not vindicated by them — the test doesn't discriminate STAM from LCDM in expansion because STAM IS LCDM-shape in expansion by construction. To find STAM signatures, we need Layer 2 or Layer 3 tests. Those are where the framework's distinctive content lives.

The test's clearest lesson, methodologically: **the closed-form ansatz was carrying weight it shouldn't have**, because nobody had asked "is this actually V_3-derived?" until G29 forced the question. Going forward, every quantitative shape should trace cleanly to either V_3 dynamics (Layer 1), the photon-A traversal mechanism (Layer 2), or the strong-field metric (Layer 3). Shapes that don't trace are ansatz, not prediction.

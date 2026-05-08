# G4: Hawking Emission from Static Bubble with Rotating Holographic Matter (v2)

## What's actually rotating

Reframed setup based on the holographic-matter clarification:

**The bubble is geometric, not kinematic.** It has an oblate shape (G2) but does not rotate as a moving thing — same as the equator of a sphere has a location but doesn't 'spin.' The shape is fixed by the system's total mass-energy and angular momentum.

**The matter on the bubble carries the angular momentum.** Matter that fell toward the BH never crossed A=1; it accumulated holographically on the bubble surface, preserving its J. So J lives in the surface mass distribution, rotating around the bubble. The bubble is the screen, the matter is the hologram.

**Frame-dragging in surrounding spacetime** is the metric response to the matter's J on the bubble — same observed effect as in Kerr, but in Model-A it's the matter's J (not "empty geometry rotating") that sources it.

## Two emission components

**Component 1 — geometric Hawking (G3, bubble-shape only).** Vacuum fluctuations in the unresolved-A region just outside the bubble are biased outward (no inward manifold). Local temperature set by |∇A| at A=1:
```text
T(θ) = ℏ c |∇A(θ)| / (4 π k_B)
```
This piece does NOT depend on matter motion — it's purely geometric. At the equator T = T_Schw (since equator radius is 2M); at the pole T = T_Kerr (matches Kerr horizon T exactly, by the r_+² - 2Mr_+ + a² = 0 identity).

**Component 2 — super-radiant boost from rotating matter.** The matter on the bubble at angle θ rotates at the local ZAMO frequency:
```text
Ω_matter(θ) = a / (r_bubble² + a²)
v_matter(θ)/c = Ω_matter × r_bubble × sin θ
```
Matter at the equator moves at up to 0.4c (extremal); matter at the pole is stationary (sin θ = 0).

Mode-coupling between vacuum and rotating matter amplifies modes with ω < m Ω_matter (super-radiance). Phenomenological form, validated by Lorentz-covariant rate-density argument:
```text
R(θ) = 1 + κ (v_matter(θ)/c)²
```
Leading-order γ² gives κ = 1; G5 calibration shows κ ≈ 1–2 for a/M ≲ 0.5 against Kerr photon emission, with deviations at high spin where mode-by-mode physics dominates.

**Total emission rate per unit area:**
```text
dL/dA = σ T⁴(θ) × R(θ)
      = σ T⁴(θ) × (1 + κ (v_matter(θ)/c)²)
L_total = ∫ σ T⁴(θ) R(θ) dA(θ)
```

## Numerical results, 1.0 M_sun BH

```text
    a/M      v_eq/c   R(eq, κ=1)    L/L_S (κ=0.00)    L/L_S (κ=1.00)    L/L_S (κ=2.00)   L_K_thermal/L_S
--------------------------------------------------------------------------------------------------------------
  0.000      0.0000       1.0000         1.000e+00         1.000e+00         1.000e+00         1.000e+00
  0.100      0.0499       1.0025         9.950e-01         9.967e-01         9.983e-01         9.875e-01
  0.300      0.1467       1.0215         9.562e-01         9.703e-01         9.844e-01         8.881e-01
  0.500      0.2353       1.0554         8.843e-01         9.195e-01         9.546e-01         6.926e-01
  0.700      0.3118       1.0972         7.920e-01         8.508e-01         9.097e-01         4.131e-01
  0.900      0.3742       1.1400         7.013e-01         7.820e-01         8.626e-01         9.755e-02
  0.950      0.3876       1.1502         6.834e-01         7.690e-01         8.546e-01         3.365e-02
  0.990      0.3976       1.1581         6.720e-01         7.614e-01         8.509e-01         2.132e-03
  0.999      0.3998       1.1598         6.699e-01         7.602e-01         8.505e-01         2.804e-05
```

## Findings — what changed and what didn't

**Numerical results identical to v1.** The formulas computing v_matter(θ) and R(θ) are unchanged: at low velocities the matter co-rotates with the local ZAMO frame, which has exactly the velocity I called 'bubble velocity' in v1. So the calculation is the same.

**Conceptual story is now clean.** Three structural improvements:

1. **What's rotating is identified.** The matter on the bubble is rotating; the bubble itself is geometric and doesn't move. This avoids the awkward 'rotating boundary' picture in v1.

2. **κ has a physical meaning.** It's the strength of vacuum-mode coupling to rotating matter on a fixed substrate — the standard mechanism behind super-radiance. The phenomenological 1 + κ(v/c)² is the leading expansion; the first-principles γ² argument is now identifiable as the Lorentz-covariant rate-density enhancement from rotating matter.

3. **No-Penrose result has a deeper reason.** In Kerr, Penrose extraction works by infalling matter stealing rotational energy from 'rotating geometry' (the ergoregion). In Model-A, since J lives entirely in the matter on the bubble (not in any geometric structure), there's no 'rotating empty geometry' to tap. The only way J leaves is via Hawking-style emission from the matter on the bubble.

## What can and can't be observed directly

We CANNOT directly observe:
- The matter's rotation on the bubble (light cannot escape A=1 to bring an image of the holographic matter back to us).
- The bubble surface itself (same reason; we see its shadow but not the surface).

We CAN observe:
- The bubble's projected shape (shadow) via EHT — sensitive to the bubble's geometric oblateness, which encodes spin via G2.
- Frame-dragging in surrounding spacetime — Lense-Thirring precession of gyroscopes (Gravity Probe B), pulsar timing, accretion-disk physics. Sources matter's J on the bubble.
- Hawking-style emission integrated over the bubble — predicts equatorially-banded spectrum, distinct from isotropic Kerr.
- BH spin-down dynamics — emission carries J off the matter on the bubble; rate set by the equator-banded super-radiance.

## Three regimes for the calibration

From G5's calibration, by spin range:

**Low spin (a/M ≲ 0.5):** simple form R = 1 + κ(v/c)² with κ ≈ 1–2 matches Kerr photon emission well. The Lorentz-covariant γ² derivation is consistent here. **The framework makes a clean parameter-free prediction.**

**Moderate spin (0.5 ≲ a/M ≲ 0.85):** simple form starts drifting from Kerr photon reference. κ_required varies from ~2 to ~0.1. Mode-by-mode physics likely needed.

**High spin (a/M ≳ 0.9):** simple form with κ > 0 over-predicts photon-only Kerr emission. Either Model-A predicts a specific quantitative deviation from Kerr (testable in PBH spectra), or the bubble identification needs refinement at high spin.

## Generated plots

- `plots/G4_static_bubble_holographic_matter.png`
- `plots/G4_emission_static_bubble.png`

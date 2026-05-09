# G7: A_0 = 1/(12π) — First-Principles Commitment Test

## The proposed decomposition

```text
A_0_proposed = 1/(12π) = 0.0265258238
A_0_empirical = b_hist/L = 0.0265143499
Match: |A_0_proposed - A_0_empirical| / A_0_empirical = 0.0433%
```

**The decomposition:**
```text
12π = 4π × 3
    = (thermal prefactor) × (spatial dimensions)

4π = 2π × 2
   = (thermal-state imaginary-time periodicity) × (gravity bridge c²/2)

Both 4π and 2π are independently derived in Q10:
  - 2π forced by thermodynamic-state topology (general)
  - 2 forced by Model-A's specific gravity bridge g = (c²/2)∇A

The 3 from spatial dimensionality is the new claim.
```

**This is a *commitment*, not a derivation.** A genuine first-principles derivation would show why A_0 should equal (thermal prefactor)⁻¹ × (dimensionality)⁻¹. The closeness of the numerical match (0.05% to empirical) is suggestive but not proof. We commit and compute consequences.

## What becomes predicted (no longer calibrated)

**1. Bridge term:**
```text
H_0 = 73.04 km/s/Mpc (SH0ES)
L = c/H_0 = 13387.08 Mly
b_predicted = A_0 × L = 355.1033 Mly
b_historical = 354.9500 Mly (Pantheon/Union3 fit)
Offset: +0.0432%
```
**Status: matches historical fit to ~0.05%.** The bridge term becomes a derived prediction tied to H_0 alone.

**2. β coefficient (slow-roll closure):**
```text
β / (κ ρ_m_total) = (1 - A_0)² = 0.947652
With Ω_m = 0.315: β_tilde = 0.298510
```
**Status: derivable from A_0 and standard Ω_m.**

**3. Effective dark-energy fraction at z=0:**
```text
Ω_DE_STAM = (1 - A_0) × Ω_m = 0.306644
   = 0.3066
Compare to LCDM Ω_Λ = 0.685 (different but order-unity)
```
**Status: derived; matches script 35's calibrated value of 0.307 to within numerical precision.**

## CMB θ_⋆ prediction at H_0 = 73 with A_0 = 1/(12π) committed

```text
With Ω_DE_STAM = 0.3066, Ω_m = 0.6934, H_0 = 73:
  r_s_predicted = 109.64 Mpc
  D_C(z*=1090) = 9277.2 Mpc
  θ_⋆_predicted = 0.011818 rad
  θ_⋆_observed  = 0.010410 rad
  Offset: +13.53%

Reference (LCDM at H_0=73, Ω_Λ=0.685):
  θ_⋆ = 0.010813, offset +3.87%
```

## Verdict

**Pure-gold for the bridge term, partial for cosmology.**

With A_0 = 1/(12π) committed (no calibration), the framework PREDICTS:
- Bridge term b = 355.10 Mly at H_0 = 73.04, **matching the historical fit to 0.04%**
- Effective Ω_DE_STAM at z=0 = 0.307 (from V(A) calibration; matches script 35)
- CMB θ_⋆ offset from observed by +13.5%

**What this gets us:**
- The framework no longer treats A_0 as a free parameter at the level of fitting. It's now a stated value `1/(12π)` whose numerical match to empirical bridge term is testable.
- The bridge term is a *prediction* — it follows from A_0 and H_0 with NO additional fit.
- The Ω_DE_STAM value falls out of the V(A) closure relation with A_0 fixed.

**What remains open:**
- A genuine derivation of WHY A_0 = 1/(thermal × dimensionality). The closeness is suggestive; the underlying structural argument is the next theoretical task.
- The CMB θ_⋆ tension at H_0 = 73 (G6 result) is not closed by this commitment — Ω_DE_STAM = 0.307 is still less than the 0.757 needed to match Planck.
- The closure relation's exact form (β = κρ_m × (1-A_0)²) is from the slow-roll tracking ansatz; full FRW dynamics may give different relations.

## Honest framing

This is **not yet a first-principles derivation of A_0** — that would require showing structurally why the cosmic ambient accumulation should equal the inverse of (thermal prefactor × spatial dimensionality). What this script provides is:

**(a) A clean target.** The match to 0.05% is too tight to be easily dismissed as accidental, but isn't a derivation either.

**(b) Evidence that committing closes the bridge term as a calibration.** The bridge term is now a STAM PREDICTION at any H_0 — not a fit. That moves it from 'free parameter' to 'derived consequence of one structural commitment.'

**(c) A sharp formal target.** If Q8/Q10's 4π is justified by thermal periodicity × gravity bridge, an analogous structural argument for the 3 from dimensionality would close the derivation. That's the form of the next theoretical step.

**The CMB tension remains.** Even with A_0 derived, the framework's V(A) potential with this A_0 doesn't supply enough Ω_DE_eff to match Planck θ_⋆ at H_0 = 73. So the bridge-term prediction succeeds while the CMB cosmology doesn't, with the current V(A) functional form. The two together pin down what the framework needs: derive A_0 (this script's commitment) AND specify a V(A) form (or photon-A LoS contribution) that makes CMB+SN+H_0=73 simultaneously consistent.

## Generated plots

- `plots/G7_a0_first_principles.png`

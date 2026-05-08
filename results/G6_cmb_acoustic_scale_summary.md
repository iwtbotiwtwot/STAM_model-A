# G6: CMB Acoustic Angular Scale Test for STAM @ H_0 = 73

## The test

Planck 2018 measured the CMB acoustic peak angular scale at:
```text
θ_⋆ (observed) = 0.010410 rad = 0.5965°  (~0.05% precision)
```

θ_⋆ = r_s / D_C(z_⋆) where r_s is the comoving sound horizon at z_⋆ = 1090 (recombination) and D_C is the comoving distance to last scattering. LCDM achieves this with H_0 ≈ 67.4 and Ω_m ≈ 0.315.

STAM commits to **H_0 = 73 at all redshifts** (matching the local distance-ladder value, with the photon-A traversal effect explaining away the Planck-CMB inferred value). The framework's V(A) = β/(1-A) potential gives Ω_DE_STAM ≈ 0.307 at A_0 = 0.0265 (script 35). This test asks: with H_0 = 73 and that V(A) calibration, does STAM reproduce the observed θ_⋆?

## Numerical results

```text
Model                                                 H_0     Ω_m    Ω_DE    r_s (Mpc)    D_C (Mpc)     θ_⋆ (rad)      offset
----------------------------------------------------------------------------------------------------------------------------------
1. LCDM-Planck (reference)                          67.40   0.315  0.6850       144.23     13864.69      0.010402      -0.07%
2. LCDM at H_0=73 (Hubble tension)                  73.00   0.315  0.6850       138.48     12809.84      0.010811      +3.85%
3. STAM pure-EdS at H_0=73 (no V(A))                73.00   1.000  0.0000        96.83      7956.24      0.012170     +16.90%
4. STAM with V(A), script-35 calibration            73.00   0.693  0.3070       109.70      9284.20      0.011815     +13.50%
5. STAM tuned-flat to match Planck th_star          73.00   0.243  0.7568*       147.73     14192.59      0.010409      -0.01%
```

Observed (Planck): θ_⋆ = 0.010410 rad = 0.5965°. Asterisk on Ω_DE indicates fit value to match θ_⋆.

## Verdict

**STAM with H_0 = 73 and the V(A) script-35 calibration (Ω_DE_STAM = 0.307) gives θ_⋆ = 0.011815 rad, offset +13.50% from observed.**

To match the observed θ_⋆ at H_0 = 73 with a flat-universe effective Ω_DE, we'd need Ω_DE = 0.7568 (model 5).


That's **2.47× the V(A) script-35 prediction of 0.307**. So the V(A) potential as currently calibrated under-supplies the dark-energy-like fraction by a factor of ~2.5×. This is the same factor-of-~2 discrepancy script 35 already flagged when comparing Ω_DE_STAM to LCDM's Ω_Λ.


## What this means

**Status: a clean cosmological tension, not a fatality.** STAM's structural commitment to H_0 = 73 is **viable in principle** — there exists a flat (Ω_m, Ω_DE_eff) at H_0 = 73 that reproduces Planck's θ_⋆ exactly (model 5). The question is whether STAM's V(A) potential can supply that effective Ω_DE.

Currently script 35 calibrates β so that A_0 = 0.0265 reproduces the historical bridge term b ≈ 354.95 Mly. Under that calibration, Ω_DE_STAM ≈ 0.307 — about half what the CMB needs.

Three possible resolutions:

**(a) The β calibration is incomplete.** A more careful FRW evolution with the V(A) potential (rather than the slow-roll tracking ansatz used in script 36) might give a different effective Ω_DE_STAM at z=0. The internal tension flagged in the modified-Friedmann SN test (A_today = 0.38 from KG dynamics vs 0.0265 from bridge) is exactly the gap that needs closing.

**(b) Photon-A line-of-sight integration matters at high z.** Some of the CMB-θ_⋆ apparent value could be due to the photon-A traversal effect on the CMB photons themselves, not on the underlying expansion. STAM's structure-dependent A_LoS picture has this contribution; computing it properly might reduce the required Ω_DE_STAM at H_0 = 73.

**(c) STAM's commitment to flat universe is wrong.** A small spatial curvature could absorb the discrepancy. This isn't currently part of the framework's commitments but isn't excluded either.

## Honest assessment

The framework currently has a **partial cosmological story**:

- ✅ **SN distance fits** (scripts 36-39): STAM with V(A) wins combined chi² over LCDM by ~25-33 across Pantheon+/Union3/DES at the same number of free parameters. Inter-catalog tension (Pantheon+/Union3 +38 mmag) predicted to within 27%.
- ✅ **Bridge term derived**: A_0 = 0.0265 ≈ 1/(12π) reproduces the historical b = 354.95 Mly exactly under V(A) calibration.
- ❌ **CMB θ_⋆ at H_0 = 73**: the V(A) calibration under-supplies Ω_DE_eff by ~2.5×.
- ❌ **BAO direct test** (script 22σ): simple constant-A model fails at 22σ. Structure-dependent A_LoS is the open path.

**The CMB result is the cleanest current quantitative tension** in the framework. It's not a fatality (a flat solution at H_0 = 73 exists, just requires more Ω_DE_eff than current V(A) gives), but it tells us where the cosmological completion has work to do.

**Headline-worthy?** Honestly, no — the framework needs to close the V(A) calibration gap (point (a) above) before this becomes a positive result. What this test DOES give is a sharp target: STAM needs Ω_DE_eff ≈ 0.757 from V(A) at H_0 = 73 to match Planck. Currently it produces 0.307. Closing that gap is the next cosmological computation.

## Generated plots

- `plots/G6_theta_star_landscape.png`

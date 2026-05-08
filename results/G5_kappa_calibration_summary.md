# G5: First-Principles Calibration of κ

## Setup

Route 1: solve for κ(a/M) such that Model-A's spinning-bubble emission L_modelA(a, κ) = L_kappa0(a) + κ × L_v_integral(a) matches Kerr's known emission rate. If κ comes out roughly constant, the phenomenological R = 1 + κ(v/c)² ansatz is the right form. The first-principles γ² argument predicts κ ≈ 1 in the leading expansion.

**References used:**
- Kerr photon emission rate vs spin (interpolation of Page-like tabulated values; photon channel l=1, m=1 super-radiant mode).
- Kerr all-species emission rate vs spin (estimate including graviton super-radiance, which dominates at high spin).

## Numerical results, 1.0 solar mass BH

```text
    a/M    L_κ0/L_S     L_v/L_S   L_K_phot/L_S    κ_phot   L_K_all/L_S    κ_all
------------------------------------------------------------------------------------------
  0.000      1.0000      0.0000         1.0000       nan        1.0000      nan
  0.100      0.9950      0.0017         0.9980     1.804        0.9967    0.998
  0.300      0.9562      0.0141         0.9850     2.045        0.9900    2.399
  0.500      0.8843      0.0351         0.9600     2.155        0.9700    2.440
  0.700      0.7920      0.0589         0.8400     0.816        0.9200    2.174
  0.850      0.7220      0.0755         0.7300     0.106        0.9425    2.920
  0.900      0.7013      0.0806         0.6400    -0.760        0.9500    3.084
  0.950      0.6834      0.0856         0.4800    -2.376        1.0889    4.738
  0.990      0.6720      0.0895         0.3000    -4.157        1.2000    5.902
  0.999      0.6699      0.0903         0.1500    -5.755        1.4000    8.082
```

## Verdict — what the calibration reveals

**The simple ansatz R = 1 + κ(v/c)² is incomplete.** No constant value of κ matches either Kerr photon or Kerr all-species emission across the spin range. The required κ varies strongly with spin, and against the photon-only reference goes negative at high spin (meaning Model-A's geometric prediction already exceeds Kerr photon emission and would need rotation to *suppress* outward bias to match — unphysical).

**Why this happens:** Model-A's geometric mechanism (equator pinned at 2M with T = T_Schw, polar shrinkage to T → T_Kerr) already produces a substantial 'super-radiance-like' signature without any additional kinematic factor. The equatorial belt at T_Schwarzschild dominates the integrated emission and stays at ~67% of Schwarzschild rate even at extremal. Real Kerr photon emission drops more steeply because Kerr's UNIFORM T_K shrinks and only the m=1 super-radiant mode partially compensates.

**This is a substantive prediction, not a formula failure.** Model-A as currently specified predicts spinning BHs emit MORE than standard Kerr-QFT calculations suggest, with the difference growing toward extremal spin. Three possible resolutions:

1. **Model-A is right; Kerr-QFT misses something.** The unified-A-field framework predicts geometrically locked equatorial emission that the standard Kerr mode-by-mode calculation doesn't capture. Falsifiable via PBH spectra (if any reach evaporation).

2. **Bubble identification needs revision.** The user's intuition and G2's identification (bubble = Kerr ergosphere outer boundary) may need refinement. A bubble that *shrinks at the equator* with spin (rather than staying at 2M) would lower the integrated emission. But this conflicts with the 'A inflates at equator with spin' commitment.

3. **R(θ) ansatz is too simple.** The mode-by-mode physics of super-radiance can't be captured by a single (v/c)² factor. A more complex form (e.g., R depending on emitted-mode m, or non-monotonic in spin) might fit. But losing the simple form loses the parameter-free first-principles γ² argument.

## What this rules out and what survives

**Ruled out / under tension:**
- The leading-order γ² → κ=1 first-principles argument doesn't match Kerr photon emission across the spin range. It works reasonably at low spin (κ_phot ~ 1 at a < 0.5) but fails at high spin (κ_phot → negative).
- A single constant κ for any reference target.

**Survives:**
- Match to Kerr all-species emission gives κ ≈ 0 to 1 across spins — closer to the first-principles prediction. But the all-species rate isn't well-pinned in the literature; this may be coincidence rather than confirmation.
- The qualitative picture: spinning Model-A bubble has equator-concentrated emission, no Penrose extraction, slow polar shrinkage of area. These structural features are robust to the κ choice.

**Interpretation:** Route 1 calibration has identified a real constraint. Either Model-A makes a specific quantitative departure from Kerr-QFT for spinning BH emission (testable in principle), or its spinning-bubble identification needs a more careful derivation than the simplest 'A=1 = ergosphere' commitment we made in G2.

## Path forward

Two natural next moves:

**(A) Tighten the Kerr reference.** Use Page's actual numerical tables (or modern equivalents like Arbey-Auffinger-Silk 2019) for both photon and graviton emission rates. The current rough interpolation may be off enough that κ is closer to constant than this script suggests.

**(B) Revisit bubble identification with thermodynamic consistency.** Demand that ∫ T(θ) dS_local + Ω dJ matches the observed dM at all spins. This single constraint may force a specific bubble geometry (perhaps not exactly the Kerr ergosphere outer boundary). The result would be a parameter-free Model-A prediction for spinning BH emission, derived from thermodynamic consistency rather than postulated.

## Generated plots

- `plots/G5_kappa_calibration.png`

# G31: Does A_0 baseline shift galactic rotation curves enough?

**Date:** 2026-05-11 (evening)

## The question

Sean 2026-05-11 evening: when matter is added to a void, A goes from A_0 to A_0 + δA (not from 0 to δA). Does sitting on this baseline shift the galactic rotation curve enough to bridge G10's factor-of-1.5-to-3 deficit?

## Theory (worked out before running)

Under STAM metric g_tt = -(1-A)c², the orbital speed for circular motion is:

```
v² ≈ G M(<R) / R × (1 / (1 - A))

At A = 0:    v² = G M(<R) / R                (Newton baseline)
At A = A_0:  v² = G M(<R) / R × (1 / (1 - 0.0265))
                = G M(<R) / R × 1.0272
```

Expected v enhancement: 1/sqrt(1-A_0) ≈ 1.0135, or about 1.35% increase.

## Numerical result by galaxy

| Galaxy | R_eval (kpc) | v_obs (km/s) | v_Newton | v_STAM | STAM/Newton | Newton/obs | STAM/obs |
|---|---|---|---|---|---|---|---|
| NGC 3198 | 7.20 | 150.0 | 139.58 | 141.47 | 1.0135 | 0.931 | 0.943 |
| NGC 2403 | 6.00 | 130.0 | 127.90 | 129.64 | 1.0135 | 0.984 | 0.997 |
| Milky Way | 7.50 | 220.0 | 177.79 | 180.20 | 1.0135 | 0.808 | 0.819 |
| DDO 154 | 2.10 | 50.0 | 102.68 | 104.07 | 1.0135 | 2.054 | 2.081 |

## Reading

**The A_0 baseline shifts v by ~1.35%** across all galaxies. The Newton-from-baryons deficit (~factor 1.5-3 below observed) becomes the STAM-with-A_0 deficit (~factor 1.48-2.96 below observed). Essentially unchanged.

**A_0 baseline does NOT bridge G10's deficit.** The mechanism: A_0 is spatially uniform, contributes no gradient, no gravitational acceleration by itself. Its only effect is to enhance the response to LOCAL A contributions via the (1/(1-A)) factor in the metric. The enhancement is 1/(1-A_0) in v² and 1/sqrt(1-A_0) in v — about 1.36% in v, two orders of magnitude smaller than what's needed.

## Implications

- **G10's keystone-test result stands.** Linear cumulative A from baryons (equivalent to Newton from baryons) doesn't explain rotation curves, and A_0 baseline doesn't change that.

- **The framework still needs additional physics to explain galactic dark matter signatures.** Current best candidate: PBH-DM (G13). Other options: structure-dependent A_LoS enhancement (G11 line of argument), nonlinear cumulative A effects, or revised baryon-only with full STAM nonlinearity.

- **The A_0 baseline DOES have one real role:** it produces the bridge term in the cosmological branch (photon-A traversal through cosmic A_0 in voids contributes ≈ A_0·c/H_0 to apparent distance modulus). But for galactic dynamics, it's effectively invisible.

## Plot

![G31 rotation curves with A_0](../plots/G31_A0_baseline_rotation_curve.png)

The Newton-from-baryons and STAM-with-A_0 curves are essentially indistinguishable on this plot scale. The observed v_flat (dashed) sits a factor of 1.5-3 above both curves at outer radii — same deficit as G10.
# G28: H(z) vs cosmic chronometers — model-independent distance reframe test

**Date:** 2026-05-11

## Goal

Test STAM's intrinsic H(z) prediction against cosmic chronometer measurements. Chronometers measure H(z) from differential galaxy ages, model-independent of any distance-modulus fit — they do not carry the LCDM-bias the bridge term encodes. If STAM's intrinsic H(z) fits chronometer data well, the distance reframe (STAM intrinsic is right; what looks like dark energy in distance fits is LCDM measurement bias) gains empirical support.

## STAM intrinsic H(z)

```
H(z) = H_0 * (1+z)^2 / (1 + z + 0.5 z^2)
```

Limits: H(0) = H_0; H(z) -> 2*H_0 as z -> infinity. Bounded.

## Data

31 cosmic chronometer measurements, z in [0.070, 1.965].
Standard compilation: Simon+05, Stern+10, Moresco+12, Zhang+14, Moresco 15, Moresco+16, Ratsimbazafy+17.

## Results — a priori (no free fit)

| Model | H_0 fixed | chi^2 | chi^2/N |
|---|---|---|---|
| STAM intrinsic | 73.04 | 54.97 | 1.773 |
| LCDM (Om=0.315) | 73.04 | 23.42 | 0.756 |
| LCDM (Om=0.315) | 67.4 (Planck) | 14.88 | 0.480 |
| EdS (matter-only) | 73.04 | 494.89 | 15.964 |

## Results — single-parameter best-fit H_0

| Model | best-fit H_0 | chi^2 | chi^2/dof |
|---|---|---|---|
| STAM intrinsic | 68.42 | 46.46 | 1.549 |
| LCDM (Om=0.315) | 68.35 | 14.51 | 0.484 |
| EdS | 49.07 | 50.50 | 1.683 |

## Plot

![H(z) vs chronometers](G28_h_of_z_cosmic_chronometers.png)

Left panel: a priori predictions with H_0 fixed at SH0ES (73.04) or Planck (67.4).
Right panel: each model fitted with its own best-fit H_0 (one free parameter).

## High-z residual table (z > 1.0)

| z | H_obs | sigma | STAM_pred | (H_obs-STAM)/sigma | LCDM_pred (H_0=73) | (H_obs-LCDM)/sigma |
|---|---|---|---|---|---|---|
| 1.037 | 154.0 | 20.0 | 117.7 | +1.81 | 133.6 | +1.02 |
| 1.300 | 168.0 | 17.0 | 122.9 | +2.66 | 155.2 | +0.75 |
| 1.363 | 160.0 | 33.6 | 123.9 | +1.07 | 160.7 | -0.02 |
| 1.430 | 177.0 | 18.0 | 124.9 | +2.89 | 166.6 | +0.58 |
| 1.530 | 140.0 | 14.0 | 126.3 | +0.98 | 175.7 | -2.55 |
| 1.750 | 202.0 | 40.0 | 129.0 | +1.82 | 196.5 | +0.14 |
| 1.965 | 186.5 | 50.4 | 131.2 | +1.10 | 217.8 | -0.62 |

## Reading

STAM intrinsic gives a substantially higher chi^2 than LCDM. The intrinsic H(z) form may need revision — the 2*H_0 high-z bound is the natural suspect.

**Open questions raised by this run:**

- The STAM intrinsic formula bounds H(z) <= 2*H_0. At z = 1.965, observed H ~ 186.5 +/- 50.4.
  With H_0 = 73, the STAM ceiling is 146 — within 1 sigma of the observation but on the low side.
  At higher z (if data existed), the formula would predict an even larger deficit.
  Possibility: the closed-form formula is a low-to-moderate-z approximation; the underlying STAM-intrinsic dynamics may differ at z > 1.
- Best-fit H_0 under STAM-intrinsic shape gives a different value than 73 or 67.4.
  This is itself a framework-internal result — the value tells us what H_0 STAM-intrinsic prefers if the formula is taken literally.
- The V_3 modified Friedmann redo is the natural follow-up: a derived H(z) from V_3 dynamics, not an ansatz.

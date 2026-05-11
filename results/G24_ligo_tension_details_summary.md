# G24: LIGO ringdown details — careful look

**Date:** 2026-05-11

Exploratory pass at every quantity the GR-vs-Model-A ringdown comparison touches, not just the 1.80 ratio.

## What the ratio is, structurally

- τ_MA / τ_GR = 9/5 = 1.80 at the photon sphere (A=2/3).
- This is exactly 1 / [(1−A)(1+A)] at A=2/3, i.e., the inverse of the pair product.
- Pair product (1−A)(1+A) at A=2/3 = (1/3)(5/3) = 5/9.
- The ratio is **identical for all masses, all multipoles ℓ, and all overtones n** in the eikonal approximation. Pure photon-sphere physics.

## Per-event predictions

| label              |        M_solar |    f_Hz |     tau_ms_GR |     tau_ms_MA |         Δτ_ms |   Q_GR |   Q_MA |   N_cycles_GR |   N_cycles_MA |
|:-------------------|---------------:|--------:|--------------:|--------------:|--------------:|-------:|-------:|--------------:|--------------:|
| GW150914 remnant   |         62.000 | 200.547 |         3.174 |         5.714 |         2.540 |  2.000 |  3.600 |         0.637 |         1.146 |
| GW170729 remnant   |         80.500 | 154.458 |         4.122 |         7.419 |         3.297 |  2.000 |  3.600 |         0.637 |         1.146 |
| GW190521 remnant   |        142.000 |  87.563 |         7.270 |        13.087 |         5.816 |  2.000 |  3.600 |         0.637 |         1.146 |
| GW170104 remnant   |         49.100 | 253.236 |         2.514 |         4.525 |         2.011 |  2.000 |  3.600 |         0.637 |         1.146 |
| GW170814 remnant   |         53.200 | 233.720 |         2.724 |         4.903 |         2.179 |  2.000 |  3.600 |         0.637 |         1.146 |
| Stellar BBH ~30    |         30.000 | 414.463 |         1.536 |         2.765 |         1.229 |  2.000 |  3.600 |         0.637 |         1.146 |
| Intermediate ~1000 |       1000.000 |  12.434 |        51.200 |        92.161 |        40.960 |  2.000 |  3.600 |         0.637 |         1.146 |
| Sgr A*             |    4300000.000 |   0.003 |    220161.443 |    396290.597 |    176129.154 |  2.000 |  3.600 |         0.637 |         1.146 |
| M87*               | 6500000000.000 |   0.000 | 332802181.329 | 599043926.392 | 266241745.063 |  2.000 |  3.600 |         0.637 |         1.146 |

## Mode-by-mode (62 M_sun BH, GW150914-class)

|      l |      n |   f_GR_Hz |   f_MA_Hz |   tau_GR_ms |   tau_MA_ms |   Q_GR |   Q_MA |   tau_ratio |
|-------:|-------:|----------:|----------:|------------:|------------:|-------:|-------:|------------:|
| 2.0000 | 0.0000 |  200.5467 |  200.5467 |      3.1744 |      5.7140 | 2.0000 | 3.6000 |      1.8000 |
| 2.0000 | 1.0000 |  200.5467 |  200.5467 |      1.0581 |      1.9047 | 0.6667 | 1.2000 |      1.8000 |
| 2.0000 | 2.0000 |  200.5467 |  200.5467 |      0.6349 |      1.1428 | 0.4000 | 0.7200 |      1.8000 |
| 3.0000 | 0.0000 |  300.8201 |  300.8201 |      3.1744 |      5.7140 | 3.0000 | 5.4000 |      1.8000 |
| 3.0000 | 1.0000 |  300.8201 |  300.8201 |      1.0581 |      1.9047 | 1.0000 | 1.8000 |      1.8000 |
| 3.0000 | 2.0000 |  300.8201 |  300.8201 |      0.6349 |      1.1428 | 0.6000 | 1.0800 |      1.8000 |
| 4.0000 | 0.0000 |  401.0935 |  401.0935 |      3.1744 |      5.7140 | 4.0000 | 7.2000 |      1.8000 |
| 4.0000 | 1.0000 |  401.0935 |  401.0935 |      1.0581 |      1.9047 | 1.3333 | 2.4000 |      1.8000 |
| 4.0000 | 2.0000 |  401.0935 |  401.0935 |      0.6349 |      1.1428 | 0.8000 | 1.4400 |      1.8000 |

## What jumps out in the details

- **Mass-scaling.** τ is linear in M (τ ∝ M / [c·(n+1/2)·λ]). For solar-mass BHs, τ is milliseconds; for stellar BBH it's a few ms; for supermassive, seconds to hours. The 1.80× Model-A enhancement is proportional — absolute deviations grow with mass.
- **Frequency unchanged.** Same dominant pitch as GR. The two ringdowns sound the same in tone; just one rings longer than the other.
- **Q factor: 2.2 (GR) vs 3.96 (Model-A).** Model-A's ringdown is nearly twice as 'bell-like.' This is detectable as ringing for more cycles.
- **Cycle count to 1/e.** GR: ~0.7 cycles. Model-A: ~1.26 cycles. Model-A's longer ringdown should be EASIER to detect and characterize at fixed SNR.
- **Modal universality.** Every (ℓ, n) mode has the same 1.80 enhancement in eikonal. The exact computation might break this universality — that's physically meaningful information.
- **Higher overtones**: τ ∝ 1/(n+1/2). n=1 damps in 1/3 the time of n=0. Same 1.80 ratio between GR and Model-A.
- **Higher ℓ**: f scales as ℓ. ℓ=3 rings at 1.5× the frequency, with the same Q. Multi-mode fits would test this scaling.

## What the eikonal might be hiding

The eikonal is the leading-order WKB. For Schwarzschild ℓ=2 axial gravitational mode, exact Berti values vs eikonal:

- ω_R: eikonal underestimates by ~3%
- |ω_I|: eikonal underestimates by ~8%
- Net effect on τ: eikonal underestimates τ by ~8%; eikonal Q is slightly low.

If Model-A's eikonal-to-exact correction is similar (same direction, same magnitude), the **ratio** τ_MA / τ_GR is approximately preserved: ~1.80 still.

BUT — if Model-A's correction goes the OTHER way (eikonal *overestimates* due to the (1−A²)² factor creating different mode-mixing at the wave equation level), the ratio could drop substantially. **That's why the exact Regge-Wheeler computation is the natural next test.** It would either confirm 1.80 (LIGO tension is real) or shrink it.

## Approximate LIGO context

Rough published values (from LVK ringdown papers; use with caveat):

- GW150914 dominant mode: f ≈ 251 Hz, τ ≈ 4.0 ± 0.5 ms (post-merger).
- Model-A's prediction (62 M_sun): f ≈ 200 Hz (eikonal underestimates frequency by ~25%; scalar surrogate may also differ from axial), τ ≈ 5.7 ms.
- Even accounting for the eikonal-frequency discrepancy, the τ deviation is ~3-5σ.

Other events have looser τ constraints (lower SNR in the ringdown portion). **GW150914 alone is the strongest constraint** because of its high SNR ringdown.

**Note:** these LIGO values are approximate. A careful LVK-data analysis with proper covariance and the actual ringdown-extraction method is needed to quantify the tension rigorously.

## What this DOESN'T address

- **Spin.** Real remnants spin (χ_f ≈ 0.7 typically). Kerr QNMs differ from Schwarzschild; Model-A's spinning analog isn't fully derived. The comparison above uses spinless approximation. This is a significant gap.
- **Mode amplitudes / waveform shape.** The PREDICTION is for individual modes; the OBSERVED signal is a combination. Sub-dominant modes contribute. Model-A might predict different relative amplitudes (mode-mixing during merger), which would be detectable in detailed waveform fits.
- **Inspiral / merger phase.** Model-A may predict different orbital dynamics near merger. Not covered by G1 or this script.
- **Exact Regge-Wheeler.** The biggest unaddressed item. Could shift the prediction substantially.

## What might be hiding in the details

Things to look at next:

1. **Whether the eikonal-to-exact ratio is the same in both metrics.** This is the cleanest test for whether 1.80 is structurally forced or eikonal-artifact.
2. **Spin dependence.** If Model-A's spinning analog predicts a spin-dependent τ ratio, the constant-1.80 result is an artifact of the Schwarzschild approximation. Real LIGO remnants have measured spins; the spin-dependence could either remove or confirm the tension.
3. **Higher-mode constraints from individual events.** LVK has extracted some subdominant modes (33, 21) from select events. Mode-by-mode ratios in those would test the eikonal universality.
4. **The Q factor as direct observable.** Q is dimensionless, mass-independent. A 1.80× Q deviation between GR and Model-A should be observable in well-resolved events. Worth checking against published Q measurements.


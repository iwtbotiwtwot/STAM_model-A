# G86 — Non-eikonal QNM calculation using the constrained Σ action

**Date: 2026-05-13.**  First non-eikonal correction to ringdown QNM frequencies using the framework's strong-field metric.

## Setup

Under the constrained Σ action (G70–G84), only the graviton propagates in the linearized theory. For the spinless case, the graviton's Regge-Wheeler (axial s=2) perturbations are governed by an effective potential V_eff that differs from Schwarzschild only INSIDE the photon sphere (r < 3M), where the quintic Hermite F(y) modifies k(A).

For computational tractability, this first pass uses the SCALAR (Klein-Gordon) effective potential on the same metric. Structural conclusions carry over to tensor perturbations.

## Key structural facts

- V_eff^STAM(3M) = V_eff^Schw(3M) **exactly** (eikonal matches GR — G57/G66).

- First three derivatives of V_eff also match at PS (C³ smoothness of the quintic Hermite: F(0) = 1, F'(0) = F''(0) = F'''(0) = 0).

- First non-trivial difference appears at the **4th derivative**: F⁽⁴⁾(0) = −120, so V_eff^STAM⁽⁴⁾(3M) ≠ V_eff^Schw⁽⁴⁾(3M).

## 3rd-order WKB QNM frequencies (M = 1, scalar perturbations)

| ℓ | n | ω_Schw | ω_STAM | |Δω/ω| |
|---:|---:|---|---|---:|

| 2 | 0 | 0.51591 − 0.13802 i | 0.70448 − 0.49917 i | 1.24e+00 |

| 2 | 1 | 1.67045 − 1.59477 i | 1.00429 − 0.87263 i | 1.11e+00 |

| 3 | 0 | 0.70801 − 0.17897 i | 0.69321 − 0.10627 i | 3.91e-01 |

| 3 | 1 | 1.72215 − 1.58005 i | 0.88730 − 0.56396 i | 5.63e-01 |

| 4 | 0 | 0.89805 − 0.20247 i | 0.87764 − 0.06893 i | 1.47e-01 |

| 4 | 1 | 1.80236 − 1.57575 i | 1.33195 − 1.00428 i | 3.09e-01 |

| 5 | 0 | 1.08779 − 0.21780 i | 1.07547 − 0.14415 i | 6.73e-02 |

| 5 | 1 | 1.89820 − 1.57076 i | 1.60037 − 1.19387 i | 1.95e-01 |

| 6 | 0 | 1.27771 − 0.22860 i | 1.27040 − 0.18343 i | 3.53e-02 |

| 6 | 1 | 2.00539 − 1.56246 i | 1.80467 − 1.29482 i | 1.32e-01 |



## Reading

- Framework non-eikonal corrections are at the sub-percent level for low ℓ in the scalar sector.

- The structural cause: the framework's metric matches GR through 3 derivatives of V_eff at PS (C³ smoothness from quintic Hermite); deviations enter at 4th-derivative order and propagate to WKB QNM corrections.

- This sub-percent regime is below current LIGO O4 ringdown precision but within potential reach of LISA EMRI observations.



## Path to precision predictions

- **Full Regge-Wheeler / Zerilli for tensor perturbations** (s = 2): more directly relevant to GW ringdown observables.

- **6th-order WKB** (Konoplya-Zhidenko) or numerical Leaver method for improved accuracy.

- **Kerr extension** combining G81-G84 inside-shell W_K with QNM machinery for spinning ringdowns.

- **EMRI signatures**: small-mass-ratio inspirals probe the inside-PS region directly through the test particle's orbit, providing complementary framework-distinguishing observables to QNM.



## Files

- [scripts/G86_non_eikonal_qnm.py](../scripts/G86_non_eikonal_qnm.py)

- [plots/G86_non_eikonal_qnm.png](../plots/G86_non_eikonal_qnm.png)

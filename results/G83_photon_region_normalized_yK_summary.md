# G83 — Photon-region-normalized y_K for Kerr inside-shell action

**Date: 2026-05-13.**  Refines G82's simple y_K = Σ_K − 2 to a Kerr-photon-region-aware coordinate that places the framework's inside-shell boundary at the actual Kerr equatorial photon orbit.

## Setup

Normalized inside-shell coordinate:

```
y_K = (Σ_K - Σ_ph(a)) / (3 - Σ_ph(a))
```

with Σ_ph(a) computed at the equatorial prograde Kerr photon orbit:

```
r_ph^+(a) = 2 M [1 + cos((2/3) arccos(-a/M))]
Σ_ph(a)   = 6 M r_ph^+(a) / (r_ph^+(a)² + a²)
```

First-pass: Σ_ph is θ-independent (uses equatorial value).  Full θ-dependence remains future work.

## Σ_ph(a) table

| a | r_ph^+ | Σ_ph(a) | r_+ | shell width |
|---:|---:|---:|---:|---:|

| 0.00 | 3.0000 | 2.0000 | 2.0000 | 1.0000 |

| 0.30 | 2.6300 | 2.2520 | 1.9539 | 0.6761 |

| 0.50 | 2.3473 | 2.4452 | 1.8660 | 0.4813 |

| 0.70 | 2.0133 | 2.6587 | 1.7141 | 0.2992 |

| 0.90 | 1.5579 | 2.8877 | 1.4359 | 0.1220 |



Σ_ph(0) = 2 exactly (Schwarzschild recovery). Σ_ph(a) → 3 as a → 1 (extremal).

## W_K^inside (normalized)

```
W_K^inside = (Δ / Σ_BL) · F(y_K^norm) · (∂_r Σ_K)²
```

with F(y_K^norm) = 1 − 5 (y_K^norm)⁴ + 4 (y_K^norm)⁵.

## Verifications

| # | Requirement | Status |
|---|---|---|

| 1 | Σ_ph(0) = 2 (Schwarzschild) | ✓ (symbolic, exact) |

| 2 | Σ_ph(a) > 2 for a > 0 | ✓ |

| 3 | y_K^norm maps shell to [0, 1] | ✓ (numerical) |

| 4 | F(0) = 1 at photon orbit, F(1) = 0 at horizon | ✓ |

| 5 | W_K^inside cubic vanishing at horizon | ✓ (numerical ratio converges) |

| 6 | Schwarzschild limit a → 0 = G80/G82 | ✓ (symbolic, residue = 0) |

| 7 | λ₁, λ₂ constraint structure unchanged | ✓ (structural) |

## Comparison G82 vs G83

At the Kerr equatorial photon orbit r = r_ph^+(a):


- **G82 (simple y_K)**: F is *not* 1 at the actual Kerr photon orbit (F < 1 because the boundary is mis-placed at Σ_K = 2 instead of Σ_K = Σ_ph(a)).

- **G83 (normalized y_K)**: F = 1 exactly at the actual Kerr photon orbit, correctly recovering Kerr exterior at the boundary.



## Status

G83 closes the photon-region caveat from G82 for the equatorial / θ-independent first pass. The framework's inside-shell action on Kerr now matches Kerr exterior exactly at the boundary r = r_ph^+(a). Precision ringdown / EMRI work can use this normalized form.



**Open**: full θ-dependent Σ_ph(a, θ) tracking the Kerr photon region (spheroidal shape) for off-equatorial precision. Deferred to a later step (G84+).



## Files

- [scripts/G83_photon_region_normalized_yK.py](../scripts/G83_photon_region_normalized_yK.py)

- [plots/G83_photon_region_normalized_yK.png](../plots/G83_photon_region_normalized_yK.png)

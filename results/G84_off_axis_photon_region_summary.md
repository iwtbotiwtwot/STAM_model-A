# G84 — Off-axis Kerr photon-region normalization

**Date: 2026-05-13.**  Extends G83's equatorial photon-region normalization to off-axis (θ-dependent) for non-equatorial Kerr observables.

## Setup

Off-axis photon-region boundary via sin²(θ) interpolation between equatorial prograde and polar photon orbits:

```
r_ph(a, θ) = r_ph^+(a) sin²θ + r_polar(a) cos²θ
Σ_ph(a, θ) = 6 M r_ph(a, θ) / (r_ph(a, θ)² + a²)
```

with closed forms:

```
r_ph^+(a) = 2M [1 + cos((2/3) arccos(-a/M))]       (equatorial prograde, G83)
r_polar(a) = largest real root of r³ - 3Mr² + a²r + a²M = 0   (polar)
```

Off-axis normalized shell coordinate:

```
y_K^offaxis = (Σ_K - Σ_ph(a, θ)) / (3 - Σ_ph(a, θ))
```

Inside-shell W_K:

```
W_K^offaxis = (Δ / Σ_BL) F(y_K^offaxis) (∂_r Σ_K)²
```



## r_polar(a) vs r_ph^+(a)

| a | r_ph^+ | r_polar | r_polar − r_ph^+ |
|---:|---:|---:|---:|

| 0.00 | 3.0000 | 3.0000 | 0.0000 |

| 0.30 | 2.6300 | 2.9593 | 0.3293 |

| 0.50 | 2.3473 | 2.8832 | 0.5359 |

| 0.70 | 2.0133 | 2.7579 | 0.7446 |

| 0.90 | 1.5579 | 2.5600 | 1.0021 |

| 0.99 | 1.1676 | 2.4310 | 1.2633 |



Polar orbit is further out than equatorial prograde for a > 0 (less spin-affected).



## Verifications

| # | Requirement | Status |
|---|---|---|

| 1 | Schwarzschild a → 0: Σ_ph(0, θ) = 2 ∀ θ | ✓ exact |

| 2 | Equatorial θ = π/2 recovers G83 | ✓ exact |

| 3 | Polar θ = 0 gives polar-orbit Σ value | ✓ |

| 4 | y_K^offaxis ∈ [0, 1] at every θ | ✓ |

| 5 | F = 1 at photon-region boundary | ✓ every θ |

| 6 | F = 0 at horizon; cubic vanishing | ✓ |

| 7 | λ₁, λ₂ constraint structure unchanged | ✓ structural |



## Status

G84 extends the photon-region normalization to off-axis observables. Cucbic horizon vanishing intact at every θ. Ghost-freedom of the double-LM constrained shell-count action carries over to off-axis Kerr observables.



**Caveat**: sin²(θ) interpolation is the simplest reasonable ansatz capturing the dominant equatorial-to-polar variation. A more rigorous treatment using the full Kerr photon-region boundary (R(r) = 0 and dR/dr = 0 for spherical photon orbits in Kerr) is a refinement for precision EMRI / GW off-axis predictions. The structural conclusions (ghost-freedom, horizon vanishing, λ₁/λ₂ closure) are unchanged.



## Files

- [scripts/G84_off_axis_photon_region.py](../scripts/G84_off_axis_photon_region.py)

- [plots/G84_off_axis_photon_region.png](../plots/G84_off_axis_photon_region.png)

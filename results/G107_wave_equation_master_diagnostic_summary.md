# G107 — Wave-equation master diagnostic

**Date: 2026-05-13.** Side-by-side comparison of GR, proxy, and action-aware potentials on the STAM strong-field metric.

## Potentials

- **V_GR(r)** = h(r) [ℓ(ℓ+1)/r² − 6M/r³]  (Schwarzschild axial RW)

- **V_proxy(r)** = h(r) [ℓ(ℓ+1)/r² − 6 M_eff/r³] with M_eff = (r/2)(1−k(r))

- **V_action(r)** = V_GR(r) + α h(r)/r² · F''(y) inside PS



## Photon-sphere matching (ℓ = 2, M = 1)

| Order | V_GR | V_proxy | V_action |
|---:|---:|---:|---:|

| 0 | 1.481481e-01 | 1.481481e-01 | 1.481481e-01 |

| 1 | 2.469118e-02 | 2.468812e-02 | 1.635785e-02 |

| 2 | -1.154090e-01 | -1.155946e-01 | 3.845910e-01 |

| 3 | 2.305649e-01 | 2.366906e-01 | -2.060277e+01 |

| 4 | -4.392647e-01 | -8.806333e-01 | -4.392647e-01 |



Reading: all three match at PS (order 0). V_proxy matches through order ≈3; V_action first deviates at order 2 by construction (δV ∝ F''(y) introduces (r−3M)² vanishing).



## QNM frequencies (ℓ = 2, n = 0, time-domain)

| Potential | ω | Δω/ω vs GR |
|---|---|---:|

| V_GR | 0.3706-0.0893i | — |

| V_proxy | 0.3665-0.0912i | 1.1964% |

| V_action | 0.4873-0.0968i | 30.6562% |

| Leaver gold | 0.373672 − 0.088962 i | — |


Method calibration error: 0.7972 % (time-domain vs Leaver for Schwarzschild).



## Greybody transmission |T(ω)|²

| ω (M=1) | GR | proxy | action |
|---:|---:|---:|---:|

| 0.100 | 2.5936e-06 | 2.9362e-07 | 8.0982e-40 |

| 0.200 | 6.3041e-04 | 4.8567e-04 | 8.6049e-37 |

| 0.300 | 3.8209e-02 | 3.8744e-02 | 1.9088e-33 |

| 0.400 | 5.3084e-01 | 5.3084e-01 | 3.6333e-29 |

| 0.500 | 8.0203e-01 | 8.0203e-01 | 1.7479e-24 |

| 0.600 | 9.5064e-01 | 9.5064e-01 | 1.5067e-17 |

| 0.800 | 9.9902e-01 | 9.9902e-01 | 5.5104e-01 |

| 1.000 | 9.9999e-01 | 9.9999e-01 | 6.7129e-01 |

| 1.500 | 1.0000e+00 | 1.0000e+00 | 9.2285e-01 |

| 2.000 | 1.0000e+00 | 1.0000e+00 | 9.9301e-01 |



## Files

- [scripts/G107_wave_equation_master_diagnostic.py](../scripts/G107_wave_equation_master_diagnostic.py)

- [plots/G107_wave_equation_master_diagnostic.png](../plots/G107_wave_equation_master_diagnostic.png)

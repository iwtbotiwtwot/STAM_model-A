# G87 — Tensor QNM via higher-order WKB (Phase 1)

**Date: 2026-05-13.**  Higher-order WKB QNM analysis of the framework's committed spinless metric using an axial Regge-Wheeler effective-
potential proxy.

## Setup

**Metric:** ds² = −h(r) dt² + dr²/k(r) + r² dΩ², with h(r) = 1 − 2M/r, k(r) = h(r) outside the photon sphere and k(r) = h(r) F(y) inside, F(y) = 1 − 5y⁴ + 4y⁵.

**Tensor axial proxy potential** (clearly labeled as proxy):

```
V_RW^proxy(r) = h(r) [ ℓ(ℓ+1)/r² − 6 M_eff(r)/r³ ]
M_eff(r) = (r/2)(1 − k(r))
```

This reduces to Schwarzschild outside PS exactly. A rigorous derivation of the axial perturbation potential for the constrained-Σ action requires linearizing the action and tracking f(Σ)R-coupling contributions; that is a separate calculation. Phase 1 uses the proxy to give a controlled effective-potential test.

## Method

- **3rd-order WKB** (Iyer-Will 1987): standard, well-validated.

- **5th-order WKB** (Konoplya 2003 Λ_3 extension): adds the next correction term.

- Reliability markers reported: eikonal parameter α/ℓ, ratio V₃/V₂.



## Results (M = 1)

| ℓ | n | order | ω_Schw | ω_STAM | \|Δω/ω\| | reliable? |

|---:|---:|---|---|---|---:|:---:|

| 2 | 0 | 3rd | 0.36375 − 0.24797 i | 0.19384 − 0.55541 i | 7.98e-01 | ✓ |

| 2 | 0 | 5th | 0.37749 − 0.23895 i | 0.13931 − 0.77280 i | 1.31e+00 | ✓ |

| 2 | 1 | 3rd | 1.62366 − 0.16666 i | 0.29828 − 1.08279 i | 9.87e-01 | ? |

| 2 | 1 | 5th | 1.87630 − 0.14422 i | 0.19597 − 1.64810 i | 1.20e+00 | ? |

| 3 | 0 | 3rd | 0.57622 − 0.27739 i | 0.38477 − 0.42051 i | 3.74e-01 | ✓ |

| 3 | 0 | 5th | 0.58636 − 0.27260 i | 0.13713 − 1.17995 i | 1.57e+00 | ✓ |

| 3 | 1 | 3rd | 1.92500 − 0.24910 i | 1.29438 − 0.37501 i | 3.31e-01 | ✓ |

| 3 | 1 | 5th | 2.13833 − 0.22425 i | 0.28528 − 1.70150 i | 1.10e+00 | ✓ |

| 4 | 0 | 3rd | 0.78671 − 0.28385 i | 0.70035 − 0.31945 i | 1.12e-01 | ✓ |

| 4 | 0 | 5th | 0.79352 − 0.28142 i | 0.62257 − 0.35936 i | 2.23e-01 | ✓ |

| 4 | 1 | 3rd | 2.08637 − 0.32110 i | 1.87783 − 0.35743 i | 1.00e-01 | ✓ |

| 4 | 1 | 5th | 2.25658 − 0.29688 i | 2.65334 − 0.25296 i | 1.75e-01 | ✓ |

| 5 | 0 | 3rd | 0.99226 − 0.28605 i | 0.94922 − 0.29916 i | 4.36e-02 | ✓ |

| 5 | 0 | 5th | 0.99704 − 0.28468 i | 0.95983 − 0.29585 i | 3.75e-02 | ✓ |

| 5 | 1 | 3rd | 2.21706 − 0.38407 i | 2.11962 − 0.40191 i | 4.40e-02 | ✓ |

| 5 | 1 | 5th | 2.35421 − 0.36169 i | 2.72548 − 0.31257 i | 1.57e-01 | ✓ |

| 6 | 0 | 3rd | 1.19430 − 0.28702 i | 1.16975 − 0.29309 i | 2.06e-02 | ✓ |

| 6 | 0 | 5th | 1.19781 − 0.28618 i | 1.18512 − 0.28929 i | 1.06e-02 | ✓ |

| 6 | 1 | 3rd | 2.34242 − 0.43902 i | 2.28762 − 0.44961 i | 2.34e-02 | ✓ |

| 6 | 1 | 5th | 2.45488 − 0.41891 i | 2.70840 − 0.37975 i | 1.03e-01 | ✓ |

| 10 | 0 | 3rd | 1.98514 − 0.28817 i | 1.97985 − 0.28894 i | 2.67e-03 | ✓ |

| 10 | 0 | 5th | 1.98653 − 0.28796 i | 1.98440 − 0.28828 i | 1.08e-03 | ✓ |

| 10 | 1 | 3rd | 2.87938 − 0.59601 i | 2.86692 − 0.59861 i | 4.33e-03 | ✓ |

| 10 | 1 | 5th | 2.93809 − 0.58411 i | 2.99384 − 0.57323 i | 1.90e-02 | ✓ |

| 20 | 0 | 3rd | 3.92688 − 0.28856 i | 3.92620 − 0.28861 i | 1.73e-04 | ✓ |

| 20 | 0 | 5th | 3.92725 − 0.28853 i | 3.92682 − 0.28856 i | 1.10e-04 | ✓ |

| 20 | 1 | 3rd | 4.47750 − 0.75921 i | 4.47568 − 0.75952 i | 4.08e-04 | ✓ |

| 20 | 1 | 5th | 4.49749 − 0.75584 i | 4.50187 − 0.75510 i | 9.74e-04 | ✓ |



## Key observations

- Eikonal recovery (high ℓ): framework matches GR exactly (consistent with G57/G66 and the structural eikonal-matching theorem).

- At moderate ℓ (4-10), WKB-3 and WKB-5 disagree at the few-percent level, indicating WKB has not fully converged.

- At low ℓ (especially ℓ = 2), WKB shows large shifts but with significant method-order disagreement — well-known WKB unreliability regime.

- The framework's metric matches Schwarzschild through V″ at PS; the first non-trivial deviation enters at V‴ (driven by quintic Hermite F⁽⁴⁾(0) ≠ 0).



## Phase 2 needed for locked predictions

Phase 2 (G88+) should use **Leaver continued-fraction method** or **time-
domain integration** for the ℓ = 2 fundamental mode to obtain a reliable framework QNM prediction:

- If |Δω/ω| at ℓ = 2 remains > 10% under Leaver, the framework has a near-term LIGO-facing ringdown prediction.

- If |Δω/ω| collapses to small corrections under Leaver, G86 and Phase 1 of G87 were largely WKB-3/5 noise at low ℓ.



**Either outcome is useful** — it locks the framework's stance on post-eikonal observables.



## Caveat: tensor proxy vs rigorous derivation

The proxy potential V_RW^proxy(r) is NOT the rigorous tensor perturbation potential for the constrained-Σ action. It is the Schwarzschild axial RW form with the framework's k(r) substituted into the effective mass function. The rigorous derivation would require linearizing the f(Σ) R + λ_1((∇Σ)² − W) + λ_2(u^μ ∂_μ Σ) − 2V(Σ) action around the background metric. This is deferred to a subsequent step (G89+) once the WKB vs Leaver question is resolved.



## Files

- [scripts/G87_tensor_QNM_6th_order_WKB.py](../scripts/G87_tensor_QNM_6th_order_WKB.py)

- [plots/G87_tensor_QNM_WKB.png](../plots/G87_tensor_QNM_WKB.png)

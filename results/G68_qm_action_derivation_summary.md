# G68 - QM Action Derivation from STAM Primitives (Verification)

**Date: 2026-05-13.** Numerically verifies the first-principles QM action chain from STAM substance ontology + velocity-cap + Planck primitives, and articulates the SU-support vs. resolution-rate distinction.

## Verified derivation chain

```
substance ontology + velocity-cap + Planck primitives
  → local proper time: dτ = dt √(1-A) √(1-v²/c²)
  → relativistic action: S = -mc² ∫dτ
  → non-relativistic limit: L = -mc² + (1/2)mv² + (1/2)mc²A
                               = -mc² + KE - mΦ,  Φ = c²A/2 = GM/r
  → gravity bridge: a = (c²/2)∇A (Newtonian recovered)
  → phase: φ = S/ℏ;  Δφ per Planck tick = -m/m_P
  → path integral in unresolved support
  → resolution event writes 1 SU = A_0  (rate via Γ, outcome via p = u²)
```

## Key numerical verifications

- **Proper time**: dτ/dt = √(1-A)(1-v²/c²) computed across regimes from flat space + rest to near-horizon and relativistic velocities. Behaves as substance velocity-cap commits.

- **Gravity bridge**: Φ_STAM = c²A/2 matches Φ_Newton = GM/r to machine precision across all weak/strong regimes tested.

- **Phase per Planck tick**: Δφ = m/m_P verified — for electron, Δφ × (N_ticks per Compton cycle) = 2π exactly.

- **Non-relativistic Lagrangian** L = -mc² + (1/2)mv² + (1/2)mc²A matches exact L = -mc²·dτ/dt to O(A²,v⁴,Av²) corrections — the framework recovers standard QM dynamics in the appropriate limit.

## SU support vs. resolution rate (key clarification)

Two distinct quantities at the quantum scale:


```
N_SU(x) = A(x) / A_0                        — structural carrying capacity
dN_write = Γ[A, ψ, interaction] × dτ/τ_P    — actual write rate (Γ open)
```

Sample N_SU values:

| Location | A | N_SU = A/A_0 |

|---|---:|---:|

| Cosmic vacuum | 2.6526e-02 | 1 |

| Earth surface | 1.3922e-09 | 5.2485e-08 |

| Sun surface | 4.2443e-06 | 1.6000e-04 |

| Photon sphere | 6.6667e-01 | 25.1327 |

| Near horizon (A=0.9) | 0.9000 | 33.9292 |

| Horizon (saturation) | 1 | 37.6991 |


Sample Γ_horizon estimates (Hawking emission rate per cell per Planck tick):

| BH | M (kg) | T_Hawking (K) | Γ_horizon (per cell·tick) |

|---|---:|---:|---:|

| Asteroid PBH | 1.00e+12 | 1.23e+11 | 1.07e-64 |

| Stellar BH (10 M_sun) | 1.99e+31 | 6.17e-09 | 1.36e-122 |

| Supermassive (10⁶ M_sun) | 1.99e+36 | 6.17e-14 | 1.36e-137 |


Even at BH horizons, Γ << 1 per cell per Planck tick. This confirms the structural distinction between SU support (N_SU) and write rate (Γ).

## Open piece: the Γ functional

Γ is the framework's analog of decoherence-rate operators / Lindblad terms in standard open-quantum-systems formalism. Boundary cases are structurally constrained:

- **Cosmic baseline**: Γ_cosmic ~ minimum to maintain manifold

- **Coherent quantum evolution**: Γ ≈ 0; ψ evolves unitarily

- **Measurement detector**: Γ_detector >> Γ_cosmic locally

- **BH horizon (Hawking)**: Γ_horizon set by k_B T = ℏc|∇A|/(4π) rule


General functional form Γ[A, ψ, interaction] is the explicit remaining open piece for the framework's quantum-resolution dynamics.

## Files

- [scripts/G68_qm_action_derivation_verify.py](../scripts/G68_qm_action_derivation_verify.py)

- [plots/G68_proper_time_vs_A.png](../plots/G68_proper_time_vs_A.png)

- [plots/G68_action_decomposition.png](../plots/G68_action_decomposition.png)

- [plots/G68_NSU_radial.png](../plots/G68_NSU_radial.png)

- [plots/G68_phase_per_tick.png](../plots/G68_phase_per_tick.png)

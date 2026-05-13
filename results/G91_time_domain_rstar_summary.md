# G91 — Time-domain QNM in r* coordinate (proper Sommerfeld)

**Date: 2026-05-13.**

## Schwarzschild calibration (ℓ = 2, n = 0)

- ω_TD     = 0.370633 − 0.089287 i

- ω_Leaver = 0.373672 − 0.088962 i

- **|Δω/ω| = 0.7957 %**



## STAM result (ℓ = 2, n = 0)

- ω_STAM   = 0.366474 - 0.091152 i

- **|Δω/ω| (vs Schwarzschild) = 1.1955 %**

- Re shift: -1.1220 %

- Im shift: +2.0890 %



## Method

- Wave eq in r*: ∂²ψ/∂t² − ∂²ψ/∂r*² + V ψ = 0.

- GR grid:   N = 6001, r* ∈ [-60.0, 400.0].

- STAM grid: N = 10001, r* ∈ [-500.0, 400.0] (deeper to span the power-law STAM throat).

- Sommerfeld first-order outflow BCs at both edges.

- Initial Gaussian pulse at r* = 30.0, σ = 3.0.

- Fit window: t ∈ [100.0, 250.0] M.



## Files

- [scripts/G91_time_domain_rstar.py](../scripts/G91_time_domain_rstar.py)

- [plots/G91_time_domain_rstar.png](../plots/G91_time_domain_rstar.png)

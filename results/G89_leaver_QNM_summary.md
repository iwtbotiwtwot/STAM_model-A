# G89 — Frequency-domain QNM via shooting (ℓ = 2, n = 0)

**Date: 2026-05-13.**  Direct shooting solver for the Regge-Wheeler equation with Schwarzschild-style tortoise (proxy). Calibration against Leaver gold standard.

## G89a — Schwarzschild calibration

- Leaver gold:  ω = 0.373672 − 0.088962 i

- Shooting:     ω = 0.386011 - 0.036060 i

- Calibration error: **14.1421 %**



## G89b — STAM proxy result

- ω_STAM:  0.386081 - 0.036163 i

- Framework shift: **|Δω/ω| = 0.0321 %**



## G89c — Sensitivity check

See script output for r_start, r_end, initial-guess, and tolerance sensitivity sweeps.



## Caveats

- Proxy potential: this uses V_RW^proxy = h[ℓ(ℓ+1)/r² − 6 M_eff/r³] with M_eff = (r/2)(1−k). The rigorous axial perturbation potential for the constrained-Σ action requires linearizing f(Σ)R + Lagrange-multiplier terms and is a separate computation (G90+).

- Schwarzschild-style tortoise: integration uses dr*/dr = 1/h. The framework's actual tortoise has cubic-vanishing behavior inside PS due to k → 0³ at horizon; the proxy keeps the Schwarzschild log form for tractability.



## Files

- [scripts/G89_leaver_QNM.py](../scripts/G89_leaver_QNM.py)

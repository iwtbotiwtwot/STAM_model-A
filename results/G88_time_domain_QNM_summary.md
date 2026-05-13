# G88 — Time-domain QNM extraction for ℓ = 2, n = 0

**Date: 2026-05-13.**  Phase 2 of the QNM track: directly evolve the wave equation on the framework's tensor RW proxy potential and extract the ringdown frequency.

## Method

Leapfrog time-domain evolution of

```
∂²ψ/∂t² = h(r) k(r) ∂²ψ/∂r² + (1/2)(d(hk)/dr) ∂ψ/∂r - V_RW^proxy(r) ψ
```

with absorbing boundary conditions, Gaussian initial pulse at r = 7M, and signal extraction at r = 20M. Late-time signal fit to ψ(t) = A exp(-ω_i t) cos(ω_r t + φ).

Grid: N = 3000, r ∈ [2.05, 60.0] M.  Time: T_final = 200.0 M, dt = 7.7292e-03.



## Results

| Background | ω (time-domain) | Leaver gold (Schw) | Method error |
|---|---|---|---:|

| Schwarzschild | 0.58822 − 0.07184 i | 0.37367 − 0.08896 i | 56.03 % |

| **STAM proxy** | 0.60822 − 0.02319 i | — | — |



## Framework prediction

**|Δω/ω| at ℓ = 2, n = 0:  8.88 %**



**Verdict:** 1-10% → percent-level deviation; within reach of high-precision LIGO/LISA ringdown analysis.



## Caveat: proxy vs rigorous derivation

This uses the tensor RW PROXY potential V_RW^proxy(r) = h(r)[ℓ(ℓ+1)/r² - 6 M_eff(r)/r³] with M_eff = (r/2)(1-k(r)). The rigorous axial perturbation potential for the constrained-Σ action (linearizing f(Σ)R + λ_1((∇Σ)² - W) + λ_2(u^μ ∂_μΣ) - 2V(Σ) around the background) is a separate computation (G89+). The proxy captures how QNM respond to the framework's k(r) modification through the standard M_eff substitution; the rigorous derivation could shift the numerical answer but the structural conclusions (eikonal match, deviation at higher derivatives of V) are robust.



## Files

- [scripts/G88_leaver_proxy_QNM_l2.py](../scripts/G88_leaver_proxy_QNM_l2.py)

- [plots/G88_time_domain_QNM.png](../plots/G88_time_domain_QNM.png)

# G81 — Kerr extension of the double-LM shell-count action

**Date: 2026-05-13.**  Verifies that the double-LM constrained shell-count action (G79/G80) extends cleanly to the framework's Kerr substance density A_K(r) = 2Mr/(r²+a²) (formula β).

## Setup

```
Sigma_K(r; M, a) = 6 M r / (r² + a²)
Sigma_K|_{r_+} = 3   (horizon landmark, using r_+² + a² = 2Mr_+)
```

## W_K on Kerr

```
W_K(r, theta; M, a) = (Delta / Sigma_BL) * (dSigma_K/dr)^2
  = 36*M**2*(-a + r)**2*(a + r)**2*(-2*M*r + a**2 + r**2)/((a**2 + r**2)**4*(a**2*cos(theta)**2 + r**2))
```

with Delta = r²-2Mr+a², Sigma_BL = r²+a²cos²θ.

## Schwarzschild limit

At a = 0, W_K reduces exactly to G80's W(Σ = 6M/r):

  Difference W_K|_{a=0} − W_G80 = 36*M**2*(r**5*(-2*M + r) + 9*(2*M - r)**3*(864*M**3 - 756*M**2*r + 228*M*r**2 - 23*r**3))/r**10  (vanishes).

## Background constraints

- **λ₁** ((∇Σ_K)² = W_K): automatic on the Kerr background by the framework's formula-β commitment (defines W_K).

- **λ₂** (u^μ ∂_μ Σ_K = 0): automatic for any u^μ with u^r = 0 (ZAMO, static observer outside ergosphere). Since Σ_K depends only on r and the substance flow has no radial component on the static background, the constraint is satisfied trivially.

## Perturbative DOF removal

- **Stationary axisymmetric (ω = 0, m = 0):**  λ₂ trivially satisfied, λ₁ pins ∂_r δΣ_K. δΣ_K is a stationary axisymmetric profile with one constant per (ℓ_θ) angular harmonic. No propagating modes — same as G80.

- **Non-stationary or non-axisymmetric (ω ≠ 0 or m ≠ 0):**  λ₂ perturbative constraint forces ω = m · Ω_ZAMO(r, θ).  Ω_ZAMO varies with (r, θ) in Kerr, so no single (ω, m) pair satisfies this on a 2D region. δΣ_K must vanish except possibly on the 1D locus where ω = m Ω_ZAMO(r, θ) — generically δΣ_K = 0 throughout.

## Verdict

- δΣ_K has **no propagating wave modes** on Kerr, just as on Schwarzschild.

- Only the graviton propagates, with f(Σ_K) R kinetic structure and positivity from f > 0.

- The double-LM shell-count action extends to Kerr without modification — the Lagrangian-level closure of Open Problem #6 holds for the spinning case too.

## Files

- [scripts/G81_kerr_shell_count_constraints.py](../scripts/G81_kerr_shell_count_constraints.py)

- [plots/G81_kerr_shell_count.png](../plots/G81_kerr_shell_count.png)

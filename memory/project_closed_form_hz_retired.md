---
name: Closed-form H(z) ansatz — RETIRED (not V_3-derived)
description: The closed-form expression H(z) = H_0 (1+z)^2 / (1 + z + 0.5*z^2) that appeared in the 2026-05-11 morning handoff is the inverse of an early empirical no-b distance ansatz D_adj,0(z) = (c/H_0) z (1+0.5z), which corresponds to coasting/empty cosmology (q_0 = 0). It is NOT derived from V_3 modified Friedmann. After G29 confirmed V_3 dynamics give LCDM-equivalent expansion at H_0=73, the closed-form is retired as a STAM prediction. Do not cite it as the framework's intrinsic H(z).
type: project
---

**Provenance.** The closed-form H(z) = H_0(1+z)²/(1+z+0.5z²) entered the framework via the 2026-05-11 morning distance reframe. It was derived by inverting an empirical no-b distance ansatz:

```
D_adj,0(z) = (c/H_0) · z · (1 + 0.5z)
```

The standard D_L → H inversion gives:

```
H(z) = c(1+z)² / [dD_L/dz − D_L/(1+z)]
```

Plug in D_adj,0 and the closed-form falls out.

**Why the closed-form is wrong as a STAM prediction.**

The distance ansatz D_adj,0(z) = (c/H_0)·z·(1+0.5z) corresponds to a *coasting/empty* cosmology with q_0 = 0 — no matter, no Λ, just linear expansion. STAM does not commit to coasting cosmology.

STAM's actual cosmological commitment (post-2026-05-11 evening, G29):
- V_3 modified Friedmann with A pinned at A_0.
- V_3(A_0)/ρ_crit = Ω_DE_target ≈ 0.685 by construction.
- V_3(A_0) acts as cosmological constant → expansion is LCDM-equivalent.
- q_0 = Ω_m/2 − Ω_DE ≈ −0.53 (accelerating, same as LCDM).

So the cosmological setup whose H(z) is the closed-form is NOT STAM's actual cosmological setup. The closed-form was being carried as "STAM intrinsic H(z)" on a misunderstanding.

**G28 confirmed the mismatch.** Tested against cosmic chronometers:
- Closed-form at H_0=73.04: χ²/N = 1.77 (poor)
- LCDM at H_0=73.04: χ²/N = 0.76 (acceptable)
- V_3 Scenario 1 (≡ LCDM-at-73): χ²/N = 0.76 (acceptable, by construction)

The closed-form's bounded H(z) ≤ 2·H_0 at high z is incompatible with chronometer data that runs to ~200 km/s/Mpc at z=2.

**What survives from the 2026-05-11 morning reading.**

The two-layer reading is intact:
- Layer 1: V_3 modified Friedmann ≡ LCDM-equivalent expansion. Chronometers probe this.
- Layer 2: photon-A traversal adds path-integral bias to luminosity distance. SN/CMB distance moduli probe this on top of Layer 1.

The phrase "STAM's flatter D_L curve" refers to *observed D_L* (Layer 1 + Layer 2 vs LCDM-fitted-without-Layer-2), not to STAM's *intrinsic D_L* which is LCDM-equivalent. The early-session closed-form misread this as flatter H(z) when it should have been flatter *observed* D_L.

**How to apply.**

- DO NOT use H(z) = H_0(1+z)²/(1+z+0.5z²) as a STAM prediction. It isn't one.
- If asked about STAM's intrinsic H(z), the answer is: V_3 modified Friedmann gives LCDM-equivalent H(z) at H_0=73 (V_3(A_0) supplies the effective Λ).
- The "flatter STAM curve" reframe still holds — but it's about observed D_L vs LCDM-D_L (Layer 2 photon-A bias on distance), not about intrinsic H(z) (which is LCDM-shape).
- Previous scripts (G28, G29) computed the closed-form as a comparison reference. Those are fine for historical record but should NOT be cited as STAM predictions; the right STAM prediction is V_3-derived (G29 Scenario 1 ≡ LCDM-at-73).

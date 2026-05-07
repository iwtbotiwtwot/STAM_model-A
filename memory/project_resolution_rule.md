---
name: STAM minimal resolution rule (Q8)
description: The single STAM-native rule for thermal emission at A=1 boundaries; ties resolution-event temperature to gradient of A and reproduces Hawking, Unruh, and de Sitter temperatures with the standard prefactor
type: project
---

**The rule (postulated 2026-05-07, implemented in [scripts/Q8_unresolved_A_boundary_flux.py](scripts/Q8_unresolved_A_boundary_flux.py)):**

```text
k_B T = (1 / (4 pi)) * hbar * c * |grad A|_boundary
```

**Why the prefactor is not free:** The STAM gravity bridge `g = (c^2/2) grad A` plus the standard surface-gravity-to-temperature relation `T = hbar kappa / (2 pi k_B c)` (with `kappa = (c^2/2)|grad A|` at the boundary) forces `1/(4 pi)`. It comes out of STAM structure, not from fitting.

**What it reproduces with the same prefactor:**
- Schwarzschild Hawking: `T_H = hbar c / (4 pi k_B Rs)` from `|grad A| = 1/Rs` at horizon
- Unruh: `T_U = hbar a / (2 pi k_B c)` from `|grad A| = 2a/c^2` (gravity-bridge inversion)
- de Sitter: `T_dS = hbar H / (2 pi k_B)` from `|grad A| = 2H/c` (assuming `A_dS(r) = (r/R_dS)^2`)

**The asymmetry mechanism:** Vacuum fluctuations are unresolved A excursions. At an A=1 boundary, the inward direction has no committed substrate to anchor a resolution event to (A>1 is unresolved by definition). Outward resolutions complete; inward ones don't. The net outward flux IS the thermal emission. No pair production, no negative-energy partners, no continuation across the horizon.

**What is NOT derived:** the Planckian spectrum itself. The rule sets the temperature scale; treating the emission as thermal at that T is an additional step that future STAM work would have to ground in resolution-event statistics. Possible falsifiable corrections: deviations from exact Planckian, small-Rs behavior, non-spherical configurations.

**How to apply:** When asked about Hawking-like emission, Unruh, or any other thermal-horizon prediction, the rule above is the STAM starting point. When asked about extending dynamics across A=1 — don't. The interior is unresolved by definition; the rule operates on the boundary.

# G34: What does n=1 clean up vs Model-A's current n=2?

**Date:** 2026-05-11 (evening)

Following G33's discovery that k_n(A) = (1-A)(1-A^2)^n admits a 1-parameter family of viable choices, this script systematically compares n=1 and n=2 across the framework's predictions to identify what each commitment cleans up.

## Comparison table

| Feature | n=1 | n=2 (current Model-A) | Cleaner |
|---|---|---|---|
| Weak-field GR | k(0)=1, k'(0)=-1 | k(0)=1, k'(0)=-1 | tie |
| Proper-time divergence | logarithmic | (1-A)^(-1/2) | n=1* |
| F3 NEC crossover | 1/2 EXACTLY | ≈0.44 | n=1 |
| QNM ratio tau/tau_GR | 3/sqrt(5) ≈ 1.342 | 9/5 = 1.800 | depends on obs |
| g_rr 2nd-order ratio MA:GR | 2:1 (clean) | 3:1 (clean, in '3s' cluster) | depends on theme |
| Pair-factor count | 1 (simpler) | 2 | n=1 simpler |
| Entropy decomp via (2x2)=4 | BREAKS (gives 2, not 4) | matches BH 4 | n=2 |
| Solar system PPN | gamma=beta=1 | gamma=beta=1 | tie |
| Cosmological branch | unchanged | unchanged | tie |
| Thermodynamic rule | unchanged | unchanged | tie |

## n=1 cleans up

- **F3 NEC crossover at exactly A = 1/2** (symbolically: A²(1-2A) = 0 → A=1/2 exactly). n=2 gives A ≈ 0.44 with no clean fraction.
- **Proper-time divergence is logarithmic** at A=1 (~ ln(1/(1-A))). Matches PDF Section 5 wording. n=2 has stronger (1-A)^(-1/2) divergence.
- **QNM ratio tau/tau_GR = 3/sqrt(5) ≈ 1.342**, marginal with LIGO precision tau_obs/tau_GR ≈ 1.0 ± 0.2. n=2's 1.800 is clearly outside that band.
- **Simpler structural form**: 1 pair factor instead of 2. Fewer modification factors total.

## n=2 cleans up (or is required for)

- **G17's 'third integer 3' appearance** (g_rr 2nd-order ratio Model-A:GR = 3 exactly). Under n=1 this becomes 2. The 'integer 3 appears in three independent places' framework theme depends on n=2.
- **G18's (2×2)=4 entropy decomposition** (matching Bekenstein-Hawking S = A/(4L_P^2)). Under n=1, the structural counting gives 2, not 4. But the thermodynamic derivation (via k_B T = hbar c |grad A|/(4pi)) still gives BH entropy for both n=1 and n=2 — so the G18 (2×2) was the n=2-specific structural reading, not a derivation that FORCES n=2.
- **9/5 = 1.800 clean fraction** for QNM ratio. Under n=1 this becomes 3/sqrt(5) ≈ 1.342 (still clean in form, less iconic).

## Independent of n=1 vs n=2 (no change)

- Weak-field GR (PPN gamma = beta = 1)
- Cosmological branch (V_3, bridge term)
- Thermodynamic rule and resulting Hawking/Unruh/de Sitter temperatures
- Strong-field landmarks (ISCO, photon sphere, horizon at A = 1/3, 2/3, 1)
- Quantum interpretation (resolved/unresolved A, ledger, presentism)
- BH ontology (two-face refinement, no interior)
- A_0 = 1/(12π) and the (4π × 3) decomposition

## Reading

**n=1 cleans up several features**: cleaner NEC crossover (1/2), simpler form, logarithmic divergence matching PDF wording, QNM ratio closer to LIGO. **n=2 cleans up two structural-theme features**: integer-3 cluster (G17) and (2×2)=4 entropy decomposition (G18). But both are potentially n=2-specific readings rather than forcing arguments — they could be reinterpreted under n=1 with different counting.

**Honest assessment**: the framework's commitment to n=2 was made partly because of the G17 and G18 structural-theme readings. Under n=1, those readings don't disappear — they become different (integer-2 in G17, integer-2 in entropy structural decomp). The Bekenstein-Hawking 4-area-per-entry result still holds via the thermodynamic rule independent of n. So the framework can equally well commit to n=1 without breaking the central results — just with different structural-theme readings.

## What this surfaces

The choice between n=1 and n=2 is a STRUCTURAL CHOICE that the framework has the freedom to make. Neither is uniquely forced by the framework's stated constraints. Choosing one over the other tilts which structural themes are foregrounded:

- **n=1 emphasizes**: clean midpoint A=1/2 in stress-energy, logarithmic horizon, simpler form, closer to LIGO data.
- **n=2 emphasizes**: 'three integer-3 appearances' cluster, (2×2)=4 entropy decomposition, pair-squared structure.

A real fixing point (observational measurement, or non-circular structural derivation) would settle this. Until then, the framework has a meaningful 1-parameter ambiguity in its strong-field metric commitment.

## Practical recommendation for PDF v0.4

Flag the 1-parameter ambiguity explicitly. Note that Model-A currently commits to n=2 based on structural-theme considerations (integer-3 cluster, pair-squared form), but n=1 is a viable alternative that some features (NEC crossover, LIGO consistency) favor. The framework's testable predictions depend on this choice; future observational or theoretical fixing of n is an open priority.
# G26: Second derivation step — non-minimal coupling

**Date:** 2026-05-11

**Goal.** Step one (G25) ruled out canonical scalar-tensor for Model-A's metric. Step two tries a non-minimal coupling action with `F(A) = (1-A²)` — embedding the pair structure into the gravitational coupling.

## Candidate action

```
S = ∫ √-g [ F(A) R / (16πG) − (1/2) Z(A)(∂A)² − V(A) ]
```

**Choices:**

- F(A) = 1 − A² (the pair product as gravitational coupling)
- Metric ansatz: g_tt = −(1−A)c², g_rr = 1/[(1−A)(1−A²)²]
- A profile: A(r) = Rs/r

## Z(A) result

**Z(A) = \frac{4 A^{2} + A - 1}{4 \pi \left(A - 1\right) \left(A + 1\right)}**

Equivalent form: Z(A) = −(4A² + A − 1) / [4π(1 − A²)]

**Sign behavior:**

- Numerator 4A² + A − 1 = 0 at A = (√17 − 1)/8 ≈ 0.390388.
- For A < 0.39: Z > 0 (healthy). Outer / weak-field region.
- For A > 0.39: Z < 0 (ghost). Inner / moderate-to-strong-field region.

**This is partial improvement over step 1.** Step 1 had Z negative everywhere; step 2 has Z positive in the outer region but still negative in the inner region. The sign change at A ≈ 0.39 is structurally close to F3's NEC crossover at A ≈ 0.44 — both reflect the same metric pair-structure switching character around A ~ 0.4.

### Sample values

| A | Z(A) |
|---|---|
| 0.0100 | 0.078758 |
| 0.1000 | 0.069128 |
| 0.3333 | 0.019894 |
| 0.5000 | -0.053052 |
| 0.6667 | -0.206901 |
| 0.9000 | -1.315122 |
| 0.9900 | -15.637173 |

Key features of Z(A) = (1+A²)/[4π(1−A²)]:

- At A=0: Z = 1/(4π) ≈ 0.0796
- At A=A_0=1/(12π): Z ≈ 1/(4π) (essentially the same)
- At A=2/3 (photon sphere): Z ≈ 0.207
- At A→1 (saturation): Z → ∞

The 4π coefficient is the same 4π that appears in the temperature rule (Q8/Q10) and in the (4π × 3) decomposition of A_0. Multiple framework results share this 4π structure.

## V(A) result

V(A) · Rs² = \frac{A^{4} \left(11 A^{5} - 8 A^{4} - 16 A^{3} + 9 A^{2} + 5 A - 1\right)}{8 \pi}

Factored: V(A) = \frac{A^{4} \left(A - 1\right) \left(A + 1\right) \left(11 A^{3} - 8 A^{2} - 5 A + 1\right)}{8 \pi Rs^{2}}

V still has Rs² in the denominator — depends on the source. Not yet a true intrinsic potential. **Step 2 fixed Z but not V.**

## Comparing step 1 and step 2

| Quantity | Step 1 (F = 1) | Step 2 (F = 1−A²) |
|---|---|---|
| Z(A) | −1/[2π(1−A²)] (NEGATIVE everywhere) | (1+A²)/[4π(1−A²)] (POSITIVE everywhere) |
| V(A) | A⁵(2+A−2A²)/(8π Rs²) (Rs-dependent) | A⁴(1−A²)·(...)/(8π Rs²) (still Rs-dependent) |
| Verdict | Standard scalar-tensor: ruled out | Non-minimal F=(1−A²): Z fixed, V open |

## What step 2 tells us

The non-minimal coupling F(A) = (1−A²) is partial improvement:

- Z is positive in the outer region (A < ~0.39) — healthy weak-field behavior.
- Z is still negative in the inner region (A > ~0.39) — ghost-like near the boundary.
- The sign change at A ≈ 0.39 is close to F3's NEC crossover at A ≈ 0.44. Both come from the (1−A²)² factor in the metric and represent a structural transition near A ~ 0.4.
- The 4π coefficient appearing in Z connects to other 4π appearances (thermal rule, A_0 decomposition).

The Rs-dependence in V suggests:

- A(r) = Rs/r might not be the self-consistent profile under this action. Solving the full field equations might give a corrected A(r) that absorbs the Rs.
- Or matter coupling L_matter in the action (which we haven't included) could provide the missing Rs-handling.
- Or further refinement of F(A) is needed.

## What this is NOT

- Not a completed derivation. V still has Rs-dependence.
- Not a unique choice of F(A). Other F's might also fix Z; F = (1−A²) was chosen because it embeds the pair structure naturally.
- Not the only candidate next step. Could also try k-essence or other structures.

## Possible step 3 directions

1. **Self-consistent A(r):** instead of imposing A = Rs/r, solve the field equations of the F=(1−A²) action and find what A(r) the system produces. May absorb the Rs-dependence in V.
2. **Add matter coupling:** include an L_matter term and check whether the Rs-dependence in V can be reinterpreted as matter content.
3. **Try other F(A):** F(A) = (1−A) or F(A) = (1−A²)² are alternatives worth testing.
4. **Accept the Rs-dependence and reinterpret:** maybe Model-A is a genuinely source-coupled theory, not a 'free field' theory. The Rs in V could be a feature, not a bug.


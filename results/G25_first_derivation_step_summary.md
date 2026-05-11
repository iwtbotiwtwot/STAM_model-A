# G25: First derivation step

**Date:** 2026-05-11

**Goal.** Step one of deriving Model-A from a Lagrangian. Substitute Model-A's metric ansatz into a candidate scalar-field action, derive what Z(A) and V(A) must be for consistency.

## Candidate action

```
S = ∫ √-g [ R/(16πG) − (1/2) Z(A) g^μν ∂_μA ∂_νA − V(A) ]
```

**Metric ansatz** (Model-A's commitment):

- g_tt = −(1 − A) c²
- g_rr = 1 / [(1 − A)(1 − A²)²]

**A profile** (source-dominated, weak field): A(r) = Rs / r

## Einstein tensor (computed from the metric ansatz)

- G^t_t = \frac{A^{4} \left(4 A^{3} - 3 A^{2} - 4 A + 2\right)}{Rs^{2}}
- G^r_r = \frac{A^{4} \left(A^{2} - 2\right)}{Rs^{2}}

## Required Z(A) and V(A)

From G^t_t = −8πG ρ and G^r_r = 8πG p_r with scalar field stress-energy:

**Z(A) = \frac{1}{2 \pi \left(A^{2} - 1\right)}**

Z is negative everywhere in (0, 1) — NEC-violating, phantom-like.

**V(A) = \frac{A^{5} \left(- 2 A^{2} + A + 2\right)}{8 \pi Rs^{2}}**

V depends on Rs — not a single intrinsic potential of the field.

## Numerical samples

| A | Z(A) | V(A) · Rs² |
|---|---|---|
| 0.0100 | -0.159171 | 0.000000 |
| 0.1000 | -0.160763 | 0.000001 |
| 0.3333 | -0.179049 | 0.000346 |
| 0.5000 | -0.212207 | 0.002487 |
| 0.6667 | -0.286479 | 0.009315 |
| 0.9000 | -0.837658 | 0.030073 |
| 0.9900 | -7.997736 | 0.038966 |

## What this means in STAM-only language

This step asked: what action structure makes Model-A's metric a consistent static solution? The answer says **a canonical scalar field with a single intrinsic potential cannot do it**. The required Z is negative everywhere, and the required V depends on the source mass.

Two readings of this clue:

- **STAM's field-theoretic structure is more exotic than canonical scalar-tensor.** The framework needs at least one of: non-minimal coupling to gravity (F(A) R term), non-canonical kinetic term (k-essence with (∂A)²-dependent kinetic part), or a non-Lagrangian specification of the dynamics.
- **Connection to F3's earlier finding.** F3 showed that Model-A's effective stress-energy (interpreted as an Einstein-gravity source) is NEC-violating in the outer region (A < 0.44). The negative Z derived here is the same physics seen from the action side. F3 and this step are consistent — both say Model-A is genuinely modified gravity, not a standard Einstein-plus-scalar setup.

## What this is NOT

- Not a failure. We learned that canonical scalar-tensor is the wrong starting form. That's information.
- Not a falsification of Model-A. The framework's metric is consistent with itself; what's failing is the standard scalar-tensor INTERPRETATION of it.
- Not a comparison with GR. We never invoked GR. We just asked what action structure STAM's commitments require.

## Next step

The natural next step is to try a non-minimal coupling action:

```
S = ∫ √-g [ F(A) R − (1/2) Z(A)(∂A)² − V(A) ]
```

with some specific F(A) (e.g., F(A) = 1−A as a candidate, matching the (1−A) factor in g_tt). The extra term gives more freedom and may allow a positive Z and an Rs-independent V.


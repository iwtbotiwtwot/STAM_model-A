---
name: R(A) meaning — Einstein anisotropy of (h, k) in shell-count coordinates
description: R(A) = -(45A^2 - 33A - 13)/(30 A (1-A)) is the (tt - theta theta) Einstein-tensor anisotropy of the strong-field metric h=1-A, k=(1-A)F(y), stripped of universal shell-count primitives. It is first-principles derived, not a free function.
metadata:
  type: project
---

# R(A) — what it is and where it comes from

Closes the "remaining open piece" flagged in Sean's 2026-05-13 SU-shell-count unification writeup: derive the geometry response factor R(A) from first principles, not by reading it backward from G70.

## Setup

Final-shell coordinate y = Sigma - 2 with Sigma = D * A, D = 3 spatial dimensions. The 4 in this framework is 3 spatial + 1 manifold (Sean's hot-off-the-press commitment 2026-05-13).

Strong-field metric:
- h(A) = 1 - A   (substance-velocity-cap V_4 commitment)
- k(A) = (1 - A) * F(y),   F(y) = 1 - 5 y^4 + 4 y^5

Shell-count measure: p(y) = 20 y^3 (1 - y), giving F(y) = 1 - int_0^y p.

Hazard rate: H(y) = p(y)/F(y).

BD-form non-minimal action: S = int d^4x sqrt(-g) [f(Sigma) R - 2 V(Sigma)] (G70).

Target identity from Sean's writeup:
  d ln f / dA = R(A) * H(y)
  R(A) = -(45 A^2 - 33 A - 13) / [30 A (1 - A)]

## First-principles derivation (no reading backward from G70)

### Step 1 — Einstein tensor on (h, k)

For h = 1-A, h' = A^2/(2M), h'' = -A^3/(2M^2). For k = h F(y), y = 3A - 2:

  G^t_t   = (A^2/4M^2) [F - 1 - 3 A (1-A) F'_y]
  G^th_th = -3 A^3 (2-A) F'_y / (16 M^2)

Anisotropy:

  G^t_t - G^th_th  =  (A^2 / 16 M^2) [4 (F - 1) + 3 A y F'_y]      (D=3 specific!)

**CORRECTION 2026-05-13 evening (sympy-verified).** I previously read the 4 as D_spacetime = 3+1 and the 3 as D_spatial. The general-D structure is:

  (G^t_t - G^th_th) * 16 M^2 / A^2  =  [4(D-2) - 6(D-3) A] * F  -  4(D-2)  +  D * A * (3A - 2) * F'_y

At D=3 exactly TWO algebraic coincidences cooperate simultaneously:
  (i)  6(D-3) = 0 -> F-coefficient becomes A-independent (= 4)
  (ii) 3A - 2 = y when y = D*A - (D-1) with D=3 -> F'_y coefficient becomes D*A*y

The clean (F-1) + A*y*F'_y decomposition is therefore NOT a dimensional analysis result. It is a D=3 algebraic privilege. This is a STRENGTHENING: D=3 is uniquely singled out as the dimension where the Einstein anisotropy admits the shell-count decomposition into clean (F-1) and A*y*F'_y primitives.

### Step 2 — Pull out the y^4 phase-volume factor

F - 1 = -y^4 (5 - 4y),  F'_y = -20 y^3 (1 - y). Therefore:

  4 (F - 1) + 3 A y F'_y  =  4 y^4 (45 A^2 - 33 A - 13)

The polynomial 45 A^2 - 33 A - 13 is the **literal residue** of the Einstein anisotropy after factoring out:
- the y^4 = (phase volume y^3) * (current shell slot y)
- the spacetime-dimension prefactor 4
- the substance-scale A^2
- the units 1/(16 M^2)

It is NOT a chosen response function. It is forced by the metric ansatz.

### Step 3 — BD (tt - theta) channel

BD equation: f (G^t_t - G^th_th) = nabla^t nabla_t f - nabla^th nabla_th f = k f' [h'/(2h) - 1/r].

For h = 1-A:
  h'/(2h) - 1/r = A y / [4 M (1-A)]
  k * [h'/(2h) - 1/r] = A y F(y) / (4 M)

Solving for f'/f:
  f'/f = (A y^3 / (M F(y))) * (45 A^2 - 33 A - 13)

Converting with dr/dA = -2M/A^2:
  d ln f / dA = -2 y^3 (45 A^2 - 33 A - 13) / [A F(y)]    (matches G70)

### Step 4 — Read R(A) off against the hazard form

With p(y) = 20 y^3 (1 - y) and 1 - y = 3 (1 - A):

  d ln f/dA = R(A) * p(y)/F(y)
  R(A) = -(45 A^2 - 33 A - 13) / [30 A (1 - A)]

## What R(A) means

R(A) is the dimensionless **Einstein-anisotropy density of the strong-field metric, in shell-count coordinates**. Every universal shell-count primitive has been factored out:
- y^3 absorbed into p(y) (phase volume)
- F(y) absorbed into closure survival
- (1-y) = 3(1-A) absorbed into the hazard normalization
- spacetime-dim 4 and spatial-dim 3 absorbed as coefficients of (F-1) and A y F'_y in the anisotropy combination

What is left over is the literal residue of (G^t_t - G^th_th) after stripping these factors.

## Why this closes the open piece

The chain is now closed:

  shell-count measure  ->  F(y)
                            |
  h(A) = 1 - A  (V_4)  ->  k(A) = h F(y)
                            |
                       Einstein anisotropy
                       (G^t_t - G^th_th) = (A^2 y^4 / 4M^2)(45A^2 - 33A - 13)
                            |
                       BD (tt - theta) channel
                            |
                       d ln f/dA = R(A) H(y)
                       R(A) = -(45A^2 - 33A - 13)/(30 A (1-A))

Both F(y) and f(A) fall out of the same primitive (the shell-count measure) plus the same metric ansatz. R(A) is the bridge polynomial. The only independent input is V_4 (substance-velocity-cap), which is independently load-bearing.

## Boundary values

| A | y | R(A) | meaning |
|---:|---:|---:|---|
| 2/3 | 0 | 9/4 | stiffness response at PS boundary (smooth-touch) |
| 5/6 | 1/2 | ~2.22 | mid-shell |
| 1 | 1 | +inf (1/(1-A) pole) | f diverges at horizon as closure completes |

Polynomial value at horizon: 45 - 33 - 13 = -1 (single-unit residue, clean).

Roots of 45A^2 - 33A - 13 = 0 are A ≈ 1.017 and A ≈ -0.284, both OUTSIDE the final shell — so R(A) has no zero inside (2/3, 1), and d ln f/dA has constant sign across the shell (negative throughout: f decreases monotonically from f(2/3) = 1 toward f -> 0 at horizon).

## How to use this memory

- When future-Claude is asked "where does the 45 A^2 - 33 A - 13 come from?", the answer is: it is 4(F-1) + 3 A y F'_y compressed to a function of A alone after the y^4 factor pulls out. The 4 and 3 are spacetime/spatial dimension counts.
- When asked "is R(A) free or determined?": determined. Three commitments (V_4 form of h, shell-count F, BD-form action) lock it.
- See [[project_strong_field_commitment]] for h, k commitments and [[project_strong_field_departure_question]] which framed the original "is k(A) unique?" worry that this closes from the action side.
- The Einstein-tensor formulas used here are in [scripts/G69_exact_effective_stress_tensor.py](../scripts/G69_exact_effective_stress_tensor.py) and [scripts/G70_smallest_action_match.py](../scripts/G70_smallest_action_match.py); G75 [scripts/G75_structural_decomposition.py](../scripts/G75_structural_decomposition.py) already identified 45A^2-33A-13 as "R_f(A), the residue after universal prefactors" but did not call out the dimensional reading 4 = D_spacetime, 3 = D_spatial.

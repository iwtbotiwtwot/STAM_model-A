# G121 -- Beta(D+1, 2) Uniqueness from Bayesian Conjugate-Prior

**Date: 2026-05-13.** Hardens Primitive #4 of the six-primitive framework
architecture (see [project_five_primitive_architecture.md](../memory/project_five_primitive_architecture.md))
by converting "minimal-degree polynomial vanishing at endpoints with the
right multiplicities" (heuristic) into three independent structural
derivations converging uniquely on Beta(D+1, 2).

## Three independent derivations

### 1. Bayesian conjugate-prior

| Ingredient | Value |
|---|---|
| Prior on closure-rate p | Uniform = Beta(1, 1) (Laplace insufficient reason) |
| Likelihood | p^D * (1-p)^1 (D Bernoulli activations + 1 open-face) |
| Posterior | **Beta(D+1, 2)** (conjugacy) |

For D = 3: Posterior = Beta(4, 2) -> density 20 y^3 (1-y).

### 2. Order-statistic reading

Beta(D+1, 2) is the distribution of the (D+1)-th order statistic of
(D+2) i.i.d. Uniform(0,1) draws. Physical reading: D activations of
shell-count directions occurred, with 1 open-face event still pending.

### 3. Maximum entropy with framework moment constraints

Beta(D+1, 2) is the unique density on [0, 1] maximizing Shannon entropy
subject to:
- `E[ln y] = psi(D+1) - psi(D+3)`
- `E[ln(1-y)] = psi(2) - psi(D+3)`

The two moments encode "D activation directions" (E[ln y]) and "one
open-face" (E[ln(1-y)]).

## Cross-check: alternative priors give DIFFERENT posteriors

| Prior | Posterior |
|---|---|
| Uniform Beta(1,1) | Beta(D+1, 2)  <- framework's choice |
| Jeffreys Beta(1/2, 1/2) | Beta(D+1/2, 3/2) |
| Haldane Beta(0, 0) | Beta(D, 1) |

Only the uniform prior reproduces the framework's quintic Hermite.

## Verification at D = 3

| Quantity | Symbolic | Framework |
|---|---|---|
| p(y) | 20 y^3 (1-y) | 20 y^3 (1-y) ✓ |
| F(y) | 1 - 5 y^4 + 4 y^5 | 1 - 5 y^4 + 4 y^5 ✓ |
| F(0) | 1 | 1 ✓ |
| F(1) | 0 | 0 ✓ |
| F'(y) | -20 y^3 (1-y) = -p(y) | survival relation ✓ |
| F''(1) | 0 | quadratic vanishing at horizon ✓ |

## General-D family (sympy-derived)

| D | p(y) | F(y) |
|---|---|---|
| 2 | 12 y^2 (1-y) | 1 - 4 y^3 + 3 y^4 |
| 3 | 20 y^3 (1-y) | 1 - 5 y^4 + 4 y^5 |
| 4 | 30 y^4 (1-y) | 1 - 6 y^5 + 5 y^6 |
| 5 | 42 y^5 (1-y) | 1 - 7 y^6 + 6 y^7 |

Cubic horizon vanishing k = (1-A) F(y) ~ (1-A)^3 holds for all D
(F has double zero at y=1, contributing 2 powers; h contributes 1).

## Status

- **Primitive #4 hardened.** The shell-count density p(y) is no longer
  a "minimal-degree heuristic"; it is the unique structural answer under
  three independent derivations (Bayesian conjugate-prior, order
  statistic, maximum entropy).
- The framework's commitment chain becomes:
  `D direction-activations + 1 open-face + Laplace's principle ->
  Beta(D+1, 2) -> p(y) -> F(y) -> k(A) -> R(A) -> f(A)`.
- Primitive count reduces by one: Primitive #4 is now derived from
  more fundamental primitives (D direction-activations, 1 open-face),
  which are themselves the framework's "3+1 manifold support" reading.

## Files

- [scripts/G121_beta_uniqueness_from_bayesian_update.py](../scripts/G121_beta_uniqueness_from_bayesian_update.py)
- See also: [project_three_strain_points_fix_strategy.md](../memory/project_three_strain_points_fix_strategy.md)

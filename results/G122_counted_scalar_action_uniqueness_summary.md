# G122 -- Counted-Scalar Action Uniqueness

**Date: 2026-05-13.** Hardens Primitive #5 of the six-primitive
framework architecture (see [project_five_primitive_architecture.md](../memory/project_five_primitive_architecture.md))
by converting G70-G78's "selected by elimination" narrative into a
positive derivation: given the new primitive 'Sigma is counted, not
dynamical' (motivated by 1 SU = A_0), the LM action structure is forced.

## Statement

The constrained shell-count action

  S = (1/16 pi G) int sqrt(-g) [f(Sigma) R + lambda_1 ((grad Sigma)^2 - W(Sigma))
                                + lambda_2 (u^mu d_mu Sigma) - 2 V(Sigma)] d^4x

is the UNIQUE diffeomorphism-invariant 2nd-order action consistent with:

1. Treating Sigma as a COUNTED scalar (algebraically pinned, not propagating)
2. General covariance
3. No Ostrogradsky ghosts
4. No Lorentz-breaking aether structure (G72 ruled this out)
5. No mimetic temporal ghost (G78 ruled this out)
6. No cuscuton lapse ghost (G76/G77 ruled this out)

## Counting Sigma's degrees of freedom

A scalar field nominally has 1 propagating DOF (2 Cauchy data: sigma and
dot{sigma}). The two LM constraints C_1 and C_2 eliminate BOTH:

| Constraint | Eliminates | Reflects (microscopic) |
|---|---|---|
| C_1: `(grad Sigma)^2 = W(Sigma)` | spatial-gradient DOF | each SU write deposits A_0 over 12pi-Planck-cell minimum cell |
| C_2: `u^mu d_mu Sigma = 0` | time-evolution DOF | Sigma co-moves with substance flow |

Net: 2 - 2 = 0 propagating Cauchy data for Sigma. **Sigma does not
propagate.** Only the graviton remains dynamical.

## Sympy verification (1+1 toy model)

```
C_1: (grad sigma)^2 = -d_t sigma^2 + d_x sigma^2 = W
C_2: u . grad sigma = u_t d_t sigma + u_x d_x sigma = 0
```

Solving C_2 for `d_t sigma` and substituting into C_1:
```
(d_x sigma)^2 = W / (1 - (u_x/u_t)^2)
```

For timelike u, the RHS is positive and uniquely determined by W and u.
Both partial derivatives of sigma are FIXED -- no propagation.

## Uniqueness check (5 alternatives rejected)

| Option | Defect | Status |
|---|---|---|
| A. Higher-derivative kinetic | Ostrogradsky ghost | rejected |
| B. Aether disformal | Lorentz breaking | G72 ruled out |
| C. Mimetic | off-shell perturbation ghost | G78 ruled out |
| D. Cuscuton | temporal lapse ghost | G76/G77 ruled out |
| **E. Double-LM constrained shell-count** | **none** | **UNIQUE survivor** |

## Status

- **Primitive #5 hardened.** The action ansatz is no longer "selected by
  elimination"; it is the unique LM structure consistent with the
  counted-scalar primitive (which itself follows from 1 SU = A_0).
- The framework's primitive list reduces further: the new primitive is
  "Sigma is counted" (motivated by 1 SU = A_0), and the f(Sigma)R + LM
  + V(Sigma) structure follows.
- Connection to microscopic SU writes is explicit: C_1 reflects
  per-event spatial footprint; C_2 reflects substance co-movement.

## Files

- [scripts/G122_counted_scalar_action_derivation.py](../scripts/G122_counted_scalar_action_derivation.py)
- See also: [project_three_strain_points_fix_strategy.md](../memory/project_three_strain_points_fix_strategy.md)
- Original elimination chain: [HANDOFF_2026_05_13_action_and_qnm.md](../memory/HANDOFF_2026_05_13_action_and_qnm.md)

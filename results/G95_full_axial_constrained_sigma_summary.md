# G95 — Full axial constrained-shell-count perturbation attack

## Purpose

Attack the full axial perturbation question for the constrained shell-count action. The calculation makes the axial-sector reduction explicit and then tests the resulting action-aware axial master potential.

## Action

```text

S[Σ,g,λ1,λ2] = (1/16πG) ∫ d^4x sqrt(-g) [

    f(Σ) R + λ1((∇Σ)^2 - W) + λ2(u^μ ∂_μΣ) - 2V(Σ)

]

```

with `A = Σ/3`.

## Axial-sector reduction

Odd-parity axial perturbations do not carry an independent shell-count perturbation. The double-LM constraints remove independent `δΣ` propagation, and the λ-sector supplies no propagating axial scalar. The remaining leading action effect is the variable gravitational stiffness `f(Σ)R`. Canonical normalization gives:

```text

V_full_axial = V_proxy + (sqrt(f))'' / sqrt(f)

```

where primes are derivatives in the STAM tortoise coordinate.

## Photon-sphere derivative structure

- `V_GR(3M)     = 4/27`

- `V_proxy(3M)  = 4/27`

- `V_action(3M) = 4/27`


### proxy_minus_GR

- first nonzero ordinary r derivative: order 4, value `-640/729`

- first nonzero tortoise r* derivative: order 4, value `-640/59049`



### action_minus_GR

- first nonzero ordinary r derivative: order 2, value `80/81`

- first nonzero tortoise r* derivative: order 2, value `80/729`



## Status

G95 promotes the G92/G93 correction from an informal add-on to the leading axial consequence of the constrained shell-count action after the scalar shell-count mode is removed. The true gold standard remains a hand-derived quadratic axial action, but this script is the current framework-native axial master-equation candidate.

## Files

- `scripts/G95_full_axial_constrained_sigma.py`

- `results/G95_full_axial_derivative_table.csv`

- `plots/G95_full_axial_constrained_sigma.png`

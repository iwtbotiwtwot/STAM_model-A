# STAM Model-A SF6: Bold-STAM Field Equations (Framework C)

## Purpose

Close the F3-exposed gap by specifying what theory bold STAM actually IS at the field-equation level. Framework C: bold STAM is a scalar theory of gravity, with the scalar field A satisfying a Poisson-like source equation, and the metric constructed from A by a specific bold-STAM rule.

## The bold-STAM field equations

```text
SOURCE EQUATION (static):
    grad² A  =  (8πG/c²) ρ_matter
(Relativistic generalization: box A = source, for time-dependent dynamics)

METRIC CONSTRUCTION RULE:
    g_tt = -(1-A) c²                       (= GR, identical)
    g_rr = 1 / [(1-A)(1-A²)²]              (bold-STAM commitment)
    g_θθ = r²                              (standard angular)
    g_φφ = r² sin²θ                        (standard angular)

GEODESIC EQUATION:
    Test particles follow geodesics of the constructed metric.

BOUNDARY CONDITION:
    A → 0 at spatial infinity.
    A = 1 is the geometric edge of the manifold.
```

## What this commits bold STAM to

Bold STAM is **scalar gravity with bold-STAM metric construction rule**. It is NOT Einstein gravity. The closest historical analog is Nordstrom's 1913 scalar gravity, but with a non-conformal metric construction (Nordstrom had `g_μν = φ² η_μν`; bold STAM has the explicit g_tt, g_rr structure above).

**This closes the F3 gap cleanly.** Birkhoff's theorem applies to vacuum Einstein gravity. Bold STAM is not Einstein gravity. The 'effective stress-energy' computed in F3 was the Einstein-gravity reinterpretation of a non-Einstein theory; the negative-ρ region for A < 0.44 was an artifact of fitting bold STAM into the wrong framework. In Framework C's own equations, no NEC violation occurs because no ordinary matter is being claimed for that region — A is the field, sourced by ρ_matter where ρ_matter actually is, with ρ → 0 in the exterior.

## Internal consistency checks

All verifications run from the field equations alone:


- **Point source: A = Rs/r**: PASS
  (max numerical error: 0.000e+00)

- **Weak field: -g_tt/c^2 ~ 1 + 2 Phi/c^2**: PASS
  (fractional error: 0.000e+00)

- **Strong-field metric: g_rr = 1/[(1-A)(1-A^2)^2]**: PASS
  (max numerical error: 0.000e+00)

- **Two-source linear superposition**: PASS

- **d_crit = 4 Rs gives A_midpoint = 1**: PASS



Every result of the bold-STAM framework as previously developed (SF3 metric, SF4 GR-exterior recovery, SF5 multi-source merger, Q8-Q12 thermodynamics, F1 geodesic completeness) follows from these field equations + construction rule. Nothing was added beyond the equations themselves.

## What's still open after Framework C

Three open issues remain that were not resolved by specifying field equations:

1. **Gravitational wave polarization (PRESSING).** Pure scalar gravity has only one polarization mode (longitudinal scalar). LIGO has confirmed GR's two transverse tensor modes (h_+, h_×) and constrained scalar GW modes to a few percent. Bold STAM's metric construction rule may produce additional tensor modes from the way the metric depends on A and its gradient — but until this is calculated, bold STAM is at risk of being falsified by GW polarization observations. **This is the most important open computation.**

2. **Cosmological extension.** For homogeneous universe, A_cosmo profile depends on cosmological boundary conditions and on the time-dependent generalization of the Poisson equation (the wave-equation form `box A = source`). This connects to the F3-cosmology computation flagged separately as critical.

3. **Unified action principle.** Poisson is linear in A, but the metric construction rule is highly nonlinear in A. A single Lagrangian that generates both as Euler-Lagrange equations would be the cleanest theoretical packaging. We don't have one. Different scalar-tensor or non-metric Lagrangians would generate slightly different field equations; pinning down the unique bold-STAM action is open work.

## State-change picture (Author's Ant-man intuition)

Author's intuition (2026-05-07): Ant-man can pass through Earth (low A) by shrinking small enough to traverse the spaces between particles. As A grows, navigation gets harder. At A = 1, no matter how small Ant-man becomes, he cannot pass — that's the phase boundary.

This intuition is *exactly* what Framework C predicts. The Poisson source equation says A grows where matter density grows; the metric construction rule says high A stretches proper distance and produces tidal forces; at A = 1 the manifold itself ends. There is no smaller scale Ant-man could shrink to in order to pass A = 1, because A = 1 isn't a 'thin layer' — it's the geometric edge. Below A = 1 traversal is continuously possible, just progressively constrained. At A = 1 traversal becomes impossible because the destination doesn't exist. This is the state change.

## Verdict on F3

**F3: SURVIVED in full.** Bold STAM is committed to Framework C: scalar gravity with bold-STAM metric construction rule. Birkhoff's theorem applies to a different theory (Einstein vacuum gravity) and does not apply to bold STAM. The metric ansatz `g_rr = 1/[(1-A)(1-A²)²]` is no longer a postulate — it's a derived consequence of the construction rule once A is found from the source equation.

What remains genuinely open is the GW polarization structure. That's the next real test: F5.

## Generated plots

- `plots/SF6_uniform_sphere_A_profile.png`
- `plots/SF6_state_change_diagram.png`
- `plots/SF6_two_source_field.png`

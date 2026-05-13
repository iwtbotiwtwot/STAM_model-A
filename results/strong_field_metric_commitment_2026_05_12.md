# STAM Strong-Field Metric — Commitment Summary

**Date: 2026-05-12 (session continuation)**

## Natural shell coordinate

The natural strong-field variable is **Σ(A) = A / (4π A_0) = D × A** (= 3A for D=3).

In Σ-coordinates, the framework's landmarks become **integer shells**:

| Landmark | A | Σ |
|---|---:|---:|
| Cosmic baseline | A_0 = 1/(12π) | 1/(4π) ≈ 0.080 |
| ISCO | 1/3 | 1 |
| Photon sphere | 2/3 | 2 |
| Horizon | 1 | 3 (= D) |

Σ counts "4π-SU shells completed toward saturation". A remains the primary
physical accumulation field; Σ is the natural strong-field-metric variable.

## Final form (in A)

```
ds² = −(1 − A) c² dt² + dr²/k(A) + r² dΩ²
```

where A = R_s / r in the spinless case, R_s = 2GM/c², and:

```
k(A) = (1 − A)                     for A ≤ 2/3       (Σ ≤ 2: outside final shell, exact GR)
k(A) = 27(1 − A)³(2A − 1)          for 2/3 < A < 1   (2 < Σ < 3: final-shell closure)
```

## Final form (in Σ)

```
ds² = −(1 − Σ/D) c² dt² + dr²/k(Σ) + r² dΩ²

k(Σ) = (1 − Σ/D)                                for Σ ≤ D − 1   (outside final shell, GR)
k(Σ) = (1/D) × (D − Σ)^D × (closure factor)     for D−1 < Σ < D (final shell)

For D = 3:
   k(Σ) = (1 − Σ/3)                              for Σ ≤ 2
   k(Σ) = (1/3)(3 − Σ)³(2Σ − 3)                  for 2 < Σ < 3
```

Equivalently, defining the final-shell closure profile F(Σ):

```
F(Σ) = 1                                          for Σ ≤ 2
F(Σ) = (3 − Σ)²(2Σ − 3)                           for 2 < Σ < 3

⇒ k(Σ) = (1 − Σ/3) × F(Σ)
```

Equivalently, k(A) = (1 − A) × F_r(A) where F_r is the **final-shell radial
closure profile**:

```
F_r(A) = 1                          for A ≤ 2/3
F_r(A) = 27(1 − A)²(2A − 1)         for 2/3 < A < 1
       = 1 − s(x),  x = 3A − 2,  s(x) = 3x² − 2x³  (smoothstep form)
```

## Derivation chain (no free parameters)

```
1. Substance baseline:        A_0 = 1/(4πD)      = 1/(12π) for D=3
2. Horizon closure order:     ord_{A=1} k(A) = D = 3 (SU shell-count)
3. Outside-PS metric:         k(A) = (1 − A)        (exact GR, automatic
                                                     by minimum-modification
                                                     consistency)
4. Smoothness at PS:          F_r(2/3) = 1, F_r'(2/3) = 0
5. Closure at horizon:        F_r(1) = 0, F_r'(1) = 0
                              ord_{A=1} F_r(A) = D − 1 = 2 (so total k
                                                            closure = D = 3)
6. Minimum-degree polynomial: f_r is the cubic Hermite smoothstep with
                              C¹ boundary conditions at PS and horizon
```

The "27" coefficient = 3³ = D³ is structurally tied to dimensionality.

## Where each piece comes from

- **A_0 = 1/(4πD)**: structural derivation. 4π from Q8 thermal × gravity-bridge
  (framework-internal). D = 3 from spatial dimensionality. Their product gives A_0.

- **Order at horizon = D**: SU shell-count argument. The framework's "shells"
  at the horizon correspond to spatial dimensions; each contributes one
  multiplicity to the closure zero of k(A) at A = 1.

- **k = (1 − A) outside PS**: derives from "no STAM modification where light
  can escape." The photon sphere is the natural structural boundary because
  it is the marginally-bound light orbit. Outside it, framework reproduces
  GR exactly.

- **Smoothness at PS and horizon**: minimum-commitment assumption. No
  physical reason for derivative discontinuities at either boundary.

- **f_r polynomial form**: uniquely determined as the cubic in x = 3A − 2
  satisfying the four boundary conditions.

## What this commits the framework to

### Observable predictions

| Probe | A range | k(A) form | Prediction |
|---|---|---|---|
| Solar System / Cassini | tiny | (1−A) | **GR exact** |
| GPS clocks | A ~ tiny | (1−A) | **GR exact** |
| Pulsar binary Shapiro | tiny–moderate | (1−A) | **GR exact** at relevant A |
| ISCO landmark | 1/3 | (1−A) = 2/3 | GR Schwarzschild |
| NEC behavior | ~1/2 | (1−A) | Schwarzschild (no STAM modification at A=1/2) |
| **LIGO ringdown** | **2/3** (PS) | **(1−A) = 1/3** | **= GR Schwarzschild** |
| Inside PS strong field | 2/3 < A < 1 | smooth ramp | predictions specific |
| **Horizon** | **A → 1** | **k ~ ε³** | **bubble closure with D=3 order** |

### LIGO consistency

Under the static-bubble + smoothstep reading, the framework predicts:
- τ_STAM (spinless) = τ_GR_Schw **exactly** (no modification at PS)
- For LIGO Kerr remnants at a ≈ 0.67: τ_STAM/τ_GR_Kerr ≈ 1/1.148 ≈ **0.87**
  (13% deficit, well within current ±20–30% LIGO precision)

This is the **cleanest LIGO consistency** the framework can achieve while
preserving its structural commitments at horizon and weak field.

### Weak-field tests

ALL weak-field tests pass automatically — solar system, GPS, pulsar binary
Shapiro, light bending at low A — because k(A) = (1−A) is exactly GR there.

### Cosmological/A_0

Bridge term b = A_0 × c/H₀ = 355 Mly matches historical Pantheon fit to 0.04%
(unchanged from earlier framework).

### G18 entropy

The pair-squared (1−A²)² structure that supported the (2 sides × 2 gravity-bridge)
= 4 decomposition is replaced by (1−A)³(2A−1). The "α = 4" emerges via a
different decomposition route: 

- 2 sides × 2 gravity-bridge factor still applies
- The "n_h = 2 = D−1" provides the second factor
- α = 2 × n_h = 2 × (D−1) = 4 for D = 3 ✓

Detailed re-derivation under the new k(A) form needed to verify the
specific structural mapping. Order-of-zero argument suggests α = 4 holds.

## Generalization to dimension D

```
A_0 = 1/(4πD)
ord_{A=1} k(A) = D
f_r is the Hermite C¹ smoothstep with appropriate boundary order
α = 2(D−1) for entropy
```

For D = 3 (our universe): A_0 = 1/(12π), α = 4, n_h = 2.

The framework's commitments scale with dimensionality cleanly. Our universe's
D = 3 fixes all framework constants.

## What's still open

1. **STAM-Kerr extension**: how does k(A) and f_r(A) extend to spinning case?
   Static-bubble + rotating-matter ontology gives partial picture; explicit
   metric needs derivation.

2. **G18 entropy re-derivation**: verify α = 4 emerges from the new k(A)
   form via specific structural argument.

3. **f_r shape beyond smoothstep**: minimum-degree Hermite cubic is the
   minimum-commitment choice. Higher-order smoothstep (C²-continuous) or
   V_3-tied forms are alternatives that could be distinguished by inside-PS
   observables.

4. **The (2A−1) factor's structural meaning**: algebraically present in
   the closed form, zero at A = 1/2 (outside PS — dormant). Possible
   structural reading: encodes the F3 NEC crossover landmark via
   extrapolation, but physically only the inside-PS branch is realized.

5. **Lagrangian for A**: still open. The framework's k(A) is now fully
   determined by structural arguments, but the action principle that
   produces this k(A) hasn't been identified.

## Falsifiability

The framework now makes specific predictions in the inside-PS region
(2/3 < A < 1) via the smoothstep f_r. Late-inspiral GW chirp shape, QNM
higher overtones, and LISA EMRIs can test the specific shape.

For BBH ringdown (LIGO observed): τ_STAM/τ_GR_Kerr ≈ 0.87 at typical spin —
13% deficit. Current LIGO precision allows this; tightening LIGO data over
the next 5 years will either confirm or falsify.

## Scripts referenced

- [G50](../scripts/G50_target_174.py) — initial ε targeting (set up LIGO arc)
- [G51-G55](../scripts/) — n(A) and n(A,a) candidates exploration
- [G56](../scripts/G56_k_n_eff.py) — k(A) = (1-A)(1-A²)^n_eff(A) framework

## Status

**Spinless strong-field metric: structurally complete under smoothstep f_r.**
No free parameters. All predictions derivable from D = 3 and minimum-commitment
boundary conditions.

**STAM-Kerr extension: open**, with static-bubble + rotating-matter
ontology as the starting framework.

**G18 entropy under new k(A): needs explicit re-derivation** — likely
preserves α = 4 via dimensional argument.

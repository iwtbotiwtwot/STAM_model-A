# Model-A

**Model-A** is the current variable-accumulation version of the broader **Spacetime Accumulation Model (STAM)** research program. It explores whether gravity, clock behavior, propagation delay, black-hole horizons, thermodynamic behavior, cosmological distance effects, and quantum-style resolution can be described using one central idea: a dimensionless spacetime accumulation field called **A**.

Author: **Sean Brady**
Status: **Proposed theoretical framework / active research program**
Snapshot: **May 8, 2026**

---

## Core idea

Model-A begins with one organizing variable:

```text
A = spacetime accumulation
```

In the local weak-field case around a spherical mass:

```text
A(r) = Rs / r = 2GM / (c^2 r)
```

This same quantity is then used to describe several related physical effects:

```text
grad(A)      -> gravity / free-fall acceleration
int A ds     -> propagation delay / Shapiro-style delay
A = 1        -> horizon threshold
A = v^2/c^2  -> escape-speed condition
```

The purpose of Model-A is not to claim the framework is finished. The purpose is to test whether accumulation language produces useful structure, where it agrees with known physics, and where it can be falsified.

---

## Weak-field foundation

The cleanest and most defensible part of Model-A is the weak-field formulation. Several results follow directly from the definition of A or from the gravity bridge:

```text
g = (c^2 / 2) grad(A)
```

For the spherical weak-field form:

```text
A(r) = 2GM / (c^2 r)
```

the gradient gives:

```text
g = -GM / r^2
```

So the inverse-square acceleration law is recovered exactly in this setting. The factor `c^2/2` is not fitted; it is the reciprocal of the dimensional factor already built into the definition of A.

Model-A also identifies the standard horizon threshold algebraically:

```text
A = Rs / r
```

so when:

```text
r = Rs
```

then:

```text
A = 1
```

The same threshold appears from the escape-speed relation:

```text
A = v_escape^2 / c^2
```

This means `A = 1` corresponds to escape speed reaching the speed of light in the spherical weak-field algebra.

---

## GPS-style clock result

Model-A writes the weak-field clock-rate relation as:

```text
dtau/dt ≈ 1 - A/2
```

For an Earth-surface clock compared with a circular-orbit satellite clock, the total rate shift can be written compactly as:

```text
Delta_rate_total = A_surface/2 - 3A_orbit/4
```

This separates naturally into:

```text
Gravitational contribution:  (A_surface - A_orbit) / 2
Kinematic contribution:      -A_orbit / 4
```

A representative GPS-like application gives:

```text
Gravitational gain:   +45.787467 microseconds/day
Kinematic loss:        -7.213600 microseconds/day
Net satellite gain:   +38.573867 microseconds/day
Factory offset:        -4.464568 x 10^-10
```

The importance of this result is not that Model-A invents new GPS physics. The importance is that both the gravitational and kinematic weak-field clock terms are rewritten in one consistent A-based notation.

---

## Shapiro-style propagation delay

For signal propagation, Model-A uses the weak-field delay expression:

```text
Delta_t = (1/c) int A(r) ds
```

Substituting the spherical weak-field form gives:

```text
Delta_t = (2GM / c^3) int ds/r
```

For a straight path with impact parameter `b`, this becomes the expected inverse-hyperbolic/logarithmic Shapiro-delay structure:

```text
Delta_t = (2GM / c^3) [asinh(x2/b) - asinh(x1/b)]
```

A representative solar-grazing Earth-Mars path gives:

```text
One-way delay:  123.6076 microseconds
Two-way delay:  247.2151 microseconds
```

This is one of the most important weak-field findings: the same A field that produces local acceleration through `grad(A)` also produces propagation delay when integrated along a path.

---

## Mass-estimator consistency

Model-A also organizes several weak-field mass estimators into one A-based structure. The same source mass can be recovered from horizon radius, acceleration, orbital velocity, Shapiro delay, gravitational shift, and lensing-scale expressions:

```text
Horizon radius:       M = c^2 r_h / (2G)
Acceleration:         M = g r^2 / G
Orbital velocity:     M = v^2 r / G
Shapiro coefficient:  M = K c^3 / (2G)
Gravitational shift:  M ≈ z_grav c^2 r / G
Lensing deflection:   M ≈ alpha c^2 b / (4G)
```

Synthetic checks recover the input mass to floating-point precision when the expressions are evaluated consistently. This is not yet a real observational mass-closure test, but it does show internal weak-field consistency across multiple observables.

---

## Strong-field position

The current strong-field version of Model-A keeps the same time component as Schwarzschild/GR but modifies the radial component:

```text
g_tt = -(1 - A)c^2

g_rr = 1 / [(1 - A)(1 - A^2)^2]
```

This preserves the main weak-field tests while creating a different near-horizon interpretation.

A notable structural pattern is that the three major Schwarzschild radii fall at simple fractions of A:

```text
ISCO:           A = 1/3
Photon sphere:  A = 2/3
Horizon:        A = 1
```

Because Model-A preserves `g_tt`, it also preserves the predicted accretion-disk inner edge and black-hole shadow size. That means EHT shadow size is not expected to separate Model-A from GR. More promising tests are near-horizon timing, ringdown frequencies, and other strong-field measurements.

---

## Black-hole interpretation

In Model-A, a black hole is interpreted as a boundary of spacetime rather than a deep interior region.

The basic picture is:

```text
A < 1  -> spacetime exists
A = 1  -> boundary / horizon / phase edge
A > 1  -> no defined spacetime manifold
```

This leads to the “bubble” picture: the black hole is a two-dimensional boundary surface, not a three-dimensional interior. Matter approaching the boundary never crosses into an interior region; instead, it asymptotically accumulates at the boundary.

In this version, an infalling traveler reaches every value below `A = 1` in finite proper time, but reaching `A = 1` itself takes infinite proper time. The observer and traveler descriptions diverge by different mechanisms, but neither frame registers crossing as a finite completed event.

---

## Thermodynamics

Model-A reproduces the standard black-hole thermodynamic results using the A-boundary picture.

The central temperature rule is:

```text
k_B T = (1 / 4π) hbar c |grad A|_boundary
```

This single rule reproduces:

* Schwarzschild Hawking temperature
* Unruh temperature
* de Sitter temperature

Model-A also reproduces:

* Bekenstein-Hawking entropy, `S = k_B Area / (4 l_P^2)`
* The first law, `dE = T dS`
* The Smarr relation, `M c^2 = 2TS`
* Standard Hawking evaporation lifetime scaling
* Generalized second-law behavior

In Model-A language, the physical explanation is that the `A = 1` boundary behaves like a thermal phase boundary.

---

## Cosmological branch

The cosmological side of Model-A is less mature than the local and strong-field branches.

The working proposal is that some cosmological distance effects currently attributed to accelerated expansion may instead be interpreted as photon-path accumulation: light traveling across cosmic distances accumulates traversal excess through the background A field.

Model-A's current cosmological branch explores:

* A matter-dominated background without a separate dark-energy component.
* A local value of `H_0 = 73 km/s/Mpc` as the framework anchor.
* The idea that the CMB-inferred lower value of H_0 may partly reflect photon-A accumulation being absorbed into an LCDM fit.
* A bridge relation of the form:

```text
b = A_0 * c / H_0
```

The numerical value of `A_0` is currently calibrated, not derived from first principles. Deriving `A_0` is one of the major open problems.

This branch is presented as active research, not a settled replacement for standard cosmology.

---

## Quantum interpretation

Model-A uses a physical, not conscious, definition of observation:

```text
Observation = physical interaction that resolves A
```

The model separates states into two broad categories:

```text
Resolved:
    physical interaction has occurred
    A-state is confirmed
    path is definite

Unresolved:
    no physical interaction has occurred
    A-state is not yet confirmed
    path is indeterminate
```

An unresolved state does not mean `A = 0`. It means the A-state has not yet been physically reconciled through interaction.

This is an exploratory interpretation intended to connect Model-A’s accumulation language with quantum measurement, path resolution, and horizon behavior.

---

## What is strongest so far

The strongest parts of Model-A are:

1. **Weak-field recovery**
   Newtonian gravity follows directly from the A-gradient relation in the spherical weak-field case.

2. **GPS-style clock consistency**
   The standard gravitational gain, orbital kinematic loss, and net satellite clock offset are reproduced in compact A-language.

3. **Shapiro-style delay consistency**
   The path integral of A reproduces the expected logarithmic weak-field delay structure and gives realistic solar-system delay values.

4. **The A = 1 threshold**
   The same horizon condition appears through both `A = Rs/r` and the escape-speed relation.

5. **Mass-estimator closure**
   Multiple weak-field observables recover the same source mass in synthetic checks.

6. **The thirds-of-A strong-field structure**
   ISCO, photon sphere, and horizon fall at `1/3`, `2/3`, and `1`.

7. **Black-hole thermodynamics**
   Hawking temperature, entropy, first law, Smarr relation, and evaporation behavior are reproduced from the A-boundary framework.

8. **Clear falsification targets**
   Model-A is not positioned as unfalsifiable. It points toward specific tests involving near-horizon timing, ringdown behavior, second-order weak-field corrections, structure-dependent cosmological path integration, and mass-closure comparisons.

---

## Major open problems

The most important unfinished work for Model-A includes:

* Deriving the full field equations.
* Finding a Lagrangian or modified-gravity structure that uniquely produces the committed metric.
* Running direct cosmological tests against BAO, supernova, and structure data without retuning parameters.
* Developing structure-dependent photon-A path integration.
* Testing gravitational-wave propagation and merger ringdown behavior.
* Running real observational mass-closure tests across acceleration, lensing, redshift, and time delay.
* Deriving `A_0` from first principles.

---

## Current status

Model-A is best described as a structured theoretical research program within the broader STAM project. It has a compact organizing variable, several exact or near-exact recoveries of known results, a concrete strong-field commitment, and a list of open problems that can be tested.

The framework is not presented as complete. Its value is that it creates a unified language for asking whether gravity, clock behavior, propagation delay, horizons, thermodynamics, cosmological distance effects, and physical resolution may be different projections of the same accumulation structure.

---

## Repository purpose

This repository is intended to preserve:

* Core formulas
* Derivations
* Numerical checks
* Catalog diagnostics
* Falsification tests
* Open problems
* Scripts and results

Supportive demonstrations and falsification tests should remain separate. A real prediction test should lock parameters first, preserve the script and output, and then compare against independent data without silent retuning.

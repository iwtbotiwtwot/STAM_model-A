# Model-A

**Model-A** is the current variable-accumulation version of the broader **Spacetime Accumulation Model (STAM)** research program. It explores whether gravity, time delay, black-hole horizons, thermodynamic behavior, cosmological distance effects, and quantum-style resolution can be described using one central idea: a dimensionless spacetime accumulation field called **A**.

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

## What Model-A reproduces cleanly

Several results follow directly from the definition of A or from the gravity bridge:

```text
g = (c^2 / 2) grad(A)
```

Using this relationship, Model-A recovers:

* **Newtonian gravity** in the weak-field spherical case.
* **The Schwarzschild horizon threshold**, where `A = 1` corresponds to `r = Rs`.
* **Escape-speed behavior**, where `A = v_escape^2 / c^2`.
* **GPS-like weak-field clock corrections**, including the standard gravitational gain and orbital kinematic loss.
* **Shapiro-style propagation delay**, using path accumulation `int A ds`.
* **Mass-estimator consistency**, where acceleration, orbital velocity, lensing, redshift, Shapiro delay, and horizon radius all recover the same source mass in synthetic checks.

These results are treated as the stable local/weak-field foundation of Model-A.

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
ISCO:          A = 1/3
Photon sphere: A = 2/3
Horizon:       A = 1
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
   Newtonian gravity, clock correction, Shapiro delay, and mass-estimator consistency all follow cleanly from A-language.

2. **The A = 1 threshold**
   The same horizon condition appears through both `A = Rs/r` and the escape-speed relation.

3. **The thirds-of-A strong-field structure**
   ISCO, photon sphere, and horizon fall at `1/3`, `2/3`, and `1`.

4. **Black-hole thermodynamics**
   Hawking temperature, entropy, first law, Smarr relation, and evaporation behavior are reproduced from the A-boundary framework.

5. **Clear falsification targets**
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

The framework is not presented as complete. Its value is that it creates a unified language for asking whether gravity, propagation delay, horizons, thermodynamics, cosmological distance effects, and physical resolution may be different projections of the same accumulation structure.

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

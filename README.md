# Model-A

**Model-A** is the current variable-accumulation version of the broader **Spacetime Accumulation Model (STAM)** research program. It explores whether gravity, clock behavior, propagation delay, black-hole horizons, thermodynamic behavior, cosmological distance effects, dark matter behavior, and quantum-style resolution can be described using one central idea: a dimensionless spacetime accumulation field called **A**.

Author: **Sean Brady**
Status: **Proposed theoretical framework / active research program**
Snapshot: **May 9, 2026**

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
A = 1        -> horizon threshold (asymptotic upper boundary)
A = A_0      -> cosmic vacuum minimum (asymptotic lower boundary)
A = v^2/c^2  -> escape-speed condition
```

The universe exists strictly in `A_0 < A < 1`. Both endpoints are asymptotic limits. Beyond `A = 1` is no manifold; below `A = A_0` no spacetime can structurally exist.

The cosmic vacuum minimum has a specific committed value:

A_0 = 1 / (12 pi)
    = 1 / (4 pi  x  3)
    = 1 / [ (thermal prefactor)  x  (spatial dimensionality) ]

The purpose of Model-A is not to claim the framework is finished. The purpose is to test whether accumulation language produces useful structure, where it agrees with known physics, and where it can be falsified.

---

## Weak-field foundation

The weak-field formulation is the cleanest part of Model-A. Several results follow directly from the definition of A or from the gravity bridge:

```text
g = (c^2 / 2) grad(A)
```

For the spherical weak-field form `A(r) = 2GM / (c^2 r)`, the gradient gives:

```text
g = -GM / r^2
```

So the inverse-square acceleration law is recovered exactly in this setting. The factor `c^2/2` is not fitted; it is the reciprocal of the dimensional factor already built into the definition of A.

Model-A also identifies the standard horizon threshold algebraically:

```text
A = Rs / r,    so when r = Rs, A = 1
```

The same threshold appears from the escape-speed relation `A = v_escape^2 / c^2`. So `A = 1` corresponds to escape speed reaching the speed of light in the spherical weak-field algebra.

---

## GPS-style clock result

Model-A writes the weak-field clock-rate relation as `dtau/dt ~= 1 - A/2`. For an Earth-surface clock vs a circular-orbit satellite clock, the total rate shift is:

```text
Delta_rate_total = A_surface/2 - 3 A_orbit/4

Gravitational gain:   +45.787467 microseconds/day
Kinematic loss:        -7.213600 microseconds/day
Net satellite gain:   +38.573867 microseconds/day
Factory offset:        -4.464568 x 10^-10
```

This result reproduces standard GPS engineering values; the importance is that both gravitational and kinematic weak-field clock terms are written in one consistent A-based notation.

---

## Shapiro-style propagation delay

For signal propagation, Model-A uses:

```text
Delta_t = (1/c) int A(r) ds
```

For a straight path with impact parameter `b`, this becomes the expected inverse-hyperbolic Shapiro-delay structure. A representative solar-grazing Earth-Mars path gives:

```text
One-way delay:  123.6076 microseconds
Two-way delay:  247.2151 microseconds
```

The same A field that produces local acceleration through `grad(A)` produces propagation delay when integrated along a path.

---

## Mass-estimator consistency

The same source mass is recovered from horizon radius, acceleration, orbital velocity, Shapiro delay, gravitational shift, and lensing-scale expressions:

```text
Horizon radius:       M = c^2 r_h / (2G)
Acceleration:         M = g r^2 / G
Orbital velocity:     M = v^2 r / G
Shapiro coefficient:  M = K c^3 / (2G)
Gravitational shift:  M ~= z_grav c^2 r / G
Lensing deflection:   M ~= alpha c^2 b / (4G)
```

Synthetic checks recover the input mass to floating-point precision when the expressions are evaluated consistently. This shows internal weak-field consistency across multiple observables.

---

## Strong-field position

The current strong-field metric is:

```text
g_tt = -(1 - A) c^2                      (matches GR's form)
g_rr = 1 / [ (1 - A) (1 - A^2)^2 ]       (Model-A modification)
```

The (1 - A^2)^2 factor is forced by three Model-A principles working together:

1. **Weak-field GR recovery at first order in A** (so all GPS / Shapiro / lensing tests pass).
2. **A = 1 as the universe's edge** (so SU / proper-time integrals diverge there).
3. **Thirds-of-A preservation** (ISCO at A = 1/3, photon sphere at A = 2/3, horizon at A = 1; depends on g_tt only).

The deviation factor (1 - A^2)^2 is within ~1% of unity for A < 0.05 (everywhere we currently measure). All weak-field GR tests pass automatically. The framework departs from GR at second order in A and at the boundary itself.

---

## Black-hole interpretation

In Model-A, a black hole is interpreted as a 2D bubble surface, not a deep interior region.

```text
A < 1  -> spacetime exists
A = 1  -> 2D phase boundary / bubble surface
A > 1  -> not part of the manifold
```

Matter that fell toward the black hole never crossed A = 1; it accumulated holographically on the bubble surface (outward-collapse picture). There is no interior; there is no singularity. Information lives on the 2D boundary surface.

For an infalling traveler, proper time to reach A = 1 is logarithmically infinite. The observer and traveler agree: there is no finite-time crossing event.

**Spinning black holes**: the bubble warps to an oblate shape, with equatorial radius preserved at 2M and polar radius shrinking to the Kerr horizon location. Matter on the bubble carries the angular momentum (the bubble is the screen, the matter is the rotating hologram). No Penrose extraction (since the rotational energy is in the matter, not in empty rotating geometry).


## A_0 void interpretation

The 4 pi is the prefactor that appears in the thermal-emission rule k_B T = hbar c |grad A| / (4 pi), independently derived in Q8/Q10 as 2 pi (thermal-state imaginary-time periodicity) x 2 (Model-A gravity bridge factor c^2/2). The 3 is the count of spatial dimensions over which A must be non-zero for a 3D manifold to exist.

Together, A_0 is proposed as the structural floor of spacetime — the minimum "amplitude per solid-angle x dimensionality" at which a spatial manifold can structurally manifest. Below this value, space cannot exist; A = 0 is not vacuum but absence of manifold.

The numerical match between 1/(12 pi) = 0.026526 and the empirical bridge term b/L = 354.95/13387 = 0.026514 is 0.04%. This is the framework's strongest evidence for treating A_0 as a derived structural constant rather than a calibrated parameter, though the first-principles derivation of the spatial-dimensionality factor of 3 remains open theoretical work.

---

## Thermodynamics

Model-A reproduces the standard black-hole thermodynamic results from a single rule:

```text
k_B T = (1 / (4 pi)) hbar c |grad A| at the boundary
```

This rule reproduces:

- Schwarzschild Hawking temperature
- Unruh temperature
- de Sitter horizon temperature
- Bekenstein-Hawking entropy `S = k_B Area / (4 ell_P^2)`
- The first law `dE = T dS`
- The Smarr relation `M c^2 = 2 T S`
- Generalized second-law behavior (1/3 surplus during evaporation)
- Standard Hawking evaporation lifetime scaling

The 4 pi prefactor factors as 2 pi (thermal-state imaginary-time periodicity) x 2 (Model-A gravity bridge factor c^2/2). Both factors are independently derived. The mechanism is phase-boundary equilibrium with asymmetric resolution at the A = 1 surface, not Schwarzschild Wick rotation.

For spinning bubbles, T(theta) is non-uniform along the bubble: equator at T_Schwarzschild for any spin, pole at T_Kerr (matching standard Kerr horizon T exactly via identity). Total Hawking emission stays substantial even near extremal, with super-radiance enhancement at the equator.

---

## Cosmological branch

The cosmological structure is now closed with one structural commitment and one calibrated parameter:



```text
V(A) = alpha / A  +  beta / (1 - A)

alpha / beta = [ A_0 / (1 - A_0) ]^2  (fixed by A_0 commitment)
beta/rho_crit calibrated to match observed cosmological dynamics
```

The potential diverges at both A = 0 and A = 1, reflecting the structural commitment that both endpoints are asymptotic boundaries.

**Bridge term (derived):**

```text
b = A_0  x  c / H_0
  = (1 / (12 pi))  x  c / H_0
  ~ 355.10 Mly at H_0 = 73.04 km/s/Mpc
```

Matches historical Pantheon/Union3 fit value (354.95 Mly) to 0.04%.

**Cosmological story:**

- Universe is matter-dominated at H_0 = 73 throughout cosmic history.
- The CMB-inferred value H_0 = 67.4 (Planck) is biased downward by an unaccounted-for photon-A line-of-sight effect.
- That LoS effect is the cumulative-A signature: photons traversing cosmic distances pick up extra apparent path through the structured cosmic web (galaxies, filaments, clusters).
- Self-consistent closure of CMB theta_star at H_0 = 73 requires LoS amplification factor ~1.5x A_0, plausibly supplied by realistic cosmic structure.

**Galactic dark matter:**

Model-A is naturally compatible with primordial-black-hole dark matter (PBH-DM). Each PBH is a small bubble with the same A = 1 boundary structure as stellar and super-massive BHs. Cumulative A from baryons + PBH-DM halo (NFW-like distribution) reproduces observed galactic rotation curves trivially.

---

## Quantum interpretation

Model-A uses a physical, not conscious, definition of observation:

```text
Observation = physical interaction that resolves A
```

The model separates states into resolved (A-state confirmed by interaction; path is definite) and unresolved (no physical interaction yet; path is indeterminate).

- Decoherence is interpreted as the dense accumulation of resolution events between a system and its environment.
- Schrodinger's cat is dead-or-alive at the moment of sealing the box, because internal interactions resolve the cat continuously. The cat was never in superposition ontologically; we simply lack epistemic access until we open the box.
- The arrow of time emerges from the irreversibility of resolution events.

The same A field plays both classical (magnitude) and quantum (resolution-status) roles. Resolution-event statistics suggest a route to deriving Born-rule probabilities, though this is open theoretical work.

---

## Falsifiable predictions

Model-A makes specific predictions in regimes LCDM does not address, each with concrete observational targets:

```text
G1: Ringdown damping factor tau / tau_GR = 1.80
    at same dominant frequency as GR.
    Testable now with LIGO/Virgo/KAGRA ringdown data.

F6: Gravitational decoherence ~0.5 s for 1 micron silica
    nanoparticle in superposition. Cavity-optomechanics frontier;
    achievable in the next decade.

G3-G4: Spinning BH Hawking emission is latitudinally banded
    (equator hot, poles cold). Total emission rate stays high
    at extremal spin, unlike standard Kerr thermal channel.
    Testable if any primordial BH evaporation is observed.

PBH-DM: requires inflationary fluctuation amplitude sigma ~ 0.05
    at the PBH-formation scale. Constrains specific inflation
    models (USR, hilltop, single-field-with-bump).
```

---

## What is strongest so far

1. **Weak-field recovery.** Newtonian gravity follows directly from `g = (c^2/2) grad A`. GPS, Shapiro, and lensing match GR weak-field predictions identically.

2. **A_0 = 1/(12 pi) committed structurally.** Bridge term `b = A_0  x  c/H_0` matches historical Pantheon/Union3 fit to 0.04%. A_0 is no longer a calibrated parameter.

3. **The thirds-of-A strong-field structure.** ISCO, photon sphere, and horizon fall at A = 1/3, 2/3, 1 — preserved exactly because g_tt is unchanged.

4. **Black-hole thermodynamics from one rule.** `k_B T = hbar c |grad A| / (4 pi)` reproduces Hawking T (Schwarzschild + Kerr + de Sitter), Bekenstein-Hawking entropy, the first law, the Smarr relation, evaporation lifetime, and the generalized second law.

5. **Strong-field metric is Model-A-derived.** The (1 - A^2)^2 modification of g_rr is forced by three Model-A principles (weak-field recovery + boundary at A = 1 + thirds-of-A preservation), not borrowed from elsewhere.

6. **Spinning bubble framework.** Static bubble + rotating holographic matter gives non-uniform Hawking T (T_pole = T_Kerr exactly via identity), super-radiance enhancement, no Penrose extraction.

7. **SN distance fits beat LCDM combined chi^2.** Model-A with V_3 modified Friedmann wins by 25-33 across Pantheon+/Union3/DES at the same number of free parameters. Predicts inter-catalog Pantheon+/Union3 tension within 27%.

8. **CMB self-consistent closure at H_0 = 73.** Cumulative-A line-of-sight amplification structurally explains the H_0 tension. Internal solution exists for realistic cosmic-structure parameters.

9. **PBH-DM compatibility for galactic dark matter.** Each PBH is a small bubble with the framework's existing thermodynamic/structural machinery. Galactic rotation curves close trivially with PBH-halo + cumulative A.

10. **Six observational regimes, one A field.** Local gravity, propagation delay, BH thermodynamics, SN distances, CMB acoustic scale, and galactic DM all from the same A field with one structural constant (A_0 = 1/(12 pi)) and one calibrated parameter (beta).

11. **Specific falsifiable predictions** (G1 ringdown, F6 decoherence, G3-G4 latitudinal Hawking, PBH-DM sigma constraint) — distinguishable from LCDM in regimes the standard framework does not address.

---

## Major open problems

1. **Lagrangian for A.** Specifying the action `S[A, g_uv]` that produces both the Poisson-like source equation and the strong-field metric construction. Currently the metric is committed by ansatz.

2. **Structural origin of the "3" factor in 1/(12 pi).** The 4 pi is independently derived (thermal prefactor); the 3 from spatial dimensionality currently rests on dimensional argument. Closing it would make A_0 fully derived.

3. **Realistic cosmic-structure modeling for f_LoS.** The CMB self-consistent solution depends on the line-of-sight A amplification factor. Computing this from realistic cosmic structure (N-body or analytic modeling) would close the cosmology branch.

4. **BAO direct test under self-consistent STAM cosmology.** Simple constant-A version failed at 22 sigma. Structure-dependent A line-of-sight test under V_3 cosmology is open.

5. **First-principles Born rule from resolution-event statistics.** The framework's quantum interpretation is consistent but not yet predictive at the level of probability rules.

6. **G1 exact computation.** The ringdown damping prediction is currently eikonal-approximation. Exact Regge-Wheeler-type calculation on the Model-A metric would tighten the prediction.

7. **A_collective for galactic dynamics if PBH-DM is rejected.** If primordial-black-hole dark matter is observationally ruled out, the framework needs to either derive an A_collective galactic-scale enhancement from first principles or accept some other DM mechanism.

---

## Speculative parking-lot ideas

These are not part of the active research line but are catalogued because they fit the framework's structure and may be worth revisiting:

- **Two-hologram horizon.** Original matter encoded on inner surface, accumulating matter on outer surface, with cosmic expansion driven by outer-layer growth.
- **Universe-as-bubble.** Our cosmos may itself be a bubble, with information encoded on its outer horizon (connects to holographic universe proposals).
- **Multiverse from non-isolated formation events.** If our universe formed inside a larger structure, similar events may have produced other universes.
- **Black-hole-as-baby-universe (Smolin-CNS-adjacent).** Time-reparameterization across the horizon dissolves the timescale-mismatch objection in principle. Mathematical correspondence with Big Bang expansion is partial but real.

These are explicitly "we cannot know" territory — speculative, unfalsifiable from inside the current framework, but structurally compatible.

---

## Current status

Model-A is a structured theoretical research program with one structurally-committed constant (A_0 = 1/(12 pi)), one calibrated cosmological parameter (beta), and a coherent unified A-field across six observational regimes. The framework reproduces all currently-tested predictions of GR identically (weak field) and provides a clean ontological story for black holes (bubble picture, holographic information, no interior), thermodynamics (one rule for T), and cosmology (V_3 modified Friedmann, structural Hubble-tension story).

Where the framework departs from LCDM and is therefore distinguishable in principle: G1 ringdown, F6 decoherence, latitudinally-banded Hawking emission for spinning BHs, the structural CMB-tension closure mechanism. None of these have been observationally settled yet.

The framework is **not** presented as complete. Significant open problems remain (Lagrangian, realistic structure modeling, BAO test, first-principles "3" derivation). Its value is that it creates a unified language with one structural constant and accountable falsification targets, rather than a collection of independent ad-hoc components.

---

## Repository purpose

This repository is intended to preserve:

* Core formulas
* Derivations
* Numerical checks
* Catalog diagnostics
* Falsification tests (passed, failed, and pending)
* Speculative parking-lot ideas (clearly labeled)
* Open problems
* Scripts and results

Supportive demonstrations and falsification tests should remain separate. A real prediction test should lock parameters first, preserve the script and output, and then compare against independent data without silent retuning.

When the framework changes direction (as it did multiple times in development — the symmetric-boundary commitment, the V_3 form selection, the PBH-DM compatibility recognition), the change is recorded in the relevant memory file rather than retroactively edited into the historical record.

---

## Acknowledgments

Development assisted by extensive conversation with Claude (Anthropic), particularly during the May 2026 push that produced the A_0 = 1/(12 pi) structural commitment, V_3 potential selection, CMB self-consistent closure, and PBH-DM compatibility analysis (G7-G14 script series).

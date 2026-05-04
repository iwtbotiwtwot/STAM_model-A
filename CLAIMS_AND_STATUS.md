# Claims and Status — STAM Model-A

This document records the current claim set for STAM Model-A and separates:

```text
claimed structure
verified identities
numerical tests already run
diagnostic/catalog results
open problems
items not yet claimed
```

It is intended as a project ledger, not as proof that the model is established physics.

---

## 1. Status labels

| Label | Meaning |
|---|---|
| **Core definition** | Part of the Model-A framework as currently defined. |
| **Algebraic identity** | Follows directly from the definitions. |
| **Numerically checked** | Verified by a script in this repository or generated test output. |
| **Catalog diagnostic** | Comparison against current observational catalog data; interpretation still open. |
| **Open** | Needed for a complete theory or stronger empirical validation. |
| **Not claimed** | Explicitly not claimed at this stage. |

---

## 2. Core Model-A definitions

### 2.1 Accumulation field

**Claim type:** Core definition

STAM Model-A uses a dimensionless accumulation field:

```text
A
```

In the local spherical weak-field limit:

```text
A(r) = Rs/r = 2GM/(c²r)
```

where:

```text
Rs = 2GM/c²
```

**Status:** Core definition retained.

---

### 2.2 Gravity bridge

**Claim type:** Core definition / algebraic weak-field recovery

Model-A uses:

```text
g⃗ = (c²/2)∇A
```

For:

```text
A(r) = 2GM/(c²r)
```

this gives the Newtonian inverse-square acceleration in the spherical weak-field limit:

```text
|g| = GM/r²
```

**Status:** Algebraically verified and numerically checked.

**Not claimed:** This alone is not a full replacement for general relativity, nor a complete strong-field metric theory.

---

### 2.3 Propagation delay

**Claim type:** Core definition

Model-A uses path accumulation for propagation delay:

```text
Δt = (1/c)∫A(r)ds
```

For the spherical weak-field form:

```text
A(r) = 2GM/(c²r)
```

this becomes:

```text
Δt = (2GM/c³)∫ds/r
```

**Status:** Numerically checked against the weak-field logarithmic Shapiro-delay structure.

---

### 2.4 Horizon threshold

**Claim type:** Core definition / algebraic identity

Model-A defines the horizon threshold as:

```text
A = 1
```

Using:

```text
A(r) = Rs/r
```

gives:

```text
A = 1 ⇔ r = Rs = 2GM/c²
```

**Status:** Algebraically verified and numerically checked.

---

## 3. Non-distance tests currently run

### 3.1 Local gravity identity

**Claim tested:**

```text
A(r)=2GM/(c²r)
```

with:

```text
g⃗=(c²/2)∇A
```

recovers:

```text
g = GM/r²
```

in the spherical weak-field limit.

**Result:** Identity recovered to numerical precision.

**Status:** Numerically checked.

**Open:** Full precision solar-system tests, non-spherical sources, and strong-field formulation.

---

### 3.2 Shapiro delay

**Claim tested:**

```text
Δt = (1/c)∫A(r)ds
```

with:

```text
A(r)=2GM/(c²r)
```

reproduces the weak-field logarithmic Shapiro-delay structure.

For a straight path with impact parameter `b_imp`:

```text
Δt = (2GM/c³)[asinh(x2/b_imp)-asinh(x1/b_imp)]
```

The script also compares this with the standard first-order logarithmic expression.

**Representative result:**

For a solar-grazing Earth-Mars-like path:

```text
one-way delay ≈ 123.6076 μs
two-way delay ≈ 247.2151 μs
```

**Status:** Numerically checked.

**Open:** Full solar-system timing model, PPN-level comparison, strong-field propagation.

---

### 3.3 GPS satellite clock adjustment

**Claim tested:**

Model-A clock-rate mapping in weak field:

```text
dτ/dt ≈ 1 - A/2
```

with:

```text
A(r)=2GM/(c²r)
```

implies a gravitational clock-rate difference between Earth geoid/surface and GPS orbit.

Using a GPS/geoid reference plus circular-orbit kinematic correction gives:

```text
STAM gravitational gain ≈ +45.787467 μs/day
kinematic loss          ≈  -7.213600 μs/day
net satellite gain      ≈ +38.573867 μs/day
```

Required factory offset:

```text
Δf/f ≈ -4.464568 × 10⁻¹⁰
```

**Status:** Numerically checked as a weak-field clock-rate test.

**Open:** Full GPS operational model, Earth multipoles, eccentricity correction, Sagnac correction, broadcast ephemeris/clock terms, atmosphere, and receiver processing.

---

### 3.4 Horizon threshold and black-hole size from mass

**Claim tested:**

```text
A(r)=Rs/r
A=1 ⇔ r=Rs
```

and:

```text
M = c²r_h/(2G)
```

The test also checks the identity:

```text
A = (v_escape/c)²
```

because:

```text
v_escape² = 2GM/r
```

Therefore:

```text
A < 1 ⇔ v_escape < c
A = 1 ⇔ v_escape = c
A > 1 ⇔ v_escape > c
```

**Representative values:**

```text
1 solar mass horizon radius  ≈ 2.953 km
10 solar mass horizon radius ≈ 29.533 km
```

**Status:** Algebraic identity and numerical check.

**Open:** Rotating black holes, charged black holes, photon sphere, observed shadow size, accretion physics, and interior/over-threshold dynamics.

---

### 3.5 Gravitational-wave propagation

**Claim tested:**

Model-A position for this test:

```text
Gravitational waves propagate locally at c.
```

Apparent delay near high accumulation is represented as traversal through accumulated spacetime, not as a lower local wave speed.

Toy propagation rule:

```text
t_obs = ∫ds/c + k∫A(s)ds/c
```

Equal-coupling case:

```text
k_GW = k_EM
```

Then gravitational waves and electromagnetic waves receive the same accumulation traversal delay through the same `A` field.

**Result:**

For equal coupling:

```text
GW propagation delay - EM propagation delay = 0
```

for the tested solar and near-horizon toy paths.

The same script checks:

```text
A = (v_escape/c)²
```

for the horizon relation.

**Status:** First-pass consistency check.

**Open:** Full gravitational-wave propagation in dynamical spacetime, waveform effects, lensing/time-delay comparison, and strong-field emission/escape behavior.

---

### 3.6 Mass-estimator consistency

**Claim tested:**

Multiple weak-field/threshold estimators infer the same source mass:

```text
Horizon:              M = c²r_h / 2G
Acceleration:         M = gr² / G
Orbital velocity:     M = v²r / G
Shapiro coefficient:  M = Kc³ / 2G
Gravitational shift:  M ≈ z_grav c²r / G
Lensing deflection:   M ≈ αc²b / 4G
```

**Result:** Synthetic tests recover input mass ratios at floating-point precision.

**Status:** Numerically checked as identity/synthetic consistency.

**Open:** Real-data mass closure across independent observations of the same object.

---

## 4. Cosmological distance structure

### 4.1 Geometric / SU layer

**Claim type:** Core cosmological structure

Earlier SU notation:

```text
SU(z) = 1231.350175 × (3.33333z + 0.50000z²)
```

is equivalent to:

```text
D_geo(z) = Lz(1 + 0.15z)
```

up to unit normalization.

**Status:** Retained as the geometric spine.

---

### 4.2 No-b Model-A distance law

**Claim type:** Current preferred forward-testing form

The no-b accumulation-adjusted Model-A distance law is:

```text
D_adj,0(z) = Lz(1 + 0.5z)
```

and:

```text
D_excess,0(z) = 0.35Lz²
```

**Status:** Current preferred physical distance form for forward tests.

---

### 4.3 Historical b bridge

**Claim type:** Catalog diagnostic, not core physics

Historical catalog comparison form:

```text
D_catalog(z) ≈ D_adj,0(z) + bz
```

Current interpretation:

```text
b
```

is a catalog-bridge diagnostic measuring the difference between current supernova catalog-inferred distances and the no-b Model-A distance relation.

Retained historical values:

```text
Pantheon/Union-style bridge: b ≈ 354.95
Original retained bridge:    b ≈ 461.3626922
DES-style bridge:            b ≈ 1335.412792
```

**Status:** Retained as diagnostic/historical record.

**Not claimed:** `b` is not currently treated as a fundamental Model-A physical constant.

---

## 5. Supernova catalog tests

### 5.1 No-b supernova run

**Claim tested:**

Run Union3, Pantheon, and DES against:

```text
D_adj,0(z)=Lz(1+0.5z)
```

without using the `b` bridge.

Common range:

```text
0.05 ≤ z ≤ 1.14418
```

Residual definition:

```text
observed catalog distance modulus - no-b Model-A prediction
```

**Results:**

| Catalog | Median residual | Mean residual | RMSE |
|---|---:|---:|---:|
| Union3 | +0.0512 mag | +0.0429 mag | 0.0539 mag |
| Pantheon | +0.0808 mag | +0.0800 mag | 0.1664 mag |
| DES | +0.2068 mag | +0.2383 mag | 0.3626 mag |

**Status:** Catalog diagnostic.

**Current interpretation:** Removing `b` does not place the catalogs on a random or incoherent curve. Union3/Pantheon sit closer to the no-b relation than DES. DES remains the larger catalog-family offset.

---

### 5.2 DES / Pantheon / Union3 discrepancy

**Claim tested:**

Check whether DES discrepancy relative to STAM is also visible as a catalog-to-catalog discrepancy relative to Pantheon/Union3.

**Previous result summary:**

```text
Pantheon/Union3 behave like one reference family.
DES behaves differently.
STAM residual against DES is similar in size and structure to the DES-vs-Pantheon/Union3 empirical discrepancy.
```

**Status:** Catalog diagnostic.

**Open:** Independent DES systematics analysis, pre-declared outlier treatment, and reproduction with final catalog files committed in `data/` and results committed or documented in `results/`.

---

### 5.3 b sensitivity

**Claim tested:**

Whether a small change in `b` can bring Union3, Pantheon, and DES into agreement.

**Result summary:**

```text
Small b adjustments do not bring all three catalogs into agreement.
A DES-style high b improves DES only by damaging the Pantheon/Union3 relation.
```

**Status:** Catalog diagnostic.

**Current interpretation:** Supports treating `b` as a catalog bridge rather than a core physical parameter.

---

## 6. BAO geometric split test

**Claim tested:**

Whether BAO transverse distances prefer the geometric layer or the adjusted layer.

Mappings tested include:

```text
D_M ≈ D_geo/(1+z)
D_M ≈ D_adj/(1+z)
```

with one fitted nuisance scale:

```text
r_d_eff
```

**Result summary:**

```text
Full all-z BAO:
    D_geo/(1+z) performed better.

Low/mid-z BAO:
    D_adj/(1+z) performed better.
```

A separate check found that the flip is not explained by `b = 354.95`; it persists when:

```text
b = 0
```

**Status:** Early independent-observable diagnostic.

**Open:** Full BAO-facing distance rule, covariance matrices, `D_H`, `D_V`, and an internally defined Model-A expansion/ruler treatment.

---

## 7. Exploratory galaxy accumulation theory

**Claim type:** Exploratory / scratch theory

This section records the current STAM-native idea for galaxy-scale behavior. It is not yet a formal Model-A claim and is not currently presented as a completed solution to galaxy rotation curves.

### 7.1 Motivation

Original STAM intuition treated accumulation as invisible and difficult to measure directly. In cosmology, redshift became the working coordinate for expressing path accumulation. For galaxies, the analogous question is whether many small, individually unobservable accumulation contributions can combine into an observable galaxy-scale effect.

The working idea is:

```text
Mass-energy produces nonzero A even at distances where each individual contribution is too small to observe directly.

In a galaxy, many stars, gas clouds, compact objects, and central mass concentrations may contribute tiny A values that sum into a coherent galaxy-scale accumulation field.
```

### 7.2 Linear cumulative A

**Claim type:** Model-A-native exploratory calculation

The linear cumulative version is:

```text
A_total(x) = Σ 2Gm_i / (c² |x - x_i|)
```

This is the direct extension of the local spherical Model-A expression:

```text
A(r) = 2GM/(c²r)
```

to many sources.

Current scratch-test result:

```text
A from one solar mass at 10 kpc      ≈ 9.571121e-18
A from 1e11 solar masses at 10 kpc   ≈ 9.571121e-7
```

The compact-source circular speed associated with:

```text
A ≈ 9.571121e-7
```

is approximately:

```text
v ≈ 207.39 km/s
```

using:

```text
v² ≈ (c²/2)A
```

**Status:** Numerically checked as a scale demonstration.

**Interpretation:** Individually tiny A contributions can accumulate into a galaxy-scale quantity. Decimal precision preserves tiny nonzero values, but the physical effect comes from cumulative summation, not from precision alone.

### 7.3 Toy visible-galaxy calculation

A simple toy visible galaxy was tested with:

```text
disk  = 6e10 solar masses
gas   = 1e10 solar masses
bulge = 1e10 solar masses
```

The linear visible cumulative-A toy produced:

```text
outer 15–35 kpc median speed ≈ 121.08 km/s
outer slope                  ≈ -2.82 km/s/kpc
```

**Status:** Scratch numerical test.

**Interpretation:** Linear cumulative A from visible toy components produces galaxy-scale orbital speeds, but the tested toy model does not by itself guarantee flat edge behavior.

### 7.4 Collective A-envelope possibility

**Claim type:** Exploratory extension, not core Model-A

The stronger idea is that a galaxy may behave as a coherent accumulation domain, not merely as a linear sum of independent point-source contributions.

This would require an additional term:

```text
A_total = A_linear + A_collective
```

where:

```text
A_collective
```

could represent a galaxy-scale accumulation envelope, coherence effect, long-lived structure effect, central-density seeding, or path-memory-like accumulation.

One exploratory toy used a small logarithmic-style envelope. Fitting a target edge speed of:

```text
220 km/s
```

near:

```text
20 kpc
```

gave:

```text
epsilon ≈ 5.775256e-7
```

At approximately 20 kpc:

```text
A_linear      ≈ 4.04e-7
A_collective  ≈ 9.10e-7
```

With this exploratory envelope:

```text
outer 15–35 kpc median speed ≈ 198.69 km/s
outer slope                  ≈ -1.71 km/s/kpc
```

**Status:** Scratch exploratory calculation only.

**Not claimed:** Model-A does not currently claim that `A_collective` is proven, derived, or required. It is a candidate direction if linear visible accumulation is insufficient.

### 7.5 Working galaxy hypothesis

Current scratch hypothesis:

```text
A galaxy can behave as a single accumulation object because many individually tiny A fields combine into a coherent galaxy-scale A field.

If observed edge-star behavior requires more than the linear sum, STAM may need a collective A-envelope term.
```

A stronger speculative form is:

```text
A_total = A_linear + A_collective
```

where the collective term may be related to central concentration, disk coherence, long-lived orbital structure, or accumulated galaxy-scale organization.

### 7.6 Connection to dark-matter interpretation

**Claim type:** Interpretive target / not established

STAM may eventually investigate whether dark-matter-like galaxy behavior is an interpretation of galaxy-scale accumulation effects.

Current careful statement:

```text
Dark-matter-like behavior may be a target for STAM reinterpretation through cumulative or collective A.

This is not yet a completed Model-A result.
```

### 7.7 Required future tests

To move this from scratch theory to a formal Model-A claim, future tests should:

```text
1. Use real galaxy rotation-curve datasets.
2. Compute A_required from observed v(r).
3. Compute A_linear from visible baryonic components.
4. Examine A_required - A_linear as an A_collective candidate.
5. Test whether A_collective has a stable shape across galaxies.
6. Test whether A_collective correlates with central concentration, disk structure, or total visible mass.
7. Pre-declare pass/fail criteria before claiming a galaxy-scale result.
```

**Status:** Open.


## 8. Not currently claimed

STAM Model-A does **not** currently claim:

```text
1. full replacement of general relativity;
2. complete strong-field metric formulation;
3. complete black-hole interior dynamics;
4. derivation of the no-b cosmological distance law from first principles;
5. complete BAO/CMB theory;
6. completed galaxy-rotation solution;
7. full GPS operational navigation model;
8. full gravitational-wave waveform theory;
9. that DES is wrong solely because STAM says so;
10. that b is fundamental physics.
```

---

## 9. Current open problems

Priority open problems:

```text
1. Derive A_path(z) or the no-b distance law from independent physical principles.
2. Define the BAO-facing distance observables in Model-A.
3. Extend local Model-A beyond spherical weak-field identities.
4. Build full clock/redshift formalism from A.
5. Test real mass closure across acceleration, lensing, redshift, and time delay for the same systems.
6. Define strong-field/horizon dynamics beyond the A=1 threshold identity.
7. Clarify gravitational-wave propagation in full dynamical settings.
8. Determine whether catalog bridge terms reveal measurement/inference structure or model incompleteness.
```

---

## 10. Current working position

Current project stance:

```text
A is the core Model-A variable.

The no-b distance law is the preferred forward-testing cosmological form.

b is retained as a catalog-bridge diagnostic and historical record.

Non-distance tests currently support the internal consistency of Model-A in:
    local gravity,
    Shapiro delay,
    GPS weak-field clocks,
    horizon threshold,
    gravitational-wave local-speed consistency,
    and mass-estimator identities.

Cosmological interpretation remains the largest open area.
```

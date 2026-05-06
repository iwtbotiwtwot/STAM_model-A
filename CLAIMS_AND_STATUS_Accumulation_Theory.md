# Claims and Status — STAM Model-A

## 0. Notable findings

This section records the strongest current findings before the longer claims ledger.


This document records the current working structure of STAM Model-A: what the model defines, what has been checked numerically, what remains diagnostic, and what still needs stronger testing.

It is a project ledger, not a declaration that the model is established physics.

---

## 1. Status labels

| Label | Meaning |
|---|---|
| **Core definition** | Part of the Model-A framework as currently defined. |
| **Algebraic identity** | Follows directly from the definitions. |
| **Numerically checked** | Verified by a script in this repository or generated test output. |
| **Catalog diagnostic** | Comparison against current observational catalog data; interpretation still open. |
| **Open** | Needed for a complete theory or stronger empirical validation. |
| **Boundary / open extension** | Active topic, but not yet completed or formalized. |

---

## Distance-interpretation principle

STAM Model-A does not assume that current catalog-inferred or model-inferred distances are identical to the underlying accumulation/geometric distance.

Current working distinction:

```text
observable / catalog-inferred distance
≠ automatically the same as
STAM geometric or accumulation distance
```

Therefore, distance discrepancies are not treated as automatic inaccuracies. They are retained as primary diagnostic targets.

The distance question is:

```text
Which observable sees which STAM layer?
```

Current STAM distance layers:

```text
D_geo(z)      = Lz(1 + 0.15z)
D_adj,0(z)    = Lz(1 + 0.5z)
D_excess,0(z) = 0.35Lz²
```

Catalog and angular-distance comparisons are kept in the record because they show how current observational distance frameworks sit relative to STAM's accumulation-distance structure.

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

Additional horizon interpretation:

```text
When an object's mass relative to its physical size produces A > 1 beyond the object's surface,
the A = 1 boundary becomes an external horizon.

Since A = (v_escape/c)², A > 1 means v_escape > c.

In STAM terms, the accumulation state beyond the surface is over-threshold:
outward escape would require faster-than-light propagation.
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

The geometric/SU layer is defined from the Hubble-scale normalization, not from an arbitrary fitted coefficient.

In Mpc:

```text
D_geo,Mpc(z) = H⁻¹ z(1 + 0.15z)
```

with:

```text
H = 0.000243635
H⁻¹ = 4104.500584 Mpc
```

So:

```text
D_geo,Mpc(z) = 4104.500584 z(1 + 0.15z)
```

In million light-years:

```text
D_geo,Mly(z) = (C/H) z(1 + 0.15z)
```

where:

```text
C = 3.261563776 Mly/Mpc
L = C/H = 13387.090426 Mly
```

so:

```text
D_geo,Mly(z) = Lz(1 + 0.15z)
```

The earlier SU notation:

```text
SU(z) = 1231.350175 × (3.33333z + 0.50000z²)
```

is the same Mpc-scale relation written in anchored form:

```text
1231.350175 = H⁻¹ / 3.33333
0.50000 / 3.33333 = 0.15
```

Therefore:

```text
1231.350175 × (3.33333z + 0.50000z²)
= H⁻¹ z(1 + 0.15z)
```

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

### 4.3 Historical supernova bridge term

**Claim type:** Catalog comparison history / diagnostic record

During earlier supernova-distance comparisons, the form:

```text
D_catalog(z) ≈ D_adj,0(z) + bz
```

was used as a bridge between current catalog-inferred distances and the no-b Model-A relation.

The useful role of `b` is historical and diagnostic:

```text
b helped show how current supernova catalogs sit relative to the no-b Model-A curve.
b helped identify that Pantheon/Union3 and DES behave like different catalog families.
```

Retained historical values:

```text
Pantheon/Union-style bridge: b ≈ 354.95
Original retained bridge:    b ≈ 461.3626922
DES-style bridge:            b ≈ 1335.412792
```

**Status:** Retained as catalog-comparison history.

Forward Model-A testing prioritizes the no-b distance relation:

```text
D_adj,0(z) = Lz(1 + 0.5z)
```

---

### Historical b derivation candidate: TE anchor linearization

**Claim type:** Diagnostic / derivation candidate

This section keeps `b` in its proper place: historical catalog bridge analysis. It is not placed in the headline findings and it is not part of the forward physical no-b distance law.

The historical bridge term:

```text
D_bridge(z) = bz
```

can be compared to Traversal Excess:

```text
TE(z) = 0.35Lz²
```

Because TE is quadratic, a low-z catalog anchor can cast a linearized bridge:

```text
b_pred = dTE/dz | z=z_anchor
       = 0.70L z_anchor
```

This gives a derivation candidate:

```text
b = 0.70L z_anchor
```

where `z_anchor` is derived from catalog/calibration structure rather than chosen to fit `b`.

Reverse-derived anchors:

```text
Pantheon/Union-style b ≈ 354.95      → z_anchor ≈ 0.037878
Original retained b ≈ 461.3626922    → z_anchor ≈ 0.049233
DES-style b ≈ 1335.412792            → z_anchor ≈ 0.142505
```

Catalog-anchor test highlights:

```text
Pantheon+Union3 weighted q25 z ≈ 0.037250
b_pred ≈ 349.07 Mly
error ≈ -1.66%

Pantheon weighted q25 z ≈ 0.037234
b_pred ≈ 348.92 Mly
error ≈ -1.70%
```

**Status:** Strong derivation candidate, not final proof.

**Current interpretation:** `b` may be the linear shadow of quadratic TE around an effective low-z catalog anchor. This keeps `b` as catalog bridge history while giving it a possible non-arbitrary origin.


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

**Current interpretation:** Removing `b` does not place the catalogs on a random or incoherent curve. Union3/Pantheon sit closer to the no-b relation than DES. DES remains the larger catalog-family offset. This is retained as a distance-interpretation result: STAM predicts current catalog-inferred distances may differ from the no-b accumulation-distance relation.

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

### 5.3 Historical bridge sensitivity

**Question tested:**

Whether small changes to the historical bridge term bring Union3, Pantheon, and DES into agreement.

**Result summary:**

```text
Small bridge adjustments do not bring all three catalogs into agreement.
A DES-style high bridge value improves DES only by damaging the Pantheon/Union3 relation.
```

**Status:** Catalog diagnostic.

**Current interpretation:** The bridge term is useful for studying catalog behavior, but it should not drive the forward Model-A distance law.

---

## 6. BAO geometric split test

**Claim type:** Independent distance-layer diagnostic

The BAO test asks whether compressed transverse BAO distances align more naturally with one STAM distance layer than another.

Mappings tested included:

```text
D_M ≈ D_geo/(1+z)
D_M ≈ D_adj,0/(1+z)
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
    D_adj,0/(1+z) performed better.
```

A separate check found that the flip is not explained by the historical supernova bridge value:

```text
b = 354.95
```

The flip persists when:

```text
b = 0
```

**Status:** Active distance-layer diagnostic.

**Interpretation:** The BAO result is retained as evidence that observable distance mappings are not yet fully settled. Because STAM separates accumulation/geometric distance from current inferred distance frameworks, this is not recorded as a simple accuracy failure. It identifies the need for a clearer BAO-facing observable-distance rule.


## 6.1 CMB angular-scale toy test

**Claim type:** High-redshift observable-distance diagnostic

The high-z stress test showed that the no-b STAM path-accumulation functions remain finite through CMB-scale redshift and approach a high-z limit:

```text
⟨A_path⟩ → 7/3
A_path,local → 7/3
D_adj,0/D_geo → 10/3
```

At approximately:

```text
z ≈ 1100
```

the no-b Model-A functions remain mathematically stable.

A separate CMB angular-scale toy test asked whether simple direct mappings of current STAM layers reproduce the observed acoustic angular scale:

```text
θ* ≈ 0.010411 rad
```

Simple tested mappings included:

```text
D_geo
D_adj,0
D_excess,0
D_geo/(1+z)
D_adj,0/(1+z)
D_geo/(1+z)²
D_adj,0/(1+z)²
```

**Result summary:**

No simple tested mapping directly reproduced the observed CMB acoustic angular scale at `z≈1090`.

**Status:** Active high-redshift distance-mapping diagnostic.

**Interpretation:** This result is retained as an open observable-mapping issue, not as a direct rejection of Model-A. STAM predicts that current standard distance interpretations may differ from accumulation/geometric distance. The CMB result therefore identifies the need for a CMB-facing rule:

```text
Which STAM distance layer maps to CMB angular observability?
What is the STAM interpretation of the sound-horizon ruler?
Is z≈1100 an accumulation marker, redshift marker, or both?
```

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


## 8. Active boundaries and open extensions

These items remain active areas of development rather than completed Model-A results:

```text
1. Full strong-field metric formulation.
2. Black-hole interior / A>1 dynamics beyond the threshold identity.
3. First-principles derivation of the no-b cosmological distance law.
4. Full BAO/CMB observable-distance mapping framework.
5. Completed galaxy-rotation reconstruction from real datasets.
6. Full GPS operational navigation model.
7. Full gravitational-wave waveform theory.
8. Independent explanation of DES/Pantheon/Union3 catalog-family behavior.
```

The historical supernova bridge term is kept as catalog-comparison history, not as a claim-boundary item.

---

---

## 11. A-conditioned observation thought experiments

**Claim type:** Thought experiment / conceptual development

This section records current STAM thought experiments related to observation, time, light-speed measurement, and cross-A comparison. These are not being presented as completed derivations. They are active conceptual tools for developing the Model-A language.

### 11.1 Local observation principle

Current thought-experiment principle:

```text
The observer is inside A too.
```

A local observer does not measure light, clocks, distance, or motion from outside the accumulation field. The observer's light path, ruler, clock, detector, neurons, cells, chemistry, and measurement process are all part of the same local A-conditioned system.

This gives the proposed STAM observation principle:

```text
Local c is preserved because both the measured light and the measuring system are governed by the same local A-condition.

Cross-A comparison reveals differences because entire physical systems translate differently across A.
```

So STAM separates:

```text
local experience
cross-A comparison
physical process rate
```

### 11.2 A-conditioning scale

A first-pass thought-experiment scaling is:

```text
S(A) = 1 + A
```

where:

```text
S(A)
```

represents local accumulation load / traversal scaling.

The toy light-clock relation is:

```text
D_A = S(A)D_geo
T_cross = D_A/c0
T_local = T_cross/S(A)
D_local = D_A/S(A)
```

Therefore:

```text
D_local/T_local = c0
```

This expresses the current STAM picture:

```text
A changes traversal, clock rate, ruler scale, and observer process rate together.

Local observation self-normalizes.

Cross-A comparison reveals the difference.
```

### 11.3 Time dilation reinterpretation

Current thought-experiment statement:

```text
Time dilation exists as a cross-A comparison.

It does not exist as a locally experienced slowing of time.
```

A high-A traveler is not sitting there experiencing slow chewing, slow thinking, or slow waving. Their chewing, neurons, light-clock, heartbeat, chemistry, and local measuring system are all in the same A-condition.

An outside observer may translate those same actions as slowed because the relation between observer and traveler passes through different accumulation conditions.

Current language:

```text
In STAM, time dilation is a cross-A translation between physical process rates.

Local experience is A-normalized.
```

### 11.4 Light speed and A-conditioned perception

The thought experiment with a light clock and a magic bell led to this STAM picture:

```text
Everyone locally measures light at c because their measuring process is A-conditioned with the light path.

A distant observer comparing across A may calculate a different traversal rate if using unconverted geometric distance.

Using accumulated traversal restores the local-c relation.
```

So STAM does not need to say:

```text
light locally slows
```

It says:

```text
local c is preserved;
cross-A traversal must be translated through A.
```

### 11.5 Traveler near a horizon

The falling-traveler picture is:

```text
The traveler chews gum and waves normally in their local A-condition.

The outside observer sees those same actions increasingly stretched because the traveler is moving through rising accumulation.

The traveler is not experiencing frozen time.

The observer is translating the traveler's motion through increasing A.
```

Near the horizon:

```text
A → 1
```

the ordinary scaling:

```text
S(A) = 1 + A
```

may not be enough. A separate near-horizon traversal factor is being considered:

```text
S_h(A) = 1/(1 - A)
```

This gives the horizon-style divergence:

```text
A → 1  ⇒  S_h(A) → ∞
```

At:

```text
A = 0.99999999
```

```text
S_h(A) = 100,000,000
```

So small local motion can translate into enormous outside comparison delay while local experience remains normal.

Current STAM-good picture:

```text
A = 0:
    clear road

A > 0:
    thicker road

A → 1:
    road becomes nearly impossible to compare across

A = 1:
    horizon threshold

A > 1:
    outward escape over-threshold
```

### 11.6 Negative A / A-deficit paths

If the same scaling is extended:

```text
S(A) = 1 + A
```

then negative relative A means reduced traversal load:

```text
-1 < A < 0
```

Current thought-experiment interpretation:

```text
A < 0 represents accumulation deficit relative to a reference baseline.

It does not imply local faster-than-light motion.

It means the path has reduced accumulated traversal.
```

Warp-style language:

```text
A warp-like effect would require an A<0 corridor where accumulated traversal is reduced below baseline,
allowing apparent >c crossing of geometric distance while local motion remains ≤c.
```

For example:

```text
A = -0.5
S(A) = 0.5
```

so the path has half the traversal load and an outside observer would infer an effective crossing speed of:

```text
c_eff = c/S(A) = 2c
```

This is not a ship locally moving faster than light. It is a reduced-A traversal path.

Boundary:

```text
A = -1  ⇒  S(A)=0
A < -1  ⇒  S(A)<0
```

Current thought-experiment stance:

```text
-1 < A < 0:
    reduced traversal / A-deficit regime

A = -1:
    zero-traversal boundary / model limit

A < -1:
    outside the current traversal interpretation unless a new rule is developed
```

### 11.7 Status of these ideas

These concepts are currently kept as thought experiments because they appear to connect several STAM pieces:

```text
local c
clock-rate changes
Shapiro delay
GPS clocks
horizon freezing
cross-A observation
negative-A traversal
```

They should be allowed to develop without forcing premature final wording. The next useful step is to keep calculating toy models and then decide which pieces deserve formal Model-A status.

## Latest diagnostic scripts

Recent b-diagnostic scripts:

```text
25_TE_clock_shapiro_b_bridge_derivation.py
26_b_from_GPS_Shapiro_scale_sniff.py
27_b_from_TE_anchor_redshift.py
```

Summary:

```text
25:
    Separates TE from historical bz bridge.
    Shows bridge/TE = b/(0.35Lz).

26:
    Tests whether GPS or Shapiro magnitudes derive b.
    Result: no. GPS/Shapiro validate categories but are many orders too small to set b.

27:
    Tests b as TE tangent at catalog anchor:
    b_pred = 0.70L z_anchor.
    Result: Pantheon/Union low-z weighted anchors closely approximate b≈354.95.
```

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

Historical bridge results are retained as catalog-comparison history.

Non-distance tests currently support the internal consistency of Model-A in:
    local gravity,
    Shapiro delay,
    GPS weak-field clocks,
    horizon threshold,
    gravitational-wave local-speed consistency,
    and mass-estimator identities.

Distance interpretation remains open by design: STAM predicts current inferred distances may not equal accumulation/geometric distances, so those discrepancies are kept as central diagnostic targets.

A-conditioned observation thought experiments are retained as conceptual development for local c, time dilation, cross-A comparison, and near-horizon traversal.
```

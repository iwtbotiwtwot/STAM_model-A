# Claims and Status

This document separates STAM Model-A results by evidentiary status.

The current working position is:

```text
A is the core Model-A variable.
b is not core physics.
b is retained as a catalog-bridge diagnostic for current supernova distance catalogs.
```

Model-A is still a proposed framework, not established physics. The results below should be read as a research status ledger: what is algebraic, what has been numerically checked, what is empirically suggestive, and what remains open.

---

## 1. Core Model-A definitions

### Local spherical accumulation

```text
A(r) = Rs / r = 2GM / (c²r)
```

where:

```text
Rs = 2GM / c²
```

### Local gravity bridge

```text
g = (c²/2) ∇A
```

For the spherical weak-field form `A(r)=Rs/r`, this reproduces the Newtonian inverse-square acceleration magnitude:

```text
g = GM / r²
```

### Propagation / traversal delay

```text
Δt = (1/c) ∫ A(r) ds
```

For `A(r)=2GM/(c²r)`, this becomes:

```text
Δt = (2GM/c³) ∫ ds/r
```

### Horizon threshold

```text
A = 1
```

For `A(r)=Rs/r`, this gives:

```text
r = Rs = 2GM/c²
```

### Cosmological distance split

The no-`b` physical Model-A distance structure is now preferred for forward testing:

```text
D_geo(z)      = Lz(1 + 0.15z)
D_adj,0(z)    = Lz(1 + 0.5z)
D_excess,0(z) = 0.35Lz²
```

The historical catalog-bridge form is retained diagnostically:

```text
D_catalog(z) ≈ D_adj,0(z) + bz
```

Current interpretation:

```text
b measures how current catalog-inferred distances differ from the no-b Model-A accumulation-distance relation.
```

---

## 2. Algebraic / definitional identities

These are important consistency checks, but by themselves they are not independent empirical validation.

| Item | Status |
|---|---|
| `A(r)=Rs/r` with `Rs=2GM/c²` | Definition in the local spherical weak-field limit |
| `g=(c²/2)∇A` gives `GM/r²` for `A=Rs/r` | Algebraic recovery of Newtonian weak-field spherical acceleration |
| `Δt=(1/c)∫A ds` gives `(2GM/c³)∫ds/r` | Algebraic substitution |
| `A=1` gives `r=Rs` | Algebraic horizon-threshold identity |
| `A=(v_escape/c)²` | Direct identity from `A=2GM/(c²r)` and `v_escape²=2GM/r` |
| `M=c²r_h/(2G)` | Direct inversion of the horizon threshold |
| `D_adj,0=D_geo+D_excess,0` | Definitional no-`b` distance decomposition |

---

## 3. Implemented non-distance Model-A tests

These are the strongest current Model-A results because they do not use supernova fitting, `b`, DES/Pantheon/Union3 adjustments, or the cosmological distance bridge.

### 3.1 Local gravity identity

**Claim tested:**

```text
A(r)=Rs/r and g=(c²/2)∇A reproduce Newtonian spherical weak-field acceleration.
```

**Status:** Passed as an algebraic/numerical identity check.

**Meaning:** Model-A’s local `A` field is normalized correctly for the weak-field spherical acceleration limit.

**Limitation:** This is not yet a full Solar-System PPN test, nor a proof of strong-field equivalence.

---

### 3.2 Shapiro delay / accumulation traversal delay

**Claim tested:**

```text
Δt = (1/c)∫A(r)ds
```

with:

```text
A(r)=2GM/(c²r)
```

reproduces the weak-field logarithmic Shapiro delay structure.

**Implemented result:**

For a solar-grazing Earth-Mars-like path:

```text
one-way delay ≈ 123.6076 μs
two-way delay ≈ 247.2151 μs
A at solar limb ≈ 4.2412 × 10⁻⁶
```

The analytic STAM accumulation integral matched the first-order logarithmic expression to:

```text
max relative error ≈ 1.7 × 10⁻¹²
```

The numerical path integral matched the analytic expression to:

```text
max relative error ≈ 1.6 × 10⁻¹¹
```

**Status:** Strong weak-field propagation structure pass.

**Meaning:** Model-A’s path-integral accumulation delay has the correct first-order Shapiro logarithmic form.

**Limitation:** This is not yet a full PPN/strong-field time-delay proof.

---

### 3.3 GPS satellite clock adjustment

**Claim tested:** Earth-surface/geoid accumulation differs from GPS-orbit accumulation, producing the gravitational clock-rate component.

Weak-field clock mapping tested:

```text
dτ/dt ≈ 1 - A/2
```

with:

```text
A(r)=2GM/(c²r)
```

Then:

```text
Δf/f_grav = (A_geoid - A_orbit)/2
```

**Implemented result at nominal GPS semimajor axis `r≈26,560 km`:**

```text
STAM gravitational gain ≈ +45.787467 μs/day
orbital kinematic loss  ≈  -7.213600 μs/day
net satellite gain      ≈ +38.573867 μs/day
```

Factory clock offset:

```text
Δf/f ≈ -4.464568 × 10⁻¹⁰
```

Reference GPS factory offset:

```text
Δf/f ≈ -4.4647 × 10⁻¹⁰
```

Difference:

```text
≈ 1.32 × 10⁻¹⁴ in fractional offset
≈ 1.14 ns/day in net timing
```

**Status:** Strong weak-field clock-rate alignment test.

**Meaning:** With the weak-field clock mapping `dτ/dt≈1-A/2`, Model-A reproduces the GPS gravitational clock correction scale and the net GPS factory offset when the standard orbital kinematic term is included.

**Limitation:** This is not a full GPS operations model. It does not include eccentricity correction, Sagnac correction, Earth multipoles, atmosphere, ephemeris processing, or satellite-specific broadcast-clock details.

---

### 3.4 Horizon threshold / black-hole size from mass

**Claim tested:**

```text
A=1
```

is both:

```text
r=Rs=2GM/c²
```

and:

```text
v_escape=c
```

because:

```text
A = (v_escape/c)²
```

**Implemented result:**

```text
max |A(r_h)-1|              = 0
max |v_escape(r_h)/c - 1|   = 0
max |M_recovered/M - 1|     = 0
max |A-(v_escape/c)²|       ≈ 4.44 × 10⁻¹⁶
```

Useful values:

```text
1 solar mass horizon radius  ≈ 2.953339 km
10 solar mass horizon radius ≈ 29.533394 km
```

**Status:** Exact spherical non-rotating threshold identity.

**Meaning:** Model-A’s horizon criterion is not arbitrary in the spherical case. `A=1` is exactly the escape-speed-equals-light-speed threshold.

**Limitation:** This does not yet model rotating Kerr horizons, charged black holes, photon spheres, observed shadow radius, accretion physics, or interior dynamics.

---

### 3.5 Gravitational-wave propagation

**Claim tested:**

```text
Gravitational waves propagate locally at c.
```

Apparent slowing near high `A` is modeled as accumulated traversal delay, not as a lower intrinsic wave speed.

Traversal rule tested:

```text
t_obs = ∫ds/c + k∫A(s)ds/c
```

Equal-coupling Model-A case:

```text
k_GW = k_EM
```

**Implemented result:**

```text
max GW-EM propagation difference, solar-path equal-coupling cases = 0 s
max GW-EM propagation difference, near-horizon equal-coupling toy cases = 0 s
```

Horizon/escape identity retained:

```text
A = (v_escape/c)²
```

so:

```text
A > 1 ⇔ outward escape would require v_escape > c
```

**GW170817-style sanity check:**

Using a simple `40 Mpc` distance scale and a `1.74 s` gamma-after-GW lag:

```text
naive fractional speed difference if all lag were propagation ≈ 4.23 × 10⁻¹⁶
```

Model-A equal-coupling propagation lag:

```text
0 s
```

**Status:** First-pass local-speed/equal-coupling consistency pass.

**Meaning:** STAM is compatible with gravitational waves moving locally at `c` if GW and EM signals share the same accumulation traversal coupling through the same `A` field.

**Limitation:** This is not a full gravitational-wave equation or dynamical spacetime calculation. Near-horizon examples are illustrative toy paths, not full strong-field waveform modeling.

---

### 3.6 Mass-inference consistency

**Claim tested:** Multiple observational routes infer the same source strength in the appropriate weak-field approximations.

Examples:

```text
Horizon radius:          M = c²r_h/(2G)
Acceleration:            M = gr²/G
Orbital velocity:        M = v²r/G
Shapiro coefficient:     M = Kc³/(2G)
Gravitational redshift:  M ≈ z_grav c²r/G
Lensing deflection:      M ≈ αc²b/(4G)
```

**Status:** Passed synthetic identity/regression tests.

**Meaning:** The formulas are internally consistent as projections of the same source strength.

**Limitation:** This is not yet an empirical lensing/dynamics closure test.

---

## 4. Cosmological / distance-status ledger

Distance remains the main nonstandard and unresolved part of Model-A.

The preferred forward physical form is no-`b`:

```text
D_adj,0(z)=Lz(1+0.5z)
```

The historical `b` term is now treated as a catalog bridge:

```text
D_catalog(z) - D_adj,0(z) ≈ bz
```

### 4.1 No-`b` supernova distance run

**Claim tested:** Removing `b` does not scatter supernovae into random or incoherent distances.

Catalogs used:

```text
Union3:   z, mb
Pantheon: zHD, MU_SH0ES
DES:      zHD, MU
```

Common range:

```text
0.05 <= z <= 1.14418
```

Residual definition:

```text
observed catalog distance modulus - no-b Model-A prediction
```

No-`b` Model-A residuals:

| Catalog | Median residual | Mean residual | RMSE | Median observed/model distance ratio |
|---|---:|---:|---:|---:|
| Union3 | +0.0512 mag | +0.0429 mag | 0.0539 mag | 1.0238 |
| Pantheon | +0.0808 mag | +0.0800 mag | 0.1664 mag | 1.0379 |
| DES | +0.2068 mag | +0.2383 mag | 0.3626 mag | 1.0999 |

**Status:** Useful no-`b` diagnostic; not a final empirical validation.

**Meaning:** The no-`b` Model-A distance curve remains structured and meaningful. Union3/Pantheon sit relatively close to it; DES remains elevated.

**Limitation:** Supernova distance-modulus catalogs are still interpreted through existing calibration pipelines. This test does not prove current supernova distances are wrong; it shows the STAM no-`b` curve is coherent and quantifies the catalog gap.

---

### 4.2 Supernova catalog discrepancies and `b` bridge

**Claim tested:** DES discrepancy against STAM may track an already-present DES-vs-Pantheon/Union3 catalog discrepancy.

Key results from the discrepancy audit:

```text
DES vs Pantheon empirical curve:
median ≈ +0.1268 mag
mean   ≈ +0.1588 mag
```

Locked historical bridge test using `b=354.95`:

```text
DES vs STAM locked-b:
median ≈ +0.1408 mag
mean   ≈ +0.1724 mag
```

Difference between the STAM DES residual and empirical DES-vs-Pantheon discrepancy:

```text
median difference ≈ +0.0160 mag
mean difference   ≈ +0.0135 mag
RMSE difference   ≈  0.0289 mag
correlation       ≈  0.9954
```

DES inferred bridge after Pantheon/Union3 reference offset:

```text
median b-like bridge ≈ 1446
mean b-like bridge   ≈ 1906
```

**Status:** Important catalog-family diagnostic.

**Meaning:** DES remains a real hurdle, but the STAM-DES miss closely resembles DES’s independent catalog discrepancy relative to Pantheon/Union3.

**Current interpretation:** Do not fit DES away with a DES-specific `b` unless an independent DES-specific catalog/systematics explanation is defined first.

---

### 4.3 `b` sensitivity

**Question tested:** Can a slight adjustment in `b` bring Union3, Pantheon, and DES closer together?

**Result:** No.

Fair reference-recalibrated scan:

```text
b=354.95 baseline balanced RMSE ≈ 0.16946 mag
best balanced RMSE near b≈631.6 ≈ 0.16939 mag
improvement ≈ 0.00007 mag
```

High DES-style `b` values can reduce DES only under an unfair fixed-offset diagnostic, while damaging Pantheon/Union3.

**Status:** `b` demoted from theory parameter to catalog-bridge diagnostic.

**Meaning:** Changing `b` is not the right lever for reconciling the catalogs. The model should move forward with no-`b` physical tests while retaining `b` as a record of catalog-distance discrepancy.

---

### 4.4 BAO geometric split test

**Claim tested:** If STAM separates geometric distance from accumulation/luminosity distance, BAO may prefer a different layer than supernovae.

Mappings tested with one nuisance scale `r_d_eff`:

```text
D_M = D_geo/(1+z)
D_M = D_adj/(1+z)
```

First-pass compressed BAO result:

```text
Combined all-z BAO:
    best mapping = D_geo/(1+z)
    chi²/dof ≈ 9.23

Combined z <= 1.6 BAO:
    best mapping = D_adj/(1+z)
    chi²/dof ≈ 2.48
```

A follow-up check found that `b=354.95` does not explain this flip. The flip persists when `b=0`.

**Status:** Mixed / under investigation.

**Meaning:** BAO is useful and nontrivial. The result does not cleanly validate the simple split, but it also does not reduce to a `b` artifact.

**Limitations:** Full covariance matrices were not used. Only transverse `D_M/r_d` points were tested. `D_H`, `D_V`, and a STAM-native `H(z)`/volume-distance rule remain open.

---

## 5. Exploratory work not yet promoted to repo-level claims

### Galaxy-scale accumulation / rotation curves

An exploratory SPARC-style reconstruction used:

```text
ε(r) = r|dA/dr| = 2v(r)²/c²
```

Preliminary scratch result:

```text
central density / central accumulation proxy correlated strongly with outer ε
Pearson r ≈ 0.789
Spearman ρ ≈ 0.834
```

**Status:** Shelved exploratory signal.

**Current wording:** Interesting pattern search only. Not a Model-A claim yet.

**Reason not promoted:** Needs better metadata, quality cuts, morphology handling, baryonic decomposition, and predefined tests.

---

## 6. Claims currently supported most strongly

The current strongest Model-A statements are:

```text
1. A(r)=Rs/r correctly normalizes weak-field spherical acceleration.

2. Δt=(1/c)∫A ds reproduces the weak-field Shapiro logarithmic delay structure.

3. With dτ/dt≈1-A/2, Earth-surface/geoid vs GPS-orbit accumulation reproduces the GPS gravitational clock component and net factory-offset scale when orbital kinematics are included.

4. A=1 is exactly the spherical threshold where v_escape=c.

5. A>1 is exactly the spherical over-threshold region where escape would require v_escape>c.

6. Gravitational waves can remain locally c-speed in Model-A if GW and EM signals share the same accumulation traversal coupling.

7. The no-b supernova curve is coherent and structured; b is best treated as a catalog bridge, not core physics.
```

---

## 7. Claims that remain open or incomplete

| Item | Status |
|---|---|
| Derive the cosmological no-`b` distance law from a deeper field equation | Open |
| Explain or eliminate catalog bridge terms from independent data | Open |
| Define a complete STAM-native BAO rule for `D_M`, `D_H`, and `D_V` | Open |
| Test full BAO covariance likelihoods | Open |
| Derive a full gravitational redshift/clock relation beyond weak-field mapping | Partly tested by GPS, still incomplete |
| Establish Solar-System precision bounds / PPN equivalence | Open |
| Model lensing normalization and photon paths beyond first-order identities | Open |
| Model observed black-hole shadow size, not just horizon radius | Open |
| Define strong-field interior dynamics for `A>=1` | Open |
| Derive gravitational-wave propagation from a wave equation on/through `A` | Open |
| Test galaxy rotation curves with predefined STAM metrics | Exploratory only |
| Test standard sirens as independent distance observables | Open |
| CMB / early-universe model | Open |

---

## 8. Hard falsification targets

Model-A becomes weaker or fails if any of the following persist after fair, predeclared tests:

1. Local acceleration, delay, redshift, and lensing require mutually inconsistent normalizations of `A`.
2. GPS/clock tests fail once full geoid/orbit corrections are modeled from STAM assumptions.
3. Shapiro delay departs from observed weak-field bounds in Solar-System regimes.
4. The `A=1` horizon threshold cannot be made dynamically or invariantly meaningful.
5. The no-`b` distance structure cannot classify independent observables into consistent geometric vs accumulation layers.
6. BAO, angular-diameter, standard-siren, and time-delay distances all require unrelated after-the-fact bridge terms.
7. DES/Pantheon/Union3-style supernova discrepancies cannot be separated from STAM residuals under predeclared catalog/systematics tests.
8. Gravitational-wave observations require GW propagation to have a different local speed or different accumulation coupling without a defined STAM reason.
9. Galaxy rotation-curve reconstructions show no stable or predictive accumulation pattern after predefined cuts.
10. The model cannot produce testable predictions beyond identities shared with standard weak-field formulas.

---

## 9. Current research priority

The next best non-distance priorities are:

```text
1. Full GPS/clock refinement:
   include geoid, eccentricity correction, and satellite-specific terms.

2. Gravitational redshift tests:
   Pound-Rebka-style and modern optical-clock height tests using dτ/dt≈1-A/2.

3. Lensing/path test:
   determine whether the same A normalization gives deflection and delay without extra factors.

4. Horizon/shadow extension:
   separate A=1 horizon radius from photon-sphere/shadow predictions.

5. BAO refinement:
   define a STAM-native BAO-facing distance rule before judging the split result.
```

The next best distance priority is:

```text
Run no-b Model-A against independent distance indicators before using b-like bridge terms.
```

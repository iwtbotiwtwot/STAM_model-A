# Claims and Status — STAM Model-A

This document records the current working structure of STAM Model-A: what the model defines, what has been checked numerically, what remains diagnostic, and what still needs stronger testing.

It is a project ledger for an active model-building effort. The current working posture is constructive: preserve what has checked out, separate what belongs in different measurement ledgers, and keep pushing the parts that remain open.

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
| **Ledger diagnostic** | A test or thought experiment that separates clock/process, path/traversal, observable-distance, or catalog-mapping effects. |

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

## 3A. Clock / path / traversal ledger

### 3A.1 Why the ledger matters

**Claim type:** Ledger diagnostic / active Model-A accounting rule

Recent Planet A+ / Planet E thought experiments clarified that Model-A should not treat every timing or distance discrepancy as one generic “A delay.”

The working ledger separates:

```text
A_clock      = local A condition of the clock/source/process being compared
A_observer   = local A condition of the observer/comparison frame
A_between    = path/traversal A along the signal route
A_traverse   = effective translation of local motion into shared observable distance
```

This distinction is now central to avoiding double-counting.

**Current rule:**

```text
GPS-style clock/process correction ≠ Shapiro-style path/traversal delay ≠ observable-distance traversal.
```

**Status:** Retained as a working Model-A accounting rule.

---

### 3A.2 GPS and Shapiro are not contradictions

**Claim type:** Ledger diagnostic

The GPS and Shapiro checks now occupy different measurement ledgers:

```text
GPS:
    A affects clock/process comparison.
    Use local A-clock mapping:
        dτ/dt ≈ 1 - A/2

Shapiro:
    A affects signal/path traversal.
    Use path integral:
        Δt = (1/c)∫A ds
```

**Current interpretation:**

```text
The same A-language can describe both effects, but they are different projections of A.
GPS corrects the clock.
Shapiro corrects the road.
```

This resolves the earlier concern that GPS clock adjustment might undermine the Shapiro/path-delay interpretation. Instead, GPS strengthens the need for a clean ledger.

**Status:** Ledger diagnostic retained.

---

### 3A.3 Magic-bell diagnostic

**Claim type:** Thought experiment / ledger diagnostic

The magic bell is not a physical device. It is a diagnostic tool that removes the signal path between two systems.

If the magic bell is used:

```text
A_between is removed from the comparison.
```

Then any remaining discrepancy must come from:

```text
A_clock ≠ A_observer
```

The diagnostic cases are:

```text
Case 1:
A_clock = A_observer
A_between = 0
→ no discrepancy

Case 2:
A_clock ≠ A_observer
A_between = 0
→ pure clock/process comparison

Case 3:
A_clock = A_observer
A_between ≠ 0
→ pure traversal/path delay

Case 4:
A_clock ≠ A_observer
A_between ≠ 0
→ mixed real-world case
```

**Status:** Retained as a core explanatory thought experiment.

---

### 3A.4 TE1 — fixed path vs clock process

**Claim tested:**

A fixed path delay shifts arrival time, while a clock/process mismatch changes tick spacing.

Using:

```text
clock/process:
dτ/dt ≈ 1 - A/2

path/traversal:
Δt_path = (1/c)∫A ds
```

TE1 tested cases with:

```text
same A / no path
magic bell / clock only
path only / same clocks
clock + path
observer higher A
```

**Representative result:**

```text
Magic bell / clock only:
A_source = 0.001
A_observer = 0
A_between = 0
source vs observer rate ≈ 0.9995
received tick interval ≈ 1.0005 s

Path only / same clocks:
A_source = 0
A_observer = 0
A_between = 0.001
one-way path delay ≈ 1.0 s
received tick interval remains ≈ 1.0000 s
```

**Interpretation:**

```text
A fixed A_between creates a one-way delay but does not by itself change repeated tick spacing.
Clock/process mismatch changes tick spacing.
Changing A_between between emissions can also change received tick spacing.
```

**Status:** Numerically checked as ledger diagnostic.

---

### 3A.5 TE2 — moving Bob / changing A-path

**Claim tested:**

If Bob moves into higher A, and the path back to Alpha+ also worsens, then Alpha+ can see received tick spacing stretch from both:

```text
1. Bob's clock/process comparison
2. changing signal/path traversal delay
```

Representative final received intervals from the toy ledger:

| Scenario | Final received interval |
|---|---:|
| Weak clock only, changing path | 19.6060 |
| Horizon-style clock only, changing path | 20.9501 |
| Horizon-style clock + traversal efficiency | 49.4106 |

**Interpretation:**

```text
TE1 showed fixed A_between shifts arrival time but not tick spacing.
TE2 shows changing A_between stretches received tick spacing.
Rising Bob-local A also stretches the clock/process comparison.
```

Barstool version retained for intuition:

```text
Bob's longer becomes Alpha+'s even longer.
```

**Status:** Numerically checked as thought-experiment ledger.

**Caution:** This is not a full black-hole geodesic calculation.

---

### 3A.6 TE3 — Bob stopwatch vs Alpha+ time

**Claim tested:**

Bob's local stopwatch can remain valid to Bob while Alpha+ translates the same trip into a different time.

Weak-field clock translation:

```text
dt_Alpha = dτ_Bob × (1 - A_Alpha/2)/(1 - A_Bob/2)
```

Toy traversal term:

```text
T(A) = 1 - A
```

Representative results:

| Scenario | Bob stopwatch | Alpha+ clock equivalent | Alpha - Bob | Bob local distance | Ledger observable distance |
|---|---:|---:|---:|---:|---:|
| Weak A | 200.0 s | 200.086 s | 0.086 s | 200.0 | 199.828 |
| Moderate A | 200.0 s | 200.812 s | 0.812 s | 200.0 | 198.386 |
| Stronger sandbox A | 200.0 s | 208.559 s | 8.559 s | 200.0 | 183.968 |

**Interpretation:**

```text
Bob can honestly report his own stopwatch time.
Alpha+ can honestly translate that same trip into a different time.
If distance inference depends on timing, source process, or received signal intervals, the relevant clock matters.
```

**Status:** Numerically checked as thought-experiment ledger.

**Caution:** Stronger-A case is sandbox only; weak-field clock formula is not a horizon model.

---

### 3A.7 Current black-hole / Bob distinction

**Claim type:** Boundary / open extension

Model-A currently distinguishes its working interpretation from the standard local-frame phrasing this way:

```text
Current relativity:
Bob locally crosses local distance normally.
Alpha+ sees Bob slow/freeze in the remote/coordinate/casual-access description.

Model-A working view:
Bob's internal clock/processes remain locally normal to Bob,
but rising A reduces Bob's effective traversal through shared observable distance.
Alpha+ then sees that already-reduced traversal through clock and path translation.
```

Clean current wording:

```text
Bob's internal clock stays normal to Bob,
but his ability to convert local motion into observable-distance progress collapses as A rises.
```

**Status:** Boundary / open extension.

**Open:** Needs a formal strong-field traversal rule and comparison against standard black-hole observables.

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

### 4.4 Updated interpretation of `b`: ledger residue

**Claim type:** Ledger diagnostic / catalog diagnostic

Recent TE-ledger tests update the interpretation of the historical `b` bridge.

Current working definition:

```text
b = unresolved ledger residue
```

More explicitly:

```text
b_total
≈ b_path/traversal_mismatch
 + b_observable-distance_translation
 + b_source/catalog_mapping
```

The historical `b` term is no longer treated as a fundamental forward physical distance term. Instead, it is treated as a diagnostic residue that appears when multiple A-ledger effects are forced into one distance bucket.

**Current status:**

```text
b is not fully derived.
b is becoming decomposable.
```

**Retained rule:**

```text
Do not put b back into the preferred forward physical distance law.
Use b as a diagnostic for missing ledger structure and catalog mapping residue.
```

---

### 4.5 Script 31 — distance ledger decomposition

**Claim tested:**

Whether TE-ledger terms reduce the old b-like catalog offset.

Models compared included:

```text
M1 path-only:
D = D_geo + λ_path D_excess

M3 path + traversal:
D = D_geo + λ_path D_excess + λ_traverse D_geo A_local
```

The broader diagnostic form was:

```text
D_model =
D_geo
+ λ_path · D_excess
+ λ_clock · D_geo · A_avg
+ λ_traverse · D_geo · A_local
+ b_mu
```

**Main result:**

Adding the traversal ledger proxy reduced the b-like residual for Pantheon, Union3, and Pantheon+Union3.

Representative shrinkage:

| Dataset | Path-only b_mu | Ledger b_mu |
|---|---:|---:|
| DES | 0.2174 | 0.0861 |
| Pantheon | 0.0941 | 0.0185 |
| Union3 | 0.1517 | -0.0225 |
| Pantheon + Union3 | 0.0992 | 0.0174 |

**Interpretation:**

```text
The historical b term appears to have been partly catching missing path + observable-traversal ledger structure.
```

**Status:** Catalog / ledger diagnostic.

---

### 4.6 Script 32 — ledger stability and transfer

**Claim tested:**

Whether the M3 path + traversal structure transfers across catalogs and redshift splits, or whether Script 31 merely overfit.

Whole-dataset M3 coefficient family for Pantheon / Union3 / Pantheon+Union3:

```text
λ_path mean/std      ≈ -3.789 / 0.327
λ_traverse mean/std  ≈  2.703 / 0.205
```

Strong transfer example:

```text
Train Pantheon → test Union3

M1 path-only test RMSE       ≈ 0.1163
M3 path + traversal test RMSE ≈ 0.0435
```

Watched split-zone result:

| Split | Median M3 - M1 test RMSE | Fraction where M3 better |
|---|---:|---:|
| z = 0.300 | -0.0075 | 60% |
| z = 1/3 | -0.0016 | 50% |
| z = 0.350 | -0.0237 | 70% |

**Interpretation:**

```text
The TE ledger is not merely an in-catalog fit.
The Pantheon/Union family shares a stable path + traversal structure.
DES remains a different catalog-family stress case.
```

**Status:** Catalog / ledger diagnostic.

---

### 4.7 Script 33 and 33.33 — b decomposition

**Claim tested:**

Whether `b` can become more than a fit by decomposing it into structured ledger pieces.

Script 33 measured:

```text
b_explained_by_ledger
=
b_path_only_required
-
b_after_path_plus_traversal
```

Script 33.33 consolidated the result.

Scoreboard:

| Dataset | Path-only b_mu | Ledger b_mu | Explained b_mu | Fraction explained | Status |
|---|---:|---:|---:|---:|---|
| DES | 0.2911 | 0.1853 | 0.1058 | 36.4% | weak |
| Pantheon | 0.0991 | 0.0243 | 0.0748 | 75.5% | strong |
| Union3 | 0.1517 | -0.0225 | 0.1742 | 85.2% | strong |
| Pantheon + Union3 | 0.1034 | 0.0231 | 0.0803 | 77.6% | strong |
| DES + Pantheon common | 0.1232 | 0.0167 | 0.1066 | 86.5% | strong |
| All three | 0.1448 | 0.0110 | 0.1337 | 92.4% | strong |

Pantheon / Union family summary:

```text
mean fraction of b explained ≈ 79.4%
mean explained b_mu          ≈ 0.1098 mag
λ_path mean/std              ≈ -3.789 / 0.327
λ_traverse mean/std          ≈  2.703 / 0.205
```

Cross-catalog examples:

| Train | Test | Path-only b needed | Ledger b needed | Fraction explained |
|---|---:|---:|---:|---:|
| Pantheon | Union3 | 0.0471 | 0.0018 | 96.2% |
| Union3 | Pantheon | 0.1615 | 0.0030 | 98.1% |
| Pantheon + Union3 | DES | 0.2765 | 0.1807 | 34.6% |
| DES | Pantheon + Union3 | 0.1138 | 0.0235 | 79.3% |

**Current conclusion:**

```text
b is not fully derived,
but it is no longer just a fit knob.
```

Best current wording:

```text
b is becoming decomposable.
It appears to be ledger residue created when path accumulation and observable-distance traversal were forced into one distance bucket.
```

**Status:** Catalog / ledger diagnostic, strengthened.

**Caution:** `b_mu` is a magnitude-offset diagnostic and is not identical in unit to the older historical distance-style `b` constants.

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


**Updated ledger interpretation:** The bridge term now appears to have been acting as an unresolved ledger residue. Recent scripts suggest that much of the Pantheon/Union-style `b` behavior is replaced by structured path + observable-distance traversal terms, while DES remains a separate catalog-family stress case.

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

## 6A. Recent BAO / SN observable-layer diagnostics

### 6A.1 BAO λ layer mapping

**Claim tested:**

BAO may see only part of the accumulation-adjusted layer rather than pure geometric or full adjusted distance.

Diagnostic mapping:

```text
D_BAO(z) = D_geo(z) + λD_excess,0(z)
D_M,BAO = D_BAO/(1+z)
```

with:

```text
λ = 0 → pure geometric
λ = 1 → full adjusted
```

**Representative result:**

```text
combined all-z:
best constant λ ≈ 0.383
χ²/dof ≈ 2.741
r_d_eff ≈ 117.887 Mpc

z ≤ 1.6:
best constant λ ≈ 0.659
χ²/dof ≈ 1.325
r_d_eff ≈ 124.062 Mpc
```

**Interpretation:**

```text
BAO appears to see a partial accumulation layer.
```

**Status:** Catalog diagnostic.

---

### 6A.2 Refined BAO no-b mapping

**Claim tested:**

Smooth λ(z) models improve BAO mapping without restoring b.

Compared:

```text
geo
adj0
constant λ
linear λ(z)=a+bz
saturating λ(z)=a+bz/(1+z)
```

**Representative all-z result:**

```text
linear λ(z):
χ²/dof ≈ 1.470
AIC ≈ 16.29
BIC ≈ 17.20
r_d_eff ≈ 131.881 Mpc
λ(z) ≈ 1.228 - 0.296z
```

Leave-one-out:

```text
saturating λ best ≈ 1.873 χ²/point
linear λ ≈ 2.021
constant λ ≈ 4.037
geo ≈ 10.572
adj0 ≈ 14.120
```

**Interpretation:**

```text
BAO improves significantly with a smooth partial-accumulation observable layer.
```

**Status:** Catalog diagnostic.

---

### 6A.3 Cross-observable λ sniff test

**Claim tested:**

Whether SN and BAO recover a shared λ curve.

Definition:

```text
λ = (D_observable - D_geo)/D_excess,0
```

**Result summary:**

```text
BAO independently recovers a smooth partial-accumulation function.

SN λ is more catalog-scale sensitive and does not yet provide one obvious shared λ curve with BAO.
```

**Interpretation:**

```text
The better next step is observable-coupling / parent-ledger testing, not forcing all observables to share one identical λ.
```

**Status:** Catalog diagnostic.

---

### 6A.4 z ≈ 0.30 to 0.35 flag zone

**Claim type:** Catalog diagnostic / active target zone

Early STAM/SU work flagged:

```text
z ≈ 0.30
```

as a crossing / anchor marker.

Recent tests updated the interpretation:

```text
z ≈ 0.30 is not supported as a simple onset point.
```

The stronger current target zone is:

```text
z ≈ 0.30–0.35
```

with special attention to:

```text
z = 0.300
z = 1/3 = 0.333333
z = 0.350
```

Current interpretation:

```text
0.30 planted the original flag.
1/3 is the cleaner mathematical suspect.
0.35 is the strongest nearby empirical split in several transfer/pivot tests.
```

This region is now treated as a possible:

```text
ledger transition / observable-mapping pivot / path-to-traversal balance zone
```

not as a proven onset.

**Status:** Active catalog diagnostic.

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


## 8. Active boundaries and open extensions

These items remain active areas of development rather than completed Model-A results:

```text
1. Full strong-field metric formulation.
2. Black-hole interior / A>1 dynamics beyond the threshold identity.
3. First-principles derivation of the no-b cosmological distance law.
4. Full BAO/CMB distance-observable framework.
5. Completed galaxy-rotation reconstruction from real datasets.
6. Full GPS operational navigation model.
7. Full gravitational-wave waveform theory.
8. Independent explanation of DES/Pantheon/Union3 catalog-family behavior.
9. Formal derivation of the TE ledger terms from first principles.
10. Determine whether the decomposed b_mu diagnostic can be converted into a unit-consistent physical b relation.
```

The historical supernova bridge term is kept as catalog-comparison history, not as a claim-boundary item.

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
9. Convert the current b-as-ledger-residue diagnostic into a derivable physical relation, if possible.
10. Define which clock/path/traversal ledger each cosmological observable actually samples.
```

---

## 10. Current working position

Current project stance:

```text
A is the core Model-A variable.

The no-b distance law remains the preferred forward-testing cosmological form.

Historical bridge results are retained as catalog-comparison history, but the interpretation has improved:
    b is now treated as unresolved ledger residue, not as a fundamental forward-model constant.

Non-distance tests currently support the internal consistency of Model-A in:
    local gravity,
    Shapiro delay,
    GPS weak-field clocks,
    horizon threshold,
    gravitational-wave local-speed consistency,
    and mass-estimator identities.

The TE ledger has strengthened the model's accounting structure by separating:
    clock/process effects,
    path/traversal effects,
    observable-distance traversal,
    and catalog/mapping residue.

Cosmological interpretation remains the largest open area, but the distance problem is now better organized:
    different observables may sample different A-ledgers rather than one universal distance bucket.
```

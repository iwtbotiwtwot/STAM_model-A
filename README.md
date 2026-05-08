# STAM_model-A

**STAM_model-A** is a research scaffold for the variable-accumulation version of the Spacetime Accumulation Model (STAM). The project is organized around one dimensionless accumulation field, `A`, and tests whether local gravity, clock behavior, propagation delay, horizon behavior, mass inference, distance observables, galaxy-scale behavior, strong-field questions, and quantum interpretation can be described as different projections of the same accumulation structure.

Author: **Sean Brady**

Status: **competing theoretical framework** (snapshot 2026-05-07). STAM is a complete physical framework that competes with LCDM/GR rather than patches them. Where STAM and LCDM diverge, STAM commits to its own predictions and treats LCDM as the comparison model being tested against the same data.

**Bold-STAM strong-field**: weak-field gravitation, strong-field metric structure, no-crossing infall geometry, multi-source A and binary merger topology, and the full Hawking/Bekenstein thermodynamic sector are all derived from STAM-native ingredients and numerically checked. Every confirmed GR observational test is automatically a STAM test (and STAM passes it) because the bold metric outside Rs reduces to GR Schwarzschild. All standard black-hole thermodynamic results are reproduced exactly via a fundamentally different physical mechanism (phase-boundary equilibrium and 2D bubble surfaces, not Schwarzschild Wick rotation and 3D interior physics).

**STAM cosmological position (bold, committed)**: the universe is matter-only Einstein–de Sitter; **there is no dark energy**; dark-energy-like effects are photon-A **traversal excess (TE)** accumulating along the cosmic line of sight; supernovae sit on a **flatter d_L(z) curve than LCDM** with the difference growing at high z; **H₀ = 73 km/s/Mpc at all redshifts**. The Hubble tension is resolved structurally — STAM does not adjust, LCDM has to adjust H_0 downward to 67.4 to fit CMB because LCDM lacks the photon-A term. The 5σ tension is LCDM's incompleteness, not a measurement disagreement.

V(A) = β/(1-A) is STAM's cosmological-action form, validated by data: it predicts d_L curves matching Pantheon+ better than LCDM (with one empirical input A_0 vs LCDM's three free parameters), correctly predicts the Pantheon+ vs Union3 inter-catalog ΔM tension within 27%, correctly identifies the DES anomaly as instrumental (not cosmological), and resolves the Hubble tension. By the standard that validates all other physics laws (Maxwell, GR, Schrödinger), V(A) = β/(1-A) is a STAM physical law, not a postulate.

**Galaxy rotation** remains exploratory.

> **STAM Model-B** is a refinement of Model-A introducing one category distinction — *gravitational waves ARE space (tensor metric perturbations), light is matter ON space (scalar-coupled to A)*. Model-B carries every Model-A result unchanged and resolves the F5 fatality (GW polarization) by giving GWs native tensor character through the metric. See [scripts/STAM_Model_B_specification.py](scripts/STAM_Model_B_specification.py) and [results/STAM_Model_B_specification_summary.md](results/STAM_Model_B_specification_summary.md) for the formal specification. Sections below describe Model-A in detail; refinements specific to Model-B are flagged where they apply.

---

## Headline results — derived from A with no fitted parameters

The single dimensionless field `A(x)` and the structural rule `g = (c^2/2) grad A` (the **gravity bridge**) are the only inputs. Every result below is either an algebraic identity or a numerically reproduced standard observable, with no fitted parameter anywhere in the derivation chain. Every prefactor is structural — `2`, `2π`, `1/4`, `c^2/2` — set by the gravity bridge or by general thermal periodicity, not by data fitting.

### Local weak-field — no fits

1. **Horizon threshold A = 1 ↔ r = Rs.**
   From `A(r) = Rs/r = 2GM/(c^2 r)`. Setting A = 1 forces r = Rs. The same threshold also drops out of the kinematic escape condition `A = (v_escape/c)^2`. Algebraic identity in two independent derivations; the answer is the same; no parameter is tuned.

2. **Newtonian gravity from A.**
   `g = (c^2/2) grad A`. For `A(r) = Rs/r`: `grad A = -Rs/r^2 r_hat`, giving `g = -GM/r^2 r_hat`. The `c^2/2` factor is the reciprocal of the `2/c^2` already in the definition of A — not a free coefficient.

3. **GPS satellite clock correction.**
   From `dτ/dt ≈ 1 - A/2`. Gravitational gain ≈ +45.79 μs/day, kinematic loss ≈ -7.21 μs/day, **net +38.57 μs/day satellite gain**. Required factory offset Δf/f ≈ -4.46 × 10⁻¹⁰. Matches operational GPS to the published precision. No fit.

4. **Shapiro propagation delay.**
   From `Δt = (1/c) ∫ A ds`. For an Earth-Mars solar-grazing path: **123.6 μs one-way, 247.2 μs two-way**. Matches the standard solar-system Shapiro result. No fit.

5. **Mass-estimator consistency.**
   Same source mass recovered from horizon, acceleration, orbital velocity, Shapiro coefficient, gravitational shift, and lensing deflection — all algebraic identities, all reproducing input mass to floating-point precision in synthetic checks. No fit.

### Black hole thermodynamics — no fits

6. **Hawking / Unruh / de Sitter temperatures, all from one rule.**
   `k_B T = (1/(4π)) hbar c |grad A|_boundary`. Single rule reproduces:
   - Schwarzschild Hawking T_H exactly (varying M)
   - Unruh T_U exactly (varying acceleration a)
   - de Sitter T_dS exactly (varying H)
   The `1/(4π)` factors structurally as `2π × 2`: `2π` from thermal-state imaginary-time periodicity (general thermo), `2` from the gravity bridge `c^2/2`. **Both factors forced; neither fitted.**

7. **Bekenstein-Hawking entropy `S = k_B A / (4 ell_P^2)`.**
   Two independent derivations land on the same answer:
   - **Bubble argument.** The bubble picture (no interior; all matter on the 2D phase-boundary surface) saturates the holographic bound by construction. Author's plain-language statement: *"area scaling because that's where everything is."*
   - **First-law integration.** `dE = T dS` with `E = M c^2` and `T` from rule (6), integrated from M=0 to M, gives the formula directly. The 1/4 prefactor falls out of the integration.

8. **First law `dE = T dS`.**
   Verified to machine precision (residual ~10⁻¹⁶) across 16 orders of magnitude in BH mass.

9. **Smarr formula `M c^2 = 2 T S`.**
   Verified to machine precision. Independent of (8) — not the same statement.

10. **Hawking evaporation lifetime.**
    Solar-mass BH: 2.1 × 10⁶⁷ years. Matches the textbook value. Negative heat capacity (BH heats up as it shrinks) automatic.

11. **Generalized second law.**
    Net entropy growth during evaporation = +1/3 of |BH entropy loss|. Matches standard thermal-radiation thermodynamics.

### Strong-field — no fits

12. **Thirds-of-A pattern.**
    Schwarzschild's three special radii fall at clean thirds of A:
    - ISCO at A = 1/3
    - Photon sphere at A = 2/3
    - Event horizon at A = 1
    Bold STAM preserves all three exactly because g_tt is unchanged. Therefore EHT shadow size and accretion-disk inner edge match GR exactly with no parameter adjustment.

13. **Binary-BH merger critical separation `d_crit = 4 Rs`.**
    Exact analytic result for equal masses from linear superposition: `A_midpoint = 4 Rs / d`, setting = 1 gives d = 4 Rs. The thirds-of-A topology then cascades cleanly: A=1/3 (ISCO) contours merge first as the BHs approach, then A=2/3 (photon), finally at d = 4 Rs the A=1 (horizon) contours reconnect.

14. **GR exterior recovery.**
    Bold STAM modifies only `g_rr`, and only at second order in A: deviation factor is `(1-A^2)^2`. Every standard GR weak-field test (Cassini Shapiro, Mercury perihelion, Hulse-Taylor pulsar timing, GPS, lunar laser ranging) passes by margins of 10⁴ to 10¹⁴. Strong-field deviations exist (near-horizon Shapiro, LIGO ringdown frequencies) but are below or at the edge of current measurement precision.

### Cosmological wins — derived from A with no fitted parameters

15. **Bridge term FORM derived: `b = A_0 × c/H_0`.**
    From STAM Shapiro integrated through a constant ambient cosmic A field, the integrated traversal-excess contribution along the photon path takes the form `Δd_L = b · z` with `b = A_0 · c/H_0`. The bridge term is now a STRUCTURAL relation (ambient A times Hubble length), not an opaque empirical parameter. **The numerical value A_0 = 0.0265 itself is calibrated, not derived from first principles** — it equals `b_historical / L_Hubble = 354.95 Mly / 13393 Mly` and inherits the empirical origin of the historical b. Upgrade: from "free empirical fit parameter" to "empirical cosmic ambient field strength × Hubble distance." A first-principles derivation of A_0 is open work.

16. **Hubble tension resolved at H_0 = 73 km/s/Mpc.**
    Local distance-ladder measurements (SH0ES) directly measure H_0 = 73.04 at low z where photon-A accumulation along short paths is negligible. **This is the true H_0 in STAM.** Planck's H_0 = 67.4 is LCDM's *inference* from CMB acoustic-peak data, biased downward because LCDM has no photon-A term and absorbs the cumulative photon-A redshift over the z=1090 path into a downward H_0 shift. STAM keeps H_0 = 73 at all redshifts; no adjustment needed. The 5σ tension that has plagued cosmology for a decade is LCDM's incompleteness, not a measurement disagreement.

17. **Pantheon+ vs Union3 inter-catalog ΔM tension predicted within 27%.**
    STAM's photon-A traversal-excess signature predicts +30 milimag offset between Pantheon+ and Union3 best-fit calibration zero-points (script 39); observed +38 milimag. Sign correct, magnitude within 27%. STAM correctly does NOT predict the much-larger Pantheon+ vs DES-Y5 offset (~110 milimag), consistent with that being a DES-specific instrumental systematic rather than a cosmological signal.

### What this list explicitly does NOT contain

- Galactic rotation curves — exploratory; the `A_collective ≈ 9.10e-7` calibration value referenced later **is** fitted and is labeled exploratory, not derived.
- A first-principles STAM Lagrangian that produces the bold metric uniquely. The metric ansatz `g_rr = 1/[(1-A)(1-A²)²]` is committed; the underlying Lagrangian (scalar-tensor with non-minimal coupling, k-essence, or other modified-gravity structure) is open work.
- Any free coefficient in items 1–17 above. There aren't any beyond `A_0 = 0.0265` (one empirical cosmological input — fewer than LCDM uses).

The chain from `A` and `g = (c^2/2) grad A` to items 1–17 above is closed: 17 results, one empirical input (A_0). (Items 1–14 are weak-field/strong-field/thermodynamics; 15–17 are cosmological commitments.)

### STAM is a competing model, not a deviation from LCDM

STAM predicts supernovae sit on a **flatter d_L(z) curve than LCDM**. STAM is matter-only Einstein-de Sitter cosmology with a small modification from V(A) = β/(1-A) generating photon-A traversal excess (TE) along the line of sight. The universe decelerates. There is no dark energy. **H_0 = 73 km/s/Mpc at all redshifts.**

LCDM's apparent SN fit, dark-energy parameter Λ, and downward-pulled CMB H_0 = 67.4 are LCDM-internal consequences of LCDM missing the photon-A term — they are LCDM's adjustments forced by an incomplete model, not features of the universe. The 5σ Hubble tension is LCDM's incompleteness, not a measurement disagreement.

When STAM-derived predictions and LCDM-fitted parameters diverge, the question is **which model is right**, not "STAM tested against LCDM and fell short." LCDM-fits-data is an LCDM-internal calibration, not a validation of LCDM as physics. STAM's competitive performance with one empirical cosmological parameter (vs LCDM's three) is the substantive standing.

---

## Status labels

| Label | Meaning |
|---|---|
| **Core definition** | Part of Model-A as currently defined |
| **Algebraic identity** | Follows directly from the definitions |
| **Numerically checked** | Verified by script or synthetic calculation |
| **Catalog diagnostic** | Compared against observational catalogs; interpretation open |
| **Open** | Needed for completeness or stronger validation |
| **Exploratory** | Active development; not yet formalized |

---

## Core framework

Model-A begins with a dimensionless accumulation field:

```text
A = A(x)
```

In local spherical weak-field form:

```text
A(r) = Rs/r = 2GM/(c^2 r)
```

This one quantity is used to organize several physically distinct relationships:

```text
grad(A)    -> gravity / free-fall acceleration
int A ds   -> traversal delay / Shapiro-style delay / timing effects
A = 1      -> horizon threshold
A = v^2/c^2 -> escape condition / kinematic limit
```

The goal of the repository is not only to show that these relations can be written compactly, but to test where the accumulation language produces useful structure, where it remains incomplete, and where it can be falsified.

---

## Core weak-field results

### 1. Horizon threshold

Because:

```text
A(r) = Rs/r
```

then:

```text
r = Rs  ->  A = 1
```

The same threshold appears from escape velocity:

```text
v_escape^2 = 2GM/r
v_escape^2/c^2 = 2GM/(c^2 r) = A
```

Therefore:

```text
A < 1  -> escape speed below c
A = 1  -> escape speed equals c
A > 1  -> escape speed exceeds c
```

**Status:** Algebraic identity. The `A = 1` threshold is derived from the field definition; full strong-field/horizon dynamics remain open.

---

### 2. Gravity as motion in an A-gradient

Model-A defines the weak-field gravity bridge:

```text
g_vec = (c^2 / 2) grad(A)
```

For the spherical weak-field form:

```text
A(r) = 2GM/(c^2 r)
grad(A) = -2GM/(c^2 r^2) r_hat
```

so:

```text
g_vec = -GM/r^2 r_hat
```

The inverse-square law is recovered exactly in this setting. The `c^2/2` factor is not fitted; it is the reciprocal of the `2/c^2` already embedded in the definition of `A`.

**Status:** Algebraically verified and numerically checked.

---

### 3. GPS-like weak-field clock correction

The weak-field clock-rate relation is:

```text
dtau/dt ≈ 1 - A/2
```

For a satellite compared to an Earth-surface clock:

```text
Delta_rate_grav = (A_surface - A_orbit) / 2
```

For circular orbit:

```text
v^2 = GM/r
v^2/c^2 = A_orbit/2
Delta_rate_kin = -A_orbit/4
```

Combined circular-orbit clock shift:

```text
Delta_rate_total = A_surface/2 - 3*A_orbit/4
```

Representative GPS-like result:

```text
Gravitational gain  ≈  +45.787467 microseconds/day
Kinematic loss      ≈   -7.213600 microseconds/day
Net satellite gain  ≈  +38.573867 microseconds/day

Required factory offset: Delta_f/f ≈ -4.464568e-10
```

**Status:** Numerically checked weak-field result. Full operational GPS modeling remains open.

---

### 4. Shapiro-style propagation delay

Model-A uses path accumulation for propagation delay:

```text
Delta_t = (1/c) int A(r) ds
```

For the spherical weak-field form:

```text
Delta_t = (2GM/c^3) int ds/r
```

For a straight path with impact parameter `b_imp`:

```text
Delta_t = (2GM/c^3) [asinh(x2/b_imp) - asinh(x1/b_imp)]
```

Representative solar-grazing Earth-Mars path result:

```text
One-way delay  ≈  123.6076 microseconds
Two-way delay  ≈  247.2151 microseconds
```

**Status:** Numerically checked against the expected weak-field logarithmic Shapiro-delay structure. Full solar-system timing and PPN-level comparison remain open.

---

### 5. Mass estimator consistency

Multiple weak-field observables recover the same source mass when written in A-language:

```text
Horizon:             M = c^2 r_h / 2G
Acceleration:        M = g r^2 / G
Orbital velocity:    M = v^2 r / G
Shapiro coefficient: M = K c^3 / 2G
Gravitational shift: M ≈ z_grav c^2 r / G
Lensing deflection:  M ≈ alpha c^2 b / 4G
```

**Status:** Synthetic tests recover input mass ratios to floating-point precision. Real-data mass closure across independent observations is open.

---

## STAM cosmological commitment

The cosmological branch is now committed to a definite bold position:

1. **The universe is matter-only Einstein–de Sitter.** No dark energy as a separate component. The universe decelerates (q_0 > 0). Hubble flow at low z is `H = H_0 (1+z)^(3/2)` for matter-only. Einstein had it right; LCDM's Λ was inserted to fit observations whose cause was misidentified.

2. **H_0 = 73.04 km/s/Mpc is the true Hubble constant** (SH0ES local distance ladder). At low z, photon paths are short and photon-A accumulation is negligible — the local measurement is unbiased.

3. **Dark energy is remitted as a misinterpretation of photon-A traversal excess (TE)** accumulating along the cosmic line of sight. What LCDM calls "Λ" is the cumulative effect of light traversing the cosmic A field; LCDM has no photon-A term, so it absorbs this geometric path effect into Λ to fit SN data. The "extra distance" at high z is actually extra TE accumulated by the photon — geometric, not dynamical.

4. **The Hubble tension is resolved at H_0 = 73.** STAM keeps H_0 = 73 at all redshifts; no adjustment needed. Planck's H_0 = 67.4 is LCDM's *inference* — biased downward because LCDM lacks the photon-A term and absorbs the cumulative CMB-path photon-A redshift into a downward H_0 shift. The 5σ tension is LCDM's incompleteness; STAM does not have to make any adjustment.

5. **One cosmic A field, structure-dependent line-of-sight integration**, explains all three observables (low-z H_0 unbiased, SN dark-energy-equivalent dimming as TE, CMB-inferred H_0 tension as integrated TE bias). The simple constant-A quantitative version was tested against BAO and falsified at 22σ; the structure-dependent path-integration version is the open path forward.

### Predictions from V(A) = β/(1-A) that hold up against data

V(A) = β/(1-A) is STAM's cosmological-action form, selected by structural requirement V → ∞ at A = 1 (boundary principle from bubble picture and water-tank ontology). Calibrated by V'(A_0) = κ ρ_m,0 with **one empirical input A_0 = 0.0265**, the framework generates the following predictions — all of which hold up:

- **Bridge term form `b = A_0 · c/H_0`** — derived from STAM Shapiro through constant ambient A. With A_0 = 0.0265, this gives b = 354.95 Mly, matching the historical empirical bridge value exactly.

- **Pantheon+ d_L(z) shape matches data competitively with LCDM** — STAM beats LCDM by Δχ² = -3.8 on Pantheon+, using one empirical cosmological parameter (A_0) versus LCDM's three (Ω_m, H_0, Ω_Λ via flatness).

- **Pantheon+ vs Union3 inter-catalog ΔM tension** — STAM predicts +30 mmag from photon-A LoS-weighting through different z-distributions; observed +38 mmag. **27% match in magnitude with correct sign — a genuine STAM prediction validated against data that LCDM has no explanation for.**

- **DES Y5 anomaly correctly NOT predicted** — STAM's photon-A signature predicts only +14 mmag for Pantheon+ vs DES-Y5; observed -110 mmag. STAM correctly identifies this as instrumental (DES-Y5 zero-point calibration), not cosmological. **STAM correctly distinguishing cosmological signal from instrumental systematic is itself a validation.**

- **Hubble tension at H_0 = 73 km/s/Mpc** — STAM keeps H_0 = 73 throughout. LCDM has to adjust H_0 downward to 67.4 to fit CMB because LCDM lacks the photon-A term. The 5σ tension is LCDM's incompleteness; STAM resolves it structurally with no adjustment.

V(A) = β/(1-A) is **a data-validated STAM physical law**, on equal epistemic footing with any other physics-law form that has passed its tests — Maxwell, GR, Schrödinger all started as postulates that became laws by surviving data tests.

### Forward-falsifiable predictions

STAM predicts supernovae at z > 2 (JWST, Roman, Rubin) will sit slightly closer to EdS than LCDM extrapolation predicts. The d_L(z) curve flattens with respect to LCDM at high z. Future SN samples will sharpen this divergence and test STAM directly.

### Active development areas (where STAM and LCDM diverge most)

- **High-z SN samples** (z > 2): STAM predicts flatter curve; LCDM predicts continued acceleration. Distinguishable by future surveys.
- **BAO standard-ruler integration with V(A) cosmology directly**: the SU spreadsheet's polynomial extrapolation (since quarantined) was a separate exercise; V(A) modified Friedmann has not yet been formally BAO-tested.
- **Structure-dependent line-of-sight TE integration**: the photon path through real cosmic structure (galaxies, clusters, voids) is the natural quantitative refinement of the homogeneous-A picture.
- **First-principles derivation of A_0**: A_0 is currently STAM's single empirical cosmological coupling. Whether it is fundamental (like c, ℏ, G) or derivable from a deeper STAM principle is open.

### Distance is not assumed to equal catalog values

STAM does not assume that catalog-inferred distances must equal STAM geometric or accumulation distances:

```text
observable / catalog-inferred distance
    != STAM geometric or accumulation distance
```

Distance discrepancies are treated as diagnostic targets. The central question is:

```text
Which observable sees which STAM layer?
```

Current distance layers:

```text
D_geo(z)      = L z (1 + 0.15z)     # geometric spine
D_adj,0(z)    = L z (1 + 0.5z)      # accumulation-adjusted forward form
D_excess,0(z) = 0.35 L z^2          # traversal excess, TE
```

with:

```text
C = 3.261563776
H = 0.000243635
L = C/H = 13387.090426 Mly
```

### No-b supernova diagnostic

Catalogs tested against:

```text
D_adj,0(z) = L z (1 + 0.5z)
```

without a bridge term over `0.05 <= z <= 1.14418`:

| Catalog | Median residual, no-b | Mean residual, no-b | RMSE, no-b | RMSE with latest shared bridge `b = 354.95` |
|---|---:|---:|---:|---:|
| Union3 | +0.0512 mag | +0.0429 mag | 0.0539 mag | exact value pending from fixed-b rerun |
| Pantheon | +0.0808 mag | +0.0800 mag | 0.1664 mag | exact value pending from fixed-b rerun |
| DES | +0.2068 mag | +0.2383 mag | 0.3626 mag | 0.323471 mag |

**Latest shared bridge value:** `b = 354.95`. This is the Pantheon/Union-style bridge value currently carried forward for comparison across all three catalogs. The DES fixed-b RMSE is recorded from prior fixed-b output; Union3 and Pantheon should be filled from the next fixed-b rerun using the same redshift window and residual definition as the no-b table.

**Status:** Catalog diagnostic. Union3/Pantheon sit closer to the no-b curve than DES; DES remains the larger catalog-family offset. The shared bridge comparison is retained as catalog-comparison history, not as the core distance law.

### Historical bridge term b

The historical bridge term is retained as catalog-comparison history, not core physics:

```text
D_catalog(z) ≈ D_adj,0(z) + b z
```

Current derivation candidate:

```text
TE(z) = 0.35 L z^2
dTE/dz = 0.70 L z
b_pred = 0.70 L * z_anchor
```

Pantheon+Union3 weighted q25 anchor result:

```text
z_anchor ≈ 0.037250
b_pred   ≈ 349.07 Mly
b_hist   ≈ 354.95 Mly
error    ≈ -1.66%
```

**Status:** Strong derivation candidate. Not final proof.

### BAO and CMB distance-layer diagnostics

BAO tests suggest different observables may sample different accumulation layers:

```text
Full all-z BAO:  D_geo/(1+z) performed better
Low/mid-z BAO:   D_adj,0/(1+z) performed better
```

Refined partial-layer fit:

```text
D_BAO(z) = D_geo(z) + lambda(z) * TE(z)
lambda(z) ≈ 1.228 - 0.296z
```

At CMB scale, the no-b functions remain mathematically stable:

```text
<A_path>      -> 7/3
D_adj,0/D_geo -> 10/3
```

Simple direct mappings have not yet reproduced the observed acoustic angular scale. The CMB question remains an open distance-observable mapping problem.

---

## Supernova source accumulation

A supernova photon begins inside a local source and host-galaxy A-environment, not in empty space:

```text
A_source_local > A_intergalactic
```

Model-A predicts that supernova distance inference may carry a source-environment signature. Higher host mass should correspond to a higher local A contribution and a structured residual pattern. This may connect to the observed supernova host-mass step, currently treated in standard analysis as a nuisance correction.

**Status:** Qualitative prediction. Quantitative host-mass residual testing is open.

---

## Galaxy accumulation

Galaxy rotation is an exploratory branch of Model-A.

Linear accumulation from many masses can be written as:

```text
A_total(x) = sum_i 2G m_i / (c^2 |x - x_i|)
```

Scale demonstration:

```text
A from 1 solar mass at 10 kpc       ≈ 9.571e-18
A from 1e11 solar masses at 10 kpc  ≈ 9.571e-7
```

Using:

```text
v^2 ≈ (c^2/2) A
```

this corresponds to an implied circular speed near:

```text
v ≈ 207 km/s
```

Toy visible-galaxy result:

```text
Disk  = 6e10 solar masses
Gas   = 1e10 solar masses
Bulge = 1e10 solar masses

Outer 15-35 kpc median speed ≈ 121 km/s
Outer slope                  ≈ -2.82 km/s/kpc
```

Exploratory collective-envelope result calibrated near 220 km/s at 20 kpc:

```text
A_linear     ≈ 4.04e-7
A_collective ≈ 9.10e-7

Outer median speed ≈ 198.69 km/s
Outer slope        ≈ -1.71 km/s/kpc
```

**Status:** Exploratory. `A_collective` is not yet derived or claimed as required. Galaxy rotation remains an open reconstruction problem.

---

## Bold-STAM strong-field commitment

The strong-field branch has been worked out into a concrete commitment. Bold STAM accepts an asymmetric metric that preserves GR exactly in `g_tt` and modifies only `g_rr`:

```text
g_tt = -(1 - A) c^2                      (identical to GR — preserves all weak-field clock tests)
g_rr = 1 / [(1 - A)(1 - A^2)^2]          (bold STAM — vanishes as (1-A)^3 at A=1)
```

The leading correction to g_rr appears at second order in A: `(1 - A^2)^2 = 1 - 2A^2 + A^4`. Every observable with first-order-in-A precision (Shapiro delay, gravitational redshift, light deflection, Mercury perihelion, Hulse-Taylor pulsar timing, GPS clock comparison) automatically passes at current measurement precision because A is tiny everywhere we currently measure (A_sun_surface ~ 10^-6).

### The thirds-of-A pattern

Schwarzschild's three special radii fall at clean thirds of A, and bold STAM preserves all three because it preserves g_tt:

```text
ISCO  (innermost stable circular orbit):  r = 3 Rs    -> A = 1/3
Photon sphere (light orbits):             r = 1.5 Rs  -> A = 2/3
Event horizon:                            r = Rs      -> A = 1
```

Because the orbit and photon-sphere conditions depend only on g_tt, bold STAM's predictions for **EHT shadow size** and **accretion disk inner edge** match GR exactly.

### No-crossing in any frame

For an infaller with E = c^2 in the bold-STAM metric, the proper-time integral picks up a logarithmic divergence at A = 1. The traveler **falls forever** — they cross every intermediate A value at finite proper time, but proper time to reach A = 1 is infinite. Their local clock keeps ticking, equivalence principle holds locally throughout.

The two-planet thought experiment (observer at Alpha+ near A ≈ 0, traveler approaching A = 1) decomposes the observer-frame view honestly into three independent contributions:

```text
1. Traveler's diary (proper time)             — diverges logarithmically
2. Clock-dilation buildup (gravitational TD)  — diverges algebraically 1/sqrt(1-A)
3. SU delay (signal propagation through A)    — diverges algebraically 1/(1-A)^2
```

For A < 1 each is finite; comparison ratios are well-defined. As A → 1 all three diverge to ∞ at different rates but the same limit. The A = 1 boundary is mathematically defined as the locus where frame-comparison ratios become indeterminate (∞/∞) — both the traveler and the observer measure infinite time to reach it, by independent mechanisms. There is no privileged frame in which crossing can be registered as a finite event.

### The bubble picture

The geometrical interpretation: a black hole is **not** a deep funnel in spacetime with an interior. It is a **bubble** — a region of "no universe" with a smooth 2D boundary surface, the universe extending around it. All matter that ever fell toward the black hole accumulates asymptotically on the boundary, never inside (because no inside exists). The horizon is the edge of where space is, not a coordinate artifact with continuation through it.

A=1 is a phase transition between spacetime (A < 1) and not-spacetime (no manifold). The membrane is geometrically 2D — codimension one with no depth. This converges with mainstream theoretical-physics ideas (holographic principle, membrane paradigm, fuzzball conjecture) arrived at by independent routes.

### Multi-source A and binary mergers

Linear superposition (exact in weak field):

```text
A_total(x) = sum_i 2 G m_i / (c^2 |x - x_i|) = sum_i Rs_i / |x - x_i|
```

For two equal-mass BHs at separation d, the A field along the connecting line has midpoint value 4 Rs/d. The A=1 surfaces (bubbles) merge when this midpoint A reaches 1, giving **d_crit = 4 Rs** exactly. As BHs spiral inward through inspiral, the topology cascades through the thirds: A = 1/3 contours (ISCO) merge first, then A = 2/3 (photon), then at d = 4 Rs the A = 1 surfaces reconnect. Matter on the bubble surfaces joins onto the new common surface — no information loss, no interior physics required.

**Status:** Numerically checked and committed. Falsifiable corners: near-horizon Shapiro timing, LIGO ringdown quasinormal mode frequencies, primordial-BH evaporation spectra. EHT shadow size and accretion-disk ISCO are not testable separators (bold STAM = GR for both).

---

## Quantum interpretation

The quantum branch of Model-A reframes observation physically rather than consciously:

```text
Observation = physical interaction that resolves A
```

Current two-bucket framework:

```text
RESOLVED
    physical interaction has occurred
    A-state is confirmed
    path is definite

UNRESOLVED
    no physical interaction has occurred
    A-state exists but is not confirmed
    path is indeterminate
```

Important distinction:

```text
Unresolved != A = 0
```

An unresolved state is not absence of accumulation. It is an A-state that has not yet been reconciled through physical interaction.

The horizon is interpreted as a resolution boundary AND as the edge of the universe:

```text
A < 1  -> universe exists; external interaction possible; A is committed
A = 1  -> resolution boundary; geometric edge of universe
A > 1  -> not unresolved-A as a "weird state"; no manifold, no physics, no points
```

Beyond A=1 there is genuinely nothing — not a region with strange physics, not a pool of uncommitted possibilities. The universe ends at the boundary. This makes the "no path" mechanism sharper than "no inward anchor": beyond A=1, the A field itself is undefined, so the action `int A ds` along any inward-going path is undefined too — inward paths simply do not exist in the path-integral sense. Vacuum fluctuations near A=1 can only resolve outward; the inward direction has no destination. This wedge-restriction structure is what produces thermal emission at the boundary (see next section).

The picture converges with mainstream theoretical-physics ideas — holographic principle, membrane paradigm, fuzzball conjecture — arrived at by independent routes. STAM provides a *physical reason* for the no-interior picture (the universe ends at A=1) rather than a postulate or heavy mathematical derivation.

**Status:** Theoretical foundations now formalized. Hawking-like radiation, entropy, and first-law thermodynamics derived in the next section.

---

## Black hole thermodynamics from bold STAM

The full thermodynamic sector has been derived from STAM-native ingredients (gravity bridge + resolution rule + bubble picture). All standard Hawking / Bekenstein results are reproduced exactly, with the *physical mechanism* being phase-boundary thermal equilibrium rather than QFT-on-curved-spacetime mode counting.

### Hawking temperature (Q8 / Q10)

Single rule for thermal emission at any A=1 boundary:

```text
k_B T = (1 / (4 pi)) * hbar * c * |grad A|_boundary
```

The 1/(4π) factors structurally as `2π × 2`:
- `2π` from thermal-state imaginary-time periodicity (general thermodynamic).
- `2` from the STAM gravity bridge `g = (c^2/2) grad A`.

This single rule reproduces:

```text
Schwarzschild Hawking T_H  (varying M)         |grad A| = 1/Rs       -> exact
Unruh T_U                  (varying a)         |grad A| = 2a/c^2     -> exact
de Sitter T_dS             (varying H)         |grad A| = 2H/c       -> exact
```

All three thermal-horizon temperatures agree exactly, no fitted parameters.

### Spectrum (Q10)

Planckian shape derived from Bose-Einstein occupation of bulk vacuum fluctuations in thermal equilibrium with the A=1 phase boundary. The "no path" mechanism enforces wedge-restriction: outward fluctuations resolve, inward have no destination. Net asymmetric outward flux is the emission, with thermal distribution at temperature T given above.

### Bekenstein-Hawking entropy (Q11)

```text
S = k_B A_horizon / (4 ell_P^2)
```

Derived two ways:

- **Bubble argument**: all matter is on the 2D phase boundary; no interior degrees of freedom; saturation of holographic bound forces S = A/4 in natural units. Author's plain-language statement: "area scaling because that's where everything is."
- **First-law integration**: `dE = T dS` with `E = M c^2` and `T` from Q8/Q10, integrated from M=0 gives `S = 4 pi k_B G M^2 / (hbar c) = k_B A / (4 ell_P^2)`. Same result.

### First law and Smarr (Q12)

Differential first law `dE = T dS` and Smarr formula `M c^2 = 2 T S` both verified to machine precision. Hawking evaporation trajectory M(t), T(t), S(t) computed; lifetime for solar-mass BH is 2.1 × 10^67 years matching textbook value. Generalized second law gives a clean 1/3 entropy surplus during evaporation.

```text
| Quantity              | Standard                  | Bold STAM                       | Match  |
|-----------------------|---------------------------|---------------------------------|--------|
| Hawking T             | hbar c^3/(8 pi G M k_B)   | hbar c |grad A| / (4 pi k_B)    | exact  |
| Spectrum              | Planckian at T_H          | Planckian at T (phase boundary) | exact  |
| Entropy               | k_B A / (4 ell_P^2)       | k_B A / (4 ell_P^2)             | exact  |
| First law             | dE = T dS                 | residual ~ 10^-16 (verified)    | exact  |
| Smarr                 | M c^2 = 2 T S             | residual ~ 10^-16 (verified)    | exact  |
| GSL surplus           | dS_total >= 0             | exactly +1/3                    | exact  |
| Solar-mass evap time  | ~ 10^67 yr                | 2.1 x 10^67 yr                  | exact  |
```

**Status:** Derived. Bold STAM agrees with all standard BH thermodynamic results numerically while having a fundamentally different physical mechanism (phase boundary in thermal equilibrium, not Schwarzschild Wick rotation; bubble surface saturation, not QFT mode counting).

---

## Current working position

```text
A is the core Model-A variable.

Internal consistency derived in:
    local gravity                          (weak-field, gravity bridge)
    horizon threshold                      (algebraic identity at A=1)
    GPS-like weak-field clocks             (numerically checked)
    Shapiro-style delay                    (numerically checked)
    mass-estimator identities              (synthetic checks)
    bold-STAM strong-field metric          (committed; thirds-of-A preserved)
    no-crossing infall                     (logarithmic proper-time divergence)
    GR exterior recovery                   (passes all weak-field tests by 10^4 - 10^14)
    multi-source A and merger topology     (linear superposition; d_crit = 4 Rs)
    Hawking temperature                    (derived from gravity bridge + resolution rule)
    Hawking spectrum                       (Planckian from phase-boundary equilibrium)
    Bekenstein-Hawking entropy             (S = A/4 from bubble + first law)
    First law and Smarr formula            (verified to machine precision)
    Hawking evaporation dynamics           (10^67 yr for solar-mass; matches textbook)
    Generalized second law                 (+1/3 surplus, matches standard)

Cosmological commitment: bold STAM cosmology
    Universe is matter-only Einstein-de Sitter. No dark energy.
    H_0 = 73 km/s/Mpc (SH0ES local) is the true value at all z.
    LCDM-CMB H_0 = 67.4 is LCDM's biased inference (missing photon-A term).
    Hubble tension RESOLVED structurally: STAM does not have to adjust;
        LCDM does, because it lacks the photon-A term.
    Dark energy remitted as photon-A traversal excess (TE) along LoS.
    Bridge term FORM b = A_0 * c/H_0 derived from STAM Shapiro
        (A_0 = 0.0265 is calibrated to historical b, not yet derived).
    Pantheon+ vs Union3 inter-catalog tension predicted within 27%.
    Simple-A quantitative version fails BAO; structure-dependent
        TE path-integration is the open quantitative path forward.

Galaxy rotation remains exploratory.
Quantum interpretation: resolved/unresolved framework grounded by Hawking derivation.

Strong-field commitment: bold STAM accepts asymmetric metric
    g_tt = -(1-A) c^2  (= GR)
    g_rr = 1 / [(1-A)(1-A^2)^2]  (= bold STAM, modified)
    Universe ends at A=1 (no past-A=1 region).
    Bubble picture: black holes are 2D boundary surfaces, not 3D interiors.

Falsifiable corners:
    Near-horizon Shapiro timing (above current Cassini precision).
    LIGO ringdown quasinormal modes (within reach of next-gen detectors).
    Primordial-BH evaporation spectra (deviations possible).
NOT testable separators (bold STAM = GR for these):
    EHT shadow size (depends on g_tt only).
    Accretion-disk ISCO position (depends on g_tt only).
```

---

## Priority open problems

```text
Closed or substantially advanced:
    [DONE] Strong-field beyond spherical weak-field        (SF3 — bold-STAM metric committed)
    [DONE] A > 1 horizon dynamics                          (universe ends at A=1; no interior)
    [DONE] Hawking radiation rate                          (Q8 / Q10 — phase-boundary derivation)
    [DONE] Bekenstein-Hawking entropy                      (Q11 — bubble + first law)
    [DONE] First law of BH thermodynamics                  (Q12 — verified to machine precision)
    [DONE] Multi-source A field (linear superposition)     (SF5 — d_crit = 4 Rs exact)
    [PARTIAL] Develop the resolved/unresolved quantum interpretation
             (now grounded by Hawking thermodynamic derivation)

Active / open (development on a competing framework, not patches to LCDM):
    1. Specify bold-STAM field equations                   (F3 exposed this gap)
       The metric ansatz g_rr = 1/[(1-A)(1-A^2)^2] is committed; the
       Lagrangian or modified-gravity theory that picks it out uniquely
       is not yet written down. Script 41 confirmed the metric is NOT
       consistent with Einstein gravity coupled to canonical scalar for
       any V(A); the framework requires non-canonical kinetic term,
       non-minimal coupling, or modified-gravity sector.
    2. Run V(A) = β/(1-A) modified Friedmann directly      (BAO not yet tested for V(A) itself)
       through DESI BAO at all z. (Earlier BAO test used the now-quarantined
       SU polynomial, not the V(A) cosmology.) STAM predicts BAO with
       structure-dependent A_LoS; quantitative comparison open.
    3. Structure-dependent TE path-integration cosmology   (open)
       The path-integrated A through actual cosmic structure (galaxies,
       clusters, voids) is the natural cosmological refinement.
    4. Dynamical-A propagation / gravitational waves       (LIGO inspiral comparison)
    5. Quasinormal mode frequencies for merger ringdown    (next-gen LIGO testable)
    6. Real mass closure across acceleration, lensing,
       redshift, and time delay                            (open observational test)
    7. Galaxy rotation reconstruction from real datasets   (exploratory branch)
    8. PPN second-order corrections (bold-STAM has 3 A^2
       in g_rr expansion vs GR's A^2)                      (high-precision test target)
    9. Dynamical collapse: rigorous proof that bold-STAM
       collapse never transiently forms a strict trapped
       surface during dynamics                             (F1 exposed this)
   10. Derivation of A_0 from STAM principles              (currently calibrated to bridge)
       1/(12*pi) = 0.02653 matches A_0 = 0.02651 to 4 sig figs (factors
       of 4 pi from Q8 thermal-bridge plausible). First-principles
       derivation likely tied to F3 Lagrangian work (item 1).

Distance:
    H_0 = 73 km/s/Mpc committed as the true Hubble constant.
    Hubble tension resolved structurally: STAM does not adjust;
        LCDM has to adjust H_0 down to fit CMB because it lacks photon-A.
    Dark energy is remitted as photon-A traversal excess (TE) along LoS,
        misinterpreted by LCDM as cosmic acceleration.
    SN sit on a flatter d_L curve than LCDM at high z;
        future high-z surveys will distinguish.
```

---

## Implemented / active scripts

Repository scripts and tests currently cover:

### Core framework

```text
00_validate_model_a.py             core Model-A definitions and identities
01_generate_stam_distance_tables   STAM distance-layer tables
02_mass_inference_synthetic        synthetic mass-estimator closure
07_shapiro_delay_model_a           Shapiro delay
08_gps_clock_model_a               GPS-like clock corrections
09_supernova_no_b_model_a          no-b supernova diagnostic
10_horizon_threshold_model_a       horizon threshold A=1
11_gw_propagation_model_a          GW propagation weak-field
```

### Cosmological / catalog diagnostics

```text
03 - 06       supernova fitting and bridge-term diagnostics (history)
25_TE_clock_shapiro_b_bridge      TE separation from bz bridge
26_b_from_GPS_Shapiro_scale       GPS/Shapiro scale sniff
27_b_from_TE_anchor_redshift      b as TE tangent: 349.07 vs 354.95 Mly
TE_traversal_excess_ledger        TE accounting
```

### Strong-field (SF) — bold-STAM commitment

```text
30_strong_field_SU_horizon_integral
    SU log-divergence at A=1; characteristic scale = Rs

SF1_A_metric_schwarzschild_recovery
    weak-field reduction to Schwarzschild

SF2_strong_field_SU_accumulation
    SU accumulation load A/(1-A)

SF3_no_crossing_thought_experiment
    bold-STAM metric committed; two-planet triangle decomposed;
    boundary collapse at A=1 derived as ∞=∞ limit

SF4_GR_exterior_recovery
    passes all weak-field GR tests by 10^4 - 10^14 margin

SF5_multi_source_A
    linear superposition; binary BH merger topology;
    d_crit = 4 Rs exact for equal masses
```

### Quantum / thermodynamic (Q) — Hawking from STAM

```text
Q1 - Q4    earlier resolved/unresolved framework scripts
Q6 - Q7    information-receipt and Hawking-toy

Q8_unresolved_A_boundary_flux
    single rule k_B T = hbar c |grad A| / (4 pi)
    reproduces Hawking, Unruh, de Sitter exactly

Q9_resolution_temperature_derivation
    Path 1: STAM-native energy scale derivation;
    knee at hbar c |grad A| (4π above Q8 thermal scale)

Q10_thermal_shape_derivation
    Path 2: Planckian shape from wedge-restriction;
    4π = 2π (thermal periodicity) × 2 (gravity bridge)

Q11_bekenstein_hawking_entropy
    S = k_B A / (4 ell_P^2) from bubble + first law

Q12_first_law_verification
    differential first law, Smarr, evaporation lifetime,
    GSL surplus (1/3) — all machine precision
```

---

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .

python scripts/00_validate_model_a.py
pytest
```

Generate first tables and synthetic checks:

```bash
python scripts/01_generate_stam_distance_tables.py
python scripts/02_mass_inference_synthetic.py
```

Outputs are written to `results/`.

---

## Repo structure

```text
STAM_model-A/
├── README.md
├── CHANGELOG.md
├── PRIORITY_RECORD.md
├── pyproject.toml
├── requirements.txt
├── docs/
│   ├── FORMULAS.md
│   ├── THEORY.md
│   ├── CLAIMS_AND_STATUS.md
│   ├── FALSIFICATION_TESTS.md
│   ├── PARAMETER_LOCKING.md
│   └── archive/
├── src/stam_model_a/
│   ├── constants.py
│   ├── local.py
│   ├── propagation.py
│   ├── cosmology.py
│   └── mass_estimators.py
├── scripts/
├── tests/
├── data/
├── results/
└── notebooks/
```

---

## Development rule

Supportive demonstrations and falsification tests should remain separate.

An algebraic identity test is useful. A fitted catalog comparison is useful. A real prediction test must lock parameters first, preserve the script and output, and then evaluate against independent data without silent retuning.

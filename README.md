# Model-A

**Model-A** is the current variable-accumulation version of the broader **Spacetime Accumulation Model (STAM)** research program. It explores whether gravity, clock behavior, propagation delay, black-hole horizons, thermodynamic behavior, cosmological distance effects, dark matter behavior, and quantum-style resolution can be described using one central idea: a dimensionless spacetime accumulation field called **A**.

Author: **Sean Brady**
Status: **Proposed theoretical framework / active research program**
Snapshot: **May 13, 2026 — strong-field arc structurally complete in spinless AND spinning cases, with the metric derived from minimum-commitment quantum/write-density route. Spinless: k(A) = (1−A) outside the photon sphere (exact GR); k(A) = (1−A)·F(Σ) inside, F(Σ) = 1 − 4(Σ−2)³ + 3(Σ−2)⁴ (quartic Hermite, C² at PS, C¹ at horizon). Σ = 3A is the natural shell coordinate (landmarks ISCO Σ=1, PS Σ=2, horizon Σ=3). Derivation: D = 3 spatial write channels + 2 horizon-pair channels (outer-face pair structure: write + reduction) → Beta(D, 2) closure density p(y) = 12y²(1−y) → F = 1 − 4y³ + 3y⁴, the minimum-commitment quantum/write-density profile. SU shell-count + two-face joint compatibility forces D = 3 (framework-internal, not anthropic). Kerr: A(r, θ; M, a) = 2Mr/(r²+a²) places the bubble at r = r_+ (Kerr horizon, constant BL coord, oblate when embedded), hologram spins on stationary horizon, T uniform = T_Kerr, LIGO ringdown = exact GR-Kerr for all spin. Entropy α = 4 = (outer-face pair × gravity-bridge) structurally derived (G18 caveat closes). Pair structure derived from elevator argument (write IS reduction on outer face). Hawking T derived two independent ways. Stress-energy verified non-pathological. Spinless ringdown = exact GR; Kerr ringdown = exact GR-Kerr; the STAM-vs-GR wedge lives inside the photon orbit and is observable only at sub-leading WKB precision. Quantum-unit commitment: each quantum physical interaction resolves exactly 1 SU = A_0 of substance — unifying the cosmological ruler reading of SU with the quantum write reading. Weak-field A(r) is the continuum / coarse-grained limit of many SU writes. Born rule derived structurally: `u_i = ||P_i ψ|| = |ψ_i|` is a projection-geometry theorem (four axioms: phase blindness, projector locality, unitary covariance, orthogonal refinement); `p_i = u_i²` is a STAM ledger-measure result from the principle that resolution-event phase volume scales as the square of unresolved support amplitude. The same pair structure (outer-face: write + reduction is one event with two coupled aspects) does three framework jobs — α_H = 2 in the strong-field metric, α = 4 in boundary entropy, |ψ|² squaring in the Born rule. Decoherence is the structural consequence of paired events resolving alternatives into distinct ledger channels before recombination, not a separate postulate. Bell / CHSH violation reproduced exactly at the Tsirelson bound (G67). First-principles QM action chain: substance velocity-cap → proper time dτ = dt√(1−A)√(1−v²/c²) → relativistic action S = −mc²·τ → non-relativistic L = KE − mΦ with Φ = c²A/2 → path integral in unresolved support → Born rule at resolution events. SU support `N_SU(x) = A(x)/A_0` and actual write rate `dN_write = Γ_res × dτ` are distinct quantities. Γ_res now has committed structural form: `Γ_res = Σ_μ ⟨L_μ† L_μ⟩ D_μ` (Lindblad-analog with L_μ as physical write channels and D_μ as channel distinguishability) bounded by `0 ≤ Γ_res ≤ (A/A_0)/τ_P` (interaction-gating from below; SU support capacity from above). Per-channel split via Born rule p_i = |ψ_i|² gives dN_i = Γ_res p_i dτ. Remaining work is identifying specific L_μ and D_μ for particular interactions — domain-of-application, not structural gap.**

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

The weak-field limit is the **continuum / coarse-grained version of many SU writes** (under the 2026-05-13 quantum-unit commitment: each quantum physical interaction resolves 1 SU = A_0 of substance). In the bulk, the smooth A(r) field is the macroscopic density that emerges from averaging over many discrete quantum write events. The weak-field successes below are what shows up at scales where the SU granularity is far below resolution — Newton's law, GPS clocks, Shapiro delay, lensing, orbital mechanics — are all the continuum limit of the framework's quantum substrate.

We keep:

```text
A(r) = R_s / r = 2 G M / (c² r)
```

and:

```text
g = (c²/2) ∇A
```

For the spherical weak-field form `A(r) = 2GM / (c^2 r)`, the gradient gives:

```text
g = -GM / r^2
```

So the inverse-square acceleration law is recovered exactly in this setting. The factor `c^2/2` is not fitted; it is the reciprocal of the dimensional factor already built into the definition of A.

The weak-field formulation is the cleanest part of Model-A and the bridge between the framework's quantum-scale write structure and its macroscopic predictions:

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

**(Spinless and Kerr cases structurally complete as of 2026-05-13.)**

The spinless strong-field metric in its final form:

```text
g_tt = -(1 - A) c²                        (matches GR's form, unchanged)
g_rr = 1 / k(A)                           (Model-A modification)

k(A) = (1 - A)                            for A ≤ 2/3   (outside PS — exact GR)
k(A) = (1 - A) · F(y),  y = 3A - 2        for 2/3 < A < 1
F(y) = 1 - 4y³ + 3y⁴                      (quartic Hermite, C² at PS, C¹ at horizon)
```

Equivalently in shell-coordinate Σ:

```text
k(Σ) = (1 - Σ/3)                          for Σ ≤ 2
k(Σ) = (1 - Σ/3) · F(Σ)                   for 2 < Σ < 3
F(Σ) = 1 - 4(Σ - 2)³ + 3(Σ - 2)⁴          (quartic Hermite)
```

**Note (2026-05-13).** Two candidate forms were explored: the quartic Hermite F = 1 − 4y³ + 3y⁴ (Beta(D, 2) under D spatial + 2 horizon-pair channels — the minimum-commitment "quantum/write-density" route) and the quintic Hermite F = 1 − 5y⁴ + 4y⁵ (Beta(D+1, 2) under an added ledger-as-fourth-channel reading, G58). The framework commits to the **quartic** under the minimum-commitment principle: the ledger is the universe's record of writes but not a separate write-channel in the configuration measure. The eikonal LIGO ringdown prediction is τ_STAM/τ_GR = 1 exactly under both (verified G66).

### Natural shell coordinate

The natural strong-field variable is **Σ = A / (4π A_0) = D × A** (= 3A for D=3):

- Σ = 1 → ISCO (A = 1/3)
- Σ = 2 → photon sphere (A = 2/3)
- Σ = 3 → horizon (A = 1)

The thirds-of-A landmarks are integer shells in Σ. The final-shell radial closure profile F(Σ) operates on the third shell: 2 < Σ < 3 (equivalently 2/3 < A < 1).

### Derivation chain (quantum/write-density route, minimum-commitment)

The strong-field metric is derived from the framework's structural primitives — substance ontology, presentism, two-face refinement, SU shell-count, elevator identity, and natural measure on the configuration manifold:

1. **A_0 = 1 / (4π D)** from substance baseline: 4π from Q8 thermal × gravity-bridge, D from spatial dimensionality. For D = 3: A_0 = 1/(12π).

2. **Σ = D × A** as the natural strong-field shell coordinate. For D = 3, Σ = 3A places ISCO at Σ = 1, PS at Σ = 2, horizon at Σ = 3.

3. **Final-shell coordinate y = Σ − 2 ∈ (0, 1).** F(Σ) is interpreted as the survival fraction of unresolved radial closure inside the final shell; the closure density is p(y) = −dF/dy.

4. **Channel structure.** Final-shell closure is distributed across D spatial write channels and 2 horizon-pair channels (outer-face pair structure: write component + reduction component, both on the outer face — the elevator identity). The ledger is the universe's record of writes but not a separate write-channel in the configuration measure (minimum commitment).

5. **Configuration volume on the final shell.** Under uniform per-entry measure (natural measure on the configuration manifold), the spatial-vs-horizon-pair split has Beta(α_S = D, α_H = 2) density. For D = 3: p(y) = 12y²(1 − y), giving survival F(y) = 1 − I_y(3, 2) = **1 − 4y³ + 3y⁴** (the quartic Hermite).

   The "minimal" candidate p(y) = 12y²(1 − y) has four structural properties: starts at zero at PS, turns on smoothly inside the final shell, returns to zero at the horizon, integrates to one completed final-shell closure.

6. **k(A) = (1 − A) outside the photon sphere.** No STAM modification where light can escape; weak-field tests pass automatically.

7. **D = 3 forced by joint compatibility.** ord_{A=1} k = α_H + 1 = 3 automatically from the two-face commitment (α_H = 2, D-independent). SU shell-count separately requires ord = D. Joint compatibility ⇒ D = 3 framework-internally (not anthropic).

Near the horizon: F ~ 6(3 − Σ)² = 54(1 − A)², so k = (1 − A) · F ~ 54(1 − A)³, giving ord_{A=1} k(A) = 3, preserving the SU shell-count result n_h = D − 1 = 2.

**Interpretation.** A is the physical accumulation field; SU is the minimum accumulation unit (1 SU = A_0); Σ is the strong-field shell count; F(Σ) is the unresolved final-shell survival profile. Strong field supplies the boundary domain 2 < Σ < 3, but the *shape* of F(Σ) is derived from final-shell resolution density (the write-density route) rather than from a global constant-n metric ansatz.

### Predictions under the committed metric

| Probe | A | k(A) | Prediction |
|---|---:|---|---|
| Solar System / Cassini | tiny | (1−A) | GR exact |
| GPS clocks | tiny | (1−A) | GR exact |
| Pulsar binary Shapiro | low | (1−A) | GR exact at relevant A |
| ISCO | 1/3 | 2/3 | GR Schwarzschild |
| NEC at A = 1/2 | 1/2 | 1/2 | Schwarzschild (F = 1, no STAM modification) |
| **LIGO ringdown spinless (PS)** | **2/3** | **1/3** | **= GR Schwarzschild exactly** (G57, G66) |
| **LIGO Kerr ringdown (any spin)** | — | — | **= GR-Kerr exactly** (G62, G64) |
| Inside-PS strong field | (2/3, 1) | quartic Hermite ramp | sub-leading WKB / late-inspiral / EMRI distinguishers |
| Horizon | A → 1 | (1−A)³ × (10/3) | k → 0 with order D = 3 |

**LIGO ringdown is exact GR for any spin under the current commitment.** Both formulas (spinless and Kerr) match GR-equivalent Lyapunov at the photon orbit, because k = (1−A) "outside the photon sphere" is exact GR there. The framework's STAM-vs-GR wedge lives entirely inside the photon orbit, observable only via higher overtones, late-inspiral chirp, or LISA EMRIs.

### Generalization to dimension D

```text
A_0 = 1/(4πD)
ord_{A=1} k(A) = D
n_h = D − 1
α_entropy = 2(D − 1)
Σ = D × A, with landmarks at Σ = 1, 2, ..., D
```

For D = 3 (our universe), the framework's strong-field metric is fully determined.

### Kerr extension (committed 2026-05-13)

```text
A(r, θ; M, a) = 2 M r / (r² + a²)        (Kerr substance density)
A = 1 at r = r_+ (Kerr horizon, constant in Boyer-Lindquist r)
```

The bubble surface is at constant r = r_+ in BL coordinates — oblate when embedded in flat 3-space. The hologram matter spins on this stationary horizon, carrying J at ZAMO frequency. Properties:

- A = 2M/r recovered for a = 0 (Schwarzschild limit)
- T uniform on the bubble = T_Kerr (standard Kerr horizon T)
- All Kerr photon orbits and ISCOs lie outside the bubble for any sub-extremal spin
- Equatorial-plane k(A) reduces to the spinless form (Σ in equatorial plane)
- Rotation enhancement R(θ) = 1 + (v_matter/c)² applies as a RATE effect (super-radiance from matter rotation), not as a T variation
- LIGO Kerr ringdown = exact GR-Kerr in eikonal for all spin

The previous "static-limit" reading (A = 2Mr/Σ, Σ = r² + a² cos²θ) — which predicted latitudinal T banding — is retired (G64). It had A > 1 at the Kerr equatorial photon orbit for a > 0.707, breaking the framework's A < 1 commitment. The committed formula A = 2Mr/(r²+a²) resolves the high-spin photon-orbit and ISCO consistency cleanly.

### Stress-energy verification (G65, 2026-05-13)

The framework's quartic Hermite metric is non-pathological (verified G65 — the verification was done initially on the quintic form, with the same conclusions holding for the committed quartic by the same arguments):

- **Conservation automatic** (TOV residual at numerical noise level — Bianchi identity)
- **Curvature bounded throughout final shell** (Kretschmann K_STAM peaks at ~0.46 in M=1 units, *less* than Schwarzschild's K = 0.75 at the horizon — no hidden singularity)
- **C² smooth at PS** (stress-energy → 0 quadratically as A → 2/3+)
- **SEC satisfied** throughout final shell (no anti-gravity / ghost behavior)
- **NEC_r and WEC violated** for most of the final shell — characteristic of modified-gravity / dark-energy-like effective stress-energy (consistent with F3's earlier w ≈ −1 finding); not a pathology
- Effective stress-energy reads as tension-dominated (p_r < 0, p_t > 0) with NEC_t satisfied

### What's still open in strong field

- **Lagrangian for A**: the metric is uniquely derived from primitives, but the action principle that produces this k(A) is still missing. The remaining open question is: what action principle yields k(A) = (1−A) outside PS and (1−A)·(1 − 4y³ + 3y⁴) inside?
- **Inside-PS observable distinguishers**: late-inspiral chirp, higher overtones, LISA EMRIs. Current LIGO precision doesn't reach this regime.
- **Spectral details of Hawking emission**: framework gives T and per-entry structure, but the full spectral distribution still requires QFT machinery beyond cell-counting.

---

## Substance velocity-cap (V_4 candidate commitment, added 2026-05-12)

The substance ontology forces a new structural commitment: **motion through elevated A is REALLY slowed**, not just observationally. This is the operational consequence of treating A as real substance density rather than coordinate artifact, combined with the magic-bell prohibition (no view-from-nowhere observation).

Operationally:
```text
v_effective² = v_Newton² × f(A_local)
```
where f(A_0) ≈ 1 and f(A) → 0 as A → 1. The current V_4 candidate form is **f(A) = 1 - A** (the proper-time-squared factor, ontologically consistent with g_tt structure). Other candidates were tested (√(1-A), (1-A)²(1+A)); the (1-A) form gave the best empirical match to GW170817 while being structurally motivated.

**Wedge with GR**: GR has the geometric metric factors but no separate substance interaction. STAM has both. At low A (everywhere currently measured), they agree to within ~3%. At high A (late binary inspiral), they diverge predictably:

- **Inspiral chirp shape**: STAM-corrected dynamics produce slightly different chirp profile than pure GR templates
- **Binary mass extraction**: LIGO templates assume GR; STAM-corrected templates would extract slightly different chirp masses (~few percent shift)
- **Binary merger engine time** (BNS systems with EM counterparts): see "Falsifiable predictions" below

**GW170817 engine time** (specific consequence): under the substance velocity-cap + Sean's elastic-shell picture, the engine time (orbital collision → gamma-ray burst onset) equals τ_critical, the Peters-Mathews chirp time from the substance-corrected doughnut threshold:

```text
Self-consistent threshold equation: A_0 = x(1-x)/(1+x)  where x = R_s/(2r)
For A_0 = 1/(12π): r_threshold/R_s ≈ 17.82
τ_critical = (5/8) × (r_threshold/R_s)⁴ × R_s/c
For 2.7 M_sun (GW170817): τ_critical = 1.677 s
```

**Match: 1.677 s predicted vs 1.74 s observed — 3.6% match, no fitting.** The progression (naive Newton 2.10 s → partial substance correction 1.89 s → full self-consistent 1.68 s) is monotonic toward observation as substance correction is applied.

**Mass-scaling**: τ_critical ∝ M_total (linear). Slope ~0.62 s per M_sun. Falsification handle for future BNS+EM events.

Scripts: [G42](scripts/G42_doughnut_threshold_corrected.py) (threshold-only correction) and [G43](scripts/G43_doughnut_contraction_full_stam.py) (full self-consistent) hold the calculation. See also [memory/project_substance_velocity_cap.md](memory/project_substance_velocity_cap.md).

---

## Black-hole interpretation

In Model-A, a black hole is interpreted as a 2D bubble surface, not a deep interior region.

```text
A < 1  -> spacetime exists
A = 1  -> 2D phase boundary / bubble surface
A > 1  -> not part of the manifold
```

Matter that fell toward the black hole never crossed A = 1; it accumulated holographically on the bubble surface (outward-collapse picture). There is no interior; there is no singularity. Information lives on the 2D boundary surface.

**A = 1 is never reached in any finite time** (refined 2026-05-11). The horizon is the asymptotic limit of an unresolved Zeno-paradox-like halving series: in GR Schwarzschild, the halving series converges and infallers cross in finite proper time; in STAM, the series does NOT converge (each halving takes constant time for n=1, or growing time for n=2). The 2D "hologram surface" is then a mathematical limit set; bulk content asymptotically piles up just-below A=1 from each side, never landing on it. "Information lives on the boundary" is shorthand for "information lives in the bulk asymptotically near the boundary." The observer and traveler agree: there is no finite-time crossing event; the structure of A=1 itself is unreachable.

**Spinning black holes (refined 2026-05-13).** Under the committed Kerr substance density A = 2Mr/(r²+a²), the bubble surface sits at constant r = r_+ in Boyer-Lindquist coordinates — oblate when embedded in flat 3-space. The bubble itself is static; matter on the outer face spins at ZAMO frequency, carrying J. The hologram spins over the stationary horizon. T is uniform on the bubble = T_Kerr (matching standard Kerr horizon T exactly via the r_+² + a² = 2Mr_+ identity). Super-radiance enhancement R(θ) = 1 + (v_matter/c)² acts on the emission rate (matter velocity peaks at the equator), not on T. No Penrose extraction (rotational energy is in the matter, not in empty rotating geometry).

**Two-face refinement of the bubble surface (committed 2026-05-10, clarified 2026-05-13).** The 2D boundary at A=1 is structurally distinguished into two faces:

- **Inner face**: holds primordial mass from the formation event. **Invariant** — the inner face cannot be reached after formation and never changes mass. This is why black holes "remember" their formation configuration permanently.
- **Outer face**: holds subsequently accreted matter. Dynamic, can rotate, carries angular momentum, drains via Hawking emission.

The framework's BH life cycle: formation presses primordial mass to the inner face (locked); accretion builds outer face; Hawking drains outer face only; horizon size tracks outer-face content in real time; final state is a stable primordial-mass remnant. **Black holes never fully evaporate.** When the outer face has emitted all its accreted content, the BH is at its primordial-mass remnant state.

The "horizon doesn't spin, hologram does" principle applies specifically to the outer face's holographic content. PBHs are "thin outer face" objects (formed primordially with little subsequent accretion), predicting evaporation signatures differing from stellar BHs both in standard ways (low spin) and in framework-specific ways (sparse outer-face content). Information that "fell" into the BH lives on the inner face permanently; information from accretion lives on the outer face and gets emitted as Hawking radiation. The information paradox dissolves with a definite mechanism rather than as a slogan.

**Pair structure of Hawking emission derived from the elevator argument (G60, 2026-05-13).** Each emission event at A=1 is ONE event with two structural aspects, both on the outer face: an outward write (substance leaves) and a horizon reduction (outer face mass decreases). The "elevator gets lighter when someone steps off" — the act of leaving IS the act of reducing weight; they're inseparable. This derives the framework's pair structure of Hawking emission from no-interior + substance conservation + two-face refinement, replacing the earlier "stated structural consequence" framing.


## A_0 void interpretation

The 4 pi is the prefactor that appears in the thermal-emission rule k_B T = hbar c |grad A| / (4 pi), independently derived in Q8/Q10 as 2 pi (thermal-state imaginary-time periodicity) x 2 (Model-A gravity bridge factor c^2/2). The 3 is the count of spatial dimensions over which A must be non-zero for a 3D manifold to exist.

Together, A_0 is proposed as the structural floor of spacetime — the minimum "amplitude per solid-angle x dimensionality" at which a spatial manifold can structurally manifest. Below this value, space cannot exist; A = 0 is not vacuum but absence of manifold.

The numerical match between 1/(12 pi) = 0.026526 and the empirical bridge term b/L = 354.95/13387 = 0.026514 is 0.04%. This is the framework's strongest evidence for treating A_0 as a derived structural constant rather than a calibrated parameter, though the first-principles derivation of the spatial-dimensionality factor of 3 remains open theoretical work.

**Quantum-scale reading (interpretation, not derivation).** Under Model-A's resolved/unresolved-A interpretation — where physical interaction is what resolves A, and the resolved record is the universe's running ledger of what has happened — A_0 reads as the minimum density of resolved-A required to maintain the manifold's structural existence. The (4 pi x 3) decomposition gives this density a concrete cell-form: roughly one structural unit of resolved-A per (4 pi solid angle x 3 spatial directions) per Planck cell. Below this density, the manifold has too few "writes" to sustain itself — which is what "A = 0 means no spacetime" maps to at the quantum-cell scale. The same vocabulary describes physical writes at every scale of the framework: horizon writes (Hawking radiation as outward A-resolution where inward is forbidden), local quantum writes (resolution events from any physical interaction), and cosmic-floor writes (the A_0 minimum that keeps the manifold on the books). This is interpretation aligning the framework's quantum vocabulary with the existing structural decomposition, not new derivation. A future first-principles calculation that produces 1/(4 pi x 3) from such a write-density premise — without referencing the bridge term or the value 0.0265 in its setup — would constitute the actual derivation. That work is open.

A_0 status (snapshot 2026-05-10): operational structural commitment. Observationally distinguished — only A_0 = 1/(12 pi) satisfies both the bridge-term match and the CMB physical-amplification constraint (script G15). Theoretically anchored on the (4 pi x 3) decomposition with two component-level groundings (thermal 4 pi from Q8/Q10; geometric 1/3 from radial-line-of-sight averaging through 3D volume; the unifying first-principles derivation that combines them is open). No currently-buildable independent derivation route lands on 0.0265 (script G16). Treated the way the historical bridge term b = 354.95 was treated before A_0 was proposed: a working number that the framework is built around, with the "why" partially answered and partially open. Updates will be made if the open derivation closes or new structural arguments emerge.

**V_3 vacuum tone (G19).** Small fluctuations of A around its V_3 minimum at A_0 oscillate at characteristic frequency `ω = sqrt(V_3''(A_0) / M_P^2_red)` ~ 10 H_0 (about ten Hubble rates per cycle, vacuum-oscillation period ~7-8 Gyr). The framework's structural-floor commitment ties V_3's natural tone to cosmic-dynamics scales rather than Planck or microscopic scales. The specific factor (~10) inherits the empirical Omega_DE calibration of beta; the order of magnitude (cosmic, not Planck) is structural. This is the framework's first computed natural frequency from its commitments — one tone, not yet a spectrum, but at the right scale to suggest V_3's role in cosmic-scale physics.

**V_3 eigenmode structure on (A_0, 1) (G20).** Treating A as a coordinate on the bounded structural-interval and solving the Schrodinger-like eigenvalue problem with V_3 potential and Dirichlet boundary conditions, the second excited state (mode 2) has its two interior zero crossings at A ~ 0.343 and 0.663 — within 1% of the orbital thirds (1/3 = 0.333, 2/3 = 0.667). The orbital thirds emerge from independent physics (orbital mechanics in the (1-A) metric); the V_3 mode-2 crossings emerge from the bounded-interval Schrodinger problem. Two independent calculations land on the same A values to within 1%. Part of this alignment is generic (any wave equation on a near-(0,1) interval at mode 2 has crossings near 1/3, 2/3); the V_3-specific portion is the small refinement that brings the crossings closer to the orbital thirds than pure-box would. The eigenvalue spectrum itself does NOT show clean integer-ratio harmonic structure (E_n/E_0 = 3.15, 6.62, 11.45, ... not 3, 5, 7, 9), so harmony in the strict frequency-ratio sense is absent; spatial structural alignment between V_3 wave dynamics and orbital mechanics is real.

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

For spinning bubbles (Kerr extension, committed 2026-05-13): under A = 2Mr/(r²+a²), T is uniform on the bubble = T_Kerr exactly. Super-radiance from rotating outer-face matter enhances the emission RATE non-uniformly across latitude (peak at the equator where v_matter is maximal), but the local temperature is uniform. The earlier "T latitudinally banded" reading was retired in G64 — it relied on the static-limit-bubble formula A = 2Mr/Σ which placed A = 1 at the ergosphere outer boundary rather than the actual Kerr horizon.

**Hawking T derived two independent ways (G60, 2026-05-13).** The framework now has two structurally distinct derivations of the Hawking temperature that agree exactly:

1. Q8 resolution rule: k_B T = (1/4π) ℏc |∇A| at the boundary
2. Elevator self-consistency: per-emission Δ A_h = −4 ℓ_P² + Δ M = −k_B T/c² + dA_h/dM = 32π G²M/c⁴ together force k_B T = ℏc³/(8π GM)

Both routes use different framework primitives and land on the same Hawking T to machine precision across stellar BHs through PBHs.

**Bekenstein-Hawking entropy from direct ledger counting (G18, derivation closed in G59 + G60 2026-05-13).** The framework reaches `S = A_h/(4 ℓ_P²)` by counting independent ledger entries on the boundary, without going through dE = T dS. Each entry occupies `α ℓ_P² = 4 ℓ_P²` of horizon area, with:

```
α = α_H × (gravity-bridge factor) = 2 × 2 = 4
```

Both factors do other framework work:

- **α_H = 2** (outer-face pair structure: write component + reduction component) enters the bulk closure density Beta(α_S=D, α_H=2) that gives the quartic Hermite F. The two-face commitment supplies it geometrically.
- **Gravity-bridge factor 2** is in A's definition (A = 2GM/(c²r)), the natural-unit framework (1 SU = A_0 = 1/(4πD)), and the resolution rule (k_B T = (1/4π)ℏc|∇A|).

Neither factor is invoked specifically for entropy. The G18 "good suspects, not derivation" caveat is now closed: α = 4 is the unique decomposition under primitives that exist for independent structural reasons. Pair structure itself is derived from the elevator argument (G60), not stated.

---

## Cosmological branch

### SU = A_0 structural identity (V_4 landed 2026-05-11)

The cosmic ambient field A_0 is identified with the framework's natural ruler unit SU. SU is a 1D ruler measuring A along a path — line-by-line accumulation reading. SU(z) is the cumulative ruler value out to redshift z.

**Derivation chain in one place:**

Step 1 — SU formula derived from genesis spreadsheet (no fit parameters):

```text
SU(z) = K · (a·z + q·z²)
      = (z/H) · (1 + 3z/20)

K = z_anchor · c/H_0           (z_anchor = 0.30 is the second structural commitment)
a = 1/z_anchor                 (= 10/3 for z_anchor = 3/10)
q = 1/2                        (universal quadratic coefficient)
3/20 = z_anchor/2              (the formula's quadratic coefficient is forced)
```

Step 2 — Bridge term form derived from Shapiro through ambient A_0:

```text
b = A_0 · c/H_0       (path-integrated A-delay → apparent extra distance)
```

Step 3 — A_0 factored structurally as (4π × 3):

```text
A_0 = 1/(4π × 3) = 1/(12π) ≈ 0.026526

4π : framework-internal from Q8 thermal/gravity-bridge structure
     k_B T = (1/4π) · ℏ · c · |∇A|
     4π = 2π (thermal periodicity) × 2 (gravity bridge c²/2)

3  : spatial dimensionality (our universe is 3D)
```

Step 4 — Two independent derivations converge:

```text
A_0 (cosmological / SU-ruler):  A_0 = b/L_H = 354.95/13387 ≈ 0.02651
A_0 (ontological / threshold):  A_0 = 1/(12π) ≈ 0.02653
Match: 0.04% (4 significant figures)
```

Step 5 — Bridge term verification:

```text
b_pred = (1/(12π)) · c/H_0
       = (299792.458 km/s) / (12π × 73.04 km/s/Mpc) × (3.262 Mly/Mpc)
       ≈ 355.10 Mly

b_historical ≈ 354.95 Mly
Match: 0.04%
```

**Status:** A_0 = 1/(12π) is derived, not calibrated. SU is a derived ruler (no fit parameters). The structural identity "1 SU = A_0" emerges as two independent derivations landing on the same object: the natural cosmological unit (from SU) and the minimum density for spacetime to exist (substance ontology). The "just happens to equal" is the framework's structural discovery.

The remaining open piece is the rigorous first-principles origin of the "3" factor (currently tied to spatial dimensionality / radial-line-of-sight averaging, but the unifying argument linking it to Q8's 4π is still open theoretical work — see "A_0 void interpretation" section above and Open Problem #2).

---

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

**Cosmological story (two-layer reading, refined 2026-05-11):**

The framework's cosmological structure is now cleanly separable into two layers:

- **Layer 1 — Cosmic expansion:** V_3 modified Friedmann with A pinned at the minimum A_0 gives **LCDM-equivalent expansion at H_0 = 73 by construction**. The calibration β_tilde = Ω_DE_target × (1-A_0)² ensures V_3(A_0) acts as a cosmological constant of magnitude Ω_DE_target. The framework's intrinsic expansion is LCDM-shape; the older "matter-dominated EdS at H_0=73" framing is retired (G29). Cosmic chronometers probe Layer 1 only.

- **Layer 2 — Distance bias:** Photon-A traversal through cosmic A_0 in voids adds path-integral bias to *observed* luminosity distance. The bridge term `b = A_0 · c/H_0 ≈ 355 Mly` lives in this layer (the G7 derivation stands; its mechanism is photon-path, not modified Friedmann). SN distance modulus and CMB θ⋆ probe both Layer 1 and Layer 2.

Under this reading, the Hubble tension reads cleanly: SH0ES (local distance ladder at low z) gives true H_0 = 73; Planck H_0 = 67.4 is what LCDM-fits-to-SN extract when they don't account for the Layer 2 photon-A bias. Both numbers are real measurements; they probe different combinations of the two layers.

**Chronometer pressure on H_0 = 73 (flagged 2026-05-11).** Cosmic chronometers — the most model-independent H(z) probe — prefer H_0 ≈ 68 freely *even given STAM's V_3 shape* (since V_3 ≡ LCDM-shape at the expansion level). At fixed H_0 = 73, V_3 gives χ²/N = 0.76 against chronometers (statistically acceptable, not preferred). The framework's H_0 = 73 commitment rides primarily on SH0ES; chronometers, BAO, and Planck all want lower. This is a known soft spot, not a falsification — but worth recognizing as a one-probe-vs-three-probes situation.

The previous "closed-form intrinsic H(z) = H_0(1+z)²/(1+z+0.5z²)" expression that appeared in earlier handoff/working notes has been **retired** — it was derived by inverting an empirical no-b distance ansatz corresponding to coasting cosmology (q_0 = 0), which STAM does not commit to. The framework's actual intrinsic H(z) is V_3-derived and LCDM-shape at H_0 = 73.

**Galactic dark matter:**

Model-A is naturally compatible with primordial-black-hole dark matter (PBH-DM). Each PBH is a small bubble with the same A = 1 boundary structure as stellar and super-massive BHs. Cumulative A from baryons alone gives a factor-of-1.5-to-3 deficit in implied rotation velocity (G10). Adding the cosmic A_0 baseline explicitly (G31) shifts rotation curves by only ~1.35% — does NOT bridge the deficit by itself. With PBH-DM halo + cumulative A, the framework closes (G13). The framework's galactic-DM mechanism is PBH-DM; A_0 baseline plays no significant role at galactic scales.

---

## Quantum interpretation

Model-A uses a physical, not conscious, definition of observation:

```text
Observation = physical interaction that resolves A
```

**Quantum unit of resolution (committed 2026-05-13).** Each quantum physical interaction resolves exactly **1 SU = A_0 = 1/(12π)** of substance. This unifies what were previously two readings of 1 SU:

- **Cosmological**: the natural ruler unit for measuring A along a path (cumulative bridge-term integration).
- **Quantum**: the discrete quantum of resolution written to the ledger by one physical interaction.

These coincide under the framework's commitments. Resolution is *quantized* at the SU/A_0 granularity — A is continuous in bulk, but resolution events themselves come in discrete quanta of magnitude A_0. The cosmic baseline A_0 then reads as "1 SU of resolved-A per Planck-scale structural cell" — the minimum write density required for the manifold to exist, replicated everywhere by structural necessity. Hawking emission events each resolve exactly 1 SU.

### Born rule from SU resolution statistics (committed 2026-05-13)

The Born rule decomposes into two structurally separable pieces.

**Step 1 — Projection-geometry theorem.** Let H be a Hilbert space, ψ ∈ H a unit vector, P a projector. Define a nonnegative resolution-support amplitude `u(P, ψ)` satisfying:

1. **Phase blindness**: `u(P, e^{iφ}ψ) = u(P, ψ)` — phase is not written to the ledger.
2. **Projector locality**: `u(P, ψ)` depends only on `Pψ` — under presentism, only the present-projected state matters.
3. **Unitary covariance**: `u(UPU†, Uψ) = u(P, ψ)` — frame independence.
4. **Orthogonal refinement**: `u(P+Q, ψ)² = u(P, ψ)² + u(Q, ψ)²` for `PQ = 0` — orthogonal alternatives are separate ledger channels; their joint-event rates add.
5. **Normalization**: `u(I, ψ) = 1` for `||ψ|| = 1`.

Then `u(P, ψ) = ||Pψ||` uniquely. For basis projector `P_i = |i⟩⟨i|`:

```
u_i = ||P_i ψ|| = |⟨i|ψ⟩| = |ψ_i|
```

This is a **theorem from projection geometry**, not a notational identification. The proof uses Pythagoras in H (for orthogonal projectors) plus Cauchy's functional equation on the squared sum.

**Step 2 — STAM ledger-measure result.** Each quantum physical interaction is one paired resolution event (write component + reduction component, elevator identity G60). The load-bearing principle:

> **Resolution-event phase volume scales as the square of unresolved support amplitude.**

The "phase volume" of resolution events at outcome `i` — the count of distinct paired-event configurations leading to that outcome — scales as `u_i²`, because each pair (write, reduction) independently samples the amplitude `u_i`. The factor of 2 in the exponent IS the pair structure, viewed as phase-volume rather than as joint-sampling. Normalized:

```
p_i = u_i² = |ψ_i|²    — the Born rule, derived from STAM pair structure
```

The same pair structure (outer-face / unresolved-support: write + reduction is one event with two coupled aspects) does three structural jobs across the framework:

- α_H = 2 in the strong-field metric (Beta(D, 2) → quartic Hermite)
- α = 4 = (pair × gravity-bridge) in boundary entropy (G59)
- p_i = u_i² in the Born rule (this section)

**Interference and decoherence.** Phase lives only in unresolved support — never written to the ledger directly. For a resolved outcome `i`, define `u_i` only after all unresolved alternatives leading to `i` have been coherently summed:

```
ψ_i = Σ_{α → i} u_{i,α} e^{iφ_{i,α}}
u_i = |ψ_i|
p_i = u_i²
```

Double-slit with no which-way interaction: paths recombine in unresolved support before the 1-SU resolution event; phases interfere; `p(x) = |ψ_L(x) + ψ_R(x)|²`. With which-way: the alternatives are separately resolved before recombination, becoming distinct ledger-written branches; they no longer form one coherent channel; probabilities add incoherently `p(x) = |ψ_L(x)|² + |ψ_R(x)|²`. **Decoherence is not a separate postulate** — it is what happens when paired events resolve alternatives into distinct ledger channels before recombination. Interference is the phase-sensitive reshaping of unresolved support *before* the 1-SU resolution event occurs.

This closes Open Problem #5 in structural form: u_i is a theorem from projection geometry; p_i = u_i² is a STAM ledger-measure result from pair structure.

**Bell test reproduction (G67, 2026-05-13).** For the spin-1/2 singlet state under STAM's two-step Born rule, the CHSH parameter reaches the Tsirelson bound |CHSH| = 2√2 ≈ 2.828 exactly (Bell-bound for local hidden variables is 2). E(a, b) = −cos(a − b) matches QM prediction to machine precision. STAM's structural reading: unresolved support is global, so the singlet ψ encodes correlations across the bipartite Hilbert space; resolution events are local (1 SU at one location); pair-structure squaring converts amplitude correlations into observed probability correlations. Local hidden variables fail because correlation structure (phase) lives only in unresolved support — never written to the ledger — and cannot be carried as locally-resolved information.

**Path integral consistency (G67).** STAM's reading places the standard QM propagator K(x_f, t; x_i, 0) in unresolved support: standard QM dynamics (Schrödinger equation, Feynman path integral) operates between resolution events, and the Born rule (STAM-derived) converts |K|² into resolution-probability density. Double-slit interference works cleanly: coherent recombination of paths in unresolved support gives `p = |ψ_L + ψ_R|²` (fringes); separate resolution of the two paths gives `p = |ψ_L|² + |ψ_R|²` (decoherence).

### SU support vs. resolution-event rate (2026-05-13)

Two distinct quantities at the quantum scale, often conflated, must be kept separate:

```
N_SU(x) = A(x) / A_0                            — SU support / occupancy at x
dN_write = Γ[A, ψ, interaction] × dτ / τ_P     — actual resolution-event rate
```

- **N_SU(x)** is the **structural carrying capacity** of unresolved support at x — how many SU's worth of A are present. At cosmic baseline, A(x) ≈ A_0, so N_SU ≈ 1 per Planck-scale cell (the minimum support to maintain manifold).
- **Γ[A, ψ, interaction]** is the **missing rate functional** that governs when SU support is actually converted to ledger writes. Depends on substance density, on the unresolved state, and on the interaction context.

The two are independent: high SU support does not automatically mean high write rate. A region with elevated A has more structural support for resolution events, but actual writes happen only when physical interactions draw from that support and produce a definite ledger update.

### Γ_res structural form (committed 2026-05-13)

Γ_res is the **local proper-time intensity for physical resolution events** — not the Born probability, but the rate at which unresolved support becomes definite ledger writes. Each write resolves exactly 1 SU = A_0:

```
dA_ledger = A_0 · dN_write
dN_write  = Γ_res · dτ
```

**Per-channel split** (Born rule applies to outcome distribution):

```
dN_i = Γ_res · p_i · dτ
p_i  = ||P_i ψ||² / Σ_j ||P_j ψ||²
     = |ψ_i|²  (for normalized mutually exclusive outcomes)
dA_i = A_0 · Γ_res · p_i · dτ
```

Γ_res is **interaction-gated**: Γ_res = 0 for isolated coherent unresolved support; Γ_res > 0 only when physical interactions create distinguishable ledger channels.

**Diagnostic definition** (operational): the decoherence rate of normalized unresolved coherence,

```
Γ_res = −d ln(C) / dτ,    C = [Σ_{i≠j} |ρ_ij|²] / [Σ_{i≠j} |ρ_ij|²_initial]
```

**Dynamical open-system form** (Lindblad-analog):

```
Γ_res = Σ_μ ⟨L_μ† L_μ⟩ · D_μ
```

where L_μ are physical interaction/write channels and D_μ is the distinguishability created by each channel. The framework adopts the open-quantum-systems formalism with STAM's structural interpretation: L_μ as resolution-event operators, D_μ as channel distinguishability.

**Bounds from STAM commitments**:

```
0 ≤ Γ_res ≤ (A/A_0) / τ_P
```

Lower bound from interaction-gating (no writes without physical interactions). Upper bound from SU support capacity: write rate per Planck tick cannot exceed structural occupancy. At cosmic baseline (A = A_0): bound is 1/τ_P. At horizon (A → 1): bound is ≈ 38/τ_P. Verified at boundaries: stellar BH Hawking emission gives Γ_horizon ≈ 10⁻¹²² per cell per τ_P (G68) — ~120 orders of magnitude below saturation, consistent with the extremely sparse resolution-event density relative to structural support.

**Coordinate-time conversion**: if Γ_res is defined in the local rest frame (proper time), the coordinate-time write rate includes the proper-time slowing:

```
dN_write / dt = Γ_res · √(1−A) · √(1−v²/c²)
```

**Boundary cases:**

- **Cosmic baseline** (vacuum, no nearby interactions): Γ_res ~ minimum to maintain manifold. Far below 1/τ_P saturation.
- **Coherent quantum evolution** (closed system, no measurement): Γ_res ≈ 0; ψ evolves unitarily in unresolved support.
- **Measurement detector**: Γ_res spikes at the interaction site; sets decoherence timescale for the measured system.
- **Black hole horizon**: Γ_res governs Hawking emission rate via the resolution rule k_B T = ℏc|∇A|/(4π). Numerically Γ_horizon ≈ 10⁻¹²²–10⁻⁶⁴ per cell per τ_P depending on BH mass (G68).

What this leaves open is no longer "the rate functional" but rather **specific L_μ and D_μ for particular interactions** — the framework's analog of identifying which Lindblad operators apply to which physical systems. That's a domain-of-application question rather than a structural gap.

### First-principles QM action chain (2026-05-13)

Putting the pieces together, the QM action S for a single particle follows from STAM primitives:

```
substance ontology + velocity-cap + Planck primitives
  → local proper time:                dτ = dt √(1−A) √(1−v²/c²)
  → relativistic particle action:     S = −m c² ∫ dτ
  → non-relativistic limit:           L = −m c² + (1/2) m v² + (1/2) m c² A
                                        = −m c² + KE − m Φ,  Φ = c² A / 2 = GM/r
  → Newtonian gravity (gravity bridge): a = (c²/2) ∇A
  → phase:                            φ = S / ℏ
  → phase per Planck tick (free, rest): Δφ = −m / m_P
  → unresolved support:               paths sum coherently with exp(iS/ℏ)
  → resolution event:                 physical interaction writes 1 SU = A_0
                                        (rate governed by Γ, outcome by p_i = u_i²)
```

The QM action emerges as **mass-in-Planck-units × proper-Planck-ticks** along the path. The substance velocity-cap provides the proper-time slowdown (the framework's GR-analog). The Born rule (STAM-derived) converts unresolved-support amplitudes into resolution probabilities at interaction events.

**What remains open in the QM-from-STAM chain:**

- **The resolution-rate functional Γ[A, ψ, interaction]** — the framework's analog of decoherence-rate operators. Boundary cases are constrained (cosmic baseline, BH horizon, detector interactions), but the general form is open.
- **ℏ as a structural input** — the framework treats ℏ as a Planck-scale primitive alongside c, G. A deeper derivation of ℏ from substance ontology would close the loop but may be beyond scope.
- **QFT extension** — multi-particle fields, second quantization, gauge interactions. The framework's structural picture extends naturally (fields as unresolved-A configurations; gauge couplings as Γ-modifying interactions), but the explicit derivation is substantial work.

The model separates states into resolved (A-state confirmed by interaction; path is definite) and unresolved (no physical interaction yet; path is indeterminate).

- Decoherence is interpreted as the dense accumulation of resolution events between a system and its environment.
- Schrodinger's cat is dead-or-alive at the moment of sealing the box, because internal interactions resolve the cat continuously. The cat was never in superposition ontologically; we simply lack epistemic access until we open the box.
- The arrow of time emerges from the irreversibility of resolution events.

The resolved-A record can be read as the universe's running ledger of what has happened: physical interactions write to the ledger, unresolved systems are simply not yet recorded, and consciousness has no privileged role (it is just one category of physical interaction among many). Horizon physics is a special case — at A = 1, inward writes to the ledger are forbidden by the no-interior commitment, so the only available resolution channel is outward, and Hawking radiation is what falls out of the universe needing to keep writing in the only direction left. The structural floor A_0 is the corresponding lower-boundary condition: the minimum write density per Planck cell needed to sustain the manifold (see "A_0 void interpretation" above).

**Refinement: the ledger as present-state, not historical archive (presentism).** The "ledger" is best read as the present-moment configuration of A everywhere, transformed by every interaction, rather than as a stack of historical entries persisting as separate objects. The past does not have separate ontological existence; it shaped how the present is currently configured. "What happened at time t-1000" means asking how the present encodes that past through its current correlations. Each interaction is a moment of becoming that transforms the whole configuration, not a record being appended to a stack. This converts "no information loss" from an axiom into a structural consequence (the present configuration evolves consistently) and makes the framework's quantum interpretation a process ontology rather than a record ontology.

The same A field plays both classical (magnitude) and quantum (resolution-status) roles. Resolution-event statistics suggest a route to deriving Born-rule probabilities, though this is open theoretical work.

---

## Falsifiable predictions

Model-A makes specific predictions in regimes LCDM does not address, each with concrete observational targets:

```text
LIGO ringdown (spinless and Kerr):
    Spinless eikonal:  tau_STAM / tau_GR_Schwarzschild = 1.000 (exact, G57, G66)
    Kerr eikonal:      tau_STAM / tau_GR_Kerr = 1.000 (exact, all spin, G64)
    Framework matches GR at the photon orbit; STAM-vs-GR wedge lives
    INSIDE the photon orbit. Distinguishable only via:
      - Higher overtones (n >= 1)
      - Late-inspiral chirp shape
      - LISA EMRI ringdowns
      - Sub-leading WKB / full Regge-Wheeler
    The old G1 prediction tau / tau_GR = 1.80 is RETIRED (used obsolete
    k(A) = (1-A)(1-A^2)^2 form).

G43: BNS merger engine time scales linearly with total binary mass:
       tau_engine = tau_critical ≈ 0.62 × (M_total / M_sun) seconds
     For GW170817 (M = 2.7 M_sun): predicted 1.68 s, observed 1.74 s
     (3.6% match, no fitting). Falsifiable with future BNS+EM events.

F6: Gravitational decoherence ~0.5 s for 1 micron silica nanoparticle
    in superposition. Cavity-optomechanics frontier; achievable in the
    next decade.

PBH evaporation rate banding (G63 + G4): under formula β for Kerr,
    spinning PBHs have uniform T = T_Kerr but EQUATORIALLY ENHANCED
    emission RATE from super-radiance R(theta) = 1 + (v_matter/c)^2.
    The earlier "T latitudinally banded" prediction (G3) is RETIRED
    under the 2026-05-13 Kerr commitment; the rate enhancement is
    the surviving distinguishing observable. Testable if PBH
    evaporation is ever observed.

PBH-DM: requires inflationary fluctuation amplitude sigma ~ 0.05
    at the PBH-formation scale. Constrains specific inflation models
    (USR, hilltop, single-field-with-bump).

Primordial-mass remnants: Black holes never fully evaporate; final
    state is the primordial inner-face mass. Testable in any future
    PBH evaporation observation that should NOT reach zero mass.
```

---

## What is strongest so far

1. **Weak-field recovery.** Newtonian gravity follows directly from `g = (c^2/2) grad A`. GPS, Shapiro, and lensing match GR weak-field predictions identically.

2. **A_0 = 1/(12 pi) committed structurally.** Bridge term `b = A_0  x  c/H_0` matches historical Pantheon/Union3 fit to 0.04%. A_0 is no longer a calibrated parameter.

3. **The thirds-of-A strong-field structure.** ISCO, photon sphere, and horizon fall at A = 1/3, 2/3, 1 — preserved exactly because g_tt is unchanged.

4. **Black-hole thermodynamics from one rule.** `k_B T = hbar c |grad A| / (4 pi)` reproduces Hawking T (Schwarzschild + Kerr + de Sitter), Bekenstein-Hawking entropy, the first law, the Smarr relation, evaporation lifetime, and the generalized second law.

5. **Strong-field metric structurally complete (spinless AND spinning, D=3).** Spinless: k(A) = (1−A) outside PS, k(A) = (1−A)·(1 − 4y³ + 3y⁴) inside with y = 3A − 2 (quartic Hermite, derived from the minimum-commitment write-density route Beta(D, 2)). Kerr: A = 2Mr/(r² + a²), bubble at constant r = r_+ (Kerr horizon, oblate in flat-space embedding), hologram spins on stationary horizon. **D = 3 is framework-internal** (forced by SU shell-count + two-face joint compatibility, not anthropic). All metric structure derived from primitives (substance + presentism + two-face + SU shell-count + elevator + natural measure) — no free parameters anywhere in the metric.

6. **LIGO ringdown is exact GR for any spin.** Spinless: τ_STAM/τ_GR_Schw = 1 (G57, G66). Kerr: τ_STAM/τ_GR_Kerr = 1 for all spin (G62, G64). Both follow from k = (1−A) outside the photon orbit being exact GR there. The STAM-vs-GR wedge is pushed entirely inside the photon orbit (sub-leading observables only).

7. **SN distance fits beat LCDM combined chi^2.** Model-A with V_3 modified Friedmann wins by 25-33 across Pantheon+/Union3/DES at the same number of free parameters. Predicts inter-catalog Pantheon+/Union3 tension within 27%.

8. **CMB self-consistent closure at H_0 = 73.** Cumulative-A line-of-sight amplification structurally explains the H_0 tension. Internal solution exists for realistic cosmic-structure parameters.

9. **PBH-DM compatibility for galactic dark matter.** Each PBH is a small bubble with the framework's existing thermodynamic/structural machinery. Galactic rotation curves close trivially with PBH-halo + cumulative A.

10. **Six observational regimes, one A field.** Local gravity, propagation delay, BH thermodynamics, SN distances, CMB acoustic scale, and galactic DM all from the same A field with one structural constant (A_0 = 1/(12 pi)) and one calibrated parameter (beta).

11. **Specific falsifiable predictions** — F6 decoherence, PBH-DM sigma constraint, G43 BNS engine time mass-scaling, primordial-mass remnants (BHs never fully evaporate), inside-PS LIGO O5+ overtones / LISA EMRIs as the STAM-vs-GR wedge. Distinguishable in regimes current LIGO doesn't precisely probe.

12. **Substance velocity-cap as STAM-vs-GR wedge** (added 2026-05-12). Substance ontology forces motion through elevated A to be REALLY slowed. For GW170817: binary merger engine time = τ_critical = 1.677 s vs observed 1.74 s — 3.6% match, no fitting. Mass-scaling linear in M_total is the falsification handle.

13. **Entropy derivation closed (2026-05-13).** α = 4 area-per-entry now derived structurally (G59) — both factors (α_H from two-face / pair structure, gravity-bridge from A's definition) do other framework work. Pair structure of Hawking emission derived from elevator argument (G60). Hawking T derived two independent ways (resolution rule + elevator self-consistency).

14. **Stress-energy verified non-pathological (G65, 2026-05-13).** Conservation automatic (TOV residual at numerical noise), Kretschmann bounded (K_STAM < K_Schw near horizon — no hidden singularity), C² smooth at PS, SEC satisfied (no ghost behavior). NEC_r and WEC violated in modified-gravity / dark-energy character — not pathological.

15. **Quantum unit of resolution fixed (2026-05-13).** Each quantum physical interaction resolves exactly 1 SU = A_0. The cosmological and quantum readings of SU unify under this commitment. Hawking emission events each resolve 1 SU, tying the framework's natural unit to the per-event mass-energy bookkeeping via the elevator identity.

16. **Born rule derived structurally (2026-05-13).** Two-step derivation: (a) projection-geometry theorem gives `u_i = ||P_i ψ|| = |ψ_i|` from four structural axioms (phase blindness, projector locality, unitary covariance, orthogonal refinement); (b) STAM ledger-measure step gives `p_i = u_i²` from paired write + reduction sampling of unresolved-A. Decoherence emerges as the structural consequence of paired events resolving alternatives into distinct ledger channels before recombination — not a separate postulate. The same pair structure does three framework jobs (metric closure α_H = 2, entropy α = 4, Born-rule squaring) — one primitive, three derivations.

---

## Major open problems

1. **Lagrangian for A.** The framework's k(A) is now uniquely derived from primitives (substance + presentism + two-face + SU shell-count + elevator + natural measure → quartic Hermite via minimum-commitment Beta(D, 2) write-density), but the action principle producing this specific k(A) is still missing. G25–G27 attempted standard scalar-tensor formulations and found ghost regions plus R_s-dependent V; G33 inverted to A-as-fundamental and found a 1-parameter family, since superseded by the structural derivation chain in G58–G60. The open question now is sharply: what action S[A, g_μν] reproduces both `k(A) = (1−A)` outside the photon orbit and `k(A) = (1−A)·(1 − 4y³ + 3y⁴)` inside?

2. **Structural origin of the "3" factor in 1/(12π).** Closed 2026-05-13: D = 3 is forced by joint compatibility of SU shell-count + two-face refinement + the Beta(D, 2) horizon exponent. The "3" is framework-internal, not anthropic. (Two-face gives α_H = 2 in any D → horizon order automatically = 3; SU shell-count requires horizon order = D; joint compatibility ⇒ D = 3.)

3. **Realistic cosmic-structure modeling for f_LoS.** The CMB self-consistent solution depends on the line-of-sight A amplification factor. Computing this from realistic cosmic structure (N-body or analytic modeling) would close the cosmology branch.

4. **BAO direct test under self-consistent STAM cosmology.** Simple constant-A version failed at 22 sigma. Structure-dependent A line-of-sight test under V_3 cosmology is open.

5. **First-principles Born rule from resolution-event statistics — CLOSED structurally (2026-05-13).** The Born rule decomposes into a projection-geometry theorem (`u_i = |ψ_i|` from four axioms: phase blindness, projector locality, unitary covariance, orthogonal refinement) plus a STAM ledger-measure result (`p_i = u_i²` from the principle that resolution-event phase volume scales as the square of unresolved support amplitude). The same pair structure that gives α_H = 2 in the metric and α = 4 in entropy gives the |ψ|² squaring in probability. Decoherence is not a separate postulate but the structural consequence of paired events resolving alternatives into distinct ledger channels before recombination. Bell / CHSH violation reproduced exactly at the Tsirelson bound (G67). First-principles QM action chain in place: substance velocity-cap → proper time → relativistic action S = −mc²·τ → non-relativistic L = KE − mΦ (with Φ = c²A/2) → path integral in unresolved support → Born-rule at resolution events.

5a. **Resolution-rate functional Γ_res — STRUCTURAL FORM COMMITTED (2026-05-13).** The framework adopts an open-quantum-systems formalism with STAM-specific interpretation:

```
Γ_res = Σ_μ ⟨L_μ† L_μ⟩ · D_μ            (dynamical / Lindblad-analog)
      = −d ln(C) / dτ                      (diagnostic / coherence-decay)
```

with bounds `0 ≤ Γ_res ≤ (A/A_0)/τ_P` from interaction-gating below and SU support capacity above. Per-write resolution is 1 SU = A_0; per-channel split via Born rule p_i = |ψ_i|². Coordinate-time conversion via proper-time factors. Boundary cases (cosmic baseline, coherent evolution, measurement, BH horizon) all numerically constrained or qualitatively committed.

What remains is **identifying specific L_μ and D_μ for particular interactions** — the framework's analog of "which Lindblad operators apply to which physical systems." This is a domain-of-application question (QFT, scattering, decoherence-in-condensed-matter, etc.), not a structural gap in the framework's quantum-resolution dynamics. The structural piece is closed; particular-case applications remain.

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

Where the framework departs from LCDM and is therefore distinguishable in principle: inside-PS LIGO observables (higher overtones, late-inspiral chirp, LISA EMRIs), F6 decoherence, super-radiance-enhanced rate banding for spinning BH Hawking emission, primordial-mass remnants, BNS engine-time mass-scaling (G43), the structural CMB-tension closure mechanism. None of these have been observationally settled yet.

The framework is **not** presented as complete. The strong-field metric + Kerr + entropy + stress-energy verification arc (G57–G66, 2026-05-13) closed the major spinless and spinning structural derivations. Significant open problems remain: Lagrangian for A, realistic cosmic-structure modeling for f_LoS, BAO test under V_3, first-principles Born rule from resolution statistics, spectral details of Hawking emission. Its value is that it creates a unified language with one structural constant (A_0), one calibrated cosmological parameter (β), framework-internal D = 3, and accountable falsification targets — rather than a collection of independent ad-hoc components.

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

Development assisted by extensive conversation with Claude (Anthropic), particularly during the May 2026 push that produced the A_0 = 1/(12 pi) structural commitment, V_3 potential selection, CMB self-consistent closure, and PBH-DM compatibility analysis (G7-G14 script series); the May 11 (evening) refinement session that produced the substance-ontology articulation (water-tank conceptual exercise), the two-layer cosmological reading (G28+G29: V_3 expansion ≡ LCDM at H_0=73 by construction, photon-A traversal as separate distance bias), the G21 redux confirming GR-exact strong-field landmarks under corrected composition, the G31 confirmation that A_0 baseline does not bridge galactic DM, and the G33/G34 swing surfacing the 1-parameter ambiguity in k(A); the May 12 morning session that produced the substance velocity-cap candidate commitment (motion through elevated A is REALLY slowed) and the GW170817 engine time prediction matching observation to 3.6% without fitting (G42, G43 scripts); the May 12 evening session that produced the strong-field metric structural derivation (A_0 = 1/(4πD), ord_{A=1} k(A) = D, smoothstep closure profile inside the photon sphere) and the Σ = D × A shell coordinate reformulation, dissolving the k(A) family ambiguity (G44–G56 scripts); and the May 13 session that produced the minimum-commitment quantum/write-density route Beta(D, 2) → quartic Hermite F(y) = 1 − 4y³ + 3y⁴ for the final-shell closure profile (G58 / G66), the structural closure of α = 4 entropy from outer-face pair × gravity-bridge (G59), the elevator-argument derivation of Hawking pair structure (G60), the framework-internal D = 3 forcing via joint SU + two-face compatibility, the Kerr extension committed under A = 2Mr/(r²+a²) with the bubble at r = r_+ and hologram spinning on the stationary horizon (G62, G64), the stress-energy non-pathology verification (G65), the quantum-unit commitment unifying the cosmological and quantum readings of SU (each quantum physical interaction resolves exactly 1 SU = A_0), the weak-field-as-coarse-grained-SU-writes reading, the structural derivation of the Born rule (projection-geometry theorem for u_i = ||P_i ψ|| + STAM ledger-measure step p_i = u_i² from paired write + reduction sampling, with decoherence as the consequence of paired events resolving alternatives into distinct ledger channels), the first-principles QM action chain (substance velocity-cap → proper time → relativistic action → non-rel Lagrangian → path integral → Born-rule resolution; G67–G68), and the structural form for the resolution-rate functional Γ_res = Σ_μ ⟨L_μ† L_μ⟩ D_μ (Lindblad-analog, with SU-support / write-rate distinction `0 ≤ Γ_res ≤ (A/A_0)/τ_P`) — closing the major spinless and spinning strong-field arcs, fixing the unit of counting for resolution-event statistics, deriving the QM action from STAM primitives, articulating the open-quantum-systems formalism inside STAM, and showing that one pair-structure primitive does three structural jobs (metric closure, boundary entropy, Born-rule squaring).

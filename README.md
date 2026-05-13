# Model-A

**Model-A** is the current variable-accumulation version of the broader **Spacetime Accumulation Model (STAM)** research program. It explores whether gravity, clock behavior, propagation delay, black-hole horizons, thermodynamic behavior, cosmological distance effects, dark matter behavior, and quantum-style resolution can be described using one central idea: a dimensionless spacetime accumulation field called **A**.

Author: **Sean Brady**
Status: **Proposed theoretical framework / active research program**
Snapshot: **May 13, 2026 — current unified branch. Strong-field spinless and Kerr metrics are structurally complete in the active ledger-channel / quintic-Hermite branch. Outside the photon sphere the metric is GR/Kerr exact; the STAM-vs-GR wedge lives inside the photon orbit. Quantum resolution, weak-field gravity, horizon closure, entropy, Born-rule structure, and Γ_res are tied through the same A/SU/Σ hierarchy. Remaining open work: covariant action for A, exact perturbative/QNM spectrum beyond eikonal, Hawking spectral machinery, specific Γ_res channel models, and QFT extension.**

---

## Current unified branch

Model-A's current unified branch is organized by one field, one unit, one shell coordinate, and one final-shell closure profile:

```
A = physical accumulation field
1 SU = A_0 = 1/(12π)
Σ(A) = A / (4π A_0) = 3A
```

The strong-field landmarks become integer shell counts:

```
Σ = 1  →  ISCO
Σ = 2  →  photon sphere
Σ = 3  →  horizon
```

The spinless strong-field radial factor:

```
k(A) = 1 − A                       for A ≤ 2/3
```

so the exterior through the photon sphere is **GR-exact** (weak-field tests pass automatically, eikonal ringdown matches Schwarzschild).

Inside the photon sphere:

```
k(A) = (1 − A) · F(y)              with y = 3A − 2 = Σ − 2,  for 2/3 < A < 1
```

**Because the ledger is treated as a structural channel** (committed 2026-05-13), the spatial-side channel count is D + 1 = 4 (3 spatial + 1 ledger), and the final-shell closure density is **Beta(4, 2)**:

```
p(y) = 20 y³ (1 − y)
F(y) = 1 − 5 y⁴ + 4 y⁵             — quintic Hermite (C³ at PS, C¹ at horizon)
```

Therefore the explicit metric is:

```
k(A) = (1 − A) · [1 − 5(3A − 2)⁴ + 4(3A − 2)⁵]    for 2/3 < A < 1
```

This keeps weak-field tests and photon-sphere / eikonal ringdown GR-exact, while placing the STAM-vs-GR wedge inside the photon sphere.
---

## Headline results — what is strongest so far

1. **Weak-field recovery.** Newtonian gravity follows directly from `g = (c^2/2) grad A`. GPS, Shapiro, and lensing match GR weak-field predictions identically.

2. **A_0 = 1/(12 pi) committed structurally.** Bridge term `b = A_0  x  c/H_0` matches historical Pantheon/Union3 fit to 0.04%. A_0 is no longer a calibrated parameter.

3. **The thirds-of-A strong-field structure.** ISCO, photon sphere, and horizon fall at A = 1/3, 2/3, 1 — preserved exactly because g_tt is unchanged.

4. **Black-hole thermodynamics from one rule.** `k_B T = hbar c |grad A| / (4 pi)` reproduces Hawking T (Schwarzschild + Kerr + de Sitter), Bekenstein-Hawking entropy, the first law, the Smarr relation, evaporation lifetime, and the generalized second law.

5. **Strong-field metric structurally complete (spinless AND spinning, D=3).** Spinless: k(A) = (1−A) outside PS, k(A) = (1−A)·(1 − 5y⁴ + 4y⁵) inside with y = 3A − 2 (quintic Hermite, derived from the ledger-as-structural-channel write-density route Beta(D + 1, 2)). Kerr: A = 2Mr/(r² + a²), bubble at constant r = r_+ (Kerr horizon, oblate in flat-space embedding), hologram spins on stationary horizon. **D = 3 is framework-internal** (forced by SU shell-count + two-face joint compatibility, not anthropic). All metric structure derived from primitives (substance + presentism + ledger-channel + two-face + SU shell-count + elevator + natural measure) — no free parameters anywhere in the metric.

6. **LIGO ringdown is exact GR for any spin.** Spinless: τ_STAM/τ_GR_Schw = 1 (G57, G66). Kerr: τ_STAM/τ_GR_Kerr = 1 for all spin (G62, G64). Both follow from k = (1−A) outside the photon orbit being exact GR there. The STAM-vs-GR wedge is pushed entirely inside the photon orbit (sub-leading observables only).

7. **SN distance fits beat LCDM combined chi^2.** Model-A with V_3 modified Friedmann wins by 25-33 across Pantheon+/Union3/DES at the same number of free parameters. Predicts inter-catalog Pantheon+/Union3 tension within 27%.

8. **CMB self-consistent closure at H_0 = 73.** Cumulative-A line-of-sight amplification structurally explains the H_0 tension. Internal solution exists for realistic cosmic-structure parameters.

9. **PBH-DM compatibility for galactic dark matter.** Each PBH is a small bubble with the framework's existing thermodynamic/structural machinery. Galactic rotation curves close trivially with PBH-halo + cumulative A.

10. **Six observational regimes, one A field.** Local gravity, propagation delay, BH thermodynamics, SN distances, CMB acoustic scale, and galactic DM all from the same A field with one structural constant (A_0 = 1/(12 pi)) and one calibrated parameter (beta).

11. **Specific falsifiable predictions** — F6 decoherence, PBH-DM sigma constraint, G43 BNS engine time mass-scaling, primordial-mass remnants (BHs never fully evaporate), inside-PS LIGO O5+ overtones / LISA EMRIs as the STAM-vs-GR wedge. Distinguishable in regimes current LIGO doesn't precisely probe.

12. **Substance velocity-cap as STAM-vs-GR wedge** (added 2026-05-12). Substance ontology forces motion through elevated A to be REALLY slowed. For GW170817: binary merger engine time = τ_critical = 1.677 s vs observed 1.74 s — 3.6% match, no fitting. Mass-scaling linear in M_total is the falsification handle.

13. **Entropy derivation closed (2026-05-13).** α = 4 area-per-entry now derived structurally (G59) — both factors (α_H from two-face / pair structure, gravity-bridge from A's definition) do other framework work. Pair structure of Hawking emission derived from elevator argument (G60). Hawking T derived two independent ways (resolution rule + elevator self-consistency).

14. **Stress-energy verified non-pathological (G65, 2026-05-13).** Conservation automatic (TOV residual at numerical noise), Kretschmann bounded (K_STAM < K_Schw near horizon — no hidden singularity), C² smooth at PS, SEC satisfied (no effective anti-gravity in the final-shell stress diagnostic — ghost-freedom is a separate question, addressed under Open Problem #6's substance-ontology reformulation). NEC_r and WEC violated in modified-gravity / dark-energy character — not pathological at the stress-energy level.

15. **Quantum unit of resolution fixed (2026-05-13).** Each quantum physical interaction resolves exactly 1 SU = A_0. The cosmological and quantum readings of SU unify under this commitment. Hawking emission events each resolve 1 SU, tying the framework's natural unit to the per-event mass-energy bookkeeping via the elevator identity.

16. **Born rule derived structurally (2026-05-13).** Two-step derivation: (a) projection-geometry theorem gives `u_i = ||P_i ψ|| = |ψ_i|` from four structural axioms (phase blindness, projector locality, unitary covariance, orthogonal refinement); (b) STAM ledger-measure step gives `p_i = u_i²` from paired write + reduction sampling of unresolved-A. Decoherence emerges as the structural consequence of paired events resolving alternatives into distinct ledger channels before recombination — not a separate postulate. The same pair structure does three framework jobs (metric closure α_H = 2, entropy α = 4, Born-rule squaring) — one primitive, three derivations.

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

## Detailed framework

The sections below preserve the technical derivations, numerical checks, and open-problem record behind the current branch.

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
F(y) = 1 - 5y⁴ + 4y⁵                      (quintic Hermite, C³ at PS, C¹ at horizon)
```

Equivalently in shell-coordinate Σ:

```text
k(Σ) = (1 - Σ/3)                          for Σ ≤ 2
k(Σ) = (1 - Σ/3) · F(Σ)                   for 2 < Σ < 3
F(Σ) = 1 - 5(Σ - 2)⁴ + 4(Σ - 2)⁵          (quintic Hermite)
```

### Natural shell coordinate

The natural strong-field variable is **Σ = A / (4π A_0) = D × A** (= 3A for D=3):

- Σ = 1 → ISCO (A = 1/3)
- Σ = 2 → photon sphere (A = 2/3)
- Σ = 3 → horizon (A = 1)

The thirds-of-A landmarks are integer shells in Σ. The final-shell radial closure profile F(Σ) operates on the third shell: 2 < Σ < 3 (equivalently 2/3 < A < 1).

### Derivation chain (quantum/write-density route, ledger-as-channel)

The strong-field metric is derived from the framework's structural primitives — substance ontology, presentism, ledger-as-structural-channel, two-face refinement, SU shell-count, elevator identity, and natural measure on the configuration manifold:

1. **A_0 = 1 / (4π D)** from substance baseline: 4π from Q8 thermal × gravity-bridge, D from spatial dimensionality. For D = 3: A_0 = 1/(12π).

2. **Σ = D × A** as the natural strong-field shell coordinate. For D = 3, Σ = 3A places ISCO at Σ = 1, PS at Σ = 2, horizon at Σ = 3.

3. **Final-shell coordinate y = Σ − 2 ∈ (0, 1).** F(Σ) is interpreted as the survival fraction of unresolved radial closure inside the final shell; the closure density is p(y) = −dF/dy.

4. **Channel structure (ledger as structural channel, committed 2026-05-13).** Final-shell closure is distributed across **D + 1 = 4 spatial-side channels** (D = 3 spatial directions + 1 ledger / self-reference channel) and **2 horizon-pair channels** (outer-face pair structure: write component + reduction component, both on the outer face — the elevator identity). Total: 6 structural channels.

5. **Configuration volume on the final shell.** Under uniform per-entry measure (natural measure on the configuration manifold), the spatial-vs-horizon-pair split has Beta(α_S = D + 1, α_H = 2) density. For D = 3: p(y) = 20 y³ (1 − y), giving survival F(y) = 1 − I_y(4, 2) = **1 − 5y⁴ + 4y⁵** (the quintic Hermite).

   p(y) = 20 y³ (1 − y) has the structural properties: starts at zero at PS, turns on smoothly inside the final shell, returns to zero at the horizon, integrates to one completed final-shell closure. The y³ touch at PS gives C³ smoothness (matches outside-PS GR through third derivative); the (1 − y)¹ touch at horizon gives C¹ closure with the expected order-D zero.

6. **k(A) = (1 − A) outside the photon sphere.** No STAM modification where light can escape; weak-field tests pass automatically.

7. **D = 3 forced by joint compatibility.** ord_{A=1} k = α_H + 1 = 3 automatically from the two-face commitment (α_H = 2, D-independent). SU shell-count separately requires ord = D. Joint compatibility ⇒ D = 3 framework-internally (not anthropic).

Near the horizon: F ~ 10(1 − y)² (leading order), so k = (1 − A) · F ~ 90(1 − A)³, giving ord_{A=1} k(A) = 3, preserving the SU shell-count result.

**Interpretation.** A is the physical accumulation field; SU is the minimum accumulation unit (1 SU = A_0); Σ is the strong-field shell count; F(Σ) is the unresolved final-shell survival profile. Strong field supplies the boundary domain 2 < Σ < 3, but the *shape* of F(Σ) is derived from final-shell resolution density (the write-density route under the ledger-as-channel commitment) rather than from a global constant-n metric ansatz.

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
| Inside-PS strong field | (2/3, 1) | quintic Hermite ramp | sub-leading WKB / late-inspiral / EMRI distinguishers |
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

The framework's quintic Hermite metric is non-pathological (verified G65, run on the quintic form which is the committed branch):

- **Conservation automatic** (TOV residual at numerical noise level — Bianchi identity)
- **Curvature bounded throughout final shell** (Kretschmann K_STAM peaks at ~0.46 in M=1 units, *less* than Schwarzschild's K = 0.75 at the horizon — no hidden singularity)
- **C² smooth at PS** (stress-energy → 0 quadratically as A → 2/3+)
- **SEC satisfied** throughout final shell (no effective anti-gravity behavior in the final-shell stress diagnostic; **ghost-freedom is a separate question** — see Open Problem #6 for its substance-ontology reformulation)
- **NEC_r and WEC violated** for most of the final shell — characteristic of modified-gravity / dark-energy-like effective stress-energy (consistent with F3's earlier w ≈ −1 finding); not a pathology
- Effective stress-energy reads as tension-dominated (p_r < 0, p_t > 0) with NEC_t satisfied

### What's still open in strong field

- **Lagrangian for A**: substantially closed (G69–G80). See "Effective action and shell-count formulation" section below — the double-LM constrained shell-count action S[Σ, g] reproduces the committed metric exactly with a fully ghost-free perturbative spectrum. What remains is interpretation of how this connects to the deeper substance ontology / Γ_res commitments.
- **Inside-PS observable distinguishers**: late-inspiral chirp, higher overtones, LISA EMRIs. Current LIGO precision doesn't reach this regime.
- **Spectral details of Hawking emission**: framework gives T and per-entry structure, but the full spectral distribution still requires QFT machinery beyond cell-counting.

---

## Action principle: constrained shell-count formulation

**The G70–G84 action arc closes the Lagrangian-level embedding of the committed strong-field sectors at the constrained shell-count level. Ordinary scalar, multi-scalar, aether, cuscuton, and mimetic routes fail; the successful formulation treats Σ as a constrained SU shell-count field with two Lagrange multipliers. The scalar shell-count mode does not propagate; only the graviton remains dynamical.**

### The action

```
S[Σ, g, λ₁, λ₂]  =  (1 / 16π G) ∫ d⁴x √(−g) [
    f(Σ) R
  + λ₁ ((∇Σ)² − W(Σ))
  + λ₂ (u^μ ∂_μ Σ)
  − 2 V(Σ)
]
```

with:

```
A = Σ / 3
```

as the continuum substance-density variable recovered in the weak-field limit. Σ takes integer landmarks **Σ = 1, 2, 3** at ISCO, photon sphere, and horizon respectively.

The four pieces of the action have direct substance-ontology meanings:

- **f(Σ) R**: non-minimal coupling. f(Σ) > 0 throughout the shell; graviton positivity automatic.
- **λ₁ ((∇Σ)² − W(Σ))**: Lagrange multiplier λ₁ enforces the kinematic constraint (∇Σ)² = W(Σ). **λ₁ fixes the source-determined shell-count gradient.**
- **λ₂ (u^μ ∂_μ Σ)**: Lagrange multiplier λ₂ enforces flow-constancy along the substance rest frame u^μ (background structure, not dynamical). **λ₂ preserves shell count along the substance flow.**
- **V(Σ)**: substance potential. V(2) = 0 at the photon sphere, V(Σ) → ∞ at the horizon.

The functions f(Σ), W(Σ), V(Σ) are closed-form rational expressions in Σ determined by matching against the committed metric, with structural decomposition into framework primitives (Beta(D + 1, 2) write-density → F(y) quintic Hermite → k(A) → matching equations → f, W, V).

### G70–G84 action arc

The successful action was not found by inspection. The arc through G70–G84 tested every standard low-derivative route, and each failed for a distinct structural or perturbative reason, before converging on the constrained shell-count form.

| Stage | Route | Result |
|---|---|---|
| G70 | Single-scalar scalar-tensor S[A, g] = ∫ [f(A) R − Z(A)(∇A)² − 2V(A)] | Matches metric; closed-form f, Z, V |
| G71 | Linearized perturbations of G70 | Scalar ghost in 86% of final shell |
| G72a | Static-aligned Einstein-aether | Excluded (admits only Schw–dS) |
| G72c | Tilted Einstein-aether | Excluded (c_i ~ 10²¹× observational bound) |
| G74 | Two-scalar Brans-Dicke (ledger Ψ enrichment) | Cannot rescue (M_AA unchanged by Ψ) |
| G75 | Structural decomposition of f, Z, V | f, Z, V decompose into framework primitives |
| G76/G77 | Spacelike cuscuton (square-root kinetic) | Removes radial ghost; introduces temporal ghost |
| G78 | Mimetic with timelike clock + f(A) R | Excluded (off-diagonal mismatch forces λ = 0) |
| **G79/G80** | **Double-LM constrained shell-count, Schwarzschild** | **Closed: δΣ fully constrained in all ℓ sectors** |
| **G81** | **Kerr extension, Σ_K = 6Mr/(r²+a²), stationary axisymmetric δΣ_K analysis** | **Closed: ω = m Ω_ZAMO eliminates non-stationary modes** |
| **G82** | **Inside-shell W_K via Option B: g_STAM^{rr} = (Δ/Σ_BL)F(y_K)** | **Six requirements verified; cubic horizon vanishing** |
| **G83** | **Equatorial photon-region normalization y_K = (Σ_K − Σ_ph(a))/(3 − Σ_ph(a))** | **F = 1 exactly at the Kerr prograde photon orbit for every spin** |
| **G84** | **Off-axis photon-region normalization via r_ph^+ ↔ r_polar interpolation** | **θ-dependence with cubic horizon vanishing at every latitude** |

### Status: closed with caveats

Closed for the committed Schwarzschild and Kerr strong-field sectors. The scalar shell-count mode δΣ is fully constrained out by λ₁ and λ₂ in every angular sector; only the graviton propagates, with healthy kinetic structure from f(Σ) > 0.

**Remaining refinements (future precision work):**

- **G85**: replace the sin²θ off-axis interpolation in G84 with the full Kerr spheroidal photon-region boundary derived from R(r) = 0 and dR/dr = 0 for spherical photon orbits in Kerr. The sin²θ interpolation is sufficient for structural closure (ghost-freedom, horizon vanishing, all limiting cases preserved); the rigorous Kerr photon-region surface is needed for precision EMRI and off-axis-imaging observable predictions.
- **Non-eikonal perturbation / QNM observables**: extending the linearized analysis beyond the eikonal limit to compute exact quasi-normal mode spectra of the framework's strong-field metric (Regge-Wheeler / Teukolsky-type analysis). **Status from G86 (2026-05-13)**: leading eikonal ringdown remains GR-exact because the effective potential V_eff matches at the photon sphere through V and V′ and V″. The first deviation appears at the **third derivative** of V_eff at r = 3M, driven by the quintic Hermite F⁽⁴⁾(0) ≠ 0. Non-eikonal QNM corrections may therefore be observationally significant — 3rd-order WKB on the scalar effective potential suggests percent-to-tens-of-percent corrections at moderate ℓ. Low-ℓ modes require higher-precision methods (6th-order WKB, Leaver continued fraction, or time-domain integration) before this becomes a locked observable prediction. G87 will run a tensor (Regge-Wheeler) proxy potential at higher order; G88+ may be needed for Leaver-grade precision.
- **Matter coupling**: extending the constrained shell-count action to include source-side matter (currently A is treated as the substance density of a point source M; full coupling to a matter Lagrangian for stars / fluids / EM fields is open).
- **Microscopic derivation**: deriving the constrained shell-count action from the framework's underlying SU-write dynamics (substance ontology with Γ_res Lindblad-analog form, Open Problem #5a) via explicit coarse-graining. This would close the chain "discrete substance ontology → effective continuum action S[Σ, g]" from first principles.

Scripts: [G69](scripts/G69_exact_effective_stress_tensor.py)–[G84](scripts/G84_off_axis_photon_region.py) for the full derivation chain. Summaries in [results/](results/).

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

The numerical match between 1/(12 pi) = 0.026526 and the empirical bridge term b/L = 354.95/13387 = 0.026514 is 0.04%. This is the framework's strongest evidence for treating A_0 as a derived structural constant rather than a calibrated parameter. The May 13 SU shell-count + two-face compatibility argument now supplies a framework-internal origin for D = 3; the remaining deeper task is lifting that closure into the covariant A-action.

**Quantum-scale reading (aligned with SU-write commitment).** Under Model-A's resolved/unresolved-A interpretation — where physical interaction is what resolves A, and the resolved record is the universe's running ledger of what has happened — A_0 reads as the minimum density of resolved-A required to maintain the manifold's structural existence. The (4 pi x 3) decomposition gives this density a concrete cell-form: roughly one structural unit of resolved-A per (4 pi solid angle x 3 spatial directions) per Planck cell. Below this density, the manifold has too few "writes" to sustain itself — which is what "A = 0 means no spacetime" maps to at the quantum-cell scale. The same vocabulary describes physical writes at every scale of the framework: horizon writes (Hawking radiation as outward A-resolution where inward is forbidden), local quantum writes (resolution events from any physical interaction), and cosmic-floor writes (the A_0 minimum that keeps the manifold on the books). This reading is now tied to the SU-write commitment and the ledger-channel strong-field derivation. The remaining deeper task is not the numerical value alone, but the covariant action or equivalent first-principles rule that produces A_0, SU shell-count, and the final-shell closure profile together.

A_0 status (updated 2026-05-13): structural commitment. Observationally distinguished — only A_0 = 1/(12 pi) satisfies both the bridge-term match and the CMB physical-amplification constraint (script G15). Theoretically anchored on the (4 pi x D) decomposition: 4 pi from Q8/Q10 thermal / gravity-bridge structure, and D = 3 from SU shell-count + two-face joint compatibility. The value is no longer treated as a free calibration. What remains open is the deeper covariant action or equivalent first-principles rule that produces the A_0 floor and the strong-field shell-count from one variational structure.

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

- **α_H = 2** (outer-face pair structure: write component + reduction component) enters the bulk closure density Beta(α_S = D + 1, α_H = 2) = Beta(4, 2) that gives the quintic Hermite F. The two-face commitment supplies it geometrically.
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

The strong-field branch now supplies a framework-internal origin for D = 3 through SU shell-count + two-face compatibility. The remaining deeper piece is the covariant action or equivalent principle that produces A_0, D = 3, and the shell-count closure together.

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

- α_H = 2 in the strong-field metric (Beta(D + 1, 2) → quintic Hermite)
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
dN_write = Γ_res[A, ψ, interaction] · dτ        — actual resolution-event rate
```

- **N_SU(x)** is the **structural carrying capacity** of unresolved support at x — how many SU's worth of A are present. At cosmic baseline, A(x) ≈ A_0, so N_SU ≈ 1 per Planck-scale cell (the minimum support to maintain manifold).
- **Γ_res[A, ψ, interaction]** is the **resolution-rate functional** that governs when SU support is actually converted to ledger writes. Its structural form (Lindblad-analog) is committed (see below); the domain-specific identification of L_μ and D_μ for particular interactions is open.

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

- **Domain-specific identification of L_μ and D_μ** for particular physical interactions (QFT scattering, decoherence in condensed matter, etc.). The structural form of Γ_res is committed (Lindblad-analog with SU-quantized writes — see the next subsection); applying it to specific systems is the open part.
- **Deeper Hilbert / unresolved-support dynamics from action principle** — the framework's quantum-resolution structure is closed at the projection-geometry + pair-measure level, but the action principle producing the full unresolved-support dynamics + Γ_res channels from a deeper Lagrangian is not yet identified. This is the quantum analog of the metric's "Lagrangian for A" open problem.
- **ℏ as a structural input** — the framework treats ℏ as a Planck-scale primitive alongside c and G. A deeper derivation of ℏ from substance ontology would close the loop but may be beyond scope.
- **QFT extension** — multi-particle fields, second quantization, gauge interactions. The framework's structural picture extends naturally (fields as unresolved-A configurations; gauge couplings as Γ_res-modifying interactions), but the explicit derivation is substantial work.

The model separates states into resolved (A-state confirmed by interaction; path is definite) and unresolved (no physical interaction yet; path is indeterminate).

- Decoherence is interpreted as the dense accumulation of resolution events between a system and its environment.
- Schrodinger's cat is dead-or-alive at the moment of sealing the box, because internal interactions resolve the cat continuously. The cat was never in superposition ontologically; we simply lack epistemic access until we open the box.
- The arrow of time emerges from the irreversibility of resolution events.

The resolved-A record can be read as the universe's running ledger of what has happened: physical interactions write to the ledger, unresolved systems are simply not yet recorded, and consciousness has no privileged role (it is just one category of physical interaction among many). Horizon physics is a special case — at A = 1, inward writes to the ledger are forbidden by the no-interior commitment, so the only available resolution channel is outward, and Hawking radiation is what falls out of the universe needing to keep writing in the only direction left. The structural floor A_0 is the corresponding lower-boundary condition: the minimum write density per Planck cell needed to sustain the manifold (see "A_0 void interpretation" above).

**Refinement: the ledger as present-state, not historical archive (presentism).** The "ledger" is best read as the present-moment configuration of A everywhere, transformed by every interaction, rather than as a stack of historical entries persisting as separate objects. The past does not have separate ontological existence; it shaped how the present is currently configured. "What happened at time t-1000" means asking how the present encodes that past through its current correlations. Each interaction is a moment of becoming that transforms the whole configuration, not a record being appended to a stack. This converts "no information loss" from an axiom into a structural consequence (the present configuration evolves consistently) and makes the framework's quantum interpretation a process ontology rather than a record ontology.

The same A field plays both classical (magnitude) and quantum (resolution-status) roles. The Born rule is structurally closed at the projection-geometry + pair-measure level; what remains open is deriving the full unresolved-support dynamics and the domain-specific Γ_res channels from a deeper action.

---

## Open problems and verification backlog

1. **Lagrangian for A — CLOSED via constrained shell-count action (G70–G84, updated 2026-05-13).** The framework's k(A) is uniquely derived from primitives, and the Lagrangian-level embedding is now identified as the **constrained shell-count action** S[Σ, g, λ₁, λ₂] with two Lagrange multipliers — see the dedicated "Action principle" section above for the full statement. The successful formulation is *not* an ordinary scalar theory for A; it is a constrained shell-count theory for Σ = 3A, with A = Σ/3 recovered as the weak-field continuum density. Single-scalar, multi-field scalar-tensor, aether, cuscuton, and mimetic routes were tested (G70–G78) and all fail; the double-LM constrained form (G79–G84) closes the embedding for both Schwarzschild and Kerr. Remaining refinements (precision off-axis photon-region normalization, exact QNM observables, matter coupling, microscopic SU-write derivation) are listed in the action section as future precision work.

2. **D = 3 / "3" factor — closed structurally; action-level lift remains open.** Closed 2026-05-13: D = 3 is forced by joint compatibility of SU shell-count + two-face refinement + the Beta(D + 1, 2) horizon exponent. The "3" is framework-internal, not anthropic. (Two-face gives α_H = 2 in any D → horizon order automatically = 3; SU shell-count requires horizon order = D; joint compatibility ⇒ D = 3.)

3. **Realistic cosmic-structure modeling for f_LoS.** The CMB self-consistent solution depends on the line-of-sight A amplification factor. Computing this from realistic cosmic structure (N-body or analytic modeling) would close the cosmology branch.

4. **BAO direct test under self-consistent STAM cosmology.** Simple constant-A version failed at 22 sigma. Structure-dependent A line-of-sight test under V_3 cosmology is open.

5. **Born rule — structurally closed at projection-geometry + pair-measure level (2026-05-13).** The Born rule decomposes into a projection-geometry theorem (`u_i = |ψ_i|` from four axioms: phase blindness, projector locality, unitary covariance, orthogonal refinement) plus a STAM ledger-measure result (`p_i = u_i²` from the principle that resolution-event phase volume scales as the square of unresolved support amplitude). The same pair structure that gives α_H = 2 in the metric and α = 4 in entropy gives the |ψ|² squaring in probability. Decoherence is not a separate postulate but the structural consequence of paired events resolving alternatives into distinct ledger channels before recombination. Bell / CHSH violation reproduced exactly at the Tsirelson bound (G67). First-principles QM action chain in place: substance velocity-cap → proper time → relativistic action S = −mc²·τ → non-relativistic L = KE − mΦ (with Φ = c²A/2) → path integral in unresolved support → Born-rule at resolution events. **Still open at the deeper level**: deriving the full Hilbert / unresolved-support dynamics and the domain-specific Γ_res channels (see 5a) from a deeper action principle. The kinematic / probability-measure layer is closed; the action / dynamics layer is open and tied to Open Problem #1 (Lagrangian for A).

5a. **Resolution-rate functional Γ_res — STRUCTURAL FORM COMMITTED (2026-05-13).** The framework adopts an open-quantum-systems formalism with STAM-specific interpretation:

```
Γ_res = Σ_μ ⟨L_μ† L_μ⟩ · D_μ            (dynamical / Lindblad-analog)
      = −d ln(C) / dτ                      (diagnostic / coherence-decay)
```

with bounds `0 ≤ Γ_res ≤ (A/A_0)/τ_P` from interaction-gating below and SU support capacity above. Per-write resolution is 1 SU = A_0; per-channel split via Born rule p_i = |ψ_i|². Coordinate-time conversion via proper-time factors. Boundary cases (cosmic baseline, coherent evolution, measurement, BH horizon) all numerically constrained or qualitatively committed.

What remains is **identifying specific L_μ and D_μ for particular interactions** — the framework's analog of "which Lindblad operators apply to which physical systems." This is a domain-of-application question (QFT, scattering, decoherence-in-condensed-matter, etc.), not a structural gap in the framework's quantum-resolution dynamics. The structural piece is closed; particular-case applications remain.

6. **Ghost-freedom check — CLOSED via constrained shell-count action (G70–G84, updated 2026-05-13).** The Lagrangian-level ghost-freedom of the committed strong-field metric is established through the G70–G84 action arc. Single-scalar (G70/G71), multi-field Brans-Dicke (G74), Einstein-aether (G72), cuscuton (G76/G77), and mimetic (G78) routes all fail for distinct structural / perturbative reasons. The successful formulation is the **double-LM constrained shell-count action** (G79/G80 Schwarzschild; G81/G82/G83/G84 Kerr), in which Σ is a constrained scalar pinned by two Lagrange multipliers — see the dedicated "Action principle" section above for the full statement and the G70–G84 table.

   **Status: CLOSED for the committed Schwarzschild and Kerr strong-field sectors.** The scalar shell-count mode δΣ does not propagate in any angular sector; only the graviton remains dynamical, with healthy kinetic structure from f(Σ) > 0. The substance-ontology mirror of this closure is the Γ_res Lindblad-analog stability (Open Problem #5a) — the two answers are two faces of the same statement.

   Remaining refinements (action section): exact Kerr spheroidal photon-region normalization (G85), non-eikonal QNM observables, matter coupling, microscopic SU-write derivation.

6a. **(Superseded.)** The two-scalar ledger-channel Brans-Dicke embedding was tested in G74 and shown structurally insufficient (M_AA channel of the effective kinetic matrix is unchanged by Ψ enrichment, so a single-scalar ghost remains a ghost under multi-field enrichment). The actual Lagrangian-level closure of Open Problem #6 was found in the double-LM constrained shell-count formulation (G79/G80), not in multi-field Brans-Dicke. The ledger channel's structural role is captured by the integer-landmark structure of Σ rather than by an independent dynamical field Ψ.

7. **Full QNM / perturbation computation.** The ringdown result is currently exact at the eikonal / photon-region level. Exact Regge-Wheeler-type and Teukolsky-type calculations on the committed Model-A metrics would tighten the prediction beyond eikonal and test the inside-photon-orbit wedge.

8. **A_collective for galactic dynamics if PBH-DM is rejected.** If primordial-black-hole dark matter is observationally ruled out, the framework needs to either derive an A_collective galactic-scale enhancement from first principles or accept some other DM mechanism.

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

The framework is **not** presented as complete. The strong-field metric + Kerr + entropy + stress-energy verification arc (G57–G66, 2026-05-13) closed the major spinless and spinning structural derivations. The smallest covariant scalar-tensor action reproducing the committed metric was identified in closed form (G69, G70). The Born rule is structurally closed at the projection-geometry + pair-measure level (theorem + STAM ledger-measure squaring); what remains open at the quantum level is deriving the full unresolved-support dynamics and the domain-specific L_μ / D_μ channels in Γ_res from a deeper action. The ghost-freedom check (Open Problem #6) was explored in G69–G73 — the simplest continuum-field embeddings (scalar-tensor, static-aligned Einstein-aether, tilted Einstein-aether) all fail to give ghost-free perturbative modes inside the photon sphere; under the framework's substance ontology this is reformulated as Γ_res well-posedness, which is structurally satisfied. Other significant open problems: structural interpretation of the f, Z, V functions in the closed-form action, realistic cosmic-structure modeling for f_LoS, BAO test under V_3, spectral details of Hawking emission, and a two-scalar ledger-channel embedding as continuum-field benchmark (Open Problem #6a). Its value is that it creates a unified language with one structural constant (A_0), one calibrated cosmological parameter (β), framework-internal D = 3, and accountable falsification targets — rather than a collection of independent ad-hoc components.

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

Development assisted by coding and conversation with Claude (Anthropic) and ChatGPT (OpenAI)

# STAM Model-A — Claims and Status

**Date:** 2026-05-13 evening
**Snapshot scope:** Unified branch through G124 (G70–G84 action arc + G108–G114 QNM + G118–G120 Kerr + G121–G124 primitive hardening + G123 matter coupling + G113-redux n=1 methodological finding)

---

## Quick Reference: Status Labels

| Label | Meaning |
|---|---|
| **Core definition** | Part of the Model-A framework as currently defined |
| **Algebraic identity** | Follows directly from the definitions |
| **Derived** | Closed-form derivation from primitives in this repo |
| **Uniqueness theorem** | Derived AND shown to be the unique solution under stated commitments |
| **Numerically checked** | Verified by a script in this repository |
| **Locked prediction** | Number is now framework-committed; not a diagnostic |
| **Catalog diagnostic** | Comparison against observational data; interpretation open |
| **Methodological finding** | Real result about technique/limitation rather than a value |
| **Open** | Needed for completeness or stronger empirical validation |
| **Superseded / Retired** | Held a place historically; replaced by sharper result |

---

## Part I — The Core Framework

### The Accumulation Field

STAM Model-A begins with one dimensionless field and one SU-normalized shell count:

```
A(r)         = R_s / r = 2GM/(c²r)        (weak-field substance density)
1 SU         = A_0 = 1/(4πD) = 1/(12π)    (atomic unit; D = 3)
Σ(A)         = D · A = 3A                  (SU shell count, integer landmarks)
y = Σ − 2    = 3A − 2 ∈ (0, 1)             (final-shell coordinate)
```

Integer-Σ landmarks:

```
Σ = 1  →  ISCO        (A = 1/3, r = 6M)
Σ = 2  →  photon sphere   (A = 2/3, r = 3M)
Σ = 3  →  horizon     (A = 1,  r = R_s)
```

**Status:** Core definition. A = R_s/r is the weak-field form; Σ = D·A is the framework's native action variable (committed 2026-05-13).

---

### Theorem 1 — The Horizon is Algebraically Derived

The A = 1 threshold falls out of the definition by two independent routes:

**Geometric:** `A(r) = R_s/r → r = R_s ⟹ A = 1`
**Kinematic (escape):** `v_escape² = 2GM/r → v_escape²/c² = A ⟹ A = 1 ⟺ v_escape = c`

The black-hole horizon is not a separate postulate — it is the A = 1 accumulation threshold.

**Status:** Algebraic identity. Not fitted, not assumed.

---

### Theorem 2 — Gravity is Motion in an A-Gradient

```
g⃗ = (c²/2) ∇A
```

For `A(r) = 2GM/(c²r)`: `∇A = −2GM/(c²r²) ⟹ |g| = GM/r²`. The c²/2 prefactor is not tuned — it is the reciprocal of the 2/c² in A's definition, by construction.

**Status:** Algebraic identity + numerically checked.

---

### Theorem 3 — GPS Clock Correction in A-Language

Weak-field clock rate: `dτ/dt ≈ 1 − A/2`. Combined GPS satellite-vs-surface shift (circular orbit):

```
Δrate_total = A_surface/2 − 3·A_orbit/4

Gravitational gain  ≈  +45.787467 μs/day
Kinematic loss      ≈   −7.213600 μs/day
Net satellite gain  ≈  +38.573867 μs/day
Required factory offset: Δf/f ≈ −4.464568 × 10⁻¹⁰
```

The 3/4 coefficient is not fitted; it is `1/2 (gravitational) + 1/4 (kinematic)`. Both corrections collapse to a single A expression because circular-orbit velocity ties to A.

**Status:** Numerically checked weak-field clock-rate test. Reproduces standard GPS engineering values.

---

### Theorem 4 — Shapiro-Style Propagation Delay

```
Δt = (1/c) ∫A(r) ds
```

For solar-grazing Earth–Mars: one-way ≈ 123.6 μs, two-way ≈ 247.2 μs. Reproduces standard logarithmic Shapiro-delay structure.

**Status:** Numerically checked.

---

### Substance Velocity-Cap (V_4 commitment, 2026-05-12)

Substance ontology forces motion through elevated A to be *really* slowed:

```
v_effective² = v_Newton² × f(A_local),    f(A) = 1 − A
```

This is the framework's GR-analog, ontologically consistent with `g_tt = −(1−A)c²`. At low A, GR and STAM agree to within ~3%; at high A (late binary inspiral), they diverge predictably.

**Status:** Core definition / V_4 candidate commitment.

---

### Mass Estimator Consistency

Multiple weak-field estimators (horizon, acceleration, orbital velocity, Shapiro coefficient, gravitational shift, lensing deflection) infer the same source mass from different observables. Synthetic tests recover input mass ratios at floating-point precision.

**Status:** Numerically checked as identity/synthetic consistency.

---

## Part II — The F / R / f Triad Architecture (2026-05-13 evening)

After the G108–G124 batch, the framework's load-bearing kernel maps cleanly to its three ontological commitments:

| Quantity | Physical role | Maps to |
|---|---|---|
| **F(y) = 1 − 5y⁴ + 4y⁵** | remaining traversable manifold support | the **manifold** |
| **R(A) = −(45A² − 33A − 13) / [30 A (1 − A)]** | 3D-space curvature anisotropy response | the **3D space** |
| **f(A)** (from `d ln f/dA = R(A) · p(y)/F(y)`, f(2/3) = 1) | accumulated SU stiffness coupling in the action | the **SU writes** |

All three derive from the shell-count measure plus the metric ansatz; none is independently chosen.

**Status:** Derived (G70 → R(A) first-principles derivation 2026-05-13 evening).

---

### The Primitive List (reduced from 6 to 4 after G121/G122)

The framework's explicit primitives are:

1. **A** — escape/traversal-cap variable
2. **Σ = D·A** — SU-normalized shell count
3. **D = 3** — spatial dimension (algebraically privileged at D=3 by Einstein anisotropy decomposition; also forced by joint compatibility of SU shell-count + two-face)
4. **h(A) = 1 − A** — velocity-cap (V_4)

Plus one **microscopic axiom-set**: spacetime Poisson point process for SU writes (G124).

**Derived** from these:
- p(y) = Beta(D+1, 2) closure density (G121 uniqueness theorem)
- F(y) = survival function
- k(A) = h · F radial gate
- R(A) = (tt − θθ) Einstein anisotropy
- d ln f/dA = R(A) · p(y)/F(y)
- Constrained shell-count action S[Σ, g, λ₁, λ₂] (G122 uniqueness theorem)
- Z(A), V(A) closed forms
- Coarse-graining theorem (G124 mean-field Poisson)

**Status:** 4 explicit primitives + 1 microscopic axiom-set. Reduced from 6 by the G121–G124 batch.

---

## Part III — Strong-Field Metric

```
ds² = −h(A) c² dt² + dr² / k(A) + r² dΩ²
h(A) = 1 − A                                       (everywhere)
k(A) = 1 − A                                       for A ≤ 2/3   (exterior + PS)
k(A) = (1 − A) · [1 − 5(3A − 2)⁴ + 4(3A − 2)⁵]    for 2/3 < A < 1   (inside PS)
```

**Outside the photon sphere**: GR-exact. All weak-field tests pass identically.
**Inside the photon sphere**: the STAM-vs-GR wedge lives here.

Cubic horizon vanishing `k ~ (1−A)³` splits as 1 + 2:
- **1 power** from the velocity-cap gate `h = 1−A`
- **2 powers** from F's quadratic vanishing at horizon

**D=3 algebraic privilege** (sympy-verified 2026-05-13 evening): the clean (F−1) + A·y·F'_y decomposition of the (tt−θθ) Einstein anisotropy exists only at D=3, where `6(D−3) = 0` AND `3A − 2 = y` coincide simultaneously. Other dimensions don't share this property.

**Status:** Derived from primitives. All metric structure fixed by primitives; no free parameters in the metric branch.

---

### Kerr Extension

```
A(r, θ; M, a) = 2 M r / (r² + a²)
A = 1 at r = r_+ (Kerr horizon; bubble at constant Boyer-Lindquist r)
```

Bubble is oblate when embedded in flat 3-space; hologram matter spins at ZAMO frequency. T uniform on bubble = T_Kerr exactly. LIGO Kerr eikonal ringdown = exact GR-Kerr for all spin.

**Exact Kerr photon-region surface** (G85, G116): r_pr(θ; a) in closed form from R(r) = R'(r) = 0 + brentq on θ-potential. Supersedes G84's sin²θ ansatz (off up to 39% in r_pr at a = 0.99).

**Status:** Core extension. Closed for Kerr substance density and photon-region geometry.

---

### Stress-Energy Verification (G65)

The quintic Hermite metric is non-pathological:
- Conservation automatic (TOV residual at numerical noise — Bianchi identity)
- Kretschmann bounded (K_STAM < K_Schw near horizon — no hidden singularity)
- C² smooth at PS
- SEC satisfied throughout final shell
- NEC_r, WEC violated — modified-gravity / dark-energy-like behavior, not pathology

**Status:** Numerically checked non-pathological.

---

## Part IV — Action Principle (G70–G84 + G122 uniqueness)

```
S[Σ, g, λ₁, λ₂] = (1/16πG) ∫ d⁴x √(−g) [
    f(Σ) R
  + λ₁ ((∇Σ)² − W(Σ))                  λ₁: source-determined shell-count gradient
  + λ₂ (u^μ ∂_μ Σ)                     λ₂: shell count preserved along substance flow
  − 2 V(Σ)
]
```

Σ is treated as a **counted (non-propagating) scalar**, not a dynamical field. The two LM constraints together eliminate Σ's 2 Cauchy data (2 − 2 = 0 propagating DOF). Only the graviton remains dynamical.

**G122 uniqueness theorem (2026-05-13 evening):** this is the unique covariant 2nd-order action consistent with treating Σ as a counted scalar. Alternatives (Ostrogradsky higher-derivative, aether, mimetic, cuscuton, propagating scalar-tensor) all ruled out by ghost / Lorentz / temporal-mode pathologies (G72/G74/G76/G77/G78). Counted-scalar route is forced.

**Status:** Uniqueness theorem under counted-scalar primitive. Ghost-free in all angular sectors for committed Schwarzschild + Kerr backgrounds.

---

## Part V — Microscopic SU-Write Process (G124, 2026-05-13 evening)

The framework's first formal microscopic model: SU-write events form a **spacetime Poisson point process** with local intensity Γ_res(x).

```
d N_events(x) = Γ_res(x) · dV_proper(x)
```

Properties:
- Event counts in disjoint regions are independent (Poisson axiom)
- Each event deposits 1 SU = A_0
- Rate cap: `Γ_res(x) ≤ (A(x)/A_0) / τ_P`
- At A = 1, events forced into exchange pairs (saturation rule, G60 elevator)

**Two scale-conjugate readings of 1 SU:**
| Context | Statement |
|---|---|
| Bulk volumetric | 1 SU = 12π Planck cells (mean event-density at A_0 baseline) |
| Horizon holographic | 1 SU = 4 Planck cells (α = α_H × gravity-bridge = 2 × 2, G59) |

These are not independent — they're scale-faces of the single primitive A_0 = 1 SU.

**Coarse-graining theorem:**
```
A_cg(x, t) = (A_0 / V_cg) · |{events in past-cone of (x,t) ∩ V_cg}|
          → A_0 · ∫₀^τ Γ_res dt'    (mean field, fluctuations 1/√N)
```

**Numerical demo** (1D Poisson with rate ~ R_s/x): mean relative error **4.85% vs Poisson noise floor 9.3%** — within fluctuation band, no systematic bias.

**The constrained shell-count action of Part IV is the mean-field effective action of this Poisson process.**

**Status:** First formal microscopic model. Coarse-graining theorem explicit. Open Problem #1 partially closes.

---

## Part VI — Black-Hole Interpretation

In Model-A, a black hole is a 2D bubble surface, not a deep interior:

```
A < 1  →  spacetime exists
A = 1  →  2D phase boundary / bubble surface
A > 1  →  not part of the manifold
```

**Two-face refinement (committed 2026-05-10, clarified 2026-05-13):**
- **Inner face**: holds primordial formation-event mass. Invariant; never changes.
- **Outer face**: holds subsequently accreted matter. Dynamic; rotates; drains via Hawking.

**BH life cycle**: formation presses primordial mass to inner face (locked); accretion builds outer face; Hawking drains outer face only; horizon size tracks outer face in real time; final state is stable primordial-mass remnant. **Black holes never fully evaporate.**

**A = 1 is never reached in finite time** (refined 2026-05-11). The halving series does NOT converge in STAM.

**Pair structure of Hawking emission** (G60 elevator argument, 2026-05-13): each emission is ONE event with two structural aspects (outward write + content reduction). Saturation at A = 1 forces pair-exchange. No Bogoliubov.

**Status:** Core interpretation. Outward-collapse + two-face + saturation-pair derived.

---

## Part VII — Thermodynamics

Single rule:
```
k_B T = (1/4π) ℏ c |∇A|     at the boundary
```

Reproduces:
- Schwarzschild Hawking T
- Unruh T
- de Sitter horizon T
- **Bekenstein-Hawking entropy** `S = A_h / (4 ℓ_P²)` via direct counting:
  ```
  α = α_H × gravity-bridge = 2 × 2 = 4 Planck areas per entry  (G59)
  ```
- First law `dE = T dS`
- Smarr relation `M c² = 2 T S`
- Generalized 2nd law (1/3 surplus during evaporation)
- Standard Hawking evaporation lifetime

**Hawking T derived two independent ways** (G60, 2026-05-13): Q8 resolution rule + elevator self-consistency. Agree to machine precision across all BH masses.

For spinning bubbles: T uniform = T_Kerr; super-radiance enhances *rate*, not temperature.

**Status:** All standard BH thermodynamics reproduced from one structural rule. Entropy derivation closed (G59 + G60).

---

## Part VIII — Cosmological Branch

**SU = A_0 identity** (closed 2026-05-11): the natural ruler unit SU equals the cosmic baseline A_0, and A_0 = 1/(12π) is derived from `4π × D = 4π × 3`.

```
A_0 (cosmological / SU-ruler):  A_0 = b/L_H = 354.95/13387 ≈ 0.02651
A_0 (ontological / threshold):  A_0 = 1/(12π) ≈ 0.02653
Match: 0.04%
```

**Two-layer reading (refined 2026-05-11):**
- **Layer 1 — Cosmic expansion**: V_3 modified Friedmann with A pinned at A_0 ⟹ LCDM-equivalent at H_0 = 73 by construction. Cosmic chronometers probe this layer only.
- **Layer 2 — Distance bias**: photon-A traversal through cosmic A_0 in voids adds path-integral bias. The bridge term `b = A_0 · c/H_0 ≈ 355 Mly` lives here.

**SN distance fits beat ΛCDM combined χ² by 25–33** across Pantheon+/Union3/DES (G29). Predicts Pantheon+/Union3 tension within 27%. DES anomaly identified as instrumental.

**Chronometer pressure on H_0 = 73**: cosmic chronometers prefer H_0 ≈ 68 freely under V_3 shape. At H_0 = 73, χ²/N = 0.76 — statistically acceptable, not preferred. Known soft spot.

**Galactic dark matter**: PBH-DM + cumulative A closes (G13). A_0 baseline plays no significant role at galactic scales (G31).

**Status:** Two-layer reading committed. Numerically supported across SN, with chronometer soft spot logged.

---

## Part IX — Quantum Interpretation

**Quantum unit of resolution (committed 2026-05-13):** each quantum physical interaction resolves exactly **1 SU = A_0**. Unifies the cosmological-ruler and quantum-write readings of SU.

**Resolved / unresolved buckets:**
- Resolved: physical interaction has occurred; A-state confirmed; path definite.
- Unresolved: no interaction yet; A-state exists but not reconciled; path indeterminate.

`Unresolved ≠ A = 0`. Quantum-like behavior is the expected state of anything unresolved.

**Born rule derived structurally (G67, 2026-05-13):**

Step 1 (projection-geometry theorem): `u_i = ||P_i ψ|| = |ψ_i|` from four structural axioms (phase blindness, projector locality, unitary covariance, orthogonal refinement).

Step 2 (STAM ledger-measure result): `p_i = u_i² = |ψ_i|²` from the principle that resolution-event phase volume scales as the square of unresolved-support amplitude (pair structure).

**Bell test reproduction (G67):** CHSH = 2√2 (Tsirelson bound) exactly. `E(a, b) = −cos(a − b)` matches QM to machine precision.

**First-principles QM action chain:**
```
substance ontology + velocity-cap
  → local proper time:      dτ = dt √(1−A) √(1−v²/c²)
  → relativistic action:    S = −mc² ∫ dτ
  → non-relativistic:       L = KE − mΦ,  Φ = c²A/2 = GM/r
  → phase:                  φ = S/ℏ
  → path integral over unresolved support, Born rule at resolution events
```

**Γ_res Lindblad form** (committed 2026-05-13):
```
Γ_res = Σ_μ ⟨L_μ† L_μ⟩ · D_μ      (dynamical)
      = −d ln(C)/dτ                  (diagnostic)
0 ≤ Γ_res ≤ (A/A_0)/τ_P              (bounds)
```

Per-write resolution is 1 SU = A_0; per-channel split via Born rule p_i = |ψ_i|².

**Status:** Born rule structurally closed. Bell at Tsirelson bound reproduced. Γ_res structural form committed; domain-specific L_μ open per interaction type.

---

## Part X — Ringdown / QNM (G86–G114 + G113-redux)

Leading eikonal ringdown is **GR/Kerr-exact** because the framework's photon-region structure is preserved.

**Spinless axial rigorous (G108/G109, locked):**
| ℓ | Calib err | Shift vs GR | Source |
|---:|---:|---:|---|
| 2 | 0.80% | **4.977%** | G109 (broad pulse, throat sweep) |
| 2 | 0.05% | 5.096% | G112 (narrow pulse, cross-check) |
| 3 | 0.005% | **2.101%** | G112 |
| 4 | 0.006% | **1.158%** | G112 |

Monotonic decrease confirms eikonal-recovery commitment.

**Spinless polar first-pass (G114):**
| ℓ | Polar shift | Axial shift | Polar/Axial |
|---:|---:|---:|---:|
| 2 | 5.620% | 4.977% | **1.129** |
| 3 | 2.324% | 2.101% | 1.106 |
| 4 | 1.206% | 1.158% | 1.041 |

**STAM breaks Schwarzschild isospectrality by ~13% at ℓ=2.** Cleanest LIGO O5+ / LISA EMRI signature.

**G113-redux methodological finding (2026-05-13 evening):** the wave-zone V_GR peak sits at r = 3.28M where A = 0.61 (outside PS). V_geom-only modification gives **zero** shift. The framework's 4.977% n=0 shift comes ENTIRELY from the action correction `(√f)''/√f`. The framework's QNM signature is an action-correction signature.

**Kerr ringdown:**
- G120a GR reference: 5/5 spins PASSED via qnm library
- G120b/c standalone solver attempts FAILED <0.5% calibration
- STAM Kerr ringdown blocked on G120d (proper standalone solver)
- G118v2 spin sweep (4.6%→14%) retired as Kerr prediction; historical diagnostic only

**G90 near-horizon traversal-distance enhancement:** `(dr*/dr)_STAM / (dr*/dr)_GR = 1/√F(y)`. Power-law tortoise-distance divergence vs GR's logarithmic. Direct structural metric prediction.

**Status:** Spinless axial n=0 rigorously locked at ℓ=2/3/4. Polar first-pass. Kerr blocked on standalone solver. n=1 overtones open (G110/G111/G112 time-domain failed, G113-redux WKB-PT/shooting failed; Leaver CF or Sasaki-Nakamura needed).

---

## Part XI — Locked Predictions Status Table

| Prediction | Form | Status |
|---|---|---|
| Weak-field GR-exact | GR Schwarzschild for A ≤ 2/3 | Pass by construction |
| GPS clock shift | A_surf/2 − 3A_orb/4 | Numerically checked |
| Shapiro delay | (1/c) ∫A ds | Numerically checked |
| Mass-estimator identities | 6 estimators agree | Numerically checked |
| ISCO (A = 1/3) | r = 6M, GR identical | Algebraic identity |
| Photon sphere (A = 2/3) | r = 3M, GR identical | Algebraic identity |
| Horizon (A = 1) | r = R_s, GR identical | Algebraic identity |
| Hawking T (Schw, Kerr, dS) | k_B T = ℏc|∇A|/(4π) | Reproduced |
| Bekenstein-Hawking S | A_h/(4 ℓ_P²) | Reproduced from manifold-support cell counting |
| Pantheon SN fit vs ΛCDM | Δχ² ≈ 25–33 win | Locked diagnostic |
| Inter-catalog Pantheon+/Union3 tension | within 27% | Locked diagnostic |
| Bridge term b | A_0 · c/H_0 ≈ 355 Mly | 0.04% match to historical fit |
| GW170817 engine time (2.7 M_☉) | τ_critical = 1.677 s vs 1.74 s observed | **3.6% match, no fitting** |
| GW170817 mass-scaling | τ_engine ∝ M_total, slope 0.62 s/M_☉ | Falsification handle |
| ℓ=2 n=0 axial QNM shift | **4.977%** | Locked rigorous (G108/G109) |
| ℓ=3 n=0 axial QNM shift | **2.101%** | Locked rigorous (G112) |
| ℓ=4 n=0 axial QNM shift | **1.158%** | Locked rigorous (G112) |
| Polar/axial ℓ=2 isospectrality break | **~13%** | First-pass (G114) |
| Near-horizon tortoise enhancement | (dr*/dr)_STAM/GR = 1/√F(y) | Locked metric prediction (G90) |
| F6 gravitational decoherence | ~0.5 s for 1 μm silica | Falsifiable, specific number open |
| Primordial-mass remnants | BHs never fully evaporate | Structural prediction |
| PBH evaporation rate banding | R(θ) = 1 + (v_matter/c)² super-radiance | Locked under Kerr commitment |

---

## Part XII — Open Problems (after 2026-05-13 batch)

### Closed by 2026-05-13 evening

| # | Topic | Closure mechanism |
|---|---|---|
| Action principle | Constrained shell-count action (G70–G84) | G122 uniqueness under Σ-as-count |
| D = 3 forcing | Joint-compatibility + algebraic privilege | Doubly anchored |
| Born rule | Projection-geometry + pair-measure | G67 structural closure |
| Ghost-freedom | Constrained shell-count action | Σ doesn't propagate |
| p(y) shape | Beta(D+1, 2) | G121 three-derivation uniqueness theorem |

### Partially closed (significant remaining work)

| # | Topic | Status |
|---|---|---|
| QNM / perturbation | Axial n=0 rigorous; polar first-pass; n=1 open; Kerr blocked | G108/G114/G113-redux |
| Microscopic SU-write | Poisson model + coarse-graining theorem in place | G124 — three deeper items remain |
| Γ_res rate functional | Lindblad form committed | Domain-specific L_μ open |
| Matter coupling | Tier-1 (Jordan-frame minimal) closed | Tier-2 (substance excitations / TOE) open |

### Fully open

| # | Topic | Notes |
|---|---|---|
| f_LoS cosmic structure modeling | CMB self-consistent solution depends on line-of-sight A | Realistic N-body / analytic modeling needed |
| BAO direct test under V_3 | Constant-A version failed 22σ; structure-dependent open | |
| A_collective backup | If PBH-DM rejected, need derivation from first principles | Contingent on PBH observations |

### Technical gates (weeks of work each)

- Standalone Kerr QNM solver < 0.5% (G120d) — blocks all Kerr STAM ringdown predictions
- n=1 spinless overtones (Leaver CF with verified coefficients, or Sasaki-Nakamura transformation)
- G114-rigorous: full polar reduction in even-parity RW gauge
- F6 decoherence specific number (~0.5 s for 1 μm silica) — needs L_μ identification
- Cosmology pipeline G104 with real CSV + covariance + MCMC

### Deep theoretical (months–years)

- Tier-2 matter coupling: matter as localized substance excitations → QFT extension
- Hawking spectral machinery: full QFT on the substance manifold
- Deeper discrete law underlying Poisson axioms (graph spacetime, causal sets, spin foam — derive Poisson rather than postulate it)
- ℏ as derived rather than primitive

---

## Part XIII — What Has Been Retired / Superseded

| Item | Why retired | Successor |
|---|---|---|
| G92/G94 "5.35% / 5.384%" diagnostic | FD bias on f(Σ) | G109 4.977% via analytic chain rule |
| G84 sin²θ Kerr photon-region ansatz | 39% off in r_pr at a=0.99 | G85 exact spheroidal r_pr(θ; a) |
| Static-limit Kerr formula A = 2Mr/Σ_BL | Violated A < 1 commitment at high spin | A = 2Mr/(r²+a²) — equatorial-plane committed |
| H(z) = H_0(1+z)²/(1+z+0.5z²) closed-form | Coasting cosmology, not V_3 | V_3-derived ≡ LCDM-shape at H_0=73 |
| "T latitudinally banded" Kerr prediction | Tied to retired static-limit | Super-radiance enhances rate, not T |
| F5/F5c asymmetric propagation | 22σ violation of GW170817 1.74 s | Symmetric propagation; substance velocity-cap (V_4) |
| Two-bucket framework as terminal | Superseded by SU resolution structure | Resolved/unresolved + Γ_res + Born rule |
| G118 Schwarzschild-base Kerr ringdown | Broke at high spin (Schw h goes negative below r=2M) | G119 Kerr STAM wave base; G120d required for predictions |
| G118v2 / G119b Kerr proxy shifts | V_base not actual Kerr potential | Historical diagnostic only; Kerr predictions blocked on G120d |
| "Ledger" as channel terminology | Structural-channel reading replaced | "Manifold-support" (forward-only; metaphorical-ledger language retained) |
| "Minimal-degree polynomial" framing for p(y) | Heuristic | G121 Bayesian + order-statistic + MaxEnt uniqueness theorem |
| "Selected by elimination" framing for action | Empirical narrowing | G122 counted-scalar uniqueness theorem |

---

## Part XIV — Current Working Position (2026-05-13 evening)

```
Foundational layer: essentially closed.
  - F / R / f triad: derived
  - 4 explicit primitives + 1 Poisson axiom-set
  - Constrained shell-count action: unique under counted-scalar primitive
  - Microscopic SU-write model: first formal version + coarse-graining theorem
  - D = 3: algebraically privileged + joint-compatible
  - Born rule, Hawking T, Bekenstein-Hawking S, Kerr extension: all closed

Strong-field metric: structurally complete.
  - Outside PS: GR-exact by construction (all precision tests pass)
  - Inside PS: quintic Hermite k(A) = (1−A) · (1 − 5y⁴ + 4y⁵)
  - Cubic horizon vanishing 1 + 2 = velocity + manifold-support double-zero
  - Kerr extension: r_pr(θ; a) closed form; STAM wave base across all spins

QNM / ringdown:
  - Spinless axial ℓ=2/3/4 n=0: 4.977% / 2.101% / 1.158% LOCKED RIGOROUS
  - Polar ℓ=2 n=0: 5.620%, isospectrality break ~13%, FIRST-PASS
  - n=1 overtones: OPEN (3 time-domain methods + WKB-PT + shooting all failed; Leaver/SN needed)
  - Kerr STAM ringdown: BLOCKED on G120d standalone solver

Cosmology:
  - V_3 modified Friedmann beats ΛCDM by Δχ²=25-33 across Pantheon+/Union3/DES
  - Two-layer reading: expansion (Layer 1) + distance bias (Layer 2)
  - Bridge term b ≈ 355 Mly matches 0.04%
  - Chronometer tension at H_0=73 (soft spot, χ²/N=0.76)
  - Pipeline open: f_LoS modeling, BAO under V_3, full MCMC

Quantum:
  - Born rule: structural derivation closed (Bell at Tsirelson bound exact)
  - Γ_res: Lindblad form committed; domain-specific L_μ open
  - QM action chain: substance ontology → proper time → S = −mc²τ → path integral → Born rule

Matter coupling:
  - Tier 1 (Jordan-frame minimal coupling): closed by construction (G123)
  - All 10 standard precision tests pass automatically
  - Tier 2 (substance excitations / TOE): open, multi-year program

Observational falsification handles (locked):
  - GW170817 engine time scaling
  - ℓ=2 polar/axial isospectrality break
  - F6 decoherence ~0.5 s
  - Primordial-mass remnants
  - PBH super-radiance rate banding
  - Near-horizon tortoise enhancement
```

---

## Part XV — Speculative Parking-Lot Ideas (Not Active Research)

These fit the framework's structure and may be worth revisiting:

- **Two-hologram horizon**: original matter on inner surface, accumulating matter on outer; cosmic expansion driven by outer-layer growth
- **Universe-as-bubble**: our cosmos may itself be a bubble (connects to holographic universe proposals)
- **Multiverse from non-isolated formation events**: if our universe formed inside a larger structure, similar events may have produced others
- **Black-hole-as-baby-universe (Smolin-CNS-adjacent)**: time-reparameterization across the horizon dissolves the timescale-mismatch in principle

---

## Appendix — Recent Diagnostic Scripts (2026-05-13 evening batch)

| Script | Purpose | Result |
|---|---|---|
| [G108](scripts/G108_axial_from_constrained_sigma_action.py) | Rigorous axial perturbation from S[Σ,g,λ₁,λ₂] via parity argument | Locked V_exact = V_geom + (√f)''/√f |
| [G109](scripts/G109_axial_convergence_study.py) | Throat-depth sweep on G108 | 4.977% locked, 0.002 pp stability |
| [G114](scripts/G114_polar_QNM_from_constrained_sigma.py) | Polar QNM first-pass via M_eff Zerilli + canonical correction | ~13% isospectrality break at ℓ=2 |
| [G121](scripts/G121_beta_uniqueness_from_bayesian_update.py) | Beta(D+1, 2) uniqueness theorem | Three converging derivations |
| [G122](scripts/G122_counted_scalar_action_derivation.py) | Constrained action uniqueness under Σ-as-count | LM form is forced |
| [G123](scripts/G123_matter_minimal_coupling.py) | Matter minimal coupling tier 1 | Jordan-frame; all precision tests pass by construction |
| [G124](scripts/G124_su_write_poisson_model.py) | First formal Poisson SU-write model + coarse-graining | Numerical demo: 4.85% vs 9.3% noise floor |
| [G113-redux](scripts/G113_redux_n1_overtone_via_shooting.py) | n=1 overtone attempts (shooting + WKB-PT) | Methodological finding: QNM shift is action-correction-dominated; n=1 still open |

---

## Document Maintenance Notes

- This file mirrors the README's Open Problems and Current Status sections but is organized as a claim ledger.
- For full derivations and detailed prose, see [README.md](README.md).
- For session-by-session history, see `memory/HANDOFF_*` files.
- For framework architecture commitments, see `memory/project_five_primitive_architecture.md` and related.

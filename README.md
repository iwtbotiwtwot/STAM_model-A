# Model-A

**Model-A** is the current variable-accumulation version of the broader **Spacetime Accumulation Model (STAM)** research program. It explores whether gravity, clock behavior, propagation delay, black-hole horizons, thermodynamic behavior, cosmological distance effects, dark matter behavior, and quantum-style resolution can be described using one central idea: a dimensionless spacetime accumulation field called **A**.

Author: **Sean Brady**
Status: **Proposed theoretical framework / active research program**
Snapshot: **May 12, 2026 — updated to add substance velocity-cap as V_4 candidate commitment (motion through elevated A is REALLY slowed, not just observed-slowed) and GW170817 engine time prediction (1.677 s predicted vs 1.74 s observed, 3.6% match without fitting). Builds on May 11 evening substance-ontology refinement, two-layer cosmological reading, and k(A) family ambiguity findings.**

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
g_rr = 1 / [ (1 - A) (1 - A^2)^2 ]       (Model-A modification, n=2 commitment)
```

**Three structural principles narrow k(A) to a 1-parameter family** (refined 2026-05-11):

1. **Weak-field GR recovery at first order in A** (so all GPS / Shapiro / lensing tests pass).
2. **A = 1 as the universe's edge** (so SU / proper-time integrals diverge there).
3. **Thirds-of-A preservation** (ISCO at A = 1/3, photon sphere at A = 2/3, horizon at A = 1; depends on g_tt only).

These three principles narrow k(A) to the family `k_n(A) = (1-A)(1-A^2)^n` for n ≥ 1, but do NOT uniquely fix n. The earlier claim that "three principles force k(A)" was overstated — they admit a family. Model-A commits to **n=2** specifically based on:

- G17's g_rr second-order coefficient ratio Model-A:GR = 3 (the "third integer 3" in the framework's recurring theme; equals n+1)
- G18's (2×2)=4 area-per-entry ledger decomposition matching Bekenstein-Hawking entropy (matches when there are two pair factors)
- The pair-squared form (1-A²)² as the framework's structural primitive

But n=1 (with k = (1-A)(1-A²)) is also viable: it gives NEC crossover at exactly A = 1/2, logarithmic horizon divergence (the original "logarithmic" wording), and QNM ratio closer to LIGO precision. Several G33/G34 framework signatures depend on this choice:

| signature | n=1 | n=2 (committed) |
|---|---|---|
| QNM eikonal τ/τ_GR | 3/√5 ≈ 1.342 | 9/5 = 1.800 |
| F3 NEC crossover | 1/2 exactly | ≈0.44 |
| Proper-time divergence | logarithmic | (1-A)^(-1/2) |
| g_rr 2nd-order ratio | 2 | 3 |

The 1-parameter ambiguity is real and is currently the framework's central open metric question (sharpening of the Lagrangian-for-A gap). A fixing point will come from either an exact LIGO QNM observation of a spinless ringdown, or a non-circular structural derivation that picks out a specific n.

The deviation factor (1 - A^2)^n is within ~1% of unity for A < 0.05 (everywhere we currently measure) regardless of n. All weak-field GR tests pass automatically. The framework departs from GR at second order in A and at the boundary itself; the size of the second-order departure (coefficient n+1 in g_rr) is what n parameterizes.

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

**Spinning black holes**: the bubble warps to an oblate shape, with equatorial radius preserved at 2M and polar radius shrinking to the Kerr horizon location. Matter on the bubble carries the angular momentum (the bubble is the screen, the matter is the rotating hologram). No Penrose extraction (since the rotational energy is in the matter, not in empty rotating geometry).

**Two-face refinement of the bubble surface (working commitment as of 2026-05-10).** The 2D boundary at A=1 is structurally distinguished into two faces: an inner face holding primordial mass from the formation event (static, doesn't rotate) and an outer face holding subsequently accreted matter (dynamic, can rotate, carries angular momentum). The "horizon doesn't spin, hologram does" principle applies specifically to the outer face's holographic content. PBHs are "thin outer face" objects (formed primordially with little subsequent accretion), which predicts evaporation signatures differing from stellar BHs both in standard ways (low spin) and in framework-specific ways (sparse outer face content). The two-face refinement gives the no-interior commitment specific physical mechanism: information that "fell" into the BH lives on the inner face permanently; information from accretion lives on the outer face and gets emitted as Hawking radiation. The information paradox dissolves with a definite mechanism rather than as a slogan.


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

For spinning bubbles, T(theta) is non-uniform along the bubble: equator at T_Schwarzschild for any spin, pole at T_Kerr (matching standard Kerr horizon T exactly via identity). Total Hawking emission stays substantial even near extremal, with super-radiance enhancement at the equator.

**Bekenstein-Hawking entropy from direct ledger counting (G18).** Beyond reproducing `S = k_B A / (4 ell_P^2)` through standard Q-series machinery, the framework reaches the same result by counting independent ledger entries on the boundary, without going through `dE = T dS`. Each independent ledger entry occupies area `alpha ell_P^2` where `alpha = (two-face structure: 2) x (gravity-bridge factor: 2) = 4`. Both factors trace to framework primitives (the two-face refinement is a working commitment; the gravity-bridge factor 2 is in A's definition). This is a parallel structural derivation route, with the honest caveat that the (2 x 2 = 4) decomposition uses framework primitives in a target-consistent way rather than being uniquely forced — a skeptical reader could note that knowing the answer 1/4 in advance shaped which framework factors got combined. What's unambiguous: the framework's no-interior + two-face + gravity-bridge commitments together produce the right entropy formula from cell-counting; the 1/4 factor receives a structural reading (4 Planck areas per ledger entry) rather than only a path-integral arithmetic origin.

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
G1: Ringdown damping factor tau / tau_GR = (9/5)^(n/2) at same dominant
    frequency as GR. For Model-A's committed n=2: tau/tau_GR = 1.80.
    For alternative n=1: tau/tau_GR ≈ 1.342. Testable now with
    LIGO/Virgo/KAGRA ringdown data (eikonal approximation; LIGO BHs
    are spinning Kerr-like, so spinless Model-A analog is the missing
    piece for direct comparison).

G43: BNS merger engine time (orbital collision → gamma-ray burst onset)
     scales linearly with total binary mass:
       τ_engine = τ_critical ≈ 0.62 × (M_total / M_sun) seconds
     For GW170817 (M = 2.7 M_sun): predicted 1.68 s, observed 1.74 s
     (3.6% match, no fitting). Falsifiable with future BNS+EM events:
     engine times should fall on a linear trend with this slope.
     If they don't, the substance velocity-cap + elastic-shell mechanism
     is refuted as the engine-time mechanism.

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

5. **Strong-field metric is Model-A-narrowed.** Three Model-A principles (weak-field recovery + boundary at A = 1 + thirds-of-A preservation) narrow k(A) to the 1-parameter family `k_n = (1-A)(1-A^2)^n` for n ≥ 1. Model-A's specific commitment to n=2 (giving the (1-A^2)^2 modification) is structurally motivated but admits a viable alternative at n=1; see "Strong-field position" section. The framework's strong-field metric is shaped by its own commitments rather than borrowed from elsewhere, but the n choice is currently a 1-parameter ambiguity.

6. **Spinning bubble framework.** Static bubble + rotating holographic matter gives non-uniform Hawking T (T_pole = T_Kerr exactly via identity), super-radiance enhancement, no Penrose extraction.

7. **SN distance fits beat LCDM combined chi^2.** Model-A with V_3 modified Friedmann wins by 25-33 across Pantheon+/Union3/DES at the same number of free parameters. Predicts inter-catalog Pantheon+/Union3 tension within 27%.

8. **CMB self-consistent closure at H_0 = 73.** Cumulative-A line-of-sight amplification structurally explains the H_0 tension. Internal solution exists for realistic cosmic-structure parameters.

9. **PBH-DM compatibility for galactic dark matter.** Each PBH is a small bubble with the framework's existing thermodynamic/structural machinery. Galactic rotation curves close trivially with PBH-halo + cumulative A.

10. **Six observational regimes, one A field.** Local gravity, propagation delay, BH thermodynamics, SN distances, CMB acoustic scale, and galactic DM all from the same A field with one structural constant (A_0 = 1/(12 pi)) and one calibrated parameter (beta).

11. **Specific falsifiable predictions** (G1 ringdown, F6 decoherence, G3-G4 latitudinal Hawking, PBH-DM sigma constraint, G43 BNS engine time mass-scaling) — distinguishable from LCDM in regimes the standard framework does not address.

12. **Substance velocity-cap as STAM-vs-GR wedge** (added 2026-05-12). Substance ontology forces motion through elevated A to be REALLY slowed (not just observationally). For GW170817, this gives binary merger engine time = τ_critical = 1.677 s vs observed 1.74 s — 3.6% match, no fitting. The progression of substance correction (naive → partial → full) is monotonic toward observation. Mass-scaling linear in M_total is the falsification handle.

---

## Major open problems

1. **Lagrangian for A; or non-Lagrangian principle that picks out k_n specifically.** Specifying the action `S[A, g_uv]` (or alternative structural principle) that produces both the Poisson-like source equation and a unique strong-field metric. G25-G27 attempted standard scalar-tensor formulations and found ghost regions plus Rs-dependent V; G33 inverted to A-as-fundamental and found the k(A) family is 1-parameter ambiguous (k_n for n ≥ 1 satisfying all stated principles). The Lagrangian gap is now the sharper question: what principle picks n out of {1, 2, 3, ...}?

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

Development assisted by extensive conversation with Claude (Anthropic), particularly during the May 2026 push that produced the A_0 = 1/(12 pi) structural commitment, V_3 potential selection, CMB self-consistent closure, and PBH-DM compatibility analysis (G7-G14 script series); the May 11 (evening) refinement session that produced the substance-ontology articulation (water-tank conceptual exercise), the two-layer cosmological reading (G28+G29: V_3 expansion ≡ LCDM at H_0=73 by construction, photon-A traversal as separate distance bias), the G21 redux confirming GR-exact strong-field landmarks under corrected composition, the G31 confirmation that A_0 baseline does not bridge galactic DM, and the G33/G34 swing surfacing the 1-parameter ambiguity in k(A); and the May 12 session that produced the substance velocity-cap candidate commitment (motion through elevated A is REALLY slowed) and the GW170817 engine time prediction matching observation to 3.6% without fitting (G42, G43 scripts).

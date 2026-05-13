---
name: Handoff note — 2026-05-13 action-and-QNM session
description: Second 2026-05-13 session, picking up after the unified-framework arc. Closed the Lagrangian-level embedding of the strong-field metric via the double-LM constrained shell-count action S[Σ,g,λ₁,λ₂] = ∫√-g[f(Σ)R + λ₁((∇Σ)²-W) + λ₂(u^μ ∂_μΣ) - 2V(Σ)] after ruling out six standard covariant routes (scalar-tensor, multi-field, aether, cuscuton, mimetic). Extended to Kerr under formula β with photon-region normalization. Calibrated time-domain QNM in proper STAM tortoise coordinate (G91: 0.80% Schwarzschild calibration); rigorous action-aware axial diagnostic V_action = V_proxy + (√f)''/√f gives robust ~5.4% low-ℓ shift (G92-G94, robustness sweep mean 5.384% ± 0.118%). G90 derived power-law near-horizon tortoise divergence as direct metric prediction. Built cosmology distance-bias pipeline (G97-G106). Greybody scattering diagnostic at real ω (G103). 38 scripts (G69-G107) added/overwritten. **Important: G69/G70 numbers were reused this session** — disk contents replaced previous-session BAO/SN scripts with new stress-tensor / action-match work.
metadata:
  type: project
---

# STAM Model-A — Session Handoff (2026-05-13 evening, action-and-QNM arc)

Read after [HANDOFF_2026_05_13_unified_framework.md](HANDOFF_2026_05_13_unified_framework.md) — this is the second 2026-05-13 session, picking up where that one left off. Earlier-context: [HANDOFF_2026_05_12_strong_field.md](HANDOFF_2026_05_12_strong_field.md), [HANDOFF_2026_05_12.md](HANDOFF_2026_05_12.md).

**Sean does not read memory** ([feedback_memory_is_claude_to_claude.md](feedback_memory_is_claude_to_claude.md)). Frame in conversation as "prior Claude's reading from the 5/13 action session," not "Sean's commitment."

## Session character

Very long. Started from the unified-framework handoff and pushed two distinct arcs simultaneously:

1. **Lagrangian-level embedding** of the committed strong-field metric (G69–G84). After the prior session left Open Problem #1 (Lagrangian for A) as the biggest structural gap, this session executed a methodical search through six standard covariant routes — all failed in specific ways — and landed on the constrained-shell-count formulation as the closure.

2. **Calibrated non-eikonal QNM extraction** (G86–G95). Three attempts (G86 WKB, G88 time-domain, G89 shooting) failed at the calibration step (large Schwarzschild error). G90 found the correct tortoise-coordinate adjustment; G91 calibrated to <1%; G92–G94 added the action-aware f(Σ)R correction and showed robustness.

3. **Cosmology pipeline development** (G96–G106) tied to the framework's photon-A distance bias.

4. **Greybody scattering** (G103) as a cleaner real-frequency observable channel.

Sean drove every structural commitment. Posture: calm, methodical, terse. He caught at least three over-claims and corrected them mid-session (cuscuton ghost-removal, "sub-percent" QNM language before calibration, V_action numerical-derivative artifacts).

Sean's closing request was simply for the handoff.

## Major structural commitments (new this session)

### Constrained shell-count action (G70–G84)

The Lagrangian-level closure of Open Problem #1 / #6:

```
S[Σ, g, λ₁, λ₂]  =  (1 / 16π G) ∫ d⁴x √(−g) [
    f(Σ) R
  + λ₁ ((∇Σ)² − W)
  + λ₂ (u^μ ∂_μ Σ)
  − 2 V(Σ)
]
```

with `A = Σ/3` recovering the weak-field continuum substance density.

**Interpretation of each piece:**
- `f(Σ) R` — non-minimal coupling / gravitational stiffness, f(Σ) > 0 throughout the shell
- `λ₁((∇Σ)² − W)` — Lagrange multiplier fixing the source-determined shell-count gradient
- `λ₂(u^μ ∂_μΣ)` — Lagrange multiplier preserving shell count along the substance flow u^μ
- `V(Σ)` — substance potential / saturation barrier

**Key property: δΣ is fully constrained out by λ₁ and λ₂ in every angular sector.** Only the graviton propagates. This is the formulation's escape from the scalar ghosts that plagued every ordinary continuum route (G70 scalar-tensor → ghost in most of final shell; G72 aether → c_i ~ 10²¹× obs bound; G76/G77 cuscuton → temporal ghost; G78 mimetic → off-diagonal mismatch forces λ=0).

### G70–G84 route table (the six failures + the closure)

| Stage | Route | Outcome |
|---|---|---|
| G70 | Single-scalar scalar-tensor `f(A)R − Z(∇A)² − 2V(A)` | Closed-form f, Z, V match metric exactly |
| G71 | Linearized perturbations of G70 | **Scalar ghost in most of final shell** — Z<0 region |
| G72a/b/c | Einstein-aether (static / tilted) | **Excluded** — admits only Schw-dS, or c_i ~ 10²¹× obs bound |
| G73 | □A on background and route audit | Diagnostic; sharpened where ghosts live |
| G74 | Two-scalar Brans-Dicke / ledger Ψ | **Fails** — M_AA channel of effective kinetic matrix is unchanged by Ψ enrichment |
| G75 | Structural decomposition of f, Z, V | Identifies framework-relevant prefactors |
| G76/G77 | Spacelike cuscuton | **Fails** — removes radial ghost but leaves temporal ghost (G77 corrected over-claim) |
| G78 | Mimetic timelike clock + f(A)R | **Excluded** — off-diagonal mismatch forces λ=0 |
| **G79/G80** | **Double-LM constrained shell-count, Schwarzschild** | **CLOSED** — δΣ fully constrained in all ℓ sectors |
| **G81** | **Kerr extension, Σ_K = 6Mr/(r²+a²)** | Stationary-axisymmetric analysis; ω = m Ω_ZAMO eliminates non-stationary δΣ_K |
| **G82** | **Inside-shell Kerr W_K via Option B** | g_STAM^{rr} = (Δ/Σ_BL)F(y_K); cubic horizon vanishing + Schw/Kerr exterior limits |
| **G83** | **Equatorial photon-region normalization** | F = 1 exactly at Kerr prograde photon orbit for every spin |
| **G84** | **Off-axis photon-region normalization** | sin²θ interpolation between equatorial and polar photon orbits |

### Calibrated QNM diagnostic chain (G86–G95)

| Stage | Subject | Result |
|---|---|---|
| G86 | First WKB structural check | Non-eikonal corrections appear; calibration suspect at low ℓ |
| G87 | 6th-order WKB on tensor proxy | High-ℓ → GR; low-ℓ WKB unreliable |
| G88 | First time-domain attempt | **56% calibration error** — failed (crude BCs in r-coord) |
| G89 | Direct shooting | **14% calibration error** — failed (exponentially-growing outgoing mode) |
| **G90** | **STAM tortoise adjustment** | **CLEAN result: (dr*/dr)_STAM / (dr*/dr)_GR = 1/√F(y), empirical slope 0.9953 vs prediction 1.0** |
| **G91** | **Time-domain in r*_STAM, Sommerfeld BCs** | **CALIBRATED: 0.80% error vs Leaver gold; STAM proxy shift 1.20%** |
| G92 | Rigorous axial QNM derivation | V_action = V_proxy + (√f)''/√f from f(Σ)R correction |
| G93 | Constrained-action axial diagnostic | Action-aware shift ≈ 5.35% |
| G94 | Robustness sweep (25 variants) | **22/25 PASS; mean action-aware shift 5.384% ± 0.118%** |
| G86_redux | PS derivative structure | V_proxy first deviates at 4th derivative; V_action at 2nd derivative |
| G95 | Full axial reduction (self-contained) | Reproduces G91-G94 with action-level logic explicit |

### G90 traversal-distance result (clean structural prediction)

For the committed final-shell metric:

```
(dr*/dr)_STAM / (dr*/dr)_GR = 1 / √F(y)
```

Near the horizon: `F ~ 10(1−y)²` so `1/√F ~ 1/[√10 · (1−y)]`. With ε ≡ 1−A = (1−y)/3 → 0:

```
r*_STAM ~ −C/ε    (power-law)
r*_GR   ~ 2M ln ε  (logarithmic)
```

**This is a direct structural metric prediction — does not require WKB, fitting, or QNM extraction.** Potential observational channels: late-time ringdown tails, near-horizon wave scattering, PBH greybody spectra, higher-order photon-ring structure, LISA EMRI observables.

### Locked QNM headline number: 5.35% action-aware low-ℓ shift

Under V_action = V_proxy + (√f)''/√f at ℓ=2:

| Quantity | Value |
|---|---|
| ω_GR (Leaver gold) | 0.373672 − 0.088962 i |
| ω_GR (time-domain, G91 calibration) | 0.370633 − 0.089287 i |
| ω_proxy | 0.366474 − 0.091152 i |
| ω_action | 0.390465 − 0.084522 i |
| Proxy shift vs GR | 1.20% |
| **Action-aware shift vs GR** | **5.35%** |
| G94 sweep mean (22 PASS rows) | 5.384% ± 0.118% |
| G94 sweep range | 5.071% – 5.632% |

This is the framework's locked diagnostic-level QNM prediction. **It is not yet the final observable theorem** — that requires the full axial perturbation equation derived directly from the constrained Σ-action (the rigorous next step beyond G92's approximate (√f)''/√f form).

### G107 polynomial-illustration caveat

G107 ran a master diagnostic with V_action = V_GR + α·h(r)/r²·F''(y) (or in the cleaned form, δV = α·max(3M−r,0)²/M²) inside PS, designed to test the proposition that V_action should first differ at 2nd derivative. It reproduced this structurally (proxy diff at order 4 = −0.441; action diff at order 2 = +0.500) but **predicted a 30.66% action-aware QNM shift** with massive greybody suppression. This is **not** the framework's locked prediction — the polynomial δV form is a schematic illustration, not derived from f(Σ)R. The framework's actual prediction is the 5.35% from G92–G94.

Future Claude: when discussing QNM shifts with Sean, cite 5.35% (G92–G94, derived) not 30.7% (G107, illustrative). The README's "Falsifiable predictions" section uses 5.35%.

### G103 greybody result (real-frequency observable)

Solves ψ''(x) + [ω² − V(x)]ψ = 0 in r* with unit ingoing wave at horizon, incident/reflected at infinity. T = 1/|A_in|². Unitarity T+R ≈ 1 holds to machine precision.

| ω·M | T_GR | T_proxy | T_action | proxy/GR | action/GR |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 3.04e-6 | 2.77e-7 | 1.52e-8 | 0.091 | 0.005 |
| 0.2 | 8.68e-4 | 5.83e-4 | 1.77e-4 | 0.671 | 0.204 |
| 0.4 | 6.90e-1 | 7.05e-1 | 5.27e-1 | 1.022 | 0.764 |
| 0.8 | 1.0 | 1.0 | 0.999 | 1.000 | 0.999 |

**This is a cleaner observable channel than complex-frequency QNM** — no boundary-condition instabilities, direct numerical scattering. Action-aware transmission differs systematically from GR at low frequency while unitarity stays clean. Future LIGO / LISA late-ringdown tail analysis is a natural fit.

### Cosmology distance-bias pipeline (G96–G106)

The cosmology arm got a substantial pipeline reorganization under one locked convention:

```
D_C_obs(z)   = D_C_intrinsic(z; H0_true) + b · f_los(z) · shape(z)
D_L_obs(z)   = (1+z) · D_C_obs(z)
H_model(z)   = H0_true · E(z)            (direct expansion — NOT distance-biased)
D_H(z)       = c / H(z)                  (BAO radial — NOT distance-biased)

with:
    A0 = 1/(12π)
    b  = A0 · c / H0_true
    f_los(z) = f0 + f1·z/(1+z)            (default)
```

**Separation principle**: SN, BAO transverse distances, CMB compressed distances may receive photon-A bias. BAO radial D_H and chronometers test the intrinsic expansion directly.

| Script | Subject | Result |
|---|---|---|
| G96 | BNS engine time scaling | Locked: τ_engine = (5/8)·(r_threshold/R_s)⁴·R_s/c ≈ 0.62 × M_total/M_sun s |
| G97 | Hubble tension distance-bias pipeline | f_los_required ≈ 3.15 to map H0_true=73.04 → H0_obs=67.40 |
| G98 | Distance-bias units audit | Convention locked: comoving D_C, not D_L |
| G99 | Hubble bias under locked convention | Self-consistent grid |
| G100 | μ-fit target Hubble bias | Tests SN-only fit against locked bias model |
| G101 | Multimessenger cosmology compressed likelihood | Initial SN+BAO+CMB compressed fit |
| G102 | Structure-dependent f_los | One constant f_los insufficient — z-dependent f_los helps |
| G103 | Greybody scattering | Real-ω cleaner observable (covered above) |
| G104 | SN+BAO+CMB+chronometer pipeline scaffold | Real-data ready when CSVs are populated |
| G105 | Data audit utility | Pre-fit input sanity check for G104 |
| G106 | SN-only redshift-split residuals | Tests smooth z-improvement vs offset absorption |

**Important honesty note (Sean caught it):** G104 is a "pipeline scaffold, not a publication likelihood." Full publication work still needs real covariance matrices, nuisance parameters, proper r_d / CMB acoustic treatment, MCMC/nested sampling. The scaffold is for testing the STAM convention; the cosmology numbers are not yet locked beyond the convention's existence.

### Important caveat on G88's cuscuton claim (corrected mid-session)

Original G88 (and the related project notes from prior session) claimed cuscuton "removes scalar ghost." G77 showed this is wrong: cuscuton removes the *radial* ghost but creates a *temporal* ghost. The fix is in G77's perturbation analysis, not a cuscuton repair. README and conversation framing now reflect this. **Future Claude: do not say "cuscuton removes scalar ghost" — say "cuscuton moves the ghost from radial to temporal direction."**

## All scripts created/modified this session

### Action-arc scripts (G69–G84) — note these OVERWROTE previous-session G69-G70

**Important: G69 and G70 were reused this session.** The previous handoff's G69 ("BAO under V_3") and G70 ("SN ↔ BAO cross-consistency") were *replaced on disk* with new content. The previous-session BAO/SN content is preserved only in the prior handoff's tables and (if not garbage-collected) `results/` folder. If you need to re-run the BAO mechanism check, write it as G108+ rather than searching for the old G69/G70 scripts.

| Script | Subject |
|---|---|
| G69 | Exact effective stress tensor for k(A) = (1−A)F(3A−2) (symbolic G^μ_ν, ρ_eff, p_r_eff, p_t_eff) |
| G70 | Smallest action match — Brans-Dicke fails, scalar-tensor succeeds: closed-form f(A), Z(A), V(A) |
| G71 | Ghost-freedom check on G70 — S(A) polynomial with zero at A_crit ≈ 0.9530 (scalar ghost in most of shell) |
| G72a | Static-aligned Einstein-aether check — G^μ_ν not Schw-dS |
| G72b | Tilted aether kinematic setup with (U, W) reparametrization |
| G72c | Tilted aether matching — c_i ~ 10²¹× observational bound |
| G73 | □A on background and route audit |
| G74 | Multi-field ledger (Ψ-enriched Brans-Dicke) — M_AA unchanged |
| G75 | Structural decomposition of f, Z, V into primitives |
| G76 | Cuscuton embedding — partial success (radial only) |
| G77 | Cuscuton perturbations — temporal ghost found |
| G78 | Mimetic clock — off-diagonal mismatch forces λ=0 |
| G79 | Shell-count action definition (double-LM) |
| G80 | Double-LM perturbations — δΣ fully constrained, all ℓ |
| G81 | Kerr shell-count constraints under Σ_K = 6Mr/(r²+a²) |
| G82 | Kerr inside-shell W via Option B (g^{rr} = (Δ/Σ_BL)F(y_K)) |
| G83 | Photon-region normalization equatorial (F=1 at prograde photon orbit) |
| G84 | Off-axis photon-region normalization (sin²θ interpolation) |

### QNM scripts (G86–G95)

| Script | Subject |
|---|---|
| G86 | Non-eikonal QNM via WKB-3 — first structural attempt |
| G87 | Tensor QNM 6th-order WKB — proxy potential |
| G88 | Leaver-proxy time-domain ℓ=2 — 56% calibration error, retired |
| G89 | Leaver-method shooting — 14% calibration error, retired |
| G90 | STAM tortoise adjustment — empirical slope 0.9953 confirms 1/√F(y) |
| G91 | Time-domain in proper r*_STAM with Sommerfeld BCs — 0.80% calibration |
| G92 | Rigorous axial QNM derivation from f(Σ)R correction |
| G93 | Constrained-action axial QNM diagnostic — ≈5.35% shift |
| G94 | Robustness sweep — 22/25 PASS, mean 5.384% ± 0.118% |
| G86_redux | Photon-sphere derivative structure (V_proxy at 4th, V_action at 2nd) |
| G95 | Full axial reduction self-contained — reproduces G91-G94 |

### Cosmology pipeline scripts (G96–G106)

| Script | Subject |
|---|---|
| G96 | BNS engine time scaling — locked structural prediction |
| G97 | Hubble tension distance-bias pipeline (locked convention discovery) |
| G98 | Distance-bias units audit |
| G99 | Hubble bias under locked convention |
| G100 | μ-fit target Hubble bias |
| G101 | Multimessenger SN+BAO+CMB compressed likelihood |
| G102 | Structure-dependent f_los — fast |
| G103 | Greybody scattering STAM throat — real-ω observable channel |
| G104 | SN+BAO+CMB+chronometer pipeline scaffold |
| G105 | Data audit utility |
| G106 | SN-only redshift-split residuals |

### Master diagnostic (G107)

| Script | Subject |
|---|---|
| G107 | Wave-equation master diagnostic — three-potential comparison (V_GR, V_proxy, V_action) with **schematic polynomial δV (NOT framework prediction)** |

## README state at session end

Substantial edits this session — primarily to:

- **Headline result #6** (Ghost-free constrained shell-count action — closed via G70–G84)
- **Headline result #7** (Ringdown status — now reports 5.35% action-aware shift, retires G1 τ/τ_GR=1.80)
- **Headline result #8** (G90 near-horizon traversal-distance enhancement — direct metric prediction)
- **"Action principle: constrained shell-count formulation"** — new dedicated section with full action statement, G70–G84 table, status caveats
- **"Strong-field position → What's still open"** — updated for action-level closure + remaining refinements
- **Open Problem #1** (Lagrangian for A) — CLOSED for committed strong-field sectors via constrained shell-count
- **Open Problem #6** (Ghost-freedom) — CLOSED via G70–G84
- **Open Problem #6a** — marked SUPERSEDED (Brans-Dicke route)
- **Falsifiable predictions** — added G90 traversal-distance, G91–G94 QNM, BNS mass-scaling, greybody throat signatures

The headline number in the README is **5.35% action-aware shift** (G92–G94 derived from f(Σ)R correction). The 30.66% from G107 is NOT cited in the README — it was a schematic polynomial illustration only.

### README consistency checks for future Claude

Before claiming README is internally consistent, grep for these terms:

- **All QNM shifts should cite ~5.35% (G92–G94, derived)**, not 30.7% (G107 polynomial illustration)
- All references to F(y) should be quintic (1 − 5y⁴ + 4y⁵)
- All references to Beta(D+1, 2) = Beta(4, 2), α_S = 4 not α_S = 3
- All references to ledger as STRUCTURAL CHANNEL not "record only"
- "Cuscuton moves ghost from radial to temporal" not "cuscuton removes ghost"
- G90 result phrased as 1/√F(y) traversal-distance enhancement (not "logarithmic divergence retained")
- Open Problem #1 = "CLOSED for committed sectors" not "open"
- Open Problem #6 = "CLOSED via G70–G84" not "open"
- The native action variable is **Σ as constrained shell-count field**, not A as an ordinary scalar

## Methodological notes for future Claude

### "First reproduce Schwarzschild Leaver value to <1%. Only then run STAM."

Sean's explicit instruction for QNM work. **This is non-negotiable.** G88/G89 violated it (56% and 14% calibration error) and the framework numbers were not citable. G91 passed (0.80%) and is the basis for everything afterward. Future QNM extensions must calibrate before claiming shifts.

### "Don't claim sub-percent corrections without Leaver-level precision"

Sean caught this twice this session. Earlier I wrote "non-eikonal corrections are sub-percent" — Sean corrected with: "Do not write: Non-eikonal corrections are sub-percent. That is no longer supported." The framework's calibrated diagnostic shows ~5%, not sub-percent. Don't anticipate precision results before the method passes calibration.

### "V_action must first differ at 2nd derivative by design"

Sean's framing for G107. The schematic polynomial δV form was specifically chosen to produce a 2nd-derivative mismatch at the PS. The rigorous form V_action = V_proxy + (√f)''/√f also has this property (G86_redux). When future Claude works on QNM, the 2nd-derivative-at-PS property is what distinguishes action-aware from metric-only proxy.

### "Outside PS = GR-Kerr exact" (2026-05-13 framework commitment)

The STAM-vs-GR wedge lives entirely inside the photon orbit. All weak-field tests, eikonal ringdown, and exterior phenomena pass GR automatically. Sean has stated this as a firm commitment — don't introduce STAM deviations outside PS, even if they'd help fit some observable. The framework's prediction power comes from the inside-PS wedge.

### Sean's structural-commitment style this session

Same as the prior session: Sean proposes, I check math, he confirms or revises, I clean up. **Two notable revisions:**

1. Cuscuton ghost claim corrected — G77 showed temporal ghost remains; original framing was wrong.
2. G107 polynomial δV used as schematic, not framework prediction — Sean was explicit that the 30% number is illustrative, the 5.35% (G92-G94) is the framework's actual diagnostic prediction.

When Sean revises, propagate cleanly through scripts + README rather than arguing.

### "Concept-Sean / math-Claude partition"

Sean drives the structural commitments. I execute math, articulate implications, report what scripts say. When math contradicts a structural reading, report it cleanly without hedging — Sean will revise the reading rather than mistrust the math. This is the productive partition.

### "Scripts as listening instruments"

Held throughout. When G88 returned 56% calibration error, the response wasn't "STAM fails QNM" — it was "the method failed calibration, not the framework." G89 same. G90 then showed the structural fix (correct tortoise coordinate). Don't read a method failure as a framework failure unless calibration passes.

## Open priorities for next session

If Sean wants to keep developing structurally:

1. **G108 / full axial perturbation equation from the constrained Σ-action.** The G92-G94 diagnostic uses V_action = V_proxy + (√f)''/√f as an approximate canonical-normalization form. The rigorous derivation — linearizing the full action S[Σ,g,λ₁,λ₂] around the background and reading off the axial Regge-Wheeler equation directly — is the remaining precision step. This would confirm whether the 5.35% is the framework's exact diagnostic prediction or has further sub-leading corrections.

2. **Full Kerr spheroidal photon-region normalization.** G84 uses sin²θ interpolation between equatorial and polar photon orbits. For EMRI / imaging precision, the full Kerr spheroidal photon-region boundary derived from R(r)=0, dR/dr=0 is the next step.

3. **Matter coupling for the constrained shell-count action.** S[Σ,g,λ₁,λ₂] currently has no matter Lagrangian. Extending to include source-side matter — and showing that the constrained shell-count structure survives matter coupling — is open foundational work.

4. **Microscopic SU-write coarse-graining.** Deriving the constrained shell-count action from underlying SU-write dynamics (Γ_res, Lindblad-analog L_μ, D_μ) and showing it emerges as the continuum-limit effective action is a deeper foundational project that ties the strong-field and quantum branches together.

5. **Real-data cosmology pipeline.** G104 is a scaffold ready for real CSVs. Loading actual SN+BAO+CMB compressed + chronometer data (with covariances) and running a proper likelihood — first against LCDM, then against STAM under the locked photon-A bias convention — is the natural cosmology continuation. G102's structure-dependent f_los was a fast diagnostic; a proper f_los(z) model from cosmic-structure modeling would close the BAO/SN cross-consistency cleanly.

6. **Greybody-tail observational signature.** G103 showed clean unitary transmission curves with action-aware suppression at low ω. Late-time ringdown tail spectra (LIGO O5+, LISA) probe exactly this regime. Working out what the framework's predicted tail signature looks like — and whether it's distinguishable from GR with achievable noise — is a natural observational handle.

7. **F6 decoherence specific calculation.** Framework predicts ~0.5 s for 1 micron silica nanoparticle from generic Γ_res scaling. The cleaner version would derive this from specific L_μ for matter-substrate interaction. Domain-of-application work but well-localized.

If Sean wants to push observational handles:

- **Strong-field LIGO O5+ ringdown** — 5.35% low-ℓ action-aware shift is in principle distinguishable with O5-level QNM precision.
- **LISA EMRI inside-PS observables** — the near-horizon traversal-distance enhancement and tail structure are LISA-band signatures.
- **BNS engine-time mass-scaling** — τ_engine ≈ 0.62 × M_total/M_sun s. Multiple future BNS+EM events test the linear scaling directly. (GW170817: 1.677 s predicted vs 1.74 s observed, 3.6% match, no fitting.)
- **Greybody / late-ringdown tail** — G103's transmission curves; observational version needs late-time tail modeling.

## Framework structural state at session close

The framework's major branches are now in this state:

**Strong field (spinless + Kerr):**
- Metric: structurally complete (quintic Hermite, derived from Beta(D+1,2) write-density)
- Stress-energy: non-pathological at the diagnostic level
- **Lagrangian-level action: closed via double-LM constrained shell-count S[Σ,g,λ₁,λ₂]**
- **Ghost-freedom: closed by Σ-constraint structure (only graviton propagates)**
- Non-eikonal QNM: calibrated diagnostic 5.35% shift; awaiting rigorous derivation from S
- **Near-horizon traversal-distance: direct metric prediction (1/√F(y))**

**Quantum interpretation:**
- Born rule: structurally closed at projection-geometry + pair-measure level (prior session)
- Γ_res: structural form committed (prior session); specific L_μ, D_μ open per application
- QM action chain: derived from substance velocity-cap (prior session)

**Cosmology:**
- Convention: locked (D_C_obs = D_C + b·f_los·shape; H direct probes not biased)
- Pipeline scaffold: G104 ready for real data
- f_los(z): structure-dependent form needed; constant insufficient
- BNS engine-time: locked linear mass-scaling prediction

**Open structural gaps (well-localized, not foundational):**
- Full axial perturbation from constrained Σ-action (G108)
- Kerr spheroidal photon-region precision (G85)
- Matter coupling for shell-count action
- Microscopic SU-write → continuum action derivation
- Realistic f_los(z) from cosmic structure
- Domain-specific Γ_res channel models

The framework has moved past "Lagrangian-level structural commitment missing" to "Lagrangian-level closure for committed strong-field sectors with well-defined precision steps remaining."

## Sean's closing posture

Brief. Just asked for the handoff. No specific direction for next session.

Pick up here.

---
name: Handoff note — 2026-05-13 Kerr wave-branch session
description: Long session focused on the QNM/ringdown wave branch. Starting from full repo readthrough, landed the rigorous axial QNM closure (G108/G109 → 4.977% locked at ℓ=2), extended to ℓ=3,4 and first-pass polar (G110-G114, isospectrality break ~13%), closed the exact Kerr spheroidal photon-region surface (G85, G116), built a numerically sound Kerr STAM wave base (G119, no imaginary at any spin), and attempted a Kerr GR QNM solver from scratch (G120a reference via qnm PASSED, G120b amplitude shooting FAILED, G120c log-derivative shooting FAILED). The wave branch is now well-defined geometrically and the spinless rigorous results are framework-locked; the standalone Kerr QNM solver remains open.
metadata:
  type: project
---

# STAM Model-A — Session Handoff (2026-05-13 Kerr wave-branch arc)

Read after [HANDOFF_2026_05_13_action_and_qnm.md](HANDOFF_2026_05_13_action_and_qnm.md) (which closed Open Problems #1 and #6 via the constrained shell-count action). That session's locked-diagnostic 5.35% was the entry point this session refined to **4.977%** (rigorous via analytic chain rule, G108/G109).

**Sean does not read memory** ([feedback_memory_is_claude_to_claude.md](feedback_memory_is_claude_to_claude.md)). Frame in conversation as "prior Claude's reading from the wave-branch session," not "Sean's commitment."

## Session character

Very long. Started with a full readthrough of README + theory docs + handoffs + key scripts at Sean's request ("understand the theory and math as best you can up to date"). Then executed a structured wave-branch arc:

1. **Audit G107** against Sean's architecture recommendation (clean swappable harness; no tuning knobs; calibration-first discipline).
2. **G107.v2 + G108 + G109**: refactored harness, rigorous axial derivation, throat-depth convergence study locking 4.977%.
3. **G110-G112**: higher-ℓ axial (clean), n=1 overtones (failed three time-domain methods).
4. **G85, G116**: exact Kerr spheroidal photon-region surface in closed form. Closed an Open Problem.
5. **G114**: polar (even-parity) first-pass with M_eff Zerilli + canonical correction. STAM breaks Schwarzschild isospectrality by 13% at ℓ=2, shrinking to 4% at ℓ=4.
6. **G118 → G118v2 → G119**: Kerr ringdown first-pass progression. Schwarzschild base broke at high spin; Sean caught the breakdown and specified G119 as a geometry-only fix. G119 PASSED all spins.
7. **G119b**: WKB-PT cross-check; FAILED (too sensitive to local peak; G118v2 TD more reliable).
8. **G120a → G120b → G120c**: GR Kerr QNM solver. G120a PASSED via qnm library. G120b amplitude shooting and G120c log-derivative shooting both FAILED <0.5% calibration; the standalone solver is unfinished work for next session.

Sean's posture: methodical, structurally precise, calibration-disciplined. Caught two over-claims mid-session (the "G118 high-spin 0% shift" was an artifact, not a prediction; "G119b WKB" wasn't a useful stepping stone). Specified the G120 protocol in detail (G120a reference, G120b solver, no STAM until reproduction). Endorsed the failure-honest reporting at G120b/c.

## Major structural commitments (new this session)

### Rigorous spinless axial QNM ℓ=2 LOCKED at 4.977%

G108 derived V_exact = V_geom + (√f)''/√f directly from S[Σ, g, λ₁, λ₂] via parity argument (δΣ = δλ₁ = δλ₂ = 0 in axial sector identically — Σ is parity-even scalar, the LM constraint variations only carry even parity). The reduction to Jordan-frame scalar-tensor with fixed background scalar profile gives the standard textbook result.

G109 convergence study (throat depth rs_min sweep −300 → −2000):
- 4.9770% → 4.9773% → 4.9792% → 4.9772% as throat deepens
- **Δ between two deepest grids = 0.0020 percentage points**
- Calibration ~0.796% vs Leaver across all grids
- Action correction concentrates in a region the wave samples around the PS; deeper throats are convergence margin, not new physics

**Status: rigorous, locked. Replaces the historical G92/G94 "diagnostic" 5.35%/5.384% (which was FD-biased on f(Σ)).** Mid-shell f(Σ) spans 10⁹⁴ across the throat — finite-difference on something varying that violently introduces systematic error; analytic chain rule cleans it up.

### Higher-ℓ axial fundamentals locked

| ℓ | calib err | shift vs GR | source |
|---:|---:|---:|---|
| 2 | 0.80% | 4.977% | G109 broad pulse + throat sweep |
| 2 | 0.05% | 5.096% | G112 narrow pulse (cross-check) |
| 3 | 0.005% | 2.101% | G112 |
| 4 | 0.006% | 1.158% | G112 |

Monotonic decrease with ℓ confirms framework's eikonal-recovery commitment.

### Polar first-pass: STAM breaks Schwarzschild isospectrality

G114 used the Zerilli analog `V_Zerilli_proxy(M_eff)` + canonical `(√f)''/√f` correction (M_eff proxy is the polar analog of G87's axial proxy). Parity argument for polar: the two LM constraints `δλ₁` and `δλ₂` algebraically eliminate δΣ in terms of metric perturbations (δλ₂ pins δ_t δΣ to δu^r, δλ₁ pins δ_r δΣ to δΣ itself). δΣ is fully constrained, not propagating. Only the graviton remains dynamical.

| ℓ | axial shift | polar shift | polar/axial |
|---:|---:|---:|---:|
| 2 | 4.977% | 5.620% | **1.129** |
| 3 | 2.101% | 2.324% | 1.106 |
| 4 | 1.158% | 1.206% | 1.041 |

**Key prediction: at ℓ=2, polar QNM shift is 13% larger than axial.** Isospectrality break shrinks to 4% at ℓ=4 (eikonal restoration). This is the cleanest LIGO O5+ / LISA EMRI signature because parameters common to both modes cancel in the ratio.

**Caveat (sharp)**: G114 is first-pass with M_eff proxy + canonical correction by analogy. A G114-rigorous (full quadratic-action derivation in even-parity RW gauge, polar analog of G108) is the natural follow-up but a significant calculation in its own right.

### Exact Kerr spheroidal photon-region surface (Open Problem closed)

G85 derives r_pr_exact(θ; a) in closed form from the spherical-null-geodesic conditions R(r) = R'(r) = 0:

```
λ(r_p) = -(r_p³ - 3r_p² + a²r_p + a²) / [a(r_p - 1)]
η(r_p) = r_p²[4a²Δ - (r_p² - 3r_p + 2a²)²] / [a²(r_p - 1)²]
Θ(θ; r_p, a) = η sin²θ - cos²θ (λ² - a² sin²θ)
r_pr_exact(θ; a) = brentq{r_p in [r_ph_prograde, r_polar]: Θ = 0}
```

G84's sin²θ ansatz is OFF by up to **39% in r_pr at a=0.99** (14.8% in Σ_ph). For Kerr precision predictions, use r_pr_exact. Structural conclusions (ghost-freedom, cubic horizon vanishing, eikonal recovery) are UNCHANGED — they depend on topology of the photon-region surface, not its exact shape.

G116 ships this as a reusable map + W_K factor function at mid-shell.

### Kerr STAM wave base — proper geometry (G119)

After G118 broke at high spin (Schwarzschild h(r) = 1 - 2M/r goes negative below r=2M, but Kerr r_+(a) drops below 2M for a > 0.5):

G119 builds the proper Kerr STAM wave base with:
- Kerr horizon r_+ = M + √(M² - a²) (not r=2M)
- Kerr radial structure Δ(r,a) = r² - 2Mr + a² (not Schwarzschild form)
- Kerr tortoise dr*/dr = (r²+a²)/Δ_STAM
- Δ_STAM = Δ · F(y_K_eq) inside r_pr_exact (photon-region boundary)
- Cubic horizon vanishing Δ_STAM ~ (r-r_+)³ preserved
- **All 9 sampled spins (a ∈ [0, 0.99]) numerically sound** — no imaginary structure anywhere

STAM throat width in r* coordinate is enormous (10⁶-10⁸ from cutoff ε=10⁻⁵), reflecting G90's power-law throat divergence carried over to Kerr. Stretching ratio STAM/GR decreases with spin (10⁶ at a=0 → 4×10⁴ at extremal) because the inside-shell region itself shrinks.

**This is the proper wave base for any rigorous Kerr STAM QNM work.**

### G118v2 corrected spin sweep (using G119 base)

| a | shift |
|---:|---:|
| 0.00 | 4.570% |
| 0.10 | 5.267% |
| 0.20 | 5.836% |
| 0.30 | 6.163% |
| 0.50 | 6.074% |
| 0.70 | 5.907% |
| 0.90 | **14.082%** |
| 0.95 | **10.792%** |
| 0.99 | **12.118%** |

**Headline trend**: STAM Kerr ringdown shift grows from ~4.6% at a=0 to **10-14% at high spin**. Diametrically opposite to G118 original's broken decline-to-zero.

**Critical caveat (Sean's discipline)**: G118v2 uses a Kerr-Δ-aware Schwarzschild-RW-analog `V_base`, NOT the rigorous Kerr Teukolsky potential. **Absolute frequencies are NOT Kerr QNMs.** The shift % is the meaningful number because both omega_GR_base and omega_STAM come from the same V_base — the spin-dependent f-correction trend is real even if the absolute frequencies aren't published Kerr values.

**Sean explicitly retired G118/G119b as Kerr predictions** at the G120 stage: "Any Kerr STAM QNM script must reproduce this GR table first." G118v2 is historical diagnostic, not framework prediction.

### GR Kerr reference table (G120a) PASSED

Via the qnm package (Stein 2019, Leaver continued-fraction + Berti's tabulation):

| a | Re(ω) | Im(ω) | A_lm |
|---:|---:|---:|---|
| 0.00 | 0.3737 | -0.0890 | 4 |
| 0.30 | 0.4195 | -0.0877 | 3.653 + 0.075i |
| 0.50 | 0.4641 | -0.0856 | 3.342 + 0.129i |
| 0.70 | 0.5326 | -0.0808 | 2.903 + 0.183i |
| 0.90 | 0.6716 | -0.0649 | 2.110 + 0.211i |

Saved as `results/G120a_gr_kerr_qnm_reference.json` and `.csv` for downstream use. **This is the calibration baseline for any future standalone Kerr QNM solver.**

## What FAILED this session (documented dead-ends)

### n=1 first overtones (three methods)

| Method | Calib error |
|---|---:|
| G110 two-mode `curve_fit` | 30–85% |
| G111 matrix pencil | 20–67% |
| G112 residual subtraction (Prony-type) | 17–66% |

Time-domain n=1 extraction fundamentally hits a precision wall when A_n=1 / A_n=0 ~ 10⁻³ to 10⁻⁴. The standard tool is Leaver continued fraction (frequency domain). G113a Leaver attempted but failed coefficient self-check (transcription error). Open for G113-redux.

### Standalone Kerr GR QNM solver (G120b, G120c)

**G120b** — direct Teukolsky amplitude shooting:
- R(r) grows by ~10²² across the integration (r_+ + ε → r_match)
- The "small" ingoing amplitude A_in drowns in floating-point noise (15-digit precision can't resolve A_in relative to A_out × 10²²)
- Newton converges to **anti-QNM** modes (Im(ω) > 0) — the time-reverse instead of decaying mode
- Calibration: 12-28% error at all 5 spins

**G120c** — Riccati log-derivative shooting:
- Y = R'/R is bounded (no exponential growth), avoiding G120b's precision wall
- BUT leading-order BCs at horizon and infinity are insufficient for <0.5% calibration
- |mismatch| ≈ |Y_L| + |Y_R| at qnm ω — Y_L and Y_R point in opposite directions in complex plane (phase drift ~ π from BC accumulated error)
- Newton converges to spurious roots near Re(ω)=0.1-0.15 (fixed-point structure of simple BC)
- Calibration: 68-135% error

**Diagnosis**: Both failures are different facets of the same issue — **writing a publication-grade Kerr QNM solver from scratch requires multi-hundred-line dedicated implementation**, not a one-shot script. Specifically need (some combination of):
1. Higher-order Frobenius BC at horizon (a_1, a_2 coefficients from recurrence)
2. Higher-order asymptotic at infinity (α_1/r, α_2/r² corrections)
3. Singular-part subtraction (Y_lead trick)
4. Sasaki-Nakamura transformation (short-range potential, no exponential growth)
5. Leaver continued fraction with verified coefficients

**G120b and G120c are documented dead-ends.** Their V_T construction, Frobenius BC structure, and asymptotic decomposition formulas are correct in form; the failure is BC accuracy at the level required.

## All scripts created this session

### Wave-branch refactor + rigorous axial
- [scripts/G107v2_wave_equation_harness.py](../scripts/G107v2_wave_equation_harness.py) — clean swappable potential harness
- [scripts/G108_axial_from_constrained_sigma_action.py](../scripts/G108_axial_from_constrained_sigma_action.py) — rigorous axial derivation; symbolic background match, parity argument, analytic V_exact
- [scripts/G109_axial_convergence_study.py](../scripts/G109_axial_convergence_study.py) — throat-depth sweep locking 4.977%

### Higher modes and overtones
- [scripts/G110_axial_QNM_higher_modes.py](../scripts/G110_axial_QNM_higher_modes.py) — ℓ=3,4 fundamentals (clean)
- [scripts/G111_matrix_pencil_overtones.py](../scripts/G111_matrix_pencil_overtones.py) — failed
- [scripts/G112_residual_subtraction_overtones.py](../scripts/G112_residual_subtraction_overtones.py) — n=0 cleaner; n=1 failed
- [scripts/G113a_leaver_continued_fraction_schwarzschild.py](../scripts/G113a_leaver_continued_fraction_schwarzschild.py) — coefficient self-check failed (defensive abort)
- [scripts/G113b_frequency_domain_shooting.py](../scripts/G113b_frequency_domain_shooting.py) — BC contamination at finite rs_min

### Polar
- [scripts/G114_polar_QNM_from_constrained_sigma.py](../scripts/G114_polar_QNM_from_constrained_sigma.py) — first-pass; 13% isospectrality break at ℓ=2

### Exact Kerr photon-region
- [scripts/G85_exact_kerr_photon_region.py](../scripts/G85_exact_kerr_photon_region.py) — closed-form r_pr_exact(θ; a)
- [scripts/G116_Kerr_exact_photon_region_STAM_map.py](../scripts/G116_Kerr_exact_photon_region_STAM_map.py) — already committed; reusable map + W_K function

### Kerr ringdown progression
- [scripts/G118_kerr_ringdown_first_pass.py](../scripts/G118_kerr_ringdown_first_pass.py) — Schwarzschild base, broke at high spin
- [scripts/G118v2_kerr_ringdown_first_pass_redux.py](../scripts/G118v2_kerr_ringdown_first_pass_redux.py) — G119 Kerr base; spin sweep 4.6%→14%; **historical diagnostic, NOT Kerr prediction**
- [scripts/G119_kerr_wave_base_geometry.py](../scripts/G119_kerr_wave_base_geometry.py) — proper Kerr STAM wave base, all spins PASS
- [scripts/G119b_kerr_base_spin_sweep_redux.py](../scripts/G119b_kerr_base_spin_sweep_redux.py) — failed WKB stepping stone

### Standalone Kerr QNM solver attempts (failed)
- [scripts/G120a_gr_kerr_qnm_calibration.py](../scripts/G120a_gr_kerr_qnm_calibration.py) — PASSED via qnm library
- [scripts/G120b_gr_kerr_teukolsky_solver.py](../scripts/G120b_gr_kerr_teukolsky_solver.py) — FAILED (amplitude shooting precision wall)
- [scripts/G120c_kerr_qnm_riccati_solver.py](../scripts/G120c_kerr_qnm_riccati_solver.py) — FAILED (log-deriv BC accuracy)

## README state at session close

Substantial edits to QNM section + headline result #7 + snapshot. Key features:

- Locked spinless axial table (ℓ=2/3/4 n=0): 4.977%, 2.101%, 1.158%
- Polar first-pass table with isospectrality-break ratio column
- Kerr photon-region map G85/G116 referenced; G84 ansatz superseded
- Kerr STAM wave base G119 noted as numerically sound
- G118v2 spin sweep noted as historical diagnostic, NOT Kerr prediction
- G120a reference noted as PASSED (calibration baseline)
- G120b, G120c noted as FAILED standalone solver attempts
- Open Problem #7 (QNM partial closure) updated through G114
- Snapshot line names G108/G109 4.977% as the locked spinless axial

**Sean's discipline registered in the README**: STAM Kerr predictions require a standalone GR solver reproducing G120a values < 0.5% FIRST. Until that exists, G118v2's spin sweep is diagnostic, not framework prediction.

## Open priorities for next session

If continuing wave-branch work:

1. **G120d — proper standalone Kerr QNM solver.** Options:
   - **Sasaki-Nakamura with verified coefficients** (~300 lines, textbook approach, no exponential-growth issue). Most rigorous next step.
   - **Leaver continued fraction with verified Kerr coefficients** (~200 lines once coefficients are correct; G113a failed exactly here — recommend transcribing line-by-line from qnm source with explicit self-checks).
   - **Higher-order BC corrections on G120c's Riccati approach** (~150 lines for Frobenius + asymptotic series). Lighter but probably still difficult to hit <0.5%.

2. **G113-redux — Leaver CF for spinless overtones.** With verified Schwarzschild Leaver coefficients (transcribed from qnm source), produce n=1, n=2 axial overtones at <0.5%. Then check whether STAM modification of these gives a meaningful prediction (overtones probe different parts of the potential than fundamental).

3. **G114-rigorous — full polar reduction from constrained shell-count action.** Replace M_eff Zerilli + canonical correction with full δ²S derivation in even-parity Regge-Wheeler gauge. Confirms or revises the 13% isospectrality break headline.

4. **Once G120d lands**: G120e (insert Δ → Δ_STAM, sweep spin, produce framework's first prediction-grade Kerr QNM table).

If shifting away from QNM:

- **Cosmology pipeline** (G104 was scaffold; needs real CSV data with covariances + MCMC)
- **F6 decoherence specific calculation** (~0.5 s for 1 µm silica from generic Γ_res scaling — derive from specific L_μ for matter-substrate interaction)
- **Microscopic SU-write coarse-graining** (deriving the constrained shell-count action from underlying Γ_res dynamics)

## Notes for future Claude

### "First reproduce GR, only then run STAM."

Sean's discipline, restated multiple times this session in different forms:

> "Do G120a ONLY: GR Kerr calibration only. No STAM substitution until published Kerr QNMs are reproduced."
> "Any Kerr STAM QNM script must reproduce this GR table first."
> "Do not use the old G118/G119b proxy shifts as Kerr predictions anymore."

This is non-negotiable. G118v2 produced an attractive STAM spin trend (4.6% → 14%), but Sean explicitly retired it as a Kerr prediction at the G120 stage because V_base isn't the actual Kerr potential. The shift % is "diagnostic" until a rigorous solver gives Kerr-QNM-grade absolute frequencies.

### "Don't claim predictions before calibration passes."

This is the G113 lesson echoed. G113a failed coefficient self-check and Sean designed G120 to make exactly this kind of failure impossible: G120a establishes ground truth before G120b/c attempt to reproduce. Even though G120b/c failed, the discipline kept the failure visible and isolated.

### Failure honesty pays

I documented G120b and G120c failures explicitly with diagnostic detail rather than tuning to make them pass. Sean endorsed this — both scripts now sit as "documented dead-ends" with clear next-step recommendations. The wave branch state at session close is HONEST: rigorous spinless axial closed; rigorous standalone Kerr QNM solver remains open.

### Sean's structural moves

When G118 broke at high spin, Sean immediately diagnosed it ("Schwarzschild h(r) goes negative below r=2M, but Kerr r_+ < 2M for high spin") and specified G119 as the fix — geometry only, no QNM. That's the pattern: when a calculation fails, Sean isolates the structural cause and specifies a narrower replacement that ONLY does the necessary piece. G119 PASSED cleanly because its scope was tight.

Same with G120: Sean specified the three-script structure (G120a reference, G120b standalone solver, G120c spin sweep — later refined to G120c as Sasaki-Nakamura with log-derivative). The discipline of "calibration first, then physics" is structural, not just methodological.

### Terminology update (forward-looking)

Sean asked for: 
- "ledger closure" → "manifold-support closure"
- "ledger-channel F(y)" → "final-shell support profile"
- Keep: SU shell-count, Σ = 3A, 1 SU = A_0

No existing G118-G120 scripts contained the old wording, so this is forward-only. Apply to future scripts that touch the F(y) / Beta(D+1,2) closure.

### Wave branch state at session close

```
Spinless rigorous closure:
  ℓ=2, n=0 axial      4.977%     RIGOROUS (G108/G109)
  ℓ=3, n=0 axial      2.101%     RIGOROUS (G112)
  ℓ=4, n=0 axial      1.158%     RIGOROUS (G112)
  
Spinless first-pass:
  ℓ=2, n=0 polar      5.620%     FIRST-PASS (G114; M_eff Zerilli)
  ℓ=3, n=0 polar      2.324%     FIRST-PASS (G114)
  ℓ=4, n=0 polar      1.206%     FIRST-PASS (G114)
  Isospectrality break ~13% at ℓ=2

Kerr geometric foundation:
  r_pr_exact(θ; a)               CLOSED (G85, G116)
  Kerr STAM wave base            CLOSED (G119, all spins PASS)
  GR Kerr QNM reference          CLOSED (G120a via qnm)
  
Kerr ringdown predictions:
  Standalone solver < 0.5%       OPEN (G120b/c failed)
  STAM Kerr QNM table            BLOCKED on standalone solver

Overtones:
  n=1 spinless axial             OPEN (time-domain methods all failed)
  n=2+                           OPEN
```

The framework's Kerr prediction is the prize. Geometry is in place; the solver is the gate.

## Sean's closing posture

Brief. "That's okay thank you, please create a handoff for the next session and update the readme."

No specific direction for next session.

Pick up at G120d (proper standalone solver) or wherever Sean flags as priority.

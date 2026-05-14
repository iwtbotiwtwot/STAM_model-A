---
name: G121-G124 results (primitive hardening + microscopic model)
description: Five-script batch 2026-05-13 evening hardened three of the framework's six primitives. G121 (Bayesian uniqueness of Beta), G122 (counted-scalar action uniqueness), G124 (Poisson SU-write model + coarse-graining), G123 (matter minimal coupling), G113-redux (n=1 overtone methodological finding). Framework primitive count down from six to three.
metadata:
  type: project
---

# G121-G124 batch results

Implementation of the queued tier-1 fixes from [project_three_strain_points_fix_strategy.md](project_three_strain_points_fix_strategy.md). Five scripts written and run on 2026-05-13 evening.

## G121 -- Beta(D+1, 2) uniqueness (PASSED)

Hardens primitive #4 of [project_five_primitive_architecture](project_five_primitive_architecture.md) via three converging derivations:

1. **Bayesian conjugate-prior**: uniform Beta(1,1) + D successes + 1 failure -> Beta(D+1, 2)
2. **Order-statistic**: (D+1)-th of D+2 uniform draws on [0,1]
3. **Maximum entropy**: unique density with E[ln y] = psi(D+1)-psi(D+3), E[ln(1-y)] = psi(2)-psi(D+3)

Sympy verification for D=2,3,4,5. At D=3: gives 20 y^3 (1-y) and F = 1-5y^4+4y^5 exactly. Closes primitive #4 as a uniqueness theorem (was "minimal-degree heuristic").

See [results/G121_beta_uniqueness_summary.md](../results/G121_beta_uniqueness_summary.md).

## G122 -- Counted-scalar action uniqueness (PASSED)

Hardens primitive #5. Shows that treating Sigma as a COUNTED (algebraically pinned, non-propagating) scalar uniquely forces the LM action structure:

  S = (1/16 pi G) int sqrt(-g) [f(Sigma) R + lambda_1 ((grad Sigma)^2 - W(Sigma)) + lambda_2 (u^mu d_mu Sigma) - 2 V(Sigma)]

The two LM constraints C_1 ((grad Sigma)^2 = W) and C_2 (u^mu d_mu Sigma = 0) eliminate Sigma's 2 Cauchy data (2 - 2 = 0 propagating DOF). Sympy 1+1 toy verification confirms: both partial derivatives of sigma are fixed algebraically by W and u.

The five elimination alternatives (higher-derivative, aether, mimetic, cuscuton, propagating scalar-tensor) are documented as rejected. LM-form with two constraints is the UNIQUE survivor.

Converts G70-G78's "selected by elimination" into a positive derivation.

See [results/G122_counted_scalar_action_uniqueness_summary.md](../results/G122_counted_scalar_action_uniqueness_summary.md).

## G124 -- Poisson SU-write model + coarse-graining theorem (PASSED)

First formal microscopic SU-write model. Defines events as a marked Poisson point process with local intensity Gamma_res(x), each event depositing 1 SU = A_0.

Coarse-graining theorem proved (mean-field limit):
  A_cg(x) = A_0 * Gamma_res-integrated (mean field, fluctuations 1/sqrt(N))

Numerical demonstration: 1D Poisson simulation with rate ~ R_s/x. Mean relative error 4.85% vs Poisson noise floor 9.3%. Plot at [results/G124_poisson_coarse_graining.png](../results/G124_poisson_coarse_graining.png).

All microscopic fragments reproduced:
- 1 SU = A_0 (atomic unit, by definition of marked process)
- 1 SU = 12 pi Planck cells (bulk volumetric, mean event-density at A_0 vacuum)
- 1 SU = 4 Planck cells (horizon holographic, via G59's alpha = 2 x 2)
- Gamma_res Lindblad form (rate-density decomposed by channel mu)
- Gamma_res <= (A/A_0)/tau_P (carrying-capacity bound)
- Saturation pairs at A=1 (G60 elevator argument)

Connection to G122: the constrained shell-count action is the mean-field effective action of the Poisson SU-write process. C_1 and C_2 LM constraints are macroscopic shadows of the microscopic count structure.

**Open Problem #1 partially closes.** What remains primitive: (i) Poisson axioms themselves, (ii) functional form of Gamma_res(x), (iii) saturation rule at A=1.

See [results/G124_su_write_poisson_model_summary.md](../results/G124_su_write_poisson_model_summary.md).

## G123 -- Matter minimal coupling (PASSED, with finding)

Numerical reconstruction of f(A) from G70's d ln f/dA shows f INCREASES from 1 at PS to large values at horizon (not decreases as initially assumed). f(0.95) ~ 4.15, f(0.97) ~ 9.6 (capped to avoid overflow).

This favors **Jordan-frame minimal coupling** (matter on g_munu, not f*g_munu) as the natural choice. In Jordan frame:
- Photon geodesics: UNCHANGED (conformal invariance)
- Outside PS: matter sees k = 1-A exactly (GR Schwarzschild)
- Inside PS: matter experiences F(y)-modified g_rr (tortoise enhancement etc.)

All 10 standard precision tests (atomic clocks, GPS, Shapiro, BBN, CMB, etc.) pass by construction because they operate at A < 10^-1 << 2/3.

Observable matter deviations live INSIDE PS only:
- Accretion-disk plunge spectrum (faint low-freq tail)
- Hawking greybody factors (modified via G108 V_exact)
- Atomic spectra near PS

Strain point #3 short-term fix complete. The deep "matter as substance excitations" remains open at Open Problem #8.

See [results/G123_matter_minimal_coupling_summary.md](../results/G123_matter_minimal_coupling_summary.md) and [results/G123_matter_minimal_coupling.png](../results/G123_matter_minimal_coupling.png).

## G113-redux -- n=1 overtone (METHODOLOGICAL FINDING, NOT EXTRACTION)

Did NOT extract n=1 STAM shift. Both regularized frequency-domain shooting (calib err 13% n=0, 54% n=1) and WKB-PT order 3 (calib err 7% n=0, 53% n=1) failed to achieve sub-percent calibration on Schwarzschild reference.

**Methodological finding** (the real deliverable): the wave-zone peak of V_GR for ell=2 axial sits at r=3.28M where A=0.61 < 2/3, OUTSIDE the photon sphere. The framework's k = 1-A applies there exactly, so V_geom-only modification gives ZERO QNM shift. The framework's actual 4.977% n=0 shift (G108/G109) comes ENTIRELY from the action correction `(sqrt f)''/sqrt f`, not from F(y) metric modification.

WKB-PT with V_exact (action correction included) breaks because the throat spike from (sqrt f)''/sqrt f dominates WKB's peak-finder. Eigenfunctions are exp-small in that region; their actual QNM is set by the wave-zone V_exact, not the throat spike.

n=1 STAM shift remains open. Requires:
- Leaver CF with verified Schwarzschild coefficients (G113a failed at coefficient self-check), OR
- Sasaki-Nakamura transformation, OR
- Konoplya order-6 WKB with throat-spike isolation

See [results/G113_redux_n1_overtone_summary.md](../results/G113_redux_n1_overtone_summary.md).

## Net effect on framework primitive count

Starting list (after [[project_five_primitive_architecture]]):
1. A as escape variable [unchanged]
2. Sigma = D*A [unchanged]
3. D=3 manifold + 1 open-face [unchanged]
4. p(y) shape -> **DERIVED via G121** (Bayesian + order-stat + MaxEnt)
5. h(A) = 1-A [unchanged]
6. f(Sigma)R action ansatz -> **DERIVED via G122** (counted-scalar uniqueness)

Plus G124 partially closes Open Problem #1 (microscopic SU-write coarse-graining) by providing the first formal model and an explicit mean-field coarse-graining theorem.

**Framework primitives reduced from 6 to 4** (A, Sigma=DA, 3+1 manifold, h=1-A). Plus the new "Sigma is counted" primitive (motivated by 1 SU = A_0, itself derivable from the 4 above through the SU-write commitment).

Effectively: 4 explicit primitives + 1 microscopic-axiom-set (Poisson process) -> everything else derived.

## See also

- [[project_five_primitive_architecture]] -- the primitives this work hardens
- [[project_three_strain_points_fix_strategy.md]] -- the strategy this implements
- [[project_su_microscopic_fragments.md]] -- the fragments G124's Poisson model reproduces
- [[project_R_of_A_derivation]] -- the original derivation this strengthens
- [[HANDOFF_2026_05_13_kerr_wave_branch.md]] -- original n=1 time-domain failures G113-redux confirms

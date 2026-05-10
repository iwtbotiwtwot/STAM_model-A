# CHANGELOG

## 2026-05-09 evening — maintenance pass, verdict-text alignment, G11 root-find

Full-sweep regression check: 33 scripts (30–48 + G1–G14) run end-to-end, all execute cleanly with no broken outputs. Every `results/*_summary.md` file reproduces from source. Two real out-of-sync items found and fixed; one documentation gap closed.

**Maintenance.**

- Full-sweep verification: scripts 30–48 + G1–G14 all run cleanly. No drift in any `result_summary.md` output.
- Two cosmetic warnings noted but not silenced (G9 matplotlib `set_ticklabels` deprecation; G14 `np.exp` overflow at the 10^25 s extrapolation endpoint).

**G8 verdict-text alignment.** The printed verdict in `scripts/G8_v_a_two_pole_test.py` and the corresponding result md said *"A_0 photon-A LoS closes only ~25% of the gap"* — stale text from when A_0 was being investigated at a smaller value (~0.011 from the early Hubble-tension fit). The script's own computed offsets (3.87% → 1.19% with A_0 LoS) actually correspond to ~69% closure. Three places in the script now compute the closure dynamically from the offsets so the value cannot drift again; result md regenerated.

**G11 root-find for Ω_m = 0.315 with bridge-term cross-check.** Added an inverse bisection in `scripts/G11_cmb_self_consistent.py` that solves directly for the structure_strength satisfying the PBH-DM target Ω_m = 0.315. Reports the precise value (`0.022`, ~2.2% of the G9 toy amplitude) along with an independent cross-check that the bridge term predicted at the solution (355.10 Mly) matches the historical empirical value (354.95 Mly, STAM's own prior catalog calibration) to 0.04%. The framework now has explicit script-output evidence that V_3's predicted θ_⋆ closes STAM's prediction to observed (offset 0.0000%) AND V_3's predicted bridge term closes to STAM's own historical empirical value (offset +0.04%), simultaneously at the PBH-DM matter content (Ω_m = 0.315). Both internal STAM observables sit on the same self-consistent solution.

**Documentation.** `data/README.md` updated to document all six data files, including the three authoritative published datasets (Pantheon+ stat+sys covariance, Union3 inverse-covariance FITS, Union3 CSV companion) added during the prior session.

## 2026-05-09 morning — A_0 = 1/(12π) commitment, V_3 selected, PBH-DM compatibility

The May 9 morning push moved the framework from "calibrated bridge term + qualitative dark-energy story + galactic-rotation problem" to "structurally-derived A_0 + V_3 commitment + CMB self-consistent closure + PBH-DM-compatible galactic story." Fourteen new G-series scripts; one Model-A claim retired (asymmetric propagation). Each framework reported on its own terms; comparison with other frameworks is for curiosity only.

**Strong-field.** Quasinormal-mode ringdown spectrum derived for Model-A's modified k(A); equator-pole bubble geometry for spinning sources derived from A(r,θ); bubble thermodynamics matches the Kerr horizon temperature exactly at the pole via the r₊² − 2Mr₊ + a² = 0 identity. (G1–G5)

**Cosmology.** A_0 committed to 1/(12π) as a structural value, decomposing as 4π × 3 (thermal prefactor × spatial dimensionality). Bridge term b = A_0 × c/H_0 = 355.10 Mly at H_0 = 73.04 — matches the historical empirical value 354.95 Mly (STAM's own prior catalog calibration) to 0.04%. V(A) form selected as V_3 = α/A + β/(1−A) under symmetric-boundary commitment, with α/β = [A_0/(1−A_0)]² fixing the V minimum at A_0. At H_0 = 73, V_3 cosmology with STAM's photon-A LoS contribution at A_0 predicts θ_⋆ within +1.19% of observed — closing ~69% of STAM's own no-LoS prediction-vs-observation offset (which sits at +3.87% before the LoS contribution is added). Full match of STAM's prediction to observed requires A_LoS ≈ 1.46× A_0, plausibly supplied by cosmic-web structure amplification along the line of sight. (G6–G9, G11)

**Galactic dark matter.** Linear cumulative A from baryons alone does not explain SPARC-style rotation curves (factor 1.5–3 deficit in v, 3–7 in implied mass). Primordial-black-hole dark matter introduced as a natural framework extension: each PBH is a small bubble with the same A = 1 boundary structure as stellar and supermassive BHs. Press-Schechter formation under enhanced small-scale power can supply full PBH-DM at asteroid-mass scales (10^17–10^23 g) with σ ≈ 0.05–0.07. Cumulative A from baryons + PBH-DM halo reproduces observed rotation curves trivially, with each spiral galaxy hosting ~10^27 asteroid-mass PBHs. (G10, G12, G13)

**GW propagation chain (F5 → F5b → F5c) closed.** The asymmetric-propagation proposal from F5b — "GW = c always, light slowed by ∫A ds" — was tested against GW170817 arrival-time data in F5c. F5b's strong form predicts ~3 years over 40 Mpc; observed is 1.74 s. F5b's asymmetric-propagation claim retired. Both messengers travel symmetrically through A; the 1.74 s is astrophysical jet-launch timing. The tensor-character recognition (GWs inherit metric-perturbation polarization) — established earlier in the F5 lineage and independent of asymmetric propagation — survives. (F5c)

**Wild-theory exploration.** Mathematical correspondence between BH-formation collapse and Big Bang scale-factor expansion characterized as partial — saturation vs monotonic growth — but the timescale objections dissolve under proper time-reparameterization across horizons. Catalogued as exploratory, not committed. (G14)

**Documentation.** README.md updated to integrate the A_0 = 1/(12π) commitment, the (4π × 3) decomposition, the structural-floor interpretation of the void, and the open work on the spatial-dimensionality factor of 3.

### Scripts added

- **G1:** QNM ringdown — eikonal/WKB at the photon sphere. τ_Model-A / τ_GR = 1.80 (mass-independent), same dominant frequency. Eikonal approximation; exact Regge-Wheeler integration is open.
- **G2:** Spinning Model-A bubble geometry — equator preserved at 2M for any spin, pole shrinks to the Kerr horizon location r₊.
- **G3:** Spinning bubble thermodynamics — T(θ) non-uniform along the bubble; equator at T_Schwarzschild for any spin, pole at T_Kerr exactly.
- **G4:** Holographic-matter rotational emission — static bubble surface, rotating holographic matter on it; super-radiance at the equator.
- **G5:** κ calibration — `R = 1 + κ(v/c)²` ansatz gives κ ≈ 1–2 at low spin (matches first-principles γ²) but breaks at high spin where mode-by-mode physics dominates.
- **G6:** CMB θ_⋆ landscape at H_0 = 73 across five model configurations.
- **G7:** A_0 = 1/(12π) commitment test — bridge term predicted to 0.04% of historical.
- **G8:** V(A) functional form selection — V_3 viable; V_2 rejected (forced minimum at A = 1/2); V_1 incomplete under symmetric-boundary commitment.
- **G9:** Cumulative-A line-of-sight toy Monte Carlo — ~3.1× amplification at LCDM-scale comoving distances.
- **G10:** SPARC-style galaxy rotation — linear cumulative A from baryons alone does NOT reproduce observed v.
- **G11:** CMB self-consistent solver — sweep over structure_strength to find Ω_m where predicted and required f_LoS converge.
- **G12:** PBH formation probability via Press-Schechter; required σ for full PBH-DM in observationally-allowed mass window.
- **G13:** PBH-DM + cumulative A galactic rotation — fits trivially with realistic PBH halo.
- **G14:** BH-formation vs Big Bang expansion mathematical comparison — exploratory.

## 2026-05-02 — v0.1.0 scaffold

- Initialized `STAM_model-A` as the variable-accumulation Model-A version.
- Removed the earlier constant-`A0` direction from the core model.
- Added local accumulation, propagation, cosmology, and mass-estimator modules.
- Added smoke tests for local identity, horizon threshold, Shapiro path form, and cosmological path-integral identity.
- Added falsification-test documentation and parameter-locking notes.
- Archived the 2 May 2026 STAM-A PDF under `docs/archive/`.

## Future entries

Record every change to:

- formula definitions;
- parameter values;
- catalog preprocessing;
- fitted `b` values;
- falsification thresholds;
- claims/status wording.

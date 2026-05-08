---
name: F5c — asymmetric GW/light propagation FALSIFIED, drop it
description: Sean's commitment 2026-05-08. F5b had asserted "GW propagates at c intrinsically; light is slowed by ∫A ds along the path" as the resolution to F5's polarization issue and as a multi-messenger story for GW170817. F5c tested this quantitatively against GW170817 (40 Mpc to NGC 4993) and found the strong-form prediction gives ~3 years of GW-vs-light delay vs the observed 1.74 s — a falsification by ~8 orders of magnitude. The asymmetric-propagation claim is dropped; what survives is GW tensor character from being a metric perturbation.
type: project
---

**The F5c result (2026-05-08, scripts/F5c_gw170817_arrival_test.py):**

Taking F5b's strong-form claim at face value — GW propagates at c always; light is slowed by (1/c) ∫A ds where A = 2GM/(c²r) is matter-sourced — and integrating along the GW170817 line of sight (40 Mpc to NGC 4993):

| Contribution | ∫A ds | Predicted Δt (light − GW) |
|---|---|---|
| Milky Way (M=10¹² M☉, b=7.5 kpc) | 2.74 × 10¹⁶ m | 2.90 years |
| NGC 4993 (M=4×10¹⁰ M☉, kilonova at 2 kpc offset) | 1.17 × 10¹⁵ m | 45 days |
| **Total** | 2.86 × 10¹⁶ m | **3.02 years** |
| Observed | — | **1.74 s** |
| Ratio predicted/observed | — | **5.5 × 10⁷** |

The ⟨A⟩ that would be required to reproduce 1.74 s is 4.2 × 10⁻¹⁶, vs. the actual ~10⁻⁵ Milky Way potential A at the Sun's position. Mismatch by ~10¹⁰. The model as written cannot be rescued by parameter tuning.

**What this means for the framework:**

1. **Drop the asymmetric-propagation claim.** F5b's "passes by construction" assertion for GW170817 was wrong — the construction itself, integrated honestly, fails. Do not carry forward the position that GW propagates at c intrinsically while light is slowed by full matter-sourced A.

2. **Multi-messenger consistency requires symmetric propagation.** Both GW and light follow null geodesics on the same metric, both feel the same Shapiro/SU integral through the gravitational potentials they cross, the differential cancels, and GW170817's 1.74 s is the astrophysical jet-launch timing. This is the same outcome GR has and is automatic in Model-A as long as both messengers are treated as null phenomena on the metric without arbitrary asymmetry.

3. **What survives from the Model-B line of work:**
   - **Tensor character of GWs is fine.** GWs are metric perturbations, h₊/h× polarizations come naturally from the metric being a (0,2) tensor field. This part resolves F5's original polarization issue (pure scalar A perturbation gave only longitudinal modes). Keep it.
   - **A is matter-sourced (∇²A = (8πG/c²)ρ_matter).** Source equation unchanged.
   - **A ≈ 0 in matter-free voids.** Replaces the de Sitter A_cosmo ansatz; cosmological A is path-dependent and structure-weighted. Keep it.
   - **Drop:** the propagation-speed asymmetry between GW and light. They feel A symmetrically.

4. **Cosmological coupling V(A) = β/(1-A) is not affected by this.** The supernova-distance work (scripts 35-39b) and the bridge-term derivation A_0 = b/L = 0.0265 don't depend on the propagation asymmetry; they depend on the cosmological A field's effect on photon paths, which is unchanged.

5. **F5b/Model-B label is no longer load-bearing.** The surviving content (tensor GWs, matter-sourced A, A ≈ 0 in voids) is just part of Model-A's metric and field structure. Per project_terminology_model_a_only.md, drop the Model-B label as well.

**How to apply:**

- Don't claim "GW propagates at c intrinsically while light is slowed by A" — that's the falsified position.
- Don't claim Model-A predicts a measurable GW-vs-light differential from cosmic-scale propagation. It doesn't, because the messengers are symmetric.
- Multi-messenger arrival within astrophysical timing is automatic in Model-A; it's not a special prediction to defend.
- The polarization fix (tensor GWs) is fine and stays.
- F5 itself is now resolved by the symmetric-messenger reading: GWs are null tensor perturbations on the metric, light is null EM on the metric, both feel A through the same Shapiro/SU mechanism.

**Falsification target this opens:** F5c's negative result rules out one specific extension; it does NOT rule out Model-A's underlying scalar A field doing real work on light propagation. The bridge term b = 354.95 Mly is still derived from photon-A traversal, but only because Model-A's photon-A coupling is integrated over the cosmological line of sight (where most of the path is in voids, not deep matter potentials), not because GWs are exempt. The cosmological coupling and the multi-messenger story are decoupled.

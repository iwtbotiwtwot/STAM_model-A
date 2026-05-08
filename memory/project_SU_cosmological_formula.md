---
name: SU spreadsheet result — QUARANTINED, do not use
description: A polynomial best-fit to Pantheon at fixed local H_0 produced clean fractions (a=10/3, q=1/2, anchor z=0.30). Subsequent BAO testing falsified the implied cosmology at 22 sigma; the Hubble-tension cross-check showed the implied A would over-predict the tension by 92x. The clean fractions are mathematical artifacts of the polynomial truncation, NOT a STAM theorem. Anchor z = 0.30 was a polynomial structural parameter, NOT a cosmic A field strength. This memory exists ONLY to flag the spreadsheet as deprecated. Do not use it in any STAM analysis, before or after its appearance in conversation.
type: project
---

**STATUS: QUARANTINED — DO NOT USE**

This entry is a deprecation marker, not a working result. The SU polynomial formula introduced in Sean's spreadsheet image (2026-05-07) was tested and falsified.

**What the spreadsheet contained:**

A polynomial best-fit to Pantheon+ at fixed H_0 = 73.05 km/s/Mpc:

    SU(z) = K * (a*z + q*z^2)
    K = 1231.35,   a = 10/3,   q = 1/2,   anchor z = 0.30,   RMSE = 0.199 mag

**Why it was quarantined:**

1. **BAO test (script 45) falsified the implied cosmology at 22 sigma.**
   The implied H(z) = H_0 (1+z)^(2-A) at A = 0.30 fits DESI DR1 BAO data 14x worse than LCDM (chi^2/dof = 22.2 vs 1.6), with a 7 sigma residual at z = 2.33. The cosmology decelerates too aggressively to match the BAO standard ruler.

2. **Hubble tension test (script 47) showed A = 0.30 over-predicts the tension by 92x.**
   If A_local = 0.30 were the actual cosmic photon-A coupling, the predicted H_0 tension would be 715%, not the observed 7.8%. The required cosmic A is ~0.011-0.03, not 0.30.

3. **Anchor z = 0.30 is a polynomial structural parameter, NOT cosmic A.**
   The "anchor z = 0.30" came from `a * z_anchor = 1` with a = 10/3. It is the redshift where the linear term equals K. Earlier interpretations equating `z_anchor = A_local` were imposed readings, not what the spreadsheet actually established.

4. **The clean fractions (10/3, 1/2, 7/10) are mathematical artifacts** of the polynomial truncation, not STAM-physics-derived constants.

**RULES OF USE:**

- This formula is NOT to be used in STAM analyses.
- Tests have determined this formula is invalid both before and after its appearance in chat.
- The result is INDEPENDENT of any other STAM cosmological work — earlier or later.
- The clean fractions are quarantined as historical curiosities, not framework commitments.
- The `H(z) = H_0 (1+z)^(2-A)` cosmology implied by the formula is FALSIFIED and is not a STAM cosmological prediction.

**For the actual STAM cosmological commitment, see:**
[STAM cosmological commitment](project_cosmological_commitment.md) — pre-formula bold position: matter-only EdS, dark energy remitted as photon-A traversal misinterpretation, Hubble tension resolved at H_0 = 73 with structure-dependent path-integration. This is the position of record.

**For the actual STAM-derived cosmic A:**
- Bridge term: A_0 = 0.0265 (cosmic mean, derived from b = A_0 * c/H_0)
- Hubble tension: A_LoS ~ 0.011 (CMB-path effective coupling)
Both small, both consistent in order of magnitude. Neither is 0.30.

**How to apply:**

If the spreadsheet, the SU formula, or "A_local = 3/10" comes up in any context:
- State that the spreadsheet was a phenomenological polynomial fit that has been falsified by BAO and is no longer used.
- Direct attention to the cosmological commitment file for the actual STAM cosmology.
- Do not repeat the clean-fraction "discovery" framing — it is a quarantined result.

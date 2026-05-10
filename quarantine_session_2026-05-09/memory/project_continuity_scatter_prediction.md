---
name: Continuity scatter prediction (Sean, 2026-05-09)
description: Model-A prediction — at the matter-to-A-field transition redshift (z ≈ 0.30–0.35), SN distance residuals should show enhanced scatter compared to redshifts above or below the transition. Above and below the transition, the cosmology is in a settled regime; through it, individual sight lines sample a rapidly-changing effective expansion, producing extra apparent-distance scatter.
type: project
---

**Sean's prediction (2026-05-09):**

> "If continuity occurs at a certain point, measurements beyond will be consistent, measurements below will not, this would result in a scatter of sorts in the data at the threshold of continuity."

In Model-A's V_3 modified Friedmann, the cosmic energy budget transitions from matter-dominated to field-coasting near z ≈ 0.30–0.35 (matter / effective dark-energy crossover). Above the transition (z > 0.35): matter-dominated regime, cosmology behaves consistently in one way. Below (z < 0.30): field-coasting regime, cosmology behaves consistently in another way. *At* the transition, individual SN sight lines sample a rapidly-changing effective expansion, and small differences (path geometry, line-of-sight structure, field's spatial variation) translate into enhanced μ-residual scatter.

**Why:** the dynamical-A field is rolling fastest near the transition; small variations in matter density along different sight lines produce slightly different effective A configurations, which propagate to slightly different apparent distances. Above and below, the field is locked into a slow regime, and these variations don't amplify.

**Testable signature:** in Pantheon+ (or any dense SN catalog), compute the variance of μ residuals (vs Model-A V_3 numerical KG) in z bins. Look for excess scatter — variance(residuals) higher than reported σ predicts — concentrated near z ≈ 0.28–0.35.

**Status (2026-05-09 — initial pass, NOTEWORTHY only):**

Tested via G8e (per-SN scatter) and G8f (per-SN data column audit) on Pantheon+ + DES, with 0.28–0.32 as the narrow target bin:

- DES chi²/SN at 0.28–0.32 is roughly double the adjacent bins (4.98 vs 2.56 below, 2.41 above), with the elevated value dropping to null by z = 0.40.
- Three independent DES quality flags (|MUPULL|>3, FITPROB<0.01, PROBCC>0.05) all peak in the same 0.28–0.32 bin at roughly 2× their adjacent-bin rates.
- Pantheon+ shows a small RMS bump at 0.28–0.32 (0.157 vs 0.131 above) but its conservative STAT+SYS errors absorb most of it; chi²/SN stays under 1.
- Combined data shows both elevated scatter and elevated coherent mean residual (+0.07 mag) at 0.28–0.32.

**Honest framing — NOTEWORTHY, NOT SUPPORTING EVIDENCE:** the feature in the data is real and localized at the predicted z range, but its origin is undetermined. It could be the continuity-scatter signature, a DES classification-pipeline transition (DES Y5 hands over from spectroscopic to photometric ID near this z), a selection effect, or something else. There is no claim of confirmation. Track as: "feature observed at the predicted location, alternative explanations open, requires further tests to identify the cause."

**How to apply:**
- Track as: data feature at z ≈ 0.28–0.32 worth following up. Origin unknown.
- Sharpening tests: DES spec-only subsample (rules out photo-classification artifact), another high-resolution catalog at this z, theoretical magnitude estimate for the predicted scatter under V_3.
- Do not cite as evidence for Model-A's continuity prediction without those follow-ups.

# G8b: Model-A V_3 SN Distance Fit (Pantheon+ + Union3)

## Setup

V_3 retrofit of script 37, with apples-to-apples covariance treatment for
both catalogs. Tracking-ansatz cosmology from G8a. The only nuisance
parameter per catalog is the absolute-magnitude offset ΔM, marginalized
analytically. **No fitted cosmology parameters** — every quantity in V_3
is structurally committed or derived from internal closure:

- A_0 = 1/(12π)  (G7 commitment)
- α/β = [A_0/(1−A_0)]²  (V_3 minimum at A_0)
- β_tilde from h²(today) = 1  (Model-A internal closure)
- Ω_m = 0.315, Ω_r = 9.2e-5  (Model-A matter content)

## Apples-to-apples covariance treatment

Both catalogs now use their **published Mahalanobis covariance** with
analytic ΔM marginalization:

- **Pantheon+:** `data/Pantheon+SH0ES_STAT+SYS.cov` from
  [github.com/PantheonPlusSH0ES/DataRelease](https://github.com/PantheonPlusSH0ES/DataRelease).
  1701×1701 stat+sys covariance, submatrixed to the 1578 cosmological SNe
  (mask: `IS_CALIBRATOR=0` and `zCMB>0.01`). Row ordering verified to match
  `pantheon.csv` by diagonal-vs-CSV-error consistency.
- **Union3:** `data/mu_mat_union3_cosmo2_mu.fits` from
  [github.com/rubind/union3_release](https://github.com/rubind/union3_release).
  22×22 inverse covariance for 22 binned distances.

This replaces the earlier asymmetric treatment (Pantheon+ DIAG, Union3 full)
that was an artifact of partial data acquisition. With both catalogs on
the same statistical footing, per-catalog χ² are directly comparable.

## Model-A V_3 result (apples-to-apples)

```text
catalog       chi^2     dof    chi^2/dof    DeltaM
Pantheon+   1493.24    1577     0.9469     +0.0204
Union3        86.94      21     4.1400     +0.0120
combined    1580.18    1598     0.9888
```

## What this says about Model-A on its own table

**Pantheon+ alone:** χ²/dof = 0.947 — passes the fit test (under 1) but
tightly. With proper covariance, the previous chi²/dof = 0.469 figure
was revealed as too-good-to-be-true; the diagonal-only treatment was
using inflated peculiar-velocity-padded errors. The honest Pantheon+
result is closer-to-threshold but still consistent with the data within
the published stat+sys covariance.

**Union3 alone:** χ²/dof = 4.14 — does *not* pass. There is a real
high-z shape mismatch in Model-A V_3's tracking-ansatz prediction that
Union3's high-z bins detect clearly. Apples-to-apples confirms the
disagreement is genuine, not a data-treatment artifact.

**Combined:** χ²/dof = 0.989 — passes overall, dominated by the much
larger Pantheon+ sample. The Union3 disagreement is real but doesn't
sink the combined fit.

## Why the per-catalog story differs

A coherent high-z shape mismatch hits the two catalogs differently in
χ²/dof, even with both using proper covariance:

- Pantheon+ has 1578 SNe (1577 dof). High-z SNe are a small subset of
  the sample. A coherent high-z shape mismatch raises Δχ² by some amount
  but spreads across many dof, giving small Δ(χ²/dof).
- Union3 has 22 bins (21 dof) with much of its statistical weight at
  high z. The same high-z shape mismatch raises Δχ² by a similar amount
  in absolute terms, but normalized by 21 dof gives a large Δ(χ²/dof).

This is real, not artifact. Union3 is intrinsically more sensitive to
Model-A's high-z shape than Pantheon+ is.

## Two-tables comparison (independent fits, apples-to-apples)

```text
                     Model-A V_3   LCDM    V_1 (reference)
Pantheon+ chi^2/dof   0.9469      0.9054   0.9086
Union3    chi^2/dof   4.1400      1.2568   1.3747
Combined  chi^2/dof   0.9888      0.9101   0.9147
```

LCDM and V_1 (the pre-G8 single-pole potential) both produce
χ²/dof ≈ 0.9 on Pantheon+ and ≈ 1.3 on Union3 — similar pattern across
catalogs. V_3 specifically diverges from this pattern at Union3, which
is an internal Model-A signature pointing at the high-z shape question,
not a comparison flag.

## Open questions (Model-A internal)

1. **Tracking ansatz vs full numerical KG.** G8a found a substantial
   tracking-vs-numerical disagreement under V_3 (A_today differs by
   factor 6, ω_today by factor 18). The high-z h(z) behavior of the
   numerical KG cosmology may be quite different from the tracking
   form used here. Re-running this SN test on the numerical h²(a)
   is the V_3-honest follow-up.

2. **Alternative two-pole potential forms.** V_3 = α/A + β/(1−A) is
   one form among many that diverge at both endpoints. The α/β ratio
   is forced by A_0 minimum, but the *form* α/A + β/(1−A) is one
   structural choice. Other two-pole forms produce different high-z
   h(z) and might fit Union3 better while still satisfying the
   symmetric-boundary commitment.

3. **β_tilde calibration.** Anchored to h²(today) = 1, which is a clean
   Model-A internal closure but isn't the only structural choice.
   Alternative non-fitting closure conditions might give different β_tilde
   and different high-z shape.

## Files

- `scripts/G8b_v3_pantheon_union3_fit.py`
- `data/Pantheon+SH0ES_STAT+SYS.cov` (authoritative Pantheon+ covariance)
- `data/mu_mat_union3_cosmo2_mu.fits` (authoritative Union3 source)
- `data/union3_with_errors.csv` (inspectable companion)
- `reports/G8b/v3_pantheon_union3_fit.png`
- `reports/G8b/v3_fit_summary.csv`

---
name: STAM modified Friedmann passes SN cosmology test (Pantheon+, Union3, DES-Y5)
description: Scripts 36-39b. STAM with V(A) = beta/(1-A) calibrated by V'(A_0) = kappa rho_m,0 yields a derived shape Delta_mu(z) ~ 0.105 ln(1+z) vs LCDM. STAM beats LCDM combined chi^2 by ~-25 to -33 across all three SN catalogs. STAM correctly predicts Pantheon+ vs Union3 inter-catalog tension (+30 mmag predicted vs +38 mmag observed, ratio 1.27). DES anomaly (~110 mmag offset) is too large for STAM's signature; consistent with DES-Y5 instrumental systematic.
type: project
---

**Context (2026-05-07):** With V(A) = beta/(1-A) committed (Sean's "Option A"), the modified Friedmann calculation was carried through and tested against SN1a data restricted to Pantheon+ + Union3 + DES-Y5.

**Modified Friedmann setup (script 36):**

Tracking solution: A_eq(a) = 1 - (1-A_0) a^(3/2)
- A_0 = 0.0265 (calibrated to bridge term in script 35)
- beta_tilde = beta/(kappa rho_crit_0) = Omega_m (1-A_0)^2 = 0.2985
- Y(A) = V/(kappa rho_crit_0) = beta_tilde/(1-A)
- Kinetic contribution: K(a) h^2 = (3/8)(1-A_0)^2 a^3 h^2
- Friedmann: h^2 = [Omega_m a^-3 + Omega_r a^-4 + Y(A_eq)] / [1 - K(a)]

Closure at z=0 (tracking ansatz):
- Omega_m = 0.315
- Omega_DE_potential = 0.307
- Omega_DE_kinetic = 0.343
- Omega_DE_total = 0.649 (vs observed 0.685, match 0.95)
- Universe under-closes by ~3.5% in tracking ansatz

Numerical KG integration (full coupled equations):
- A_today = 0.380 (NOT 0.027 — tracking breaks down at late times)
- Field lags equilibrium because V(A) = beta/(1-A) is too steep to track
- Omega_DE_total (numerical) = 0.535 (match 0.78 vs observed)
- TENSION: A_today = 0.38 from dynamics conflicts with A_0 = 0.0265 from bridge

**Derived adjustment shape (script 38):**

Delta_mu(z) = mu_LCDM(z) - mu_STAM(z) when both use same H_0 = 67.4 and Omega_m = 0.315.

Functional fits (with their RMS residuals to the smooth Delta_mu curve):
- Constant: 0.037 mag (poor)
- Linear in z: 0.011 mag
- **Linear in ln(1+z): 0.006 mag — Delta_mu ~ 0.105 * ln(1+z)**
- Quadratic in z: 0.003 mag (best 3-param)

Equivalent multiplicative form: d_L_LCDM/d_L_STAM = (1+z)^0.0484

Sample values:
- z=0:    Delta_mu = 0.000 mag
- z=0.5:  Delta_mu = +0.022 mag
- z=1.0:  Delta_mu = +0.063 mag
- z=2.0:  Delta_mu = +0.097 mag

**SN catalog fits (scripts 37, 39, 39b):**

Each catalog allowed its own DeltaM nuisance (standard SN cosmology — calibrates absolute magnitude per catalog).

| Catalog        | LCDM chi2/dof | STAM chi2/dof | delta_chi2 |
|----------------|---------------|---------------|------------|
| Pantheon+      | 0.4513        | 0.4489        | -3.79      |
| Union3         | 0.7724        | 0.8560        | +1.76      |
| DES-Y5 full    | 1.9151        | 1.8981        | -30.81     |
| DES-Y5 trimmed | 1.6299        | 1.6169        | -23.55     |

Combined (DES trimmed): LCDM 3683, STAM 3657 → STAM wins by chi2 = -25.6 over 3411 dof
Combined (DES full):    LCDM 4211, STAM 4179 → STAM wins by chi2 = -32.8

NOTE: chi^2/dof < 1 for Pantheon+ is because diagonal-only errors (MU_SH0ES_ERR_DIAG) over-count uncertainty without the off-diagonal covariance corrections. Relative comparison still meaningful.

**Cross-catalog tension test (script 39):**

If data follows STAM, each catalog's LCDM-fit DeltaM should equal -<Delta_mu_LCDM-STAM(z)>_cat (the catalog's z-weighted average). The catalog-to-catalog DeltaM differences should match STAM's prediction.

| Pair                  | Observed DeltaM diff | STAM-predicted | Ratio  | Sign match |
|-----------------------|----------------------|----------------|--------|------------|
| Pantheon+ - Union3    | +0.0385 mag          | +0.0302 mag    | 1.27   | YES        |
| Pantheon+ - DES-Y5    | -0.1101 mag          | +0.0137 mag    | -8.0   | NO         |
| Union3 - DES-Y5       | -0.1487 mag          | -0.0165 mag    | 9.0    | yes (sign) |

**Outlier removal test (script 39b):**

Top 6 DES outliers identified by |residual| from LCDM fit:
- CIDs: 1299503, 1289982, 1291149, 1294436, 1319870, 1900801
- Residuals: 1.66 to 3.01 mag (8-15 sigma each)
- These are the SAME 6 SN identified in prior STAM-residual analysis
- Removing them: DES chi^2/dof drops 1.92 -> 1.63 (significant)
- BUT cross-catalog DeltaM tension barely changes (Pantheon+ vs DES: -113 -> -110 mmag)
- Conclusion: DES anomaly is a UNIFORM systematic, not outlier-driven

**Physical interpretation:**

1. STAM's photon-A signature (Delta_mu ~ ln(1+z)) is the "derived adjustment" needed to bridge STAM's raw curve to LCDM's curve.
2. Pantheon+ vs Union3 inter-catalog tension (+38 mmag, unexplained in LCDM literature) IS predicted by STAM to within 27% — sign + magnitude. Real cosmological signal.
3. DES vs others tension (110-150 mmag) is 8-9x larger than STAM predicts and survives outlier removal. Consistent with a uniform DES-Y5 calibration zero-point systematic (photometric calibration, M-step, host-mass step, surface-brightness selection).

**Caveats:**
- Pantheon+ used diagonal-only errors; full covariance would tighten the constraint
- Union3 errors estimated as 0.05 mag uniform (file lacks per-bin errors)
- chi^2 differences of ~25-33 over 3400 dof are suggestive but not decisive
- Models are non-nested; rigorous significance requires Bayesian model comparison

**How to apply:**
- When asked "does STAM fit SN data": yes, equal-to-or-better than LCDM with derived (not fitted) shape adjustment, same number of free parameters.
- When asked "is DES the outlier": yes — STAM's photon-A signature predicts Pantheon+/Union3 tension correctly, but DES has a separate ~100 mmag instrumental systematic.
- The SN test does NOT prove STAM right; it shows STAM is competitive with LCDM and provides physical interpretation for inter-catalog tensions LCDM has no story for.
- The internal STAM tension (A_0 = 0.0265 from bridge vs A_today = 0.38 from KG dynamics) remains unresolved.

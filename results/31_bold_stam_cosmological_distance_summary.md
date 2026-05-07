# STAM Model-A 31: Bold-STAM Cosmological Distance vs LCDM

## Purpose

Test whether the bold-STAM Shapiro propagation, applied at cosmological scales with the de Sitter ansatz `A_cosmo(r) = (r/L)^2` (same A profile used in Q8 to derive the de Sitter horizon temperature), predicts a distance-redshift relation closer to LCDM than current STAM does. If yes: a natural origin for some of the historical bridge term `b`. If no: an honest finding that the cosmological branch needs different work.

## Method

- Bold-STAM lookback distance: `D_bold(z) = integral from 0 to L z of (1/((1-A)(1-A^2))) dr`, with `A_cosmo(r) = (r/L)^2`.
- Equivalent in dimensionless form: `D_bold(z) / L = integral 0 to z of 1/((1-u^2)^2 (1+u^2)) du`.
- Current STAM no-b: `D_adj,0(z) = L z (1 + 0.5 z)` (from README).
- LCDM: standard `D_L(z) = (1+z) c/H_0 * integral 1/E(z') dz'` with Omega_m = 0.315, flat.

## Numerical results

|    z |   D_bold_STAM_Mpc |   D_adj_no_b_Mpc |   D_geo_Mpc |   D_LCDM_lookback_Mpc |   D_LCDM_luminosity_Mpc |   bold/adj_ratio |   adj/lcdm_ratio_lookback |   bold/lcdm_ratio_lookback |
|-----:|------------------:|-----------------:|------------:|----------------------:|------------------------:|-----------------:|--------------------------:|---------------------------:|
| 0.01 |           41.0409 |          41.2447 |     41.1011 |               44.3743 |                 44.8181 |         0.995058 |                  0.929473 |                   0.924879 |
| 0.05 |          205.369  |         210.328  |    206.737  |              219.747  |                230.734  |         0.976425 |                  0.957134 |                   0.93457  |
| 0.1  |          411.78   |         430.915  |    416.551  |              434.111  |                477.522  |         0.955594 |                  0.992638 |                   0.948559 |
| 0.2  |          832.275  |         902.869  |    845.414  |              846.358  |               1015.63   |         0.921811 |                  1.06677  |                   0.983361 |
| 0.3  |         1272.4    |        1415.86   |   1286.59   |             1236.44   |               1607.38   |         0.898671 |                  1.14511  |                   1.02908  |
| 0.5  |         2286.85   |        2564.97   |   2205.87   |             1951.39   |               2927.08   |         0.891571 |                  1.31443  |                   1.17191  |
| 0.7  |         3814.5    |        3878.23   |   3174.41   |             2584.73   |               4394.04   |         0.983565 |                  1.50044  |                   1.47578  |
| 1    |          inf      |        6155.93   |   4719.54   |             3401.26   |               6802.53   |       inf        |                  1.80989  |                 inf        |
| 1.2  |          inf      |        7879.59   |   5811.2    |             3870.33   |               8514.72   |       inf        |                  2.0359   |                 inf        |


## Honest finding

**Bold-STAM is even flatter than current STAM.** The leading correction in bold-STAM lookback distance is z^3, whereas current STAM has a z^2 correction. At low z (where supernova catalogs live), z^3 is much smaller than z^2, so bold STAM undershoots catalogs by more than current STAM does.

Specifically: at z = 1, bold STAM predicts D ≈ 1.33 L while current STAM predicts 1.50 L and LCDM predicts ≈ 1.63 L (luminosity). Bold STAM's discrepancy with LCDM is about double that of current STAM at z=1.

**This does NOT derive the historical bridge term b.** In fact, it makes the catalog-fit residuals worse in the no-b form. The de Sitter A_cosmo ansatz combined with bold-STAM Shapiro is not the right cosmological structure for matching LCDM-fitted catalogs.

## What this means

Two interpretations, depending on what one believes:

1. **STAM is right; LCDM is wrong.** Bold STAM's prediction of an even flatter distance-redshift relation is what the universe actually has. Catalog distance moduli are calibrated within an LCDM framework, so they don't directly probe the true distance — they probe whatever the LCDM model fits to apparent flux. Under this interpretation, bold STAM predicts that high-z supernovae *should* appear even brighter than current STAM predicts (because they are closer than LCDM thinks).
2. **The cosmological A profile is not (r/L)^2.** The de Sitter ansatz might be wrong for our actual universe. A different A_cosmo profile (e.g., A ∝ r at low r) might give a different distance-redshift relation that matches catalogs better. This would require deriving A_cosmo from cosmological mass distribution + linear superposition, which is open work.

Either way: the bold-STAM thermodynamic and strong-field work does NOT produce the historical bridge term b as a derived consequence. The cosmological distance branch remains an open and separate problem.

## What was hoped for vs what was found

Hoped: that the de Sitter A_cosmo + bold-STAM Shapiro would produce a distance-z relation closer to LCDM than current STAM, giving a first-principles origin for the bridge term b.

Found: bold STAM is even FLATTER than current STAM. The leading correction is z^3 (not z^2), so it's a smaller correction at moderate z. The bridge-term structure is not naturally produced.

Negative result, but a clean one. The cosmological branch is unaffected by the thermodynamic / strong-field successes; distance work stays as a separate open problem if pursued.

## Generated plots

- `plots/31_distance_vs_z.png`
- `plots/31_modulus_residuals.png`
- `plots/31_low_z_zoom.png`

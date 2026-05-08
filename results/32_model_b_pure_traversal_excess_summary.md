# STAM Model-B Cosmological Distance — Pure Traversal Excess Only

## Purpose

Compute the cosmological distance prediction in STAM Model-B using ONLY the path-traversal excess of light through the cosmic A field. Clock dilation effects are kept SEPARATE (per author's 2026-05-07 sharpening). The full Shapiro effect lumps the two together; here we peel them apart and use only the path piece.

## Setup

- L = c/H_0 = 13387 Mly  (STAM Hubble length)
- Historical bridge term: b = 354.95 Mly (fitted, retained for comparison)
- Cosmic A profile (Model-B): `A_cosmo(r) ≈ (r/L)²` from critical-density Poisson
- Pure traversal excess: `TE(z) = ∫_0^Lz A_cosmo(r') dr' = L z³/3`
- Apparent distance: `d_Model-B(z) = L z + L z³/3`

## Why A_cosmo ≈ (r/L)² is a derived consequence, not a postulate

For a homogeneous critical-density universe, Poisson's equation
```text
grad² A = (8π G / c²) ρ_crit
```
applied to a sphere of radius r centered on observer gives
```text
A(r) ~ (8π G / 3c²) ρ_crit r² = r²/L²
```
using `ρ_crit = 3H₀²/(8πG)` and `L = c/H₀`. This is NOT the de Sitter ansatz from Q8 (which was postulated). It is a *derived* consequence of Model-B's matter-sourced framework under the observed critical-density assumption.

## Numerical results

|    z |   TE_analytical_Mly |   TE_numerical_Mly |   d_Model_B_Mly |   d_Model_A_no_b_Mly |   d_Model_A_bridge_Mly |   d_LCDM_DL_Mly |   Model_B_over_no_b |   Model_B_over_LCDM |
|-----:|--------------------:|-------------------:|----------------:|---------------------:|-----------------------:|----------------:|--------------------:|--------------------:|
| 0.01 |          0.00446236 |         0.00446237 |         133.875 |              134.54  |                138.09  |         146.196 |            0.995058 |            0.915722 |
| 0.05 |          0.557795   |         0.557796   |         669.912 |              686.088 |                703.836 |         752.656 |            0.976423 |            0.890065 |
| 0.1  |          4.46236    |         4.46237    |        1343.17  |             1405.64  |               1441.14  |        1557.68  |            0.955556 |            0.862291 |
| 0.2  |         35.6989     |        35.6989     |        2713.12  |             2945.16  |               3016.15  |        3312.98  |            0.921212 |            0.818934 |
| 0.3  |        120.484      |       120.484      |        4136.61  |             4618.55  |               4725.03  |        5243.27  |            0.895652 |            0.788938 |
| 0.5  |        557.795      |       557.796      |        7251.34  |             8366.93  |               8544.41  |        9548.14  |            0.866667 |            0.759451 |
| 0.7  |       1530.59       |      1530.59       |       10901.6   |            12650.8   |              12899.3   |       14333.4   |            0.861728 |            0.760571 |
| 1    |       4462.36       |      4462.37       |       17849.5   |            20080.6   |              20435.6   |       22189.8   |            0.888889 |            0.804398 |
| 1.2  |       7710.96       |      7710.97       |       23775.5   |            25703.2   |              26129.2   |       27775     |            0.925    |            0.856002 |


## Findings

- **Model-B pure TE = L z³/3.** Cubic in z, with no z² or z-linear term.
- **Model-A's TE = 0.35 L z² is a quadratic form.** Differs in functional shape.
- **Bridge term b z is z-linear.** Different again.
- All three differ in z-dependence, but at z ≈ 1 they happen to give comparable magnitudes (Model-B: 0.33L, Model-A: 0.35L, b: 0.027L).
- Model-B is FLATTER than Model-A no-b at low z (z³ < z² for z<1).
- All three are flatter than LCDM D_L. None reproduces LCDM.

## Honest verdict

**Pure traversal excess alone does not derive the historical bridge term b.** Model-B's matter-sourced A field at cosmological scales gives a z³ correction to the Hubble linear distance, not the z-linear b z form. The bridge term remains a fitting artifact when comparing STAM (any version) to LCDM-fitted catalogs.

What Model-B's pure TE does give us:
- A clean, derived prediction `D = Lz + Lz³/3` for cosmological distance under matter-sourced A and pure path-excess only.
- Confirmation that Model-B's distance prediction is FLATTER than Model-A's no-b form, which is in turn flatter than LCDM. The discrepancy with LCDM is real and widening, not narrowing.

If STAM is correct, this prediction needs to be tested against catalogs honestly:
- Either the catalogs (calibrated within LCDM) systematically overestimate supernova distances at high z, in which case Model-B's flatter prediction is closer to the truth.
- Or STAM's distance prediction is wrong, and the framework needs additional structure beyond pure path traversal excess.

**The F3-extended cosmology calculation (modified Friedmann with bold-STAM effective stress-energy) remains the next critical follow-up.** It might introduce additional dark-energy-like structure that moves the prediction closer to LCDM. This script does NOT include that effect.

## Generated plots

- `plots/32_model_b_distance_comparison.png`
- `plots/32_traversal_excess_only.png`
- `plots/32_distance_residuals_vs_LCDM.png`

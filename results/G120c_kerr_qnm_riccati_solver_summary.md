# G120c -- Riccati / log-derivative GR Kerr QNM solver

Log-derivative matching on the Teukolsky radial equation, with the log-derivative Y = R'/R bounded throughout integration.

## Settings
- eps_h = 1e-05
- r_max = 400.0
- r_match = 25.0
- horizon BC: leading order (2 - i σ_+) / eps_h
- infinity BC: leading order i ω + (3 + 2 i ω M) / r

## Stage 1: mismatch at omega_qnm
| a | omega_qnm | |Y_L - Y_R| | |Y_L| | |Y_R| |
|---:|---|---:|---:|---:|
| 0.00 | (0.21949044937597292+0.39789895122296853j) | 8.7529e-01 | 4.5442e-01 | 4.2487e-01 |
| 0.30 | (0.21769572066444065+0.4490764505797835j) | 9.6817e-01 | 4.9906e-01 | 4.7279e-01 |
| 0.50 | (0.21503998409375064+0.4985289426014869j) | 1.0590e+00 | 5.4293e-01 | 5.1944e-01 |
| 0.70 | (0.20929675651417265+0.5740690549396271j) | 1.1992e+00 | 6.1103e-01 | 5.9122e-01 |
| 0.90 | (0.19133632939375486+0.7265949992118796j) | 1.4864e+00 | 7.5137e-01 | 7.3758e-01 |

## Stage 2: Newton calibration
| a | omega_qnm | omega_solver | rel err | status |
|---:|---|---|---:|---|
| 0.00 | (0.37367168441804177-0.08896231568893546j) | (-0.11532019656417102+0.00884220809767389j) | 129.8247% | FAIL |
| 0.30 | (0.4195266817638516-0.08772927189431394j) | (-0.1531635819130917+0.009439909841916429j) | 135.5281% | FAIL |
| 0.50 | (0.46412302597593885-0.08563883498806352j) | (0.15330087158637865+0.00974899631914406j) | 68.8895% | FAIL |
| 0.70 | (0.5326002435510185-0.08079287315500769j) | (0.15299770137960056+0.009878762575035895j) | 72.4496% | FAIL |
| 0.90 | (0.6716142721321626-0.0648692358757954j) | (0.10623632495074346+0.009235824990416976j) | 84.5087% | FAIL |

## Verdict

FAIL  --  5 spins above 2%.
Check Riccati equation, BC implementation, r_max, r_match.

## Files
- [scripts/G120c_kerr_qnm_riccati_solver.py](../scripts/G120c_kerr_qnm_riccati_solver.py)
- CSV: `results/G120c_kerr_qnm_riccati_solver_table.csv`
- Plot: `plots/G120c_kerr_qnm_riccati_solver.png`

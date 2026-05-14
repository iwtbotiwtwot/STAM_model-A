# G120b -- standalone GR Kerr Teukolsky shooting solver

Direct shooting on the Teukolsky radial equation, calibrated
against G120a's qnm reference table.  No STAM.

## Settings
- eps_horizon = 1e-06
- r_match = 200.0
- Frobenius BC at horizon: leading order only (a_1 = 0)
- Asymptotic at infinity: leading order only (α = 0)

## Stage 1: sanity (|A_in/A_out| at omega_qnm)
| a | omega_qnm | |A_in| | |A_out| | |A_in/A_out| |
|---:|---|---:|---:|---:|
| 0.00 | (0.37367168441804177-0.08896231568893546j) | 2.6282e+21 | 6.9757e-01 | 3.7676e+21 |
| 0.30 | (0.4195266817638516-0.08772927189431394j) | 8.7062e+20 | 5.9597e-01 | 1.4608e+21 |
| 0.50 | (0.46412302597593885-0.08563883498806352j) | 1.9548e+20 | 4.9707e-01 | 3.9327e+20 |
| 0.70 | (0.5326002435510185-0.08079287315500769j) | 8.9922e+18 | 3.6125e-01 | 2.4892e+19 |
| 0.90 | (0.6716142721321626-0.0648692358757954j) | 3.5068e+15 | 1.5773e-01 | 2.2232e+16 |

## Stage 2: Newton calibration
| a | omega_qnm | omega_solver | rel err | status |
|---:|---|---|---:|---|
| 0.00 | (0.37367168441804177-0.08896231568893546j) | (0.37322944416764764+0.017546518462939593j) | 27.7286% | FAIL |
| 0.30 | (0.4195266817638516-0.08772927189431394j) | (0.4173771375446208+0.018092272542225403j) | 24.6951% | FAIL |
| 0.50 | (0.46412302597593885-0.08563883498806352j) | (0.47587676232367787+0.01800373714726506j) | 22.1009% | FAIL |
| 0.70 | (0.5326002435510185-0.08079287315500769j) | (0.53421965392042+0.01757838674647772j) | 18.2636% | FAIL |
| 0.90 | (0.6716142721321626-0.0648692358757954j) | (0.672274159081819+0.019834449788004877j) | 12.5539% | FAIL |

## Verdict

FAIL  --  5 spins above 2% error.
Iterate: increase r_match, decrease eps_h, add first-order
Frobenius BC and first-order 1/r asymptotic corrections.

## Files
- [scripts/G120b_gr_kerr_teukolsky_solver.py](../scripts/G120b_gr_kerr_teukolsky_solver.py)
- CSV: `results/G120b_gr_kerr_teukolsky_solver_table.csv`
- Plot: `plots/G120b_gr_kerr_teukolsky_solver.png`

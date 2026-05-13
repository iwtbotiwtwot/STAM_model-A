# G109 - Axial QNM convergence study

Convergence sweep of the analytic V_exact axial QNM with deepening throat.

## Sweep
| rs_min | N | eps_min | f_max | calib | omega_exact | shift |
|---:|---:|---:|---:|---:|---|---:|
| -300 | 8001 | 7.22e-04 | 9.46e+15 | 0.796% | (0.38837731849782914-0.08257234484597935j) | 4.9770% |
| -500 | 11001 | 4.29e-04 | 2.41e+25 | 0.796% | (0.38837952490553196-0.08257606280052045j) | 4.9773% |
| -1000 | 18001 | 2.13e-04 | 4.07e+48 | 0.798% | (0.3883912319297978-0.08259489894508167j) | 4.9792% |
| -2000 | 32001 | 1.06e-04 | 5.30e+94 | 0.796% | (0.3883785586974569-0.08257468007113407j) | 4.9772% |

## Verdict
CONVERGED.  Framework's rigorous axial QNM shift = 4.977%  (at rs_min = -2000).

## Reference values
- G92 / G93 reported (G92 finite-difference): 5.35%
- G94 robustness sweep mean (FD): 5.384% +/- 0.118%
- G108 analytic (single grid, rs_min = -300): 4.977%

## Files
- [scripts/G109_axial_convergence_study.py](../scripts/G109_axial_convergence_study.py)
- [plots/G109_axial_convergence.png](../plots/G109_axial_convergence.png)

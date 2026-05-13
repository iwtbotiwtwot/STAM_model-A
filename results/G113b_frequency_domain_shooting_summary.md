# G113b -- frequency-domain shooting QNM solver

Schwarzschild calibration of the numerical-shooting method.

Best grid: rs_min=-50.0, rs_max=100.0, N=6001

## Results
| ell | n | omega_solved | omega_Leaver | rel error | iter | status |
|---:|---:|---|---|---:|---:|---|
| 2 | 0 | (0.39320305540273215-0.03346559215208037j) | (0.373672-0.088962j) | 1.5316e-01 | 60 | FAIL |
| 2 | 1 | (0.3652723979051539-0.03156876606045974j) | (0.346711-0.273915j) | 5.5008e-01 | 60 | FAIL |
| 3 | 0 | (0.6226758077672265-0.03444124890248701j) | (0.599443-0.092703j) | 1.0341e-01 | 60 | FAIL |

## Status
PARTIAL  --  shooting method needs further tuning; see per-mode errors above.

## Files
- [scripts/G113b_frequency_domain_shooting.py](../scripts/G113b_frequency_domain_shooting.py)

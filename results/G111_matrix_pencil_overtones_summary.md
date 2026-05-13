# G111 - Higher-mode axial QNM via matrix pencil

Replacement for G110's two-mode curve_fit n=1 extraction (which failed).
Matrix pencil method extracts all dominant poles simultaneously via SVD
+ generalized eigenvalue; no local-minimum trap.

## Results
| ell | n | calib | omega_GR | Leaver gold | omega_exact | shift |
|---:|---:|---:|---|---|---|---:|
| 2 | 0 | 37.499% | (0.5084800204735982-0.038222400899226835j) | (0.373672-0.088962j) | (0.5118347050833365-0.03990663398186503j) | 0.7362% |
| 2 | 1 | 65.732% | (0.21537496976915396-0.014866461306359395j) | (0.346711-0.273915j) | (0.21442779553144184-0.014555598008008982j) | 0.4618% |
| 3 | 0 | 20.359% | (0.5169822454732316-0.0007776764952397302j) | (0.599443-0.092703j) | (0.5171052673334757-0.00039864133233840764j) | 0.0771% |
| 3 | 1 | 49.794% | (0.3934781629843617-0.020520832878413584j) | (0.582644-0.281298j) | (0.3962327978213393-0.021057205340535637j) | 0.7123% |
| 4 | 0 | 49.557% | (0.41304410624350996-0.016298558394051057j) | (0.809178-0.094164j) | (0.41268915196245726-0.016245793865918436j) | 0.0868% |
| 4 | 1 | 66.601% | (0.2908822346990312-0.036179083324371065j) | (0.796632-0.284334j) | (0.2896484225777107-0.03643695417743104j) | 0.4300% |

## Files
- [scripts/G111_matrix_pencil_overtones.py](../scripts/G111_matrix_pencil_overtones.py)
- [plots/G111_matrix_pencil_overtones.png](../plots/G111_matrix_pencil_overtones.png)

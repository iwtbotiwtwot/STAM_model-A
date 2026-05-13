# G114 -- Polar (even-parity) QNM from constrained shell-count action

First-pass polar prediction.  Uses V_Zerilli_STAM_proxy(M_eff) + (sqrt f)'' / sqrt f, by analogy with axial G108 result.

## Setup
- Two LM constraints algebraically eliminate delta Sigma in polar
  (delta lambda_1 pins partial_r delta Sigma; delta lambda_2 pins
  partial_t delta Sigma to delta u^r).  Scalar mode does not propagate.
- Effective polar master equation is Zerilli-type with f-modification.
- First-pass:  V_polar = V_Zerilli_STAM_proxy + (sqrt f)'' / sqrt f.

## Results
| ell | calib | omega_GR (Zerilli) | omega_polar STAM | shift |
|---:|---:|---|---|---:|
| 2 | 0.783% | (0.37081708376785244-0.08990951649618498j) | (0.39214947465758404-0.08773749652290519j) | 5.6197% |
| 3 | 1.306% | (0.5934773926056195-0.09791080709677492j) | (0.6067846229162843-0.0936249132412091j) | 2.3243% |
| 4 | 1.241% | (0.8072471642586169-0.10408733371323621j) | (0.816014585269741-0.09967247414908367j) | 1.2060% |

## Polar vs axial
| ell | axial shift | polar shift | polar / axial |
|---:|---:|---:|---:|
| 2 | 4.9770% | 5.6197% | 1.129 |
| 3 | 2.1010% | 2.3243% | 1.106 |
| 4 | 1.1580% | 1.2060% | 1.041 |

## Caveats
- M_eff substitution in Zerilli (polar analog of G87 axial proxy);
  rigorous polar reduction with k != h inside PS would replace this.
- (sqrt f)'' / sqrt f correction by analogy with axial; rigorous polar
  reduction may augment with constraint-feedback terms.

Follow-up:  G114-rigorous would derive V_polar directly by varying
S[Sigma,g,lambda_1,lambda_2] in even-parity Regge-Wheeler gauge,
eliminating delta Sigma via the two LM equations, and reading off
the Sturm-Liouville form.  This is the polar analog of G108.

## Files
- [scripts/G114_polar_QNM_from_constrained_sigma.py](../scripts/G114_polar_QNM_from_constrained_sigma.py)
- [plots/G114_polar_QNM.png](../plots/G114_polar_QNM.png)

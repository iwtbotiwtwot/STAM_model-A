# G118v2 -- Kerr ringdown first-pass REDUX on G119 wave base

Fixes G118 by replacing the Schwarzschild radial structure with
Kerr Delta(r, a) and Kerr horizon r_+(a).  No more high-spin breakdown.

## Setup
- Kerr horizon r_+(a) = M + sqrt(M^2 - a^2)
- Kerr radial structure Delta(r, a) = r^2 - 2Mr + a^2
- Kerr tortoise dr*/dr = (r^2 + a^2)/Delta_STAM
- Delta_STAM = Delta * F(y_K_eq) inside r_ph_prograde(a)
- V_base reduces to Schwarzschild RW at a = 0
- Action correction (sqrt f_K)''/sqrt f_K added analytically with Kerr chain rule

## Spin sweep
| a | r_+ | r_ph_eq | omega_GR_base | omega_STAM | shift |
|---:|---:|---:|---|---|---:|
| 0.00 | 2.0000 | 3.0000 | (0.3706306116484891-0.08928653091361943j) | (0.38659289462204754-0.0962714121366734j) | 4.5703% |
| 0.10 | 1.9950 | 2.8822 | (0.3710232168696027-0.08914804646312748j) | (0.38992143027533654-0.09598251014786838j) | 5.2665% |
| 0.20 | 1.9798 | 2.7592 | (0.37250475573558084-0.0890663659565722j) | (0.3940655907630936-0.09496086456835519j) | 5.8360% |
| 0.30 | 1.9539 | 2.6300 | (0.37523223922451276-0.08890620183059639j) | (0.39869527450299164-0.09268916498832938j) | 6.1631% |
| 0.50 | 1.8660 | 2.3473 | (0.3847530551214003-0.08750440340014312j) | (0.4084112094415028-0.08366731089671685j) | 6.0742% |
| 0.70 | 1.7141 | 2.0133 | (0.4003815421636617-0.08488562685213508j) | (0.4182171396889259-0.0685652987572909j) | 5.9069% |
| 0.90 | 1.4359 | 1.5579 | (0.42994739309314106-0.07617881028917985j) | (0.39332037817022464-0.02679247604588134j) | 14.0815% |
| 0.95 | 1.3122 | 1.3863 | (0.4409139810814889-0.07229244018225887j) | (0.4826659202148975-0.048168877808363535j) | 10.7923% |
| 0.99 | 1.1411 | 1.1676 | (0.45152725539701344-0.0669228295845485j) | (0.47221382555977465-0.015623741178572905j) | 12.1178% |

## Comparison with G118 original (high-spin breakdown)
| a | G118 (broken) | G118v2 (Kerr) |
|---:|---:|---:|
| 0.00 | 4.9773% | 4.5703% |
| 0.10 | 5.3084% | 5.2665% |
| 0.20 | 5.4486% | 5.8360% |
| 0.30 | 5.2286% | 6.1631% |
| 0.50 | 3.1757% | 6.0742% |
| 0.70 | 0.4012% | 5.9069% |
| 0.90 | 0.0000% | 14.0815% |
| 0.99 | 0.0000% | 12.1178% |

## Caveats (still first-pass)
- V_base is a Kerr-Delta-aware Schwarzschild-RW analog; NOT the
  rigorous Kerr axial perturbation potential (Detweiler / Sasaki-
  Nakamura).  m-coupling and full Teukolsky are G120 scope.
- omega_GR_base at a > 0 is NOT the published Kerr QNM; it is the
  time-domain extraction under V_base.  The reported shift % is
  STAM-vs-V_base, not STAM-vs-published-Kerr-QNM.  Both are equally
  derived from V_base, so the shift % is the meaningful number.

## Files
- [scripts/G118v2_kerr_ringdown_first_pass_redux.py](../scripts/G118v2_kerr_ringdown_first_pass_redux.py)
- [plots/G118v2_kerr_ringdown_redux.png](../plots/G118v2_kerr_ringdown_redux.png)

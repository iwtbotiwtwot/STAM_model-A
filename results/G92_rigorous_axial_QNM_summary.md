# G92 - Rigorous axial perturbation equation from constrained shell-count action

**Date: 2026-05-13.**

## Derivation result

Odd parity does not perturb the shell-count field: `delta Sigma = 0`. The two LM constraints remove the scalar shell-count mode, so the propagating axial sector is the tensor odd mode of `f(Sigma) R` on the committed background.

The canonical axial equation used here is:

```text

d_t^2 Psi - d_*^2 Psi + V_ax Psi = 0

d/dr* = sqrt(h k) d/dr

V_ax = V_geom + 1/2 d_*P + 1/4 P^2

P = d_* ln f

V_geom = h [ ell(ell+1)/r^2 - 2(1-k)/r^2 - k'/(2r) - k h'/(2hr) ]

```

For Schwarzschild, `f = 1` and `k = h`, so this reduces exactly to the standard Regge-Wheeler potential.


## Schwarzschild calibration

- omega_TD     = 0.370633 - 0.089287 i

- omega_Leaver = 0.373672 - 0.088962 i

- **|Delta omega / omega| = 0.7957 %**


## STAM rigorous axial result

- omega_STAM = 0.388389 - 0.082590 i

- **|Delta omega / omega| vs calibrated Schwarzschild = 4.9778 %**

- Re shift: +4.7908 %

- Im shift: -7.5003 %


## Method

- GR grid: N = 6001, r* in [-60.0, 400.0].

- STAM grid: N = 10001, r* in [-500.0, 400.0], eps_min = 4.2886e-04.

- Same Sommerfeld time-domain pipeline as G91.

- Fit window: t in [100.0, 250.0] M.


## Files

- [scripts/G92_rigorous_axial_QNM.py](../scripts/G92_rigorous_axial_QNM.py)

- [plots/G92_rigorous_axial_QNM.png](../plots/G92_rigorous_axial_QNM.png)

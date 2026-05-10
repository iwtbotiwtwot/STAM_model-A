# G20: V_3 Eigenmode Spectrum on (A_0, 1)

## Setup

Schrödinger-like wave equation for A as a coordinate on the bounded interval (A_0, 1) with V_3 potential and Dirichlet boundary conditions.

```text
-(1/2) ψ''(A) + V_3(A) ψ(A) = E ψ(A)
ψ(A_0) = ψ(1) = 0     (hard walls at structural boundaries)
V_3(A) = α/A + β/(1-A) with α/β = [A_0/(1-A_0)]²
A_0 = 1/(12π), β = 1 (sets scale)
```

Numerical solution: discretize on N=1200 grid points, build tridiagonal Hamiltonian, diagonalize.

## Eigenvalue spectrum (lowest 10)

```text
  n             E_n   E_n / E_0
  0        7.614888      1.0000
  1       23.972658      3.1481
  2       50.443528      6.6243
  3       87.200357     11.4513
  4      134.301223     17.6367
  5      191.771946     25.1838
  6      259.626004     34.0945
  7      337.871155     44.3698
  8      426.512131     56.0103
  9      525.551902     69.0164
```

## Zero crossings

Interior zero crossings of low-order modes:


```text
  Mode 0: no interior crossings
  Mode 1: 0.4938
  Mode 2: 0.3432, 0.6625
  Mode 3: 0.2660, 0.5062, 0.7479
  Mode 4: 0.2190, 0.4118, 0.6049, 0.7991
  Mode 5: 0.1874, 0.3484, 0.5096, 0.6710, 0.8332
  Mode 6: 0.1647, 0.3029, 0.4412, 0.5796, 0.7182, 0.8574
  Mode 7: 0.1476, 0.2687, 0.3898, 0.5110, 0.6322, 0.7537, 0.8755
```

Compare to orbital thirds: A = 1/3 ≈ 0.3333, A = 2/3 ≈ 0.6667.


**Mode 3 crossings:** A = 0.2660 and A = 0.5062


- vs orbital ISCO (A=1/3): +20.19% deviation

- vs orbital photon sphere (A=2/3): +24.07% deviation


## Generated plot

- `plots/G20_v3_eigenmode_spectrum.png`

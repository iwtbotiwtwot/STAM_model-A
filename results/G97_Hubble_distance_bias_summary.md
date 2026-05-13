# G97 — Hubble tension as STAM photon-A distance bias

## Interpretation

This script tests the idea that the Hubble tension may be partly a **distance-inference tension** rather than a pure expansion-rate disagreement.

- Layer 1: intrinsic expansion \(H(z)=H_0E(z)\).
- Layer 2: photon-A traversal bias adds apparent distance to photon observables.
- Direct \(H(z)\) probes such as cosmic chronometers should respond mainly to Layer 1.

## Locked constants

- \(A_0=1/(12\pi)=0.026525823849\)
- \(H_0^{true}=73.0400\) km/s/Mpc
- \(H_0^{target}=67.4000\) km/s/Mpc
- Bias model: `linear`
- Bridge term \(b=A_0c/H_0=108.875163\) Mpc = `355.103289` Mly

## Low-z analytic estimate

For a linear low-z bias, \(D_{obs}\approx(c/H_0)z(1+f_{los}A_0)\), so a no-bias fit infers:

\[
H_0^{inferred}\approx\frac{H_0^{true}}{1+f_{los}A_0}.
\]

Solving for \(f_{los}\):

\[
f_{los}=\frac{H_0^{true}/H_0^{target}-1}{A_0}.
\]

- Required \(f_{los}\) = `3.154644`
- Nearest grid \(f_{los}\) = `1.000`
- H0 inferred from low-z slope = `67.4207`
- H0 inferred from full distance-modulus fit = `71.7392`

## Caveat

This is a scaffold, not a publication-grade cosmology pipeline. A real test must include SN covariance, BAO likelihoods, CMB acoustic angle, chronometers, structure growth, nuisance parameters, and locked priors. The point of G97 is to separate Layer-1 expansion from Layer-2 photon-distance bias and quantify the required scale.

## Files

- Grid CSV: `/mnt/data/G97_run/results/G97_Hubble_distance_bias_grid.csv`
- Plot: `/mnt/data/G97_run/plots/G97_Hubble_inferred_vs_flos.png`
- Plot: `/mnt/data/G97_run/plots/G97_distance_modulus_residuals.png`
- Plot: `/mnt/data/G97_run/plots/G97_Hz_direct_probe_unchanged.png`

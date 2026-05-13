#!/usr/bin/env python3
"""
G98_distance_bias_units_audit.py

Fast audit of STAM / Model-A Hubble-tension distance-bias bookkeeping.

Why this exists
===============
G97 exposed a bookkeeping mismatch:

    analytic low-z estimate: f_los_required ≈ 3.1546
    G97 printed nearest grid: f_los ≈ 1.0

Those cannot both be right under the same distance convention.

G98 separates the conventions:
    1. bias added directly to luminosity distance D_L
    2. bias added to comoving/traversal distance D_C and then multiplied by (1+z)
    3. pure low-z Hubble-law distance
    4. SU-shaped versions of 1 and 2

Main conclusion
===============
The analytic relation

    H0_inferred = H0_true / (1 + f_los * A0)

is the correct low-z result when the bias term f_los*b*z is added to the same
distance slope being used to infer H0.

For STAM, the recommended convention is physically:
    D_C_obs(z) = D_C_intrinsic(z) + f_los * b * shape(z)
    D_L_obs(z) = (1+z) * D_C_obs(z)

That makes photon-A traversal a path/comoving correction first, then applies
the usual luminosity-distance factor.

Outputs
=======
results/G98_distance_bias_units_audit_summary.md
results/G98_distance_bias_units_grid.csv
plots/G98_H0_inference_conventions.png
plots/G98_lowz_residuals_conventions.png
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

c_km_s = 299_792.458
MPC_TO_MLY = 3.261563776
A0 = 1.0 / (12.0 * math.pi)


def E_flat_lcdm(z, omega_m):
    return np.sqrt(omega_m * (1.0 + z)**3 + (1.0 - omega_m))


def cumulative_comoving_distance(z_grid, H0, omega_m):
    """Fast trapezoid cumulative integral D_C(z)=∫ c/H dz over sorted z_grid."""
    z = np.asarray(z_grid, dtype=float)
    # Include z=0 for integration anchor.
    z_full = np.concatenate([[0.0], z])
    invH = c_km_s / (H0 * E_flat_lcdm(z_full, omega_m))
    dz = np.diff(z_full)
    trap = 0.5 * (invH[1:] + invH[:-1]) * dz
    return np.cumsum(trap)


def bridge_b_mpc(H0):
    return A0 * c_km_s / H0


def bias_shape(z, model):
    z = np.asarray(z, dtype=float)
    if model == "linear":
        return z
    if model == "stam_su":
        return z * (1.0 + 3.0*z/20.0)
    if model == "saturating":
        return z / (1.0 + z)
    raise ValueError(model)


def infer_H0_slope(z, D):
    slope = np.sum(z * D) / np.sum(z*z)
    return c_km_s / slope


def D_obs_from_convention(z, D_C, H0, f_los, convention):
    b = bridge_b_mpc(H0)
    D_L = (1.0 + z) * D_C
    if convention == "bias_to_luminosity":
        return D_L + f_los * b * bias_shape(z, "linear")
    if convention == "bias_to_comoving_then_luminosity":
        return (1.0 + z) * (D_C + f_los * b * bias_shape(z, "linear"))
    if convention == "bias_to_hubble_law_distance":
        return (c_km_s / H0) * z + f_los * b * bias_shape(z, "linear")
    if convention == "stam_su_to_luminosity":
        return D_L + f_los * b * bias_shape(z, "stam_su")
    if convention == "stam_su_to_comoving_then_luminosity":
        return (1.0 + z) * (D_C + f_los * b * bias_shape(z, "stam_su"))
    raise ValueError(convention)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="/mnt/data/G98_run")
    parser.add_argument("--H0-true", type=float, default=73.04)
    parser.add_argument("--H0-target", type=float, default=67.40)
    parser.add_argument("--omega-m", type=float, default=0.30)
    args = parser.parse_args()

    root = Path(args.outdir).resolve()
    results = root / "results"
    plots = root / "plots"
    results.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)

    H0_true = args.H0_true
    H0_target = args.H0_target
    omega_m = args.omega_m
    b = bridge_b_mpc(H0_true)
    f_required = (H0_true/H0_target - 1.0) / A0

    conventions = [
        "bias_to_luminosity",
        "bias_to_comoving_then_luminosity",
        "bias_to_hubble_law_distance",
        "stam_su_to_luminosity",
        "stam_su_to_comoving_then_luminosity",
    ]

    grids = {
        "ultra_low_z_0p001_0p01": np.linspace(0.001, 0.010, 80),
        "low_z_0p005_0p10": np.linspace(0.005, 0.100, 120),
        "sne_like_0p01_2p0": np.linspace(0.010, 2.000, 300),
    }

    f_values = np.linspace(0.0, 6.0, 301)
    rows = []

    for grid_name, z in grids.items():
        D_C = cumulative_comoving_distance(z, H0_true, omega_m)
        for conv in conventions:
            for f in [0.0, 1.0, f_required]:
                D = D_obs_from_convention(z, D_C, H0_true, f, conv)
                rows.append({
                    "grid": grid_name,
                    "convention": conv,
                    "f_los": f,
                    "H0_slope": infer_H0_slope(z, D),
                    "note": "zero" if f == 0 else ("f=1" if abs(f-1)<1e-12 else "analytic_required")
                })

            # Numerical f required by slope
            Hs = []
            for f in f_values:
                D = D_obs_from_convention(z, D_C, H0_true, f, conv)
                Hs.append(infer_H0_slope(z, D))
            Hs = np.array(Hs)
            idx = int(np.argmin(np.abs(Hs - H0_target)))
            rows.append({
                "grid": grid_name,
                "convention": conv,
                "f_los": f_values[idx],
                "H0_slope": Hs[idx],
                "note": "numerical_required_slope"
            })

    csv_path = results / "G98_distance_bias_units_grid.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["grid", "convention", "f_los", "H0_slope", "note"])
        writer.writeheader()
        writer.writerows(rows)

    # Console compact table
    print("="*96)
    print("G98: distance-bias units / definition audit")
    print("="*96)
    print()
    print(f"A0 = {A0:.12f}")
    print(f"H0_true   = {H0_true:.4f}")
    print(f"H0_target = {H0_target:.4f}")
    print(f"b = A0*c/H0_true = {b:.6f} Mpc = {b*MPC_TO_MLY:.6f} Mly")
    print(f"Analytic f_los_required = {f_required:.6f}")
    print(f"Analytic H0(f=1) = {H0_true/(1+A0):.6f}")
    print(f"Analytic H0(f_required) = {H0_true/(1+f_required*A0):.6f}")
    print()
    print("Compact low-z audit:")
    print(f"{'convention':<40}{'H0 f=1':>12}{'H0 f_req':>12}{'num f_req':>12}")
    print("-"*76)
    for conv in conventions:
        subset = [r for r in rows if r["grid"]=="low_z_0p005_0p10" and r["convention"]==conv]
        f1 = [r for r in subset if r["note"]=="f=1"][0]
        fr = [r for r in subset if r["note"]=="analytic_required"][0]
        nr = [r for r in subset if r["note"]=="numerical_required_slope"][0]
        print(f"{conv:<40}{f1['H0_slope']:>12.4f}{fr['H0_slope']:>12.4f}{nr['f_los']:>12.3f}")

    # Plots
    z = grids["low_z_0p005_0p10"]
    D_C = cumulative_comoving_distance(z, H0_true, omega_m)

    plt.figure(figsize=(10,6))
    for conv in conventions:
        Hs = []
        for f in f_values:
            D = D_obs_from_convention(z, D_C, H0_true, f, conv)
            Hs.append(infer_H0_slope(z, D))
        plt.plot(f_values, Hs, label=conv)
    plt.axhline(H0_target, color="black", linestyle=":", label="target H0")
    plt.axhline(H0_true, color="gray", linestyle="--", label="true H0")
    plt.axvline(f_required, color="black", linestyle="-.", label=f"analytic f={f_required:.2f}")
    plt.xlabel("f_los")
    plt.ylabel("H0 inferred from low-z distance slope")
    plt.title("G98: H0 inference depends on distance-bias convention")
    plt.grid(alpha=0.3)
    plt.legend(fontsize=8)
    p1 = plots / "G98_H0_inference_conventions.png"
    plt.tight_layout()
    plt.savefig(p1, dpi=200)
    plt.close()

    z2 = np.linspace(0.001, 0.20, 250)
    D_C2 = cumulative_comoving_distance(z2, H0_true, omega_m)
    D_base = (1+z2)*D_C2
    plt.figure(figsize=(10,6))
    for conv in conventions:
        D = D_obs_from_convention(z2, D_C2, H0_true, f_required, conv)
        plt.plot(z2, (D-D_base)/z2, label=conv)
    plt.axhline(f_required*b, color="black", linestyle=":", label="expected f*b low-z slope")
    plt.xlabel("z")
    plt.ylabel("(D_obs - D_base)/z [Mpc]")
    plt.title("G98: low-z bias residual slope by convention")
    plt.grid(alpha=0.3)
    plt.legend(fontsize=8)
    p2 = plots / "G98_lowz_residuals_conventions.png"
    plt.tight_layout()
    plt.savefig(p2, dpi=200)
    plt.close()

    summary = results / "G98_distance_bias_units_audit_summary.md"
    md = f"""# G98 — Distance-bias units / definition audit

## Why G98 exists

G97 exposed a bookkeeping mismatch: the analytic low-z estimate gave `f_los_required ≈ 3.15`, while the numerical printout appeared to give `f_los ≈ 1`.

## Locked constants

- `A0 = {A0:.12f}`
- `H0_true = {H0_true:.4f}` km/s/Mpc
- `H0_target = {H0_target:.4f}` km/s/Mpc
- `b = A0*c/H0_true = {b:.6f}` Mpc = `{b*MPC_TO_MLY:.6f}` Mly

## Correct low-z relation

If the photon-A bias adds to the same low-z fitted distance slope,

```text
D_obs ≈ (c/H0_true) z + f_los * (A0*c/H0_true) z
      ≈ (c/H0_true) z * (1 + f_los*A0)

H0_inferred ≈ H0_true / (1 + f_los*A0)
f_los_required = (H0_true/H0_target - 1) / A0
```

- `f_los_required = {f_required:.6f}`
- `H0(f=1) = {H0_true/(1+A0):.6f}` km/s/Mpc
- `H0(f_required) = {H0_true/(1+f_required*A0):.6f}` km/s/Mpc

## Main conclusion

The analytic value `f_los ≈ 3.15` is the correct low-z requirement if the bias term is `f_los*b*z` applied to the low-z fitted distance slope. If a script appears to return `f_los≈1` for the same target, it is mixing conventions or reporting the wrong row.

## Recommended convention for STAM going forward

Treat photon-A traversal as a path/comoving correction first:

```text
Intrinsic expansion distance:
    D_C(z) = ∫ c/H(z) dz

Photon-A traversal bias:
    ΔD_C(z) = f_los * b * shape(z)

Observed luminosity distance:
    D_L_obs(z) = (1+z) * [D_C(z) + ΔD_C(z)]
```

This makes the distance-bias layer explicit and keeps it separate from the expansion layer.

## Files

- CSV: `{csv_path}`
- Plot: `{p1}`
- Plot: `{p2}`
"""
    summary.write_text(md, encoding="utf-8")

    print()
    print(f"CSV written: {csv_path}")
    print(f"Summary written: {summary}")
    print(f"Plot written: {p1}")
    print(f"Plot written: {p2}")


if __name__ == "__main__":
    main()

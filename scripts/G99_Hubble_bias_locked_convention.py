#!/usr/bin/env python3
"""
G99_Hubble_bias_locked_convention.py

STAM / Model-A Hubble-tension distance-bias test using the G98-approved
distance convention.

G98 result
==========
Photon-A traversal should be treated as a path / comoving-distance correction:

    D_C,obs(z) = D_C,intrinsic(z) + ΔD_C(z)

with:

    ΔD_C(z) = f_los * b * shape(z)

and then converted to luminosity distance in the standard way:

    D_L,obs(z) = (1+z) * D_C,obs(z)

where:

    A0 = 1/(12π)
    b  = A0 * c / H0_true

Purpose
=======
G99 does NOT claim to solve the Hubble tension. It is a locked-convention
scaffold that asks:

    Can one A0-based photon-A distance-bias layer make distance probes
    infer a lower H0 while direct H(z) probes remain tied to the intrinsic
    expansion?

It separates:

    Layer 1: intrinsic expansion H(z)
    Layer 2: photon-A traversal / distance bias

Default model
=============
- H0_true = 73.04 km/s/Mpc
- H0_target = 67.40 km/s/Mpc
- flat LCDM-like E(z), Omega_m = 0.30
- A0 = 1/(12π)
- b = A0*c/H0_true ≈ 108.875 Mpc ≈ 355.103 Mly

Bias shapes
===========
linear:
    shape(z)=z

stam_su:
    shape(z)=z*(1+3z/20)

saturating:
    shape(z)=z/(1+z)

Outputs
=======
results/G99_Hubble_bias_locked_summary.md
results/G99_Hubble_bias_locked_grid.csv
plots/G99_inferred_H0_distance_vs_flos.png
plots/G99_mu_residuals_locked_convention.png
plots/G99_probe_separation.png
plots/G99_bias_shapes.png

Optional data
=============
SN-like CSV:
    --sn-csv path.csv
Columns:
    z,mu
optional:
    sigma_mu

Chronometer-like CSV:
    --hz-csv path.csv
Columns:
    z,H
optional:
    sigma_H

These are simple fits only. A publication-grade cosmology pipeline still needs
covariance matrices, BAO likelihoods, CMB theta*, nuisance parameters, and
locked priors.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar


# Constants
c_km_s = 299_792.458
MPC_TO_MLY = 3.261563776
A0 = 1.0 / (12.0 * math.pi)


def E_flat_lcdm(z, omega_m):
    z = np.asarray(z, dtype=float)
    return np.sqrt(omega_m * (1.0 + z)**3 + (1.0 - omega_m))


def H_z(z, H0, omega_m):
    return H0 * E_flat_lcdm(z, omega_m)


def cumulative_comoving_distance(z_grid, H0, omega_m):
    """Fast trapezoid D_C(z) on sorted z_grid."""
    z = np.asarray(z_grid, dtype=float)
    z_full = np.concatenate([[0.0], z])
    H = H_z(z_full, H0, omega_m)
    integrand = c_km_s / H
    dz = np.diff(z_full)
    trap = 0.5 * (integrand[1:] + integrand[:-1]) * dz
    return np.cumsum(trap)


def bridge_b_mpc(H0_true):
    return A0 * c_km_s / H0_true


def bias_shape(z, model):
    z = np.asarray(z, dtype=float)
    if model == "linear":
        return z
    if model == "stam_su":
        return z * (1.0 + 3.0*z/20.0)
    if model == "saturating":
        return z / (1.0 + z)
    raise ValueError(f"Unknown bias model: {model}")


def DC_obs_locked(z, DC_intrinsic, H0_true, f_los, shape_model):
    b = bridge_b_mpc(H0_true)
    return DC_intrinsic + f_los * b * bias_shape(z, shape_model)


def DL_from_DC(z, DC):
    return (1.0 + z) * DC


def distance_modulus(DL_mpc):
    DL_mpc = np.asarray(DL_mpc, dtype=float)
    return 5.0 * np.log10(DL_mpc) + 25.0


def infer_H0_lowz_slope(z, DL_obs):
    """At low z, D_L ≈ (c/H0)z. Fit through origin."""
    slope = np.sum(z * DL_obs) / np.sum(z*z)
    return c_km_s / slope


def infer_H0_mu_fit(z, mu_obs, omega_m, H0_bounds=(50, 90)):
    """Fit no-bias flat LCDM H0 to observed distance moduli at fixed omega_m."""
    z = np.asarray(z, dtype=float)
    mu_obs = np.asarray(mu_obs, dtype=float)

    def loss(H0):
        DC = cumulative_comoving_distance(z, H0, omega_m)
        mu_pred = distance_modulus(DL_from_DC(z, DC))
        return float(np.mean((mu_pred - mu_obs)**2))

    res = minimize_scalar(loss, bounds=H0_bounds, method="bounded", options={"xatol": 1e-9})
    return float(res.x), float(res.fun)


def read_csv_rows(path):
    with Path(path).open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fit_sn_csv(path, H0_true, omega_m, shape_model):
    rows = read_csv_rows(path)
    z = np.array([float(r["z"]) for r in rows])
    mu = np.array([float(r["mu"]) for r in rows])
    if "sigma_mu" in rows[0] and rows[0]["sigma_mu"] != "":
        sig = np.array([float(r.get("sigma_mu", 1.0) or 1.0) for r in rows])
    else:
        sig = np.ones_like(mu)

    DC_true = cumulative_comoving_distance(z, H0_true, omega_m)

    def chi2_f(f_los):
        DC_o = DC_obs_locked(z, DC_true, H0_true, f_los, shape_model)
        mu_pred = distance_modulus(DL_from_DC(z, DC_o))
        return float(np.sum(((mu - mu_pred) / sig)**2))

    res = minimize_scalar(chi2_f, bounds=(-20, 20), method="bounded")
    return {"f_los_best": float(res.x), "chi2": float(res.fun), "N": len(z)}


def fit_hz_csv(path, omega_m):
    rows = read_csv_rows(path)
    z = np.array([float(r["z"]) for r in rows])
    H = np.array([float(r["H"]) for r in rows])
    if "sigma_H" in rows[0] and rows[0]["sigma_H"] != "":
        sig = np.array([float(r.get("sigma_H", 1.0) or 1.0) for r in rows])
    else:
        sig = np.ones_like(H)

    def chi2_H0(H0):
        pred = H_z(z, H0, omega_m)
        return float(np.sum(((H - pred) / sig)**2))

    res = minimize_scalar(chi2_H0, bounds=(50, 90), method="bounded")
    return {"H0_best": float(res.x), "chi2": float(res.fun), "N": len(z)}


def main():
    parser = argparse.ArgumentParser(description="G99 locked-convention STAM Hubble distance-bias scaffold")
    parser.add_argument("--outdir", default="/mnt/data/G99_run")
    parser.add_argument("--H0-true", type=float, default=73.04)
    parser.add_argument("--H0-target", type=float, default=67.40)
    parser.add_argument("--omega-m", type=float, default=0.30)
    parser.add_argument("--shape", choices=["linear", "stam_su", "saturating"], default="linear")
    parser.add_argument("--zmax", type=float, default=2.0)
    parser.add_argument("--sn-csv", default=None)
    parser.add_argument("--hz-csv", default=None)
    args = parser.parse_args()

    root = Path(args.outdir).resolve()
    results = root / "results"
    plots = root / "plots"
    results.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)

    H0_true = args.H0_true
    H0_target = args.H0_target
    omega_m = args.omega_m
    shape_model = args.shape
    b = bridge_b_mpc(H0_true)

    # Analytic low-z required f_los for linear low-z behavior.
    # Note all three shapes are linear in z at very low z, so same leading estimate.
    f_required_analytic = (H0_true / H0_target - 1.0) / A0

    # Grids
    z_low = np.linspace(0.005, 0.10, 120)
    z_full = np.linspace(0.01, args.zmax, 300)
    z_plot = np.linspace(0.001, args.zmax, 400)
    f_values = np.linspace(0.0, 6.0, 301)

    rows = []
    for shape in ["linear", "stam_su", "saturating"]:
        DC_low = cumulative_comoving_distance(z_low, H0_true, omega_m)
        DC_full = cumulative_comoving_distance(z_full, H0_true, omega_m)

        H_low_values = []
        H_mu_values = []
        for f in f_values:
            DC_low_obs = DC_obs_locked(z_low, DC_low, H0_true, f, shape)
            DL_low_obs = DL_from_DC(z_low, DC_low_obs)
            H_low = infer_H0_lowz_slope(z_low, DL_low_obs)

            DC_full_obs = DC_obs_locked(z_full, DC_full, H0_true, f, shape)
            DL_full_obs = DL_from_DC(z_full, DC_full_obs)
            mu_full_obs = distance_modulus(DL_full_obs)
            H_mu, mu_loss = infer_H0_mu_fit(z_full, mu_full_obs, omega_m)

            H_low_values.append(H_low)
            H_mu_values.append(H_mu)
            rows.append({
                "shape": shape,
                "f_los": f,
                "H0_inferred_lowz_slope": H_low,
                f"H0_inferred_mu_fit_zmax_{args.zmax}": H_mu,
                "mu_loss": mu_loss,
            })

        H_low_values = np.array(H_low_values)
        H_mu_values = np.array(H_mu_values)
        idx_low = int(np.argmin(np.abs(H_low_values - H0_target)))
        idx_mu = int(np.argmin(np.abs(H_mu_values - H0_target)))

        # special rows
        for idx, note in [(idx_low, "closest_lowz_to_target"), (idx_mu, "closest_mufit_to_target")]:
            rows.append({
                "shape": shape,
                "f_los": f_values[idx],
                "H0_inferred_lowz_slope": H_low_values[idx],
                f"H0_inferred_mu_fit_zmax_{args.zmax}": H_mu_values[idx],
                "mu_loss": "",
                "note": note,
            })

    # CSV
    csv_path = results / "G99_Hubble_bias_locked_grid.csv"
    all_fields = sorted(set().union(*(r.keys() for r in rows)))
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=all_fields)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    # Specific model values for selected shape
    DC_low_sel = cumulative_comoving_distance(z_low, H0_true, omega_m)
    DC_full_sel = cumulative_comoving_distance(z_full, H0_true, omega_m)
    selected = {}
    for f in [0.0, 1.0, f_required_analytic]:
        DL_low = DL_from_DC(z_low, DC_obs_locked(z_low, DC_low_sel, H0_true, f, shape_model))
        H_low = infer_H0_lowz_slope(z_low, DL_low)
        DL_full = DL_from_DC(z_full, DC_obs_locked(z_full, DC_full_sel, H0_true, f, shape_model))
        mu_full = distance_modulus(DL_full)
        H_mu, loss = infer_H0_mu_fit(z_full, mu_full, omega_m)
        selected[f] = {"H_low": H_low, "H_mu": H_mu, "loss": loss}

    # Optional CSV fits
    sn_fit = fit_sn_csv(args.sn_csv, H0_true, omega_m, shape_model) if args.sn_csv else None
    hz_fit = fit_hz_csv(args.hz_csv, omega_m) if args.hz_csv else None

    # Print report
    print("="*100)
    print("G99: locked-convention STAM Hubble distance-bias scaffold")
    print("="*100)
    print()
    print("Convention:")
    print("  D_C_obs(z) = D_C_intrinsic(z) + f_los*b*shape(z)")
    print("  D_L_obs(z) = (1+z)*D_C_obs(z)")
    print()
    print(f"A0 = {A0:.12f}")
    print(f"H0_true   = {H0_true:.4f}")
    print(f"H0_target = {H0_target:.4f}")
    print(f"Omega_m   = {omega_m:.4f}")
    print(f"shape     = {shape_model}")
    print(f"b = {b:.6f} Mpc = {b*MPC_TO_MLY:.6f} Mly")
    print(f"analytic low-z f_los_required = {f_required_analytic:.6f}")
    print()
    print(f"{'f_los':>12}{'H0 low-z slope':>18}{'H0 full mu fit':>18}")
    print("-"*52)
    for f, vals in selected.items():
        print(f"{f:>12.6f}{vals['H_low']:>18.4f}{vals['H_mu']:>18.4f}")
    print()

    # plots
    plt.figure(figsize=(10, 6))
    for shape in ["linear", "stam_su", "saturating"]:
        subset = [r for r in rows if r.get("shape")==shape and "note" not in r]
        farr = np.array([float(r["f_los"]) for r in subset])
        Harr = np.array([float(r["H0_inferred_lowz_slope"]) for r in subset])
        plt.plot(farr, Harr, label=f"{shape}: low-z slope")
    plt.axhline(H0_true, color="gray", linestyle="--", label="true H0")
    plt.axhline(H0_target, color="black", linestyle=":", label="target H0")
    plt.axvline(f_required_analytic, color="black", linestyle="-.", label=f"analytic f={f_required_analytic:.2f}")
    plt.xlabel("f_los")
    plt.ylabel("H0 inferred by no-bias low-z distance fit")
    plt.title("G99 locked convention: distance-inferred H0 vs f_los")
    plt.grid(alpha=0.3)
    plt.legend(fontsize=8)
    p1 = plots / "G99_inferred_H0_distance_vs_flos.png"
    plt.tight_layout()
    plt.savefig(p1, dpi=200)
    plt.close()

    # residuals for selected shape at f_required
    DC_plot = cumulative_comoving_distance(z_plot, H0_true, omega_m)
    DL_base = DL_from_DC(z_plot, DC_plot)
    DL_bias = DL_from_DC(z_plot, DC_obs_locked(z_plot, DC_plot, H0_true, f_required_analytic, shape_model))
    mu_base = distance_modulus(DL_base)
    mu_bias = distance_modulus(DL_bias)

    plt.figure(figsize=(10, 6))
    plt.plot(z_plot, mu_bias - mu_base, label=f"bias residual, {shape_model}, f={f_required_analytic:.2f}")
    plt.xlabel("z")
    plt.ylabel("Δμ relative to intrinsic H0 distance")
    plt.title("G99 STAM photon-A distance-modulus residual")
    plt.grid(alpha=0.3)
    plt.legend()
    p2 = plots / "G99_mu_residuals_locked_convention.png"
    plt.tight_layout()
    plt.savefig(p2, dpi=200)
    plt.close()

    # Probe separation H(z)
    zH = np.linspace(0, 2, 250)
    plt.figure(figsize=(10, 6))
    plt.plot(zH, H_z(zH, H0_true, omega_m), label=f"Layer 1 intrinsic H(z), H0={H0_true:.2f}")
    plt.plot(zH, H_z(zH, H0_target, omega_m), label=f"comparison lower H0={H0_target:.2f}")
    plt.xlabel("z")
    plt.ylabel("H(z) km/s/Mpc")
    plt.title("G99 probe separation: direct H(z) not changed by photon-distance bias")
    plt.grid(alpha=0.3)
    plt.legend()
    p3 = plots / "G99_probe_separation.png"
    plt.tight_layout()
    plt.savefig(p3, dpi=200)
    plt.close()

    # Bias shapes
    plt.figure(figsize=(10, 6))
    for shape in ["linear", "stam_su", "saturating"]:
        plt.plot(z_plot, bias_shape(z_plot, shape), label=shape)
    plt.xlabel("z")
    plt.ylabel("shape(z)")
    plt.title("G99 photon-A bias shapes")
    plt.grid(alpha=0.3)
    plt.legend()
    p4 = plots / "G99_bias_shapes.png"
    plt.tight_layout()
    plt.savefig(p4, dpi=200)
    plt.close()

    # summary
    summary = results / "G99_Hubble_bias_locked_summary.md"
    md = f"""# G99 — Hubble tension as STAM distance bias, locked convention

## Convention

G99 uses the G98-approved convention:

```text
D_C_obs(z) = D_C_intrinsic(z) + f_los * b * shape(z)
D_L_obs(z) = (1+z) * D_C_obs(z)
```

This treats photon-A traversal as a path/comoving correction first, then applies the standard luminosity-distance factor.

## Constants

- `A0 = {A0:.12f}`
- `H0_true = {H0_true:.4f}` km/s/Mpc
- `H0_target = {H0_target:.4f}` km/s/Mpc
- `Omega_m = {omega_m:.4f}`
- `b = A0*c/H0_true = {b:.6f}` Mpc = `{b*MPC_TO_MLY:.6f}` Mly
- selected shape = `{shape_model}`

## Low-z estimate

All supported shapes are linear at very low z, so the leading low-z estimate is:

```text
H0_inferred ≈ H0_true / (1 + f_los*A0)
f_los_required = (H0_true/H0_target - 1)/A0
```

- `f_los_required = {f_required_analytic:.6f}`

## Selected-shape results

| f_los | H0 inferred low-z slope | H0 inferred full mu-fit z<={args.zmax} |
|---:|---:|---:|
"""
    for f, vals in selected.items():
        md += f"| {f:.6f} | {vals['H_low']:.4f} | {vals['H_mu']:.4f} |\n"
    md += """
## Interpretation

The distance-bias layer can make a true local `H0_true` appear lower in distance-only fits. Direct H(z) probes are not changed by the photon-distance bias layer. This is the core STAM interpretation of the Hubble tension as a possible distance-inference tension.

This is still a scaffold, not a publication-grade cosmology likelihood. The next serious test must include SN covariance, BAO likelihoods, CMB acoustic angle, chronometers, structure growth, nuisance parameters, and locked priors.

"""
    if sn_fit:
        md += f"""## Optional SN CSV fit

- `f_los_best = {sn_fit['f_los_best']:.6f}`
- `chi2 = {sn_fit['chi2']:.3f}`
- `N = {sn_fit['N']}`

"""
    if hz_fit:
        md += f"""## Optional H(z) CSV fit

- `H0_best = {hz_fit['H0_best']:.6f}`
- `chi2 = {hz_fit['chi2']:.3f}`
- `N = {hz_fit['N']}`

"""
    md += f"""## Files

- Grid CSV: `{csv_path}`
- Plot: `{p1}`
- Plot: `{p2}`
- Plot: `{p3}`
- Plot: `{p4}`
"""
    summary.write_text(md, encoding="utf-8")

    print(f"Grid CSV written: {csv_path}")
    print(f"Summary written: {summary}")
    print(f"Plot written: {p1}")
    print(f"Plot written: {p2}")
    print(f"Plot written: {p3}")
    print(f"Plot written: {p4}")


if __name__ == "__main__":
    main()

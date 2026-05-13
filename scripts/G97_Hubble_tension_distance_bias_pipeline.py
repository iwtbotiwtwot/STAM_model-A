#!/usr/bin/env python3
"""
G97_Hubble_tension_distance_bias_pipeline.py

STAM / Model-A cosmology test scaffold:
Hubble tension as distance-bias, not pure expansion-rate disagreement.

Purpose
=======
This script separates two layers:

Layer 1 — intrinsic expansion:
    H(z) = H0_true * E(z)
    Default: flat LCDM-like expansion with H0_true = 73.04 km/s/Mpc.

Layer 2 — photon-A distance bias:
    Observed photon distances receive an additive traversal / accumulated-A
    bias controlled by:

        A0 = 1/(12*pi)
        b = A0 * c/H0_true

    In Mpc:
        b_Mpc = A0 * c/H0_true

    In Mly:
        b_Mly = b_Mpc * 3.261563776

The goal is NOT to declare the Hubble tension solved.
The goal is to test whether one locked A0-based photon-distance-bias model can
make distance probes infer a lower H0 while direct H(z) probes remain tied to
the intrinsic expansion.

Core idea
=========
If observed photon distances are larger than the intrinsic expansion distance,
then a model that ignores the bias will infer a lower H0.

At low redshift, with a simple linear bias:

    D_obs(z) ≈ (c/H0_true) z + f_los * b_Mpc * z
             = (c/H0_true) z * (1 + f_los * A0)

A no-bias fit would infer:

    H0_inferred ≈ H0_true / (1 + f_los * A0)

Therefore the f_los required to make H0_true appear as H0_target is:

    f_los_required = (H0_true / H0_target - 1) / A0

Defaults:
    H0_true   = 73.04
    H0_target = 67.40

This gives f_los_required ≈ 3.15.

Bias models
===========
--bias-model linear:
    delta_D = f_los * b_Mpc * z

--bias-model saturating:
    delta_D = f_los * b_Mpc * z/(1+z)

--bias-model stam_su:
    delta_D = f_los * b_Mpc * z * (1 + 3z/20)
    This mirrors the current SU(z) shape factor without claiming it is final.

Outputs
=======
- results/G97_Hubble_distance_bias_summary.md
- results/G97_Hubble_distance_bias_grid.csv
- plots/G97_Hubble_inferred_vs_flos.png
- plots/G97_distance_modulus_residuals.png
- plots/G97_Hz_direct_probe_unchanged.png

Optional data input
===================
This script can also fit simple external CSVs.

SN-like CSV:
    --sn-csv data.csv
Required columns:
    z,mu
Optional:
    sigma_mu

Chronometer-like CSV:
    --hz-csv data.csv
Required columns:
    z,H
Optional:
    sigma_H

The external data support is intentionally simple. A real publication-grade
cosmology pipeline still needs covariance matrices, BAO likelihoods, CMB
theta*, selection effects, nuisance parameters, and locked priors.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Callable, Dict, Tuple, List

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize_scalar


# ---------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------

c_km_s = 299_792.458
MPC_TO_MLY = 3.261563776
A0 = 1.0 / (12.0 * math.pi)


# ---------------------------------------------------------------------
# Cosmology functions
# ---------------------------------------------------------------------

def E_flat_lcdm(z: float, omega_m: float) -> float:
    omega_l = 1.0 - omega_m
    return math.sqrt(omega_m * (1.0 + z)**3 + omega_l)


def H_of_z(z: float, H0: float, omega_m: float) -> float:
    return H0 * E_flat_lcdm(z, omega_m)


def comoving_distance_mpc(z: float, H0: float, omega_m: float) -> float:
    val, _ = quad(lambda zz: c_km_s / H_of_z(zz, H0, omega_m), 0.0, z,
                  epsabs=1e-9, epsrel=1e-9, limit=200)
    return val


def luminosity_distance_mpc(z: float, H0: float, omega_m: float) -> float:
    return (1.0 + z) * comoving_distance_mpc(z, H0, omega_m)


def distance_modulus(d_l_mpc: float) -> float:
    return 5.0 * math.log10(d_l_mpc) + 25.0


def mu_model(z: float, H0: float, omega_m: float) -> float:
    return distance_modulus(luminosity_distance_mpc(z, H0, omega_m))


def bridge_b_mpc(H0_true: float) -> float:
    return A0 * c_km_s / H0_true


def bias_shape(z: np.ndarray | float, model: str):
    z_arr = np.asarray(z, dtype=float)
    if model == "linear":
        out = z_arr
    elif model == "saturating":
        out = z_arr / (1.0 + z_arr)
    elif model == "stam_su":
        out = z_arr * (1.0 + 3.0 * z_arr / 20.0)
    else:
        raise ValueError(f"Unknown bias model: {model}")
    if np.ndim(z) == 0:
        return float(out)
    return out


def observed_distance_with_bias_mpc(z: float, H0_true: float, omega_m: float,
                                    f_los: float, bias_model: str) -> float:
    d_intrinsic = luminosity_distance_mpc(z, H0_true, omega_m)
    b = bridge_b_mpc(H0_true)
    delta = f_los * b * bias_shape(z, bias_model)
    return d_intrinsic + delta


def mu_observed_with_bias(z: float, H0_true: float, omega_m: float,
                          f_los: float, bias_model: str) -> float:
    return distance_modulus(observed_distance_with_bias_mpc(z, H0_true, omega_m, f_los, bias_model))


def fit_H0_no_bias_to_biased_distances(z_grid: np.ndarray, mu_obs: np.ndarray,
                                       omega_m: float, H0_bounds=(50.0, 90.0)) -> float:
    """Fit H0 using no-bias LCDM distances to biased distance moduli."""
    def chi2(H0):
        pred = np.array([mu_model(float(z), H0, omega_m) for z in z_grid])
        return float(np.mean((pred - mu_obs)**2))

    res = minimize_scalar(chi2, bounds=H0_bounds, method="bounded",
                          options={"xatol": 1e-8})
    return float(res.x)


def fit_H0_from_lowz_slope(z_grid: np.ndarray, d_obs_mpc: np.ndarray) -> float:
    """Fit D ~= (c/H0) z through origin at low-z."""
    # d = slope*z; slope = sum(z*d)/sum(z^2)
    slope = float(np.sum(z_grid * d_obs_mpc) / np.sum(z_grid**2))
    return c_km_s / slope


# ---------------------------------------------------------------------
# Optional CSV fitting
# ---------------------------------------------------------------------

def read_csv_cols(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fit_sn_csv(path: Path, H0_true: float, omega_m: float, bias_model: str):
    rows = read_csv_cols(path)
    z = np.array([float(r["z"]) for r in rows])
    mu = np.array([float(r["mu"]) for r in rows])
    if "sigma_mu" in rows[0] and rows[0]["sigma_mu"] != "":
        sig = np.array([float(r.get("sigma_mu", 1.0) or 1.0) for r in rows])
    else:
        sig = np.ones_like(mu)

    def chi2_for(f_los):
        pred = np.array([mu_observed_with_bias(float(zz), H0_true, omega_m, f_los, bias_model)
                         for zz in z])
        return float(np.sum(((mu - pred) / sig)**2))

    res = minimize_scalar(chi2_for, bounds=(-20, 20), method="bounded")
    return {"f_los_best": float(res.x), "chi2": float(res.fun), "N": len(z)}


def fit_hz_csv(path: Path, H0_guess: float, omega_m: float):
    rows = read_csv_cols(path)
    z = np.array([float(r["z"]) for r in rows])
    H = np.array([float(r["H"]) for r in rows])
    if "sigma_H" in rows[0] and rows[0]["sigma_H"] != "":
        sig = np.array([float(r.get("sigma_H", 1.0) or 1.0) for r in rows])
    else:
        sig = np.ones_like(H)

    def chi2_for(H0):
        pred = np.array([H_of_z(float(zz), H0, omega_m) for zz in z])
        return float(np.sum(((H - pred) / sig)**2))

    res = minimize_scalar(chi2_for, bounds=(50, 90), method="bounded")
    return {"H0_best": float(res.x), "chi2": float(res.fun), "N": len(z)}


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="G97 Hubble tension as STAM distance-bias pipeline")
    parser.add_argument("--outdir", default=".", help="Output root directory")
    parser.add_argument("--H0-true", type=float, default=73.04)
    parser.add_argument("--H0-target", type=float, default=67.40)
    parser.add_argument("--omega-m", type=float, default=0.30)
    parser.add_argument("--bias-model", choices=["linear", "saturating", "stam_su"], default="linear")
    parser.add_argument("--zmax", type=float, default=2.0)
    parser.add_argument("--sn-csv", default=None)
    parser.add_argument("--hz-csv", default=None)
    args = parser.parse_args()

    root = Path(args.outdir).resolve()
    results_dir = root / "results"
    plots_dir = root / "plots"
    results_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    H0_true = args.H0_true
    H0_target = args.H0_target
    omega_m = args.omega_m
    bias_model = args.bias_model

    b_mpc = bridge_b_mpc(H0_true)
    b_mly = b_mpc * MPC_TO_MLY

    f_required_lowz = (H0_true / H0_target - 1.0) / A0

    # Grids
    z_grid = np.linspace(0.01, args.zmax, 220)
    z_low = np.linspace(0.005, 0.10, 80)

    f_values = np.linspace(0.0, 6.0, 121)
    grid_rows = []

    for f_los in f_values:
        # Full-range no-bias fit to biased mu
        mu_obs = np.array([mu_observed_with_bias(float(z), H0_true, omega_m, f_los, bias_model)
                           for z in z_grid])
        H0_fit_full = fit_H0_no_bias_to_biased_distances(z_grid, mu_obs, omega_m)

        # Low-z slope fit
        d_obs_low = np.array([observed_distance_with_bias_mpc(float(z), H0_true, omega_m, f_los, bias_model)
                              for z in z_low])
        H0_fit_lowz = fit_H0_from_lowz_slope(z_low, d_obs_low)

        grid_rows.append({
            "f_los": f_los,
            "H0_inferred_full_mu_fit": H0_fit_full,
            "H0_inferred_lowz_slope": H0_fit_lowz,
        })

    # Required f from numerical low-z target nearest
    best_row = min(grid_rows, key=lambda r: abs(r["H0_inferred_lowz_slope"] - H0_target))

    # Bias residual curves using f_required
    f_demo = f_required_lowz
    mu_intrinsic = np.array([mu_model(float(z), H0_true, omega_m) for z in z_grid])
    mu_biased = np.array([mu_observed_with_bias(float(z), H0_true, omega_m, f_demo, bias_model)
                          for z in z_grid])
    mu_target_lcdm = np.array([mu_model(float(z), H0_target, omega_m) for z in z_grid])

    # H(z) curves
    z_h = np.linspace(0.0, 2.0, 200)
    H_true_curve = np.array([H_of_z(float(z), H0_true, omega_m) for z in z_h])
    H_target_curve = np.array([H_of_z(float(z), H0_target, omega_m) for z in z_h])

    print("=" * 96)
    print("G97: Hubble tension as STAM photon-A distance bias")
    print("=" * 96)
    print()
    print(f"A0 = 1/(12π) = {A0:.12f}")
    print(f"H0_true   = {H0_true:.4f} km/s/Mpc")
    print(f"H0_target = {H0_target:.4f} km/s/Mpc")
    print(f"Omega_m   = {omega_m:.4f}")
    print(f"bias model = {bias_model}")
    print()
    print(f"Bridge term b = A0*c/H0_true = {b_mpc:.6f} Mpc = {b_mly:.6f} Mly")
    print(f"Low-z analytic f_los required to make {H0_true:.2f} look like {H0_target:.2f}:")
    print(f"  f_los_required = {f_required_lowz:.6f}")
    print()
    print("Nearest grid result:")
    print(f"  f_los = {best_row['f_los']:.3f}")
    print(f"  H0 inferred low-z slope = {best_row['H0_inferred_lowz_slope']:.4f}")
    print(f"  H0 inferred full mu-fit = {best_row['H0_inferred_full_mu_fit']:.4f}")
    print()

    # Optional CSV fits
    sn_fit = None
    hz_fit = None
    if args.sn_csv:
        sn_fit = fit_sn_csv(Path(args.sn_csv), H0_true, omega_m, bias_model)
        print("SN CSV fit:")
        print(sn_fit)
    if args.hz_csv:
        hz_fit = fit_hz_csv(Path(args.hz_csv), H0_true, omega_m)
        print("H(z) CSV fit:")
        print(hz_fit)

    # Write grid CSV
    grid_csv = results_dir / "G97_Hubble_distance_bias_grid.csv"
    with grid_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["f_los", "H0_inferred_full_mu_fit", "H0_inferred_lowz_slope"])
        writer.writeheader()
        for r in grid_rows:
            writer.writerow(r)

    # Plots
    plt.figure(figsize=(9, 6))
    plt.plot([r["f_los"] for r in grid_rows],
             [r["H0_inferred_lowz_slope"] for r in grid_rows],
             label="No-bias H0 inferred from low-z slope")
    plt.plot([r["f_los"] for r in grid_rows],
             [r["H0_inferred_full_mu_fit"] for r in grid_rows],
             label=f"No-bias H0 inferred from mu fit, z<={args.zmax}")
    plt.axhline(H0_true, linestyle="--", label="true local H0")
    plt.axhline(H0_target, linestyle=":", label="target lower inferred H0")
    plt.axvline(f_required_lowz, linestyle="-.", label=f"analytic f_los ≈ {f_required_lowz:.2f}")
    plt.xlabel("Line-of-sight amplification f_los")
    plt.ylabel("H0 inferred by model that ignores photon-A bias")
    plt.title("STAM distance-bias effect on inferred H0")
    plt.grid(alpha=0.3)
    plt.legend()
    p1 = plots_dir / "G97_Hubble_inferred_vs_flos.png"
    plt.tight_layout()
    plt.savefig(p1, dpi=200)
    plt.close()

    plt.figure(figsize=(9, 6))
    plt.plot(z_grid, mu_biased - mu_intrinsic, label="STAM bias residual vs true H0 distance")
    plt.plot(z_grid, mu_biased - mu_target_lcdm, label="STAM-biased H0=73 minus no-bias H0=67.4")
    plt.xlabel("z")
    plt.ylabel("distance modulus residual Δμ")
    plt.title(f"Photon-A distance bias residuals ({bias_model}, f_los={f_demo:.3f})")
    plt.grid(alpha=0.3)
    plt.legend()
    p2 = plots_dir / "G97_distance_modulus_residuals.png"
    plt.tight_layout()
    plt.savefig(p2, dpi=200)
    plt.close()

    plt.figure(figsize=(9, 6))
    plt.plot(z_h, H_true_curve, label=f"Intrinsic H(z), H0={H0_true:.2f}")
    plt.plot(z_h, H_target_curve, label=f"Lower-H0 comparison, H0={H0_target:.2f}")
    plt.xlabel("z")
    plt.ylabel("H(z) km/s/Mpc")
    plt.title("Direct H(z) probes are Layer-1, not photon-distance bias")
    plt.grid(alpha=0.3)
    plt.legend()
    p3 = plots_dir / "G97_Hz_direct_probe_unchanged.png"
    plt.tight_layout()
    plt.savefig(p3, dpi=200)
    plt.close()

    # Summary
    summary = results_dir / "G97_Hubble_distance_bias_summary.md"
    md = []
    md.append("# G97 — Hubble tension as STAM photon-A distance bias\n\n")
    md.append("## Interpretation\n\n")
    md.append("This script tests the idea that the Hubble tension may be partly a **distance-inference tension** rather than a pure expansion-rate disagreement.\n\n")
    md.append("- Layer 1: intrinsic expansion \(H(z)=H_0E(z)\).\n")
    md.append("- Layer 2: photon-A traversal bias adds apparent distance to photon observables.\n")
    md.append("- Direct \(H(z)\) probes such as cosmic chronometers should respond mainly to Layer 1.\n\n")
    md.append("## Locked constants\n\n")
    md.append(f"- \(A_0=1/(12\\pi)={A0:.12f}\)\n")
    md.append(f"- \(H_0^{{true}}={H0_true:.4f}\\) km/s/Mpc\n")
    md.append(f"- \(H_0^{{target}}={H0_target:.4f}\\) km/s/Mpc\n")
    md.append(f"- Bias model: `{bias_model}`\n")
    md.append(f"- Bridge term \(b=A_0c/H_0={b_mpc:.6f}\\) Mpc = `{b_mly:.6f}` Mly\n\n")
    md.append("## Low-z analytic estimate\n\n")
    md.append("For a linear low-z bias, \(D_{obs}\\approx(c/H_0)z(1+f_{los}A_0)\), so a no-bias fit infers:\n\n")
    md.append("\\[\nH_0^{inferred}\\approx\\frac{H_0^{true}}{1+f_{los}A_0}.\n\\]\n\n")
    md.append("Solving for \(f_{los}\):\n\n")
    md.append("\\[\nf_{los}=\\frac{H_0^{true}/H_0^{target}-1}{A_0}.\n\\]\n\n")
    md.append(f"- Required \(f_{{los}}\) = `{f_required_lowz:.6f}`\n")
    md.append(f"- Nearest grid \(f_{{los}}\) = `{best_row['f_los']:.3f}`\n")
    md.append(f"- H0 inferred from low-z slope = `{best_row['H0_inferred_lowz_slope']:.4f}`\n")
    md.append(f"- H0 inferred from full distance-modulus fit = `{best_row['H0_inferred_full_mu_fit']:.4f}`\n\n")
    md.append("## Caveat\n\n")
    md.append("This is a scaffold, not a publication-grade cosmology pipeline. A real test must include SN covariance, BAO likelihoods, CMB acoustic angle, chronometers, structure growth, nuisance parameters, and locked priors. The point of G97 is to separate Layer-1 expansion from Layer-2 photon-distance bias and quantify the required scale.\n\n")
    if sn_fit:
        md.append("## SN CSV fit\n\n")
        md.append(f"- best f_los = `{sn_fit['f_los_best']:.6f}`\n")
        md.append(f"- chi2 = `{sn_fit['chi2']:.3f}` for N = `{sn_fit['N']}`\n\n")
    if hz_fit:
        md.append("## H(z) CSV fit\n\n")
        md.append(f"- best H0 = `{hz_fit['H0_best']:.6f}`\n")
        md.append(f"- chi2 = `{hz_fit['chi2']:.3f}` for N = `{hz_fit['N']}`\n\n")
    md.append("## Files\n\n")
    md.append(f"- Grid CSV: `{grid_csv}`\n")
    md.append(f"- Plot: `{p1}`\n")
    md.append(f"- Plot: `{p2}`\n")
    md.append(f"- Plot: `{p3}`\n")
    summary.write_text("".join(md), encoding="utf-8")

    print(f"Grid CSV written: {grid_csv}")
    print(f"Summary written: {summary}")
    print(f"Plot written: {p1}")
    print(f"Plot written: {p2}")
    print(f"Plot written: {p3}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
G100_Hubble_bias_mu_fit_target.py

STAM / Model-A Hubble-tension distance-bias test using the G98/G99 locked
comoving-distance convention, but fitting the FULL distance-modulus curve
rather than using a crude low-z slope.

Why G100 exists
===============
G99 showed that low-z slope fits can be misleading because even f_los=0
does not recover H0_true if the slope is fit over z=0.005–0.10 with a
linear Hubble-law approximation. The full distance-modulus fit correctly
recovers H0_true when f_los=0.

Therefore G100 asks the better question:

    For a true intrinsic H0 = 73.04,
    what f_los makes a biased photon-distance curve look like
    a no-bias LCDM curve with H0 = 67.4
    over a chosen redshift range?

Locked convention
=================
Photon-A traversal is a comoving/path correction first:

    D_C_obs(z) = D_C_intrinsic(z; H0_true) + f_los * b * shape(z)

Then luminosity distance:

    D_L_obs(z) = (1+z) * D_C_obs(z)

where:

    A0 = 1/(12π)
    b  = A0*c/H0_true

Bias shapes
===========
linear:
    shape(z)=z

stam_su:
    shape(z)=z*(1+3z/20)

saturating:
    shape(z)=z/(1+z)

Deliverables
============
For each redshift range and shape, solve for f_los such that the no-bias
mu-fit inferred H0 is as close as possible to H0_target.

Outputs:
    results/G100_Hubble_bias_mu_fit_target_summary.md
    results/G100_Hubble_bias_mu_fit_target_grid.csv
    plots/G100_H0_mu_fit_vs_flos.png
    plots/G100_best_fit_residuals.png
    plots/G100_required_flos_by_zmax.png

Optional:
    --sn-csv path.csv
      columns: z,mu[,sigma_mu]
      Fits f_los directly to provided SN-like data under the locked convention.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar


# ---------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------

c_km_s = 299_792.458
MPC_TO_MLY = 3.261563776
A0 = 1.0 / (12.0 * math.pi)


def E_flat_lcdm(z, omega_m):
    z = np.asarray(z, dtype=float)
    return np.sqrt(omega_m * (1.0 + z)**3 + (1.0 - omega_m))


def H_z(z, H0, omega_m):
    return H0 * E_flat_lcdm(z, omega_m)


def cumulative_comoving_distance(z_grid, H0, omega_m):
    z = np.asarray(z_grid, dtype=float)
    z_full = np.concatenate([[0.0], z])
    integrand = c_km_s / H_z(z_full, H0, omega_m)
    dz = np.diff(z_full)
    trap = 0.5 * (integrand[1:] + integrand[:-1]) * dz
    return np.cumsum(trap)


def DL_from_DC(z, DC):
    return (1.0 + np.asarray(z, dtype=float)) * np.asarray(DC, dtype=float)


def mu_from_DL(DL):
    DL = np.asarray(DL, dtype=float)
    return 5.0 * np.log10(DL) + 25.0


def bridge_b_mpc(H0_true):
    return A0 * c_km_s / H0_true


def shape_z(z, model):
    z = np.asarray(z, dtype=float)
    if model == "linear":
        return z
    if model == "stam_su":
        return z * (1.0 + 3.0*z/20.0)
    if model == "saturating":
        return z / (1.0 + z)
    raise ValueError(model)


def mu_no_bias(z, H0, omega_m):
    DC = cumulative_comoving_distance(z, H0, omega_m)
    return mu_from_DL(DL_from_DC(z, DC))


def mu_biased(z, H0_true, omega_m, f_los, shape):
    DC = cumulative_comoving_distance(z, H0_true, omega_m)
    DC_obs = DC + f_los * bridge_b_mpc(H0_true) * shape_z(z, shape)
    return mu_from_DL(DL_from_DC(z, DC_obs))


def infer_H0_from_mu(z, mu_obs, omega_m, bounds=(50, 90)):
    z = np.asarray(z, dtype=float)
    mu_obs = np.asarray(mu_obs, dtype=float)

    def loss(H0):
        mu_pred = mu_no_bias(z, H0, omega_m)
        return float(np.mean((mu_pred - mu_obs)**2))

    res = minimize_scalar(loss, bounds=bounds, method="bounded", options={"xatol": 1e-10})
    return float(res.x), float(res.fun)


def solve_f_for_target_H0(z, H0_true, H0_target, omega_m, shape, f_bounds=(-10, 20)):
    def objective(f_los):
        mu_obs = mu_biased(z, H0_true, omega_m, f_los, shape)
        H0_fit, loss = infer_H0_from_mu(z, mu_obs, omega_m)
        return (H0_fit - H0_target)**2

    res = minimize_scalar(objective, bounds=f_bounds, method="bounded", options={"xatol": 1e-8})
    f_best = float(res.x)
    mu_obs = mu_biased(z, H0_true, omega_m, f_best, shape)
    H0_fit, loss = infer_H0_from_mu(z, mu_obs, omega_m)
    return f_best, H0_fit, loss


def read_sn_csv(path):
    with Path(path).open("r", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    z = np.array([float(r["z"]) for r in rows])
    mu = np.array([float(r["mu"]) for r in rows])
    if rows and "sigma_mu" in rows[0] and rows[0]["sigma_mu"] != "":
        sig = np.array([float(r.get("sigma_mu", 1.0) or 1.0) for r in rows])
    else:
        sig = np.ones_like(mu)
    return z, mu, sig


def fit_f_to_sn(z, mu_data, sig_mu, H0_true, omega_m, shape, f_bounds=(-20, 20)):
    def chi2(f_los):
        pred = mu_biased(z, H0_true, omega_m, f_los, shape)
        return float(np.sum(((mu_data - pred) / sig_mu)**2))

    res = minimize_scalar(chi2, bounds=f_bounds, method="bounded", options={"xatol": 1e-8})
    return float(res.x), float(res.fun)


def main():
    parser = argparse.ArgumentParser(description="G100 Hubble bias full mu-fit target")
    parser.add_argument("--outdir", default="/mnt/data/G100_run")
    parser.add_argument("--H0-true", type=float, default=73.04)
    parser.add_argument("--H0-target", type=float, default=67.40)
    parser.add_argument("--omega-m", type=float, default=0.30)
    parser.add_argument("--sn-csv", default=None)
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

    shapes = ["linear", "stam_su", "saturating"]
    z_ranges = [
        ("low_SN_0p01_0p15", 0.01, 0.15, 180),
        ("pantheon_like_0p01_2p0", 0.01, 2.00, 420),
        ("low_mid_0p01_0p8", 0.01, 0.80, 300),
        ("high_z_0p1_2p0", 0.10, 2.00, 360),
    ]

    rows = []
    grids = {}

    for name, zmin, zmax, n in z_ranges:
        z = np.linspace(zmin, zmax, n)
        grids[name] = z
        for shape in shapes:
            # baseline check f=0
            mu0 = mu_biased(z, H0_true, omega_m, 0.0, shape)
            H0_f0, loss0 = infer_H0_from_mu(z, mu0, omega_m)

            # f=1 check
            mu1 = mu_biased(z, H0_true, omega_m, 1.0, shape)
            H0_f1, loss1 = infer_H0_from_mu(z, mu1, omega_m)

            # solve f target
            f_best, H0_best, loss_best = solve_f_for_target_H0(
                z, H0_true, H0_target, omega_m, shape
            )

            rows.append({
                "range": name,
                "zmin": zmin,
                "zmax": zmax,
                "shape": shape,
                "H0_fit_f0": H0_f0,
                "H0_fit_f1": H0_f1,
                "f_los_target": f_best,
                "H0_fit_target": H0_best,
                "loss_target": loss_best,
            })

    # optional SN fits
    sn_fit_rows = []
    if args.sn_csv:
        z_data, mu_data, sig_mu = read_sn_csv(args.sn_csv)
        for shape in shapes:
            f_best, chi2 = fit_f_to_sn(z_data, mu_data, sig_mu, H0_true, omega_m, shape)
            sn_fit_rows.append({
                "shape": shape,
                "f_los_best": f_best,
                "chi2": chi2,
                "N": len(z_data),
                "chi2_per_N": chi2 / len(z_data) if len(z_data) else float("nan"),
            })

    # Write CSV
    csv_path = results / "G100_Hubble_bias_mu_fit_target_grid.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(rows[0].keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    if sn_fit_rows:
        sn_csv = results / "G100_Hubble_bias_sn_fit.csv"
        with sn_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(sn_fit_rows[0].keys()))
            writer.writeheader()
            writer.writerows(sn_fit_rows)
    else:
        sn_csv = None

    # Console output
    print("=" * 100)
    print("G100: Hubble distance-bias full mu-fit target")
    print("=" * 100)
    print()
    print("Convention:")
    print("  D_C_obs = D_C_intrinsic + f_los*b*shape(z)")
    print("  D_L_obs = (1+z)*D_C_obs")
    print()
    print(f"A0 = {A0:.12f}")
    print(f"H0_true   = {H0_true:.4f}")
    print(f"H0_target = {H0_target:.4f}")
    print(f"Omega_m   = {omega_m:.4f}")
    print(f"b = {b:.6f} Mpc = {b*3.261563776:.6f} Mly")
    print()
    print(f"{'range':<24}{'shape':<12}{'H0 f=0':>10}{'H0 f=1':>10}{'f_target':>12}{'H0 target fit':>15}")
    print("-"*86)
    for r in rows:
        print(f"{r['range']:<24}{r['shape']:<12}{r['H0_fit_f0']:>10.3f}{r['H0_fit_f1']:>10.3f}"
              f"{r['f_los_target']:>12.3f}{r['H0_fit_target']:>15.3f}")
    print()

    # Plots: f_target by range/shape
    plt.figure(figsize=(11, 6))
    x = np.arange(len(z_ranges))
    width = 0.25
    for i, shape in enumerate(shapes):
        vals = [r["f_los_target"] for r in rows if r["shape"] == shape]
        plt.bar(x + (i-1)*width, vals, width, label=shape)
    plt.xticks(x, [r[0] for r in z_ranges], rotation=20, ha="right")
    plt.ylabel("f_los required for no-bias mu-fit to infer H0_target")
    plt.title("G100: required f_los by redshift range and bias shape")
    plt.grid(axis="y", alpha=0.3)
    plt.legend()
    p1 = plots / "G100_required_flos_by_zrange.png"
    plt.tight_layout()
    plt.savefig(p1, dpi=200)
    plt.close()

    # Plot H0 fit vs f for pantheon-like range
    range_name = "pantheon_like_0p01_2p0"
    z = grids[range_name]
    f_grid = np.linspace(0, 8, 250)

    plt.figure(figsize=(10, 6))
    for shape in shapes:
        Hs = []
        for f in f_grid:
            mu = mu_biased(z, H0_true, omega_m, f, shape)
            H_fit, loss = infer_H0_from_mu(z, mu, omega_m)
            Hs.append(H_fit)
        plt.plot(f_grid, Hs, label=shape)
    plt.axhline(H0_true, color="gray", linestyle="--", label="H0 true")
    plt.axhline(H0_target, color="black", linestyle=":", label="H0 target")
    plt.xlabel("f_los")
    plt.ylabel("H0 inferred by no-bias full mu-fit")
    plt.title("G100: full distance-modulus inferred H0 vs f_los")
    plt.grid(alpha=0.3)
    plt.legend()
    p2 = plots / "G100_H0_mu_fit_vs_flos.png"
    plt.tight_layout()
    plt.savefig(p2, dpi=200)
    plt.close()

    # Best residuals for pantheon-like range
    plt.figure(figsize=(10, 6))
    for shape in shapes:
        rr = [r for r in rows if r["range"] == range_name and r["shape"] == shape][0]
        f = rr["f_los_target"]
        mu_obs = mu_biased(z, H0_true, omega_m, f, shape)
        mu_target = mu_no_bias(z, H0_target, omega_m)
        plt.plot(z, mu_obs - mu_target, label=f"{shape}, f={f:.2f}")
    plt.xlabel("z")
    plt.ylabel("mu_biased(H0=73) - mu_no_bias(H0=67.4)")
    plt.title("G100: residuals after matching H0 target by full mu-fit")
    plt.grid(alpha=0.3)
    plt.legend()
    p3 = plots / "G100_best_fit_residuals.png"
    plt.tight_layout()
    plt.savefig(p3, dpi=200)
    plt.close()

    # Summary
    summary = results / "G100_Hubble_bias_mu_fit_target_summary.md"
    md = []
    md.append("# G100 — Hubble bias full distance-modulus target\n\n")
    md.append("## Purpose\n\n")
    md.append("G100 uses the G98/G99 locked convention but targets the full distance-modulus fit rather than a crude low-z slope. It asks: for true intrinsic `H0=73.04`, what `f_los` makes a biased photon-distance curve look like a no-bias `H0=67.4` curve over chosen redshift ranges?\n\n")
    md.append("## Locked convention\n\n")
    md.append("```text\n")
    md.append("D_C_obs(z) = D_C_intrinsic(z) + f_los * b * shape(z)\n")
    md.append("D_L_obs(z) = (1+z) * D_C_obs(z)\n")
    md.append("```\n\n")
    md.append("## Constants\n\n")
    md.append(f"- `A0 = {A0:.12f}`\n")
    md.append(f"- `H0_true = {H0_true:.4f}` km/s/Mpc\n")
    md.append(f"- `H0_target = {H0_target:.4f}` km/s/Mpc\n")
    md.append(f"- `Omega_m = {omega_m:.4f}`\n")
    md.append(f"- `b = {b:.6f}` Mpc = `{b*MPC_TO_MLY:.6f}` Mly\n\n")
    md.append("## Results\n\n")
    md.append("| Redshift range | Shape | H0 fit at f=0 | H0 fit at f=1 | f_los target | H0 target fit |\n")
    md.append("|---|---|---:|---:|---:|---:|\n")
    for r in rows:
        md.append(f"| {r['range']} | {r['shape']} | {r['H0_fit_f0']:.3f} | {r['H0_fit_f1']:.3f} | {r['f_los_target']:.3f} | {r['H0_fit_target']:.3f} |\n")
    md.append("\n## Interpretation\n\n")
    md.append("The full distance-modulus fit correctly recovers `H0_true` when `f_los=0`. Nonzero photon-A bias pushes the inferred no-bias H0 lower. The required `f_los` depends on redshift range and bias-shape model, so a real test must use the actual SN/BAO/CMB/chronometer likelihoods rather than a single low-z estimate.\n\n")
    md.append("## Status\n\n")
    md.append("G100 is still a scaffold, not a final cosmology result. It locks the distance convention and identifies the line-of-sight amplification scale required to turn local `H0≈73` into distance-inferred `H0≈67.4` under different bias shapes and redshift ranges. The next real test is to plug in actual SN and BAO/CMB data with covariance and fixed priors.\n\n")
    if sn_fit_rows:
        md.append("## Optional SN CSV fit\n\n")
        md.append("| Shape | f_los best | chi2 | N | chi2/N |\n")
        md.append("|---|---:|---:|---:|---:|\n")
        for r in sn_fit_rows:
            md.append(f"| {r['shape']} | {r['f_los_best']:.4f} | {r['chi2']:.3f} | {r['N']} | {r['chi2_per_N']:.3f} |\n")
        md.append("\n")
    md.append("## Files\n\n")
    md.append(f"- Grid CSV: `{csv_path}`\n")
    if sn_csv:
        md.append(f"- SN fit CSV: `{sn_csv}`\n")
    md.append(f"- Plot: `{p1}`\n")
    md.append(f"- Plot: `{p2}`\n")
    md.append(f"- Plot: `{p3}`\n")
    summary.write_text("".join(md), encoding="utf-8")

    print(f"Grid CSV written: {csv_path}")
    if sn_csv:
        print(f"SN CSV written: {sn_csv}")
    print(f"Summary written: {summary}")
    print(f"Plot written: {p1}")
    print(f"Plot written: {p2}")
    print(f"Plot written: {p3}")


if __name__ == "__main__":
    main()

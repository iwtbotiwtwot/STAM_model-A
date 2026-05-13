#!/usr/bin/env python3
"""
G101_multimessenger_cosmo_compressed.py

Fast compressed STAM / Model-A multi-probe cosmology stress test:
SN + BAO + CMB-distance + chronometer-like H(z)

This is NOT a publication-grade likelihood. It is a diagnostic scaffold.

Core question
=============
Can a single STAM photon-A distance-bias layer make an intrinsic
H0_true = 73.04 look like a lower no-bias H0_ref = 67.40 in distance probes,
while direct H(z) probes remain tied to the intrinsic expansion?

Locked convention from G98/G99
==============================
D_C_obs(z) = D_C_intrinsic(z; H0_true) + f_los * b * shape(z)
D_L_obs(z) = (1+z) * D_C_obs(z)

where:
    A0 = 1/(12π)
    b  = A0*c/H0_true

Direct H(z) probes are NOT distance-biased.

Important method note
=====================
For fixed Omega_m, flat LCDM distances scale exactly as 1/H0:
    D_C(z; H0) = (H0_ref/H0) D_C(z; H0_ref)

So the best no-bias H0 fit to a biased distance-modulus curve can be computed
analytically from the mean distance-modulus offset, avoiding slow nested
optimizers.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

c_km_s = 299_792.458
MPC_TO_MLY = 3.261563776
A0 = 1.0 / (12.0 * math.pi)


def E_flat_lcdm(z, omega_m):
    z = np.asarray(z, dtype=float)
    return np.sqrt(omega_m*(1+z)**3 + (1-omega_m))


def H_z(z, H0, omega_m):
    return H0 * E_flat_lcdm(z, omega_m)


def DC_grid(z, H0, omega_m):
    z = np.asarray(z, dtype=float)
    z_full = np.concatenate([[0.0], z])
    integ = c_km_s / H_z(z_full, H0, omega_m)
    dz = np.diff(z_full)
    trap = 0.5*(integ[1:]+integ[:-1])*dz
    return np.cumsum(trap)


def DL(z, DC):
    return (1+np.asarray(z))*np.asarray(DC)


def mu(DL_mpc):
    return 5*np.log10(np.asarray(DL_mpc)) + 25


def bridge_b(H0_true):
    return A0 * c_km_s / H0_true


def shape_z(z, shape):
    z = np.asarray(z, dtype=float)
    if shape == "linear":
        return z
    if shape == "stam_su":
        return z*(1+3*z/20)
    if shape == "saturating":
        return z/(1+z)
    if shape == "log":
        return np.log1p(z)
    raise ValueError(shape)


def biased_DC(z, H0_true, omega_m, f_los, shape):
    dc = DC_grid(z, H0_true, omega_m)
    return dc + f_los * bridge_b(H0_true) * shape_z(z, shape)


def best_H0_from_mu_offset(z, mu_obs, H0_ref, omega_m):
    """For fixed Omega_m, no-bias mu(H0)=mu(H0_ref)-5log10(H0/H0_ref).
    Min mean-square residual over one scalar offset."""
    mu_ref = mu(DL(z, DC_grid(z, H0_ref, omega_m)))
    delta = float(np.mean(mu_obs - mu_ref))
    # mu_obs ≈ mu_ref - 5 log10(H0_fit/H0_ref)
    return H0_ref * 10**(-delta/5)


def chi2_components(H0_true, H0_ref, omega_m, f_los, shape, setup):
    out = {}

    # SN-like distance modulus
    z_sn = setup["z_sn"]
    mu_ref = mu(DL(z_sn, DC_grid(z_sn, H0_ref, omega_m)))
    mu_model = mu(DL(z_sn, biased_DC(z_sn, H0_true, omega_m, f_los, shape)))
    sig_mu = setup["sigma_mu"]
    out["SN"] = (float(np.sum(((mu_model-mu_ref)/sig_mu)**2)), len(z_sn))

    # BAO transverse D_M
    z_bao = setup["z_bao"]
    DM_ref = DC_grid(z_bao, H0_ref, omega_m)
    DM_model = biased_DC(z_bao, H0_true, omega_m, f_los, shape)
    frac_DM = setup["sigma_bao_DM_frac"]
    out["BAO_DM"] = (float(np.sum(((DM_model-DM_ref)/(frac_DM*DM_ref))**2)), len(z_bao))

    # BAO radial D_H = c/H(z) no bias
    DH_ref = c_km_s / H_z(z_bao, H0_ref, omega_m)
    DH_model = c_km_s / H_z(z_bao, H0_true, omega_m)
    frac_DH = setup["sigma_bao_DH_frac"]
    out["BAO_DH"] = (float(np.sum(((DH_model-DH_ref)/(frac_DH*DH_ref))**2)), len(z_bao))

    # CMB compressed distance
    zstar = np.array([setup["zstar"]])
    DM_ref_cmb = DC_grid(zstar, H0_ref, omega_m)[0]
    DM_model_cmb = biased_DC(zstar, H0_true, omega_m, f_los, shape)[0]
    frac_cmb = setup["sigma_cmb_D_frac"]
    out["CMB_D"] = (float(((DM_model_cmb-DM_ref_cmb)/(frac_cmb*DM_ref_cmb))**2), 1)

    # Chronometer H(z) no bias
    z_h = setup["z_hz"]
    H_ref = H_z(z_h, H0_ref, omega_m)
    H_model = H_z(z_h, H0_true, omega_m)
    frac_H = setup["sigma_H_frac"]
    out["H_z"] = (float(np.sum(((H_model-H_ref)/(frac_H*H_ref))**2)), len(z_h))

    total = sum(v[0] for v in out.values())
    N = sum(v[1] for v in out.values())
    out["TOTAL"] = (float(total), N)
    return out


def fit_f(H0_true, H0_ref, omega_m, shape, setup, include):
    def obj(f):
        comps = chi2_components(H0_true, H0_ref, omega_m, f, shape, setup)
        return sum(comps[k][0] for k in include)
    res = minimize_scalar(obj, bounds=(-5,20), method="bounded", options={"xatol":1e-6})
    f = float(res.x)
    return f, chi2_components(H0_true, H0_ref, omega_m, f, shape, setup)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="/mnt/data/G101_run")
    p.add_argument("--H0-true", type=float, default=73.04)
    p.add_argument("--H0-ref", type=float, default=67.40)
    p.add_argument("--omega-m", type=float, default=0.30)
    args = p.parse_args()

    root = Path(args.outdir)
    results = root/"results"; plots = root/"plots"
    results.mkdir(parents=True, exist_ok=True); plots.mkdir(parents=True, exist_ok=True)

    H0_true, H0_ref, omega_m = args.H0_true, args.H0_ref, args.omega_m

    setup = {
        "z_sn": np.linspace(0.01,2.0,80),
        "sigma_mu": 0.10,
        "z_bao": np.array([0.38,0.51,0.61,1.48,2.33]),
        "sigma_bao_DM_frac": 0.02,
        "sigma_bao_DH_frac": 0.03,
        "zstar": 1089.0,
        "sigma_cmb_D_frac": 0.005,
        "z_hz": np.array([0.1,0.2,0.35,0.5,0.7,1.0,1.5,2.0]),
        "sigma_H_frac": 0.07,
    }

    shapes = ["linear","stam_su","saturating","log"]
    fit_sets = {
        "SN_only": ["SN"],
        "CMB_only": ["CMB_D"],
        "distance_only_SN_BAO_CMB": ["SN","BAO_DM","CMB_D"],
        "distance_plus_radial_BAO": ["SN","BAO_DM","CMB_D","BAO_DH"],
        "all_SN_BAO_CMB_Hz": ["SN","BAO_DM","CMB_D","BAO_DH","H_z"],
    }

    rows = []
    for shape in shapes:
        for fit_name, include in fit_sets.items():
            fbest, comps = fit_f(H0_true,H0_ref,omega_m,shape,setup,include)
            row = {"shape":shape, "fit_set":fit_name, "include":"+".join(include), "f_los_best":fbest}
            for k,(chi2,N) in comps.items():
                row[f"chi2_{k}"]=chi2
                row[f"N_{k}"]=N
                row[f"chi2_per_{k}"]=chi2/N if N else np.nan
            rows.append(row)

    csv_path = results/"G101_multimessenger_cosmo_grid.csv"
    fields = sorted(set().union(*(r.keys() for r in rows)))
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)

    print("="*110)
    print("G101: compressed SN+BAO+CMB+chronometer stress test")
    print("="*110)
    print()
    print(f"A0 = {A0:.12f}")
    print(f"H0_true = {H0_true:.4f}")
    print(f"H0_ref  = {H0_ref:.4f}")
    print(f"b = {bridge_b(H0_true):.6f} Mpc = {bridge_b(H0_true)*MPC_TO_MLY:.6f} Mly")
    print()
    print(f"{'shape':<12}{'fit_set':<28}{'f_best':>9}{'chi2/N total':>14}{'SN':>10}{'BAO_DM':>10}{'BAO_DH':>10}{'CMB':>10}{'Hz':>10}")
    print("-"*120)
    for r in rows:
        print(f"{r['shape']:<12}{r['fit_set']:<28}{r['f_los_best']:>9.3f}{r['chi2_per_TOTAL']:>14.2f}"
              f"{r['chi2_per_SN']:>10.2f}{r['chi2_per_BAO_DM']:>10.2f}"
              f"{r['chi2_per_BAO_DH']:>10.2f}{r['chi2_per_CMB_D']:>10.2f}{r['chi2_per_H_z']:>10.2f}")

    # Plots
    all_rows = [r for r in rows if r["fit_set"]=="all_SN_BAO_CMB_Hz"]
    probes = ["SN","BAO_DM","BAO_DH","CMB_D","H_z"]
    x = np.arange(len(shapes)); width=0.15
    plt.figure(figsize=(12,6))
    for i,probe in enumerate(probes):
        vals=[r[f"chi2_per_{probe}"] for r in all_rows]
        plt.bar(x+(i-2)*width, vals, width, label=probe)
    plt.xticks(x, shapes); plt.yscale("log")
    plt.ylabel("chi2 / N (log scale)")
    plt.title("G101: compressed residuals at best combined f_los")
    plt.grid(axis="y", alpha=0.3); plt.legend()
    p1 = plots/"G101_chi2_by_probe.png"; plt.tight_layout(); plt.savefig(p1,dpi=200); plt.close()

    # Residual plot distance-only fits.
    z = np.linspace(0.01,2.0,300)
    mu_ref = mu(DL(z, DC_grid(z,H0_ref,omega_m)))
    plt.figure(figsize=(10,6))
    for shape in shapes:
        rbest = [r for r in rows if r["shape"]==shape and r["fit_set"]=="distance_only_SN_BAO_CMB"][0]
        fbest = rbest["f_los_best"]
        mu_model = mu(DL(z, biased_DC(z,H0_true,omega_m,fbest,shape)))
        plt.plot(z, mu_model-mu_ref, label=f"{shape}, f={fbest:.2f}")
    plt.axhline(0,color="black",lw=0.8)
    plt.xlabel("z"); plt.ylabel("mu_STAM_biased - mu_ref_lowerH0")
    plt.title("G101: distance residuals after distance-only compressed fit")
    plt.grid(alpha=0.3); plt.legend(fontsize=8)
    p2=plots/"G101_best_residuals.png"; plt.tight_layout(); plt.savefig(p2,dpi=200); plt.close()

    # Shape comparison.
    zshape=np.linspace(0,5,300)
    plt.figure(figsize=(10,6))
    for shape in shapes:
        plt.plot(zshape, shape_z(zshape,shape), label=shape)
    plt.xlabel("z"); plt.ylabel("shape(z)")
    plt.title("G101: tested photon-A bias shapes")
    plt.grid(alpha=0.3); plt.legend()
    p3=plots/"G101_shape_comparison.png"; plt.tight_layout(); plt.savefig(p3,dpi=200); plt.close()

    summary = results/"G101_multimessenger_cosmo_summary.md"
    md=[]
    md.append("# G101 — compressed SN+BAO+CMB+chronometer stress test\n\n")
    md.append("This is a compressed toy likelihood, not a publication-grade cosmology result. It uses synthetic reference data generated from no-bias LCDM with lower `H0_ref`, then tests whether a STAM model with intrinsic `H0_true` plus photon-A distance bias can mimic it across probes.\n\n")
    md.append("## Locked convention\n\n```text\nD_C_obs(z) = D_C_intrinsic(z; H0_true) + f_los * b * shape(z)\nD_L_obs(z) = (1+z) * D_C_obs(z)\nDirect H(z) probes are not distance-biased.\n```\n\n")
    md.append(f"- `A0 = {A0:.12f}`\n")
    md.append(f"- `H0_true = {H0_true:.4f}`\n")
    md.append(f"- `H0_ref = {H0_ref:.4f}`\n")
    md.append(f"- `b = {bridge_b(H0_true):.6f}` Mpc = `{bridge_b(H0_true)*MPC_TO_MLY:.6f}` Mly\n\n")
    md.append("## Results\n\n")
    md.append("| Shape | Fit set | f_los best | chi2/N total | SN | BAO_DM | BAO_DH | CMB_D | H_z |\n")
    md.append("|---|---|---:|---:|---:|---:|---:|---:|---:|\n")
    for r in rows:
        md.append(f"| {r['shape']} | {r['fit_set']} | {r['f_los_best']:.3f} | {r['chi2_per_TOTAL']:.2f} | {r['chi2_per_SN']:.2f} | {r['chi2_per_BAO_DM']:.2f} | {r['chi2_per_BAO_DH']:.2f} | {r['chi2_per_CMB_D']:.2f} | {r['chi2_per_H_z']:.2f} |\n")
    md.append("\n## Interpretation\n\n")
    md.append("Distance-only probes can be shifted by `f_los`; direct expansion probes cannot. If a distance-bias shape improves SN/CMB transverse distances while BAO radial and chronometer H(z) remain badly off, the model is not yet a full Hubble-tension solution. It is a distance-bias mechanism that still needs a realistic multi-probe likelihood and likely a structure-dependent `f_los(z, sightline)` rather than one constant.\n\n")
    md.append("## Files\n\n")
    md.append(f"- CSV: `{csv_path}`\n- Plot: `{p1}`\n- Plot: `{p2}`\n- Plot: `{p3}`\n")
    summary.write_text("".join(md), encoding="utf-8")

    print(f"CSV written: {csv_path}")
    print(f"Summary written: {summary}")
    print(f"Plot written: {p1}")
    print(f"Plot written: {p2}")
    print(f"Plot written: {p3}")


if __name__ == "__main__":
    main()

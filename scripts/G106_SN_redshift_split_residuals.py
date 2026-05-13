#!/usr/bin/env python3
"""
G106_SN_redshift_split_residuals.py

STAM / Model-A SN-only redshift-split residual test.

Purpose
=======
Isolate data/sn.csv and test whether STAM improves the Hubble diagram smoothly
across redshift bins, or merely absorbs a constant offset.

Input
=====
data/sn.csv with columns:

    z,mu,sigma_mu

Model comparison
================
Baseline LCDM:
    H0, Omega_m, Mcal

STAM distance-bias:
    H0, Omega_m, f0, f1, Mcal

Locked STAM convention:
    D_C_obs(z) = D_C_intrinsic(z) + b * f_los(z) * shape(z)
    D_L_obs(z) = (1+z) * D_C_obs(z)

where:
    A0 = 1/(12π)
    b  = A0*c/H0
    f_los(z) = f0 + f1*z/(1+z)       default
    shape(z) = linear                default

Outputs
=======
results/G106_SN_redshift_split_summary.md
results/G106_SN_best_fit_params.csv
results/G106_SN_bin_stats.csv
results/G106_SN_row_residuals.csv
plots/G106_SN_residuals_vs_z.png
plots/G106_SN_bin_chi2.png
plots/G106_SN_flos_shape.png
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, minimize

c_km_s = 299_792.458
A0 = 1.0 / (12.0 * math.pi)

def E_flat_lcdm(z, omega_m):
    z = np.asarray(z, dtype=float)
    return np.sqrt(omega_m * (1.0 + z)**3 + (1.0 - omega_m))

def H_z(z, H0, omega_m):
    return float(H0) * E_flat_lcdm(z, omega_m)

def comoving_distance_interp(z_values, H0, omega_m):
    z_values = np.asarray(z_values, dtype=float)
    if len(z_values) == 0:
        return np.array([])
    zmax = float(np.max(z_values))
    grid = np.linspace(0.0, zmax, max(4000, int(2500 * max(zmax, 1.0))))
    integrand = c_km_s / H_z(grid, H0, omega_m)
    dz = np.diff(grid)
    dc_grid = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * dz)])
    return np.interp(z_values, grid, dc_grid)

def distance_modulus(D_L_mpc):
    D_L_mpc = np.asarray(D_L_mpc, dtype=float)
    return 5.0 * np.log10(np.maximum(D_L_mpc, 1e-300)) + 25.0

def bridge_b_mpc(H0):
    return A0 * c_km_s / float(H0)

def bias_shape(z, model):
    z = np.asarray(z, dtype=float)
    if model == "linear":
        return z
    if model == "stam_su":
        return z * (1.0 + 3.0 * z / 20.0)
    if model == "saturating":
        return z / (1.0 + z)
    if model == "log":
        return np.log1p(z)
    raise ValueError(model)

def flos_z(z, f0, f1, model):
    z = np.asarray(z, dtype=float)
    if model == "constant":
        return f0 + 0.0 * z
    if model == "one_plus_z":
        return f0 + f1 * z / (1.0 + z)
    if model == "log":
        return f0 + f1 * np.log1p(z)
    raise ValueError(model)

def model_mu(z, params, model_kind, bias_shape_name, flos_model):
    H0 = params["H0"]
    om = params["omega_m"]
    D_C = comoving_distance_interp(z, H0, om)
    if model_kind == "stam":
        D_C = D_C + bridge_b_mpc(H0) * flos_z(z, params.get("f0",0.0), params.get("f1",0.0), flos_model) * bias_shape(z, bias_shape_name)
    D_L = (1.0 + np.asarray(z)) * D_C
    return distance_modulus(D_L) + params.get("Mcal", 0.0)

def read_sn_csv(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing SN file: {path}")
    rows = []
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        required = {"z", "mu", "sigma_mu"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"SN CSV missing columns: {sorted(missing)}")
        for r in reader:
            try:
                z = float(r["z"])
                mu = float(r["mu"])
                sig = float(r["sigma_mu"])
                if np.isfinite(z) and np.isfinite(mu) and np.isfinite(sig) and sig > 0 and z > 0:
                    rows.append({"z": z, "mu": mu, "sigma_mu": sig})
            except Exception:
                continue
    rows.sort(key=lambda r: r["z"])
    return rows

def unpack_params(x, model_kind, fit_Mcal):
    if model_kind == "lcdm":
        H0, om = x[0], x[1]
        idx = 2
        f0 = 0.0; f1 = 0.0
    else:
        H0, om, f0, f1 = x[0], x[1], x[2], x[3]
        idx = 4
    Mcal = x[idx] if fit_Mcal else 0.0
    return {"H0": float(H0), "omega_m": float(om), "f0": float(f0), "f1": float(f1), "Mcal": float(Mcal)}

def bounds_for(model_kind, fit_Mcal):
    if model_kind == "lcdm":
        b = [(50.0, 90.0), (0.05, 0.60)]
    else:
        b = [(50.0, 90.0), (0.05, 0.60), (-10.0, 20.0), (-30.0, 30.0)]
    if fit_Mcal:
        b.append((-5.0, 5.0))
    return b

def chi2_model(x, model_kind, z, mu_obs, sig, bias_shape_name, flos_model, fit_Mcal):
    params = unpack_params(x, model_kind, fit_Mcal)
    if params["H0"] <= 0 or params["omega_m"] <= 0 or params["omega_m"] >= 1:
        return 1e99
    if model_kind == "stam":
        zcheck = np.array([0.0, 0.01, 0.1, 0.5, 1.0, 2.5])
        if np.min(flos_z(zcheck, params["f0"], params["f1"], flos_model)) < -1.0:
            return 1e99
    pred = model_mu(z, params, model_kind, bias_shape_name, flos_model)
    return float(np.sum(((pred - mu_obs) / sig)**2))

def fit_model(model_kind, z, mu_obs, sig, bias_shape_name, flos_model, fit_Mcal):
    bounds = bounds_for(model_kind, fit_Mcal)
    def obj(x):
        return chi2_model(x, model_kind, z, mu_obs, sig, bias_shape_name, flos_model, fit_Mcal)
    de = differential_evolution(obj, bounds=bounds, seed=13, polish=False, tol=1e-7)
    local = minimize(obj, de.x, method="Nelder-Mead", options={"maxiter": 10000, "xatol": 1e-8, "fatol": 1e-8})
    xbest = local.x if local.fun <= de.fun else de.x
    params = unpack_params(xbest, model_kind, fit_Mcal)
    chi2 = chi2_model(xbest, model_kind, z, mu_obs, sig, bias_shape_name, flos_model, fit_Mcal)
    k = len(xbest); N = len(z)
    return params, chi2, chi2 + 2*k, chi2 + k * math.log(max(N,1))

def bin_edges_default():
    return [
        (0.000, 0.010, "z < 0.01"),
        (0.010, 0.050, "0.01 ≤ z < 0.05"),
        (0.050, 0.150, "0.05 ≤ z < 0.15"),
        (0.150, 0.500, "0.15 ≤ z < 0.5"),
        (0.500, 1.000, "0.5 ≤ z < 1.0"),
        (1.000, 10.000, "z ≥ 1.0"),
    ]

def compute_bin_stats(residual_rows):
    stats = []
    for zlo, zhi, label in bin_edges_default():
        subset = [r for r in residual_rows if zlo <= r["z"] < zhi]
        if not subset:
            continue
        pulls_l = np.array([r["pull_lcdm"] for r in subset])
        pulls_s = np.array([r["pull_stam"] for r in subset])
        resid_l = np.array([r["resid_lcdm"] for r in subset])
        resid_s = np.array([r["resid_stam"] for r in subset])
        chi_l = float(np.sum(pulls_l**2)); chi_s = float(np.sum(pulls_s**2)); N = len(subset)
        stats.append({
            "bin": label, "z_min": zlo, "z_max": zhi, "N": N,
            "chi2_lcdm": chi_l, "chi2_stam": chi_s,
            "chi2_per_lcdm": chi_l/N, "chi2_per_stam": chi_s/N,
            "mean_resid_lcdm": float(np.mean(resid_l)), "mean_resid_stam": float(np.mean(resid_s)),
            "rms_resid_lcdm": float(np.sqrt(np.mean(resid_l**2))), "rms_resid_stam": float(np.sqrt(np.mean(resid_s**2))),
        })
    return stats

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/mnt/data/G104_run/data")
    ap.add_argument("--outdir", default="/mnt/data/G106_run")
    ap.add_argument("--bias-shape", choices=["linear","stam_su","saturating","log"], default="linear")
    ap.add_argument("--flos-model", choices=["constant","one_plus_z","log"], default="one_plus_z")
    ap.add_argument("--no-fit-Mcal", action="store_true")
    args = ap.parse_args()
    data_dir = Path(args.data_dir); outdir = Path(args.outdir)
    results = outdir / "results"; plots = outdir / "plots"
    results.mkdir(parents=True, exist_ok=True); plots.mkdir(parents=True, exist_ok=True)
    sn_path = data_dir / "sn.csv"
    sn_rows = read_sn_csv(sn_path)
    z = np.array([r["z"] for r in sn_rows]); mu_obs = np.array([r["mu"] for r in sn_rows]); sig = np.array([r["sigma_mu"] for r in sn_rows])
    fit_Mcal = not args.no_fit_Mcal
    params_lcdm, chi_lcdm, aic_lcdm, bic_lcdm = fit_model("lcdm", z, mu_obs, sig, args.bias_shape, args.flos_model, fit_Mcal)
    params_stam, chi_stam, aic_stam, bic_stam = fit_model("stam", z, mu_obs, sig, args.bias_shape, args.flos_model, fit_Mcal)
    pred_lcdm = model_mu(z, params_lcdm, "lcdm", args.bias_shape, args.flos_model)
    pred_stam = model_mu(z, params_stam, "stam", args.bias_shape, args.flos_model)
    residual_rows = []
    for i in range(len(z)):
        residual_rows.append({
            "index": i, "z": float(z[i]), "mu_obs": float(mu_obs[i]), "sigma_mu": float(sig[i]),
            "mu_lcdm": float(pred_lcdm[i]), "mu_stam": float(pred_stam[i]),
            "resid_lcdm": float(mu_obs[i] - pred_lcdm[i]), "resid_stam": float(mu_obs[i] - pred_stam[i]),
            "pull_lcdm": float((pred_lcdm[i] - mu_obs[i]) / sig[i]), "pull_stam": float((pred_stam[i] - mu_obs[i]) / sig[i]),
        })
    bin_stats = compute_bin_stats(residual_rows)
    # CSVs
    params_csv = results / "G106_SN_best_fit_params.csv"
    with params_csv.open("w", newline="", encoding="utf-8") as f:
        fields = ["model","H0","omega_m","f0","f1","Mcal","chi2","N","chi2_per_N","AIC","BIC"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for name,p,chi,aic,bic in [("lcdm",params_lcdm,chi_lcdm,aic_lcdm,bic_lcdm),("stam",params_stam,chi_stam,aic_stam,bic_stam)]:
            w.writerow({"model":name,"H0":p["H0"],"omega_m":p["omega_m"],"f0":p.get("f0",0),"f1":p.get("f1",0),"Mcal":p.get("Mcal",0),"chi2":chi,"N":len(z),"chi2_per_N":chi/len(z),"AIC":aic,"BIC":bic})
    resid_csv = results / "G106_SN_row_residuals.csv"
    with resid_csv.open("w", newline="", encoding="utf-8") as f:
        fields = list(residual_rows[0].keys()); w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(residual_rows)
    bin_csv = results / "G106_SN_bin_stats.csv"
    with bin_csv.open("w", newline="", encoding="utf-8") as f:
        fields = list(bin_stats[0].keys()); w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(bin_stats)
    # plots
    order = np.argsort(z); z_sorted = z[order]
    res_l = np.array([r["resid_lcdm"] for r in residual_rows])[order]; res_s = np.array([r["resid_stam"] for r in residual_rows])[order]
    plt.figure(figsize=(11,6)); plt.scatter(z_sorted,res_l,s=8,alpha=0.35,label="SN - LCDM"); plt.scatter(z_sorted,res_s,s=8,alpha=0.35,label="SN - STAM")
    plt.axhline(0,color="black",lw=0.8); plt.xlabel("z"); plt.ylabel("distance modulus residual"); plt.title("G106 SN residuals vs redshift"); plt.grid(alpha=0.3); plt.legend()
    p1 = plots / "G106_SN_residuals_vs_z.png"; plt.tight_layout(); plt.savefig(p1,dpi=200); plt.close()
    labels = [b["bin"] for b in bin_stats]; x = np.arange(len(labels))
    lcdm_vals = [b["chi2_per_lcdm"] for b in bin_stats]; stam_vals = [b["chi2_per_stam"] for b in bin_stats]
    width=0.35; plt.figure(figsize=(11,6)); plt.bar(x-width/2,lcdm_vals,width,label="LCDM"); plt.bar(x+width/2,stam_vals,width,label="STAM")
    plt.xticks(x,labels,rotation=25,ha="right"); plt.ylabel("chi2 / N"); plt.title("G106 SN redshift-bin chi2/N"); plt.grid(axis="y",alpha=0.3); plt.legend()
    p2=plots/"G106_SN_bin_chi2.png"; plt.tight_layout(); plt.savefig(p2,dpi=200); plt.close()
    zgrid=np.linspace(0,max(2.5,float(np.max(z))),300); plt.figure(figsize=(10,6)); plt.plot(zgrid,flos_z(zgrid,params_stam["f0"],params_stam["f1"],args.flos_model)); plt.axhline(0,color="black",lw=0.8)
    plt.xlabel("z"); plt.ylabel("f_los(z)"); plt.title("G106 best-fit STAM line-of-sight amplification"); plt.grid(alpha=0.3)
    p3=plots/"G106_SN_flos_shape.png"; plt.tight_layout(); plt.savefig(p3,dpi=200); plt.close()
    # console
    print("="*100); print("G106: SN-only redshift-split residual test"); print("="*100)
    print(f"SN file: {sn_path}"); print(f"N = {len(z)}, z range = {np.min(z):.5f} .. {np.max(z):.5f}\n")
    print(f"LCDM: H0={params_lcdm['H0']:.4f}, Om={params_lcdm['omega_m']:.4f}, Mcal={params_lcdm.get('Mcal',0):.4f}, chi2={chi_lcdm:.3f}, chi2/N={chi_lcdm/len(z):.3f}, AIC={aic_lcdm:.3f}, BIC={bic_lcdm:.3f}")
    print(f"STAM: H0={params_stam['H0']:.4f}, Om={params_stam['omega_m']:.4f}, f0={params_stam['f0']:.4f}, f1={params_stam['f1']:.4f}, Mcal={params_stam.get('Mcal',0):.4f}, chi2={chi_stam:.3f}, chi2/N={chi_stam/len(z):.3f}, AIC={aic_stam:.3f}, BIC={bic_stam:.3f}\n")
    print(f"Delta chi2 LCDM-STAM = {chi_lcdm-chi_stam:.3f}"); print(f"Delta AIC  LCDM-STAM = {aic_lcdm-aic_stam:.3f}"); print(f"Delta BIC  LCDM-STAM = {bic_lcdm-bic_stam:.3f}\n")
    print(f"{'bin':<18}{'N':>8}{'LCDM chi2/N':>16}{'STAM chi2/N':>16}{'Δ chi2/N':>12}"); print("-"*74)
    for b in bin_stats:
        delta = b["chi2_per_lcdm"] - b["chi2_per_stam"]
        print(f"{b['bin']:<18}{b['N']:>8}{b['chi2_per_lcdm']:>16.3f}{b['chi2_per_stam']:>16.3f}{delta:>12.3f}")
    # summary
    summary=results/"G106_SN_redshift_split_summary.md"
    md=f"""# G106 — SN-only redshift-split residual test

## Inputs

- SN file: `{sn_path}`
- N = `{len(z)}`
- z range = `{np.min(z):.5f}` to `{np.max(z):.5f}`
- bias shape = `{args.bias_shape}`
- f_los model = `{args.flos_model}`
- fit Mcal = `{fit_Mcal}`

## Best fits

| Model | H0 | Omega_m | f0 | f1 | Mcal | chi2 | chi2/N | AIC | BIC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LCDM | {params_lcdm['H0']:.4f} | {params_lcdm['omega_m']:.4f} | 0 | 0 | {params_lcdm.get('Mcal',0):.4f} | {chi_lcdm:.3f} | {chi_lcdm/len(z):.3f} | {aic_lcdm:.3f} | {bic_lcdm:.3f} |
| STAM | {params_stam['H0']:.4f} | {params_stam['omega_m']:.4f} | {params_stam['f0']:.4f} | {params_stam['f1']:.4f} | {params_stam.get('Mcal',0):.4f} | {chi_stam:.3f} | {chi_stam/len(z):.3f} | {aic_stam:.3f} | {bic_stam:.3f} |

- Δχ² LCDM−STAM = `{chi_lcdm-chi_stam:.3f}`
- ΔAIC LCDM−STAM = `{aic_lcdm-aic_stam:.3f}`
- ΔBIC LCDM−STAM = `{bic_lcdm-bic_stam:.3f}`

## Redshift-bin stats

| Bin | N | LCDM chi2/N | STAM chi2/N | Δ chi2/N |
|---|---:|---:|---:|---:|
"""
    for b in bin_stats:
        delta=b["chi2_per_lcdm"]-b["chi2_per_stam"]
        md += f"| {b['bin']} | {b['N']} | {b['chi2_per_lcdm']:.3f} | {b['chi2_per_stam']:.3f} | {delta:.3f} |\n"
    md += f"""
## Interpretation guide

- Smooth improvement across several bins supports a redshift-structured distance-bias interpretation.
- Improvement only in one bin suggests a localized artifact or nuisance absorption.
- If STAM improves raw chi2 but loses AIC/BIC, it is competitive but not statistically preferred under this parameterization.

## Files

- Best-fit params: `{params_csv}`
- Row residuals: `{resid_csv}`
- Bin stats: `{bin_csv}`
- Plot: `{p1}`
- Plot: `{p2}`
- Plot: `{p3}`
"""
    summary.write_text(md,encoding="utf-8")
    print(f"\nSummary written: {summary}"); print(f"Params CSV: {params_csv}"); print(f"Bin stats CSV: {bin_csv}")

if __name__ == "__main__":
    main()

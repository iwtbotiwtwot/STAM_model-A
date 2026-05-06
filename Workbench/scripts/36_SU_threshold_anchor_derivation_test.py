#!/usr/bin/env python3
"""
Script 36 — SU threshold-anchor derivation test

Accumulation Theory / STAM Model-A

Purpose:
    Test whether the SU distance structure is naturally tied to the original
    z≈0.300 threshold/flag anchor rather than being an arbitrary polynomial.

Clean SU form:
    SU(z) = (z/H)(1 + 3z/20)
          = H^-1 z(1 + 0.15z)

Spreadsheet/factored form:
    SU(z) = K(a z + q z^2)

with:
    a = 3.33333 = 1/0.300
    q = 0.5
    K = H^-1/a

Therefore:
    beta = q/a = q*z_anchor
    for q=0.5 and z_anchor=0.300:
        beta = 0.15

Core question:
    Do the catalogs prefer a curvature beta near 0.15?
    If q is fixed at 0.5, does the implied z_anchor = beta/q land near 0.300?
    If z_anchor is fixed at 0.300, does q land near 0.5?

Important:
    With an additive magnitude offset b_mu, the catalog primarily constrains
    the shape beta = q/a, not the absolute SU scale. This script reports both
    with-b and no-b diagnostics.
"""

from pathlib import Path
import argparse, json, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from scipy.optimize import minimize_scalar
    SCIPY_OK = True
except Exception:
    SCIPY_OK = False

H = 0.000243635
L_MPC = 1.0 / H

def mu_from_mpc(D):
    D = np.maximum(np.asarray(D, dtype=float), 1e-12)
    return 5*np.log10(D) + 25.0

def D_SU_beta(z, beta):
    return L_MPC * z * (1.0 + beta*z)

def D_SU_anchor_q(z, z_anchor, q=0.5):
    beta = q*z_anchor
    return D_SU_beta(z, beta)

def find_col(cols, names):
    for target in names:
        for c in cols:
            if c.lower() == target.lower():
                return c
    for target in names:
        for c in cols:
            if target.lower() in c.lower():
                return c
    return None

def infer_name(path):
    n = Path(path).stem.lower()
    if "union" in n: return "Union3"
    if "des" in n: return "DES"
    if "pantheon" in n: return "Pantheon"
    return Path(path).stem

def load_catalog(path):
    path = Path(path)
    df = pd.read_csv(path)
    cols = list(df.columns)
    zcol = find_col(cols, ["z", "zHD", "zCMB", "zcmb", "redshift"])
    mucol = find_col(cols, ["MU", "MU_SH0ES", "mu", "distance_modulus", "distmod", "mb"])
    if zcol is None or mucol is None:
        return None
    out = pd.DataFrame({
        "catalog": infer_name(path),
        "z": pd.to_numeric(df[zcol], errors="coerce"),
        "mu": pd.to_numeric(df[mucol], errors="coerce"),
        "source_file": path.name,
    })
    out = out.replace([np.inf,-np.inf], np.nan).dropna(subset=["z","mu"])
    out = out[out["z"] > 0].copy()
    return out

def fit_beta(df, include_b=True, beta_bounds=(-1.0, 1.0)):
    z = df["z"].values.astype(float)
    y = df["mu"].values.astype(float)

    def eval_beta(beta):
        mu = mu_from_mpc(D_SU_beta(z, beta))
        b = float(np.mean(y-mu)) if include_b else 0.0
        pred = mu + b
        sse = float(np.sum((y-pred)**2))
        return sse, b, pred

    if SCIPY_OK:
        opt = minimize_scalar(lambda b: eval_beta(b)[0], bounds=beta_bounds, method="bounded", options={"xatol":1e-12})
        beta = float(opt.x)
    else:
        grid = np.linspace(beta_bounds[0], beta_bounds[1], 20001)
        vals = [eval_beta(b)[0] for b in grid]
        beta = float(grid[int(np.argmin(vals))])
    sse, bmu, pred = eval_beta(beta)
    n = len(y); k = 1 + (1 if include_b else 0)
    rmse = float(np.sqrt(np.mean((y-pred)**2)))
    aic = n*np.log(max(sse/n, 1e-30)) + 2*k
    return {"beta": beta, "b_mu": bmu, "sse": sse, "rmse": rmse, "aic": float(aic), "n": n, "k": k}

def eval_fixed_beta(df, beta, include_b=True):
    z = df["z"].values.astype(float)
    y = df["mu"].values.astype(float)
    mu = mu_from_mpc(D_SU_beta(z, beta))
    b = float(np.mean(y-mu)) if include_b else 0.0
    pred = mu + b
    sse = float(np.sum((y-pred)**2))
    n = len(y); k = (1 if include_b else 0)
    rmse = float(np.sqrt(np.mean((y-pred)**2)))
    aic = n*np.log(max(sse/n,1e-30)) + 2*k
    return {"beta": beta, "b_mu": b, "sse": sse, "rmse": rmse, "aic": float(aic), "n": n, "k": k}

def run(inputs, outdir, zmin=0.05):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dfs = [load_catalog(p) for p in inputs]
    dfs = [d for d in dfs if d is not None]
    all_df = pd.concat(dfs, ignore_index=True)
    all_df = all_df[all_df["z"] >= zmin].copy()
    all_df.to_csv(outdir/"script36_loaded_catalog_points.csv", index=False)

    datasets = {c: all_df[all_df["catalog"]==c].copy() for c in sorted(all_df["catalog"].unique())}
    if "Pantheon" in datasets and "Union3" in datasets:
        datasets["Pantheon_plus_Union3"] = pd.concat([datasets["Pantheon"], datasets["Union3"]], ignore_index=True)
    if "DES" in datasets and "Pantheon" in datasets:
        zlo = max(datasets["DES"]["z"].min(), datasets["Pantheon"]["z"].min(), zmin)
        zhi = min(datasets["DES"]["z"].max(), datasets["Pantheon"]["z"].max())
        datasets["DES_plus_Pantheon_common"] = pd.concat([
            datasets["DES"][(datasets["DES"]["z"]>=zlo)&(datasets["DES"]["z"]<=zhi)],
            datasets["Pantheon"][(datasets["Pantheon"]["z"]>=zlo)&(datasets["Pantheon"]["z"]<=zhi)]
        ], ignore_index=True)
    datasets["All_three"] = all_df.copy()

    rows = []
    fixed_betas = {
        "anchored_z0_0p300_q_0p5_beta_0p150": 0.15,
        "z0_one_third_q_0p5_beta_0p166667": 1/6,
        "z0_0p350_q_0p5_beta_0p175": 0.175,
        "no_quadratic_beta_0": 0.0,
    }

    for name, df in datasets.items():
        if len(df) < 8: 
            continue

        for include_b in [False, True]:
            free = fit_beta(df, include_b=include_b, beta_bounds=(-1.0, 1.0))
            beta = free["beta"]
            implied_anchor_q05 = beta/0.5
            implied_q_anchor03 = beta/0.300
            rows.append({
                "dataset": name,
                "model": "free_beta",
                "include_b_mu": include_b,
                "beta": beta,
                "z_anchor_if_q_0p5": implied_anchor_q05,
                "q_if_z_anchor_0p300": implied_q_anchor03,
                "distance_from_beta_0p15": beta-0.15,
                "distance_from_z_anchor_0p300": implied_anchor_q05-0.300,
                **{k:v for k,v in free.items() if k not in ["beta"]}
            })

            for label, fixed_beta in fixed_betas.items():
                rec = eval_fixed_beta(df, fixed_beta, include_b=include_b)
                rows.append({
                    "dataset": name,
                    "model": label,
                    "include_b_mu": include_b,
                    "beta": fixed_beta,
                    "z_anchor_if_q_0p5": fixed_beta/0.5,
                    "q_if_z_anchor_0p300": fixed_beta/0.300,
                    "distance_from_beta_0p15": fixed_beta-0.15,
                    "distance_from_z_anchor_0p300": fixed_beta/0.5-0.300,
                    **{k:v for k,v in rec.items() if k not in ["beta"]}
                })

    results = pd.DataFrame(rows)
    results.to_csv(outdir/"script36_SU_anchor_fit_results.csv", index=False)

    # Best/free vs anchored comparison
    comp_rows = []
    for (dname, incb), g in results.groupby(["dataset","include_b_mu"]):
        free = g[g["model"]=="free_beta"].iloc[0]
        anchor = g[g["model"]=="anchored_z0_0p300_q_0p5_beta_0p150"].iloc[0]
        z0333 = g[g["model"]=="z0_one_third_q_0p5_beta_0p166667"].iloc[0]
        z035 = g[g["model"]=="z0_0p350_q_0p5_beta_0p175"].iloc[0]
        comp_rows.append({
            "dataset": dname,
            "include_b_mu": incb,
            "free_beta": free["beta"],
            "free_implied_z_anchor_q0p5": free["z_anchor_if_q_0p5"],
            "free_implied_q_z0p300": free["q_if_z_anchor_0p300"],
            "free_rmse": free["rmse"],
            "anchor_0p300_rmse": anchor["rmse"],
            "anchor_1third_rmse": z0333["rmse"],
            "anchor_0p350_rmse": z035["rmse"],
            "anchor_0p300_minus_free_rmse": anchor["rmse"]-free["rmse"],
            "anchor_1third_minus_free_rmse": z0333["rmse"]-free["rmse"],
            "anchor_0p350_minus_free_rmse": z035["rmse"]-free["rmse"],
            "anchor_0p300_delta_aic_vs_free": anchor["aic"]-free["aic"],
            "anchor_1third_delta_aic_vs_free": z0333["aic"]-free["aic"],
            "anchor_0p350_delta_aic_vs_free": z035["aic"]-free["aic"],
        })
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(outdir/"script36_anchor_vs_free_summary.csv", index=False)

    # Plots: implied anchor
    plot_df = comp[comp["include_b_mu"]==True].copy()
    plt.figure(figsize=(10,5))
    plt.bar(plot_df["dataset"], plot_df["free_implied_z_anchor_q0p5"])
    plt.axhline(0.300, linestyle="--", label="0.300 original flag")
    plt.axhline(1/3, linestyle="--", label="1/3")
    plt.axhline(0.350, linestyle="--", label="0.350")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("implied z_anchor if q=0.5")
    plt.title("Script 36: free shape implies threshold anchor")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir/"script36_implied_z_anchor_with_b.png", dpi=180)
    plt.close()

    plt.figure(figsize=(10,5))
    x = np.arange(len(plot_df))
    width = 0.25
    plt.bar(x-width, plot_df["anchor_0p300_minus_free_rmse"], width, label="0.300 anchor")
    plt.bar(x, plot_df["anchor_1third_minus_free_rmse"], width, label="1/3 anchor")
    plt.bar(x+width, plot_df["anchor_0p350_minus_free_rmse"], width, label="0.350 anchor")
    plt.axhline(0, linewidth=1)
    plt.xticks(x, plot_df["dataset"], rotation=45, ha="right")
    plt.ylabel("RMSE penalty vs free beta")
    plt.title("Script 36: penalty for fixed threshold anchors")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir/"script36_fixed_anchor_rmse_penalty.png", dpi=180)
    plt.close()

    # For P+U family summary
    fam_names = ["Pantheon", "Union3", "Pantheon_plus_Union3"]
    fam = comp[(comp["include_b_mu"]==True)&(comp["dataset"].isin(fam_names))]
    family_summary = {
        "family": fam_names,
        "mean_free_beta": float(fam["free_beta"].mean()) if len(fam) else None,
        "std_free_beta": float(fam["free_beta"].std()) if len(fam)>1 else None,
        "mean_implied_z_anchor_q0p5": float(fam["free_implied_z_anchor_q0p5"].mean()) if len(fam) else None,
        "std_implied_z_anchor_q0p5": float(fam["free_implied_z_anchor_q0p5"].std()) if len(fam)>1 else None,
        "mean_implied_q_z0p300": float(fam["free_implied_q_z0p300"].mean()) if len(fam) else None,
        "std_implied_q_z0p300": float(fam["free_implied_q_z0p300"].std()) if len(fam)>1 else None,
        "mean_anchor_0p300_rmse_penalty": float(fam["anchor_0p300_minus_free_rmse"].mean()) if len(fam) else None,
        "mean_anchor_1third_rmse_penalty": float(fam["anchor_1third_minus_free_rmse"].mean()) if len(fam) else None,
        "mean_anchor_0p350_rmse_penalty": float(fam["anchor_0p350_minus_free_rmse"].mean()) if len(fam) else None,
    }

    verdict = {
        "script": "36_SU_threshold_anchor_derivation_test.py",
        "purpose": "Test whether SU(z)=(z/H)(1+3z/20) is naturally tied to the z≈0.300 threshold anchor.",
        "core_relations": {
            "SU_clean": "SU(z)=(z/H)(1+3z/20)=H^-1 z(1+0.15z)",
            "factored": "SU(z)=K(a z + q z^2)",
            "anchor": "a=1/z_anchor",
            "beta": "beta=q/a=q*z_anchor",
            "original": "z_anchor=0.300 and q=0.5 gives beta=0.15"
        },
        "important_limitation": "With additive b_mu, SN catalogs mainly constrain beta=q/a, not absolute SU scale.",
        "family_summary_Pantheon_Union": family_summary,
        "interpretation_rules": [
            "If free beta implies z_anchor near 0.300 for q=0.5, the original anchor is supported.",
            "If anchored beta=0.15 performs close to free beta, the original SU form is robust.",
            "If 1/3 or 0.35 anchors perform better, the original 0.300 flag may be a rough marker rather than exact."
        ]
    }
    with open(outdir/"script36_verdict.json","w") as f:
        json.dump(verdict, f, indent=2)

    print("Script 36 SU threshold-anchor derivation test complete.")
    print(f"Output directory: {outdir}")
    print()
    print("Anchor vs free summary:")
    print(comp.to_string(index=False))
    print()
    print("Pantheon/Union family summary:")
    print(json.dumps(family_summary, indent=2))
    print()
    print("Main read:")
    print("  The test asks whether the data-preferred beta implies z_anchor≈0.300 when q=0.5.")
    print("  It also measures how much RMSE penalty comes from forcing the original beta=0.15 anchor.")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--inputs", nargs="+", required=True)
    p.add_argument("--outdir", default="results/script36_SU_threshold_anchor")
    p.add_argument("--zmin", type=float, default=0.05)
    args = p.parse_args()
    run(args.inputs, args.outdir, zmin=args.zmin)

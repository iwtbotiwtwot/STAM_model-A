#!/usr/bin/env python3
"""
G102_structure_dependent_flos_fast.py

Fast STAM / Model-A structure-dependent line-of-sight amplification test.

G101 result:
    one constant f_los is too crude.

G102 test:
    replace constant f_los with simple redshift-dependent forms and check whether
    the compressed SN + BAO + CMB distance fit improves.

Locked convention:
    D_C_obs(z) = D_C_intrinsic(z; H0_true) + b * f_los(z) * shape(z)
    D_L_obs(z) = (1+z) * D_C_obs(z)

Direct H(z) and BAO radial D_H are NOT distance-biased.

This is a compressed toy likelihood, not a publication-grade cosmology fit.
"""

from __future__ import annotations

import argparse, csv, math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

c_km_s = 299_792.458
MPC_TO_MLY = 3.261563776
A0 = 1.0 / (12.0 * math.pi)


def E(z, om):
    z = np.asarray(z, dtype=float)
    return np.sqrt(om*(1+z)**3 + (1-om))


def H(z, H0, om):
    return H0 * E(z, om)


def DC(z, H0, om):
    """Fast trapezoid integration over sorted z array."""
    z = np.asarray(z, dtype=float)
    z_full = np.concatenate([[0.0], z])
    integ = c_km_s / H(z_full, H0, om)
    dz = np.diff(z_full)
    return np.cumsum(0.5*(integ[1:]+integ[:-1])*dz)


def DL(z, dc):
    return (1+np.asarray(z))*dc


def mu(dl):
    return 5*np.log10(dl) + 25


def b_mpc(H0):
    return A0*c_km_s/H0


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


def flos_z(z, model, p):
    z = np.asarray(z, dtype=float)
    if model == "constant":
        return p[0] + 0*z
    if model == "one_plus_z":
        return p[0] + p[1]*z/(1+z)
    if model == "log":
        return p[0] + p[1]*np.log1p(z)
    raise ValueError(model)


def biased_DC(z, H0_true, om, shape, flos_model, p):
    return DC(z, H0_true, om) + b_mpc(H0_true)*flos_z(z, flos_model, p)*shape_z(z, shape)


def components(H0_true, H0_ref, om, shape, flos_model, p, setup):
    out = {}
    # SN
    z = setup["z_sn"]
    ref = mu(DL(z, DC(z, H0_ref, om)))
    mod = mu(DL(z, biased_DC(z, H0_true, om, shape, flos_model, p)))
    out["SN"] = (float(np.sum(((mod-ref)/setup["sigma_mu"])**2)), len(z))
    # BAO transverse
    z = setup["z_bao"]
    ref = DC(z, H0_ref, om)
    mod = biased_DC(z, H0_true, om, shape, flos_model, p)
    out["BAO_DM"] = (float(np.sum(((mod-ref)/(setup["sigma_bao_DM_frac"]*ref))**2)), len(z))
    # BAO radial
    ref = c_km_s/H(z, H0_ref, om)
    mod = c_km_s/H(z, H0_true, om)
    out["BAO_DH"] = (float(np.sum(((mod-ref)/(setup["sigma_bao_DH_frac"]*ref))**2)), len(z))
    # CMB distance
    z = np.array([setup["zstar"]])
    ref = DC(z, H0_ref, om)[0]
    mod = biased_DC(z, H0_true, om, shape, flos_model, p)[0]
    out["CMB_D"] = (float(((mod-ref)/(setup["sigma_cmb_D_frac"]*ref))**2), 1)
    # H(z)
    z = setup["z_hz"]
    ref = H(z, H0_ref, om)
    mod = H(z, H0_true, om)
    out["H_z"] = (float(np.sum(((mod-ref)/(setup["sigma_H_frac"]*ref))**2)), len(z))
    total = sum(v[0] for v in out.values())
    N = sum(v[1] for v in out.values())
    out["TOTAL"] = (total, N)
    return out


def score(H0_true, H0_ref, om, shape, flos_model, p, setup, include):
    c = components(H0_true, H0_ref, om, shape, flos_model, p, setup)
    return sum(c[k][0] for k in include)


def grid_fit(H0_true, H0_ref, om, shape, flos_model, setup, include):
    if flos_model == "constant":
        best = None
        for f0 in np.linspace(-2, 8, 501):
            p = np.array([f0])
            s = score(H0_true,H0_ref,om,shape,flos_model,p,setup,include)
            if best is None or s < best[0]:
                best = (s,p)
        return best[1], components(H0_true,H0_ref,om,shape,flos_model,best[1],setup)
    # coarse-to-fine 2D grid for f0,f1
    best = None
    for span, n in [(((-2,8),(-12,12)), 81), (None, 81)]:
        if span is None:
            f0c, f1c = best[1]
            span = ((f0c-1.0, f0c+1.0), (f1c-2.0, f1c+2.0))
        f0s = np.linspace(span[0][0], span[0][1], n)
        f1s = np.linspace(span[1][0], span[1][1], n)
        for f0 in f0s:
            for f1 in f1s:
                p = np.array([f0,f1])
                # mild physical penalty if f_los goes badly negative over z=[0,1100]
                zcheck = np.array([0.0, 0.5, 2.0, 10.0, 1100.0])
                if np.min(flos_z(zcheck, flos_model, p)) < -0.5:
                    continue
                s = score(H0_true,H0_ref,om,shape,flos_model,p,setup,include)
                if best is None or s < best[0]:
                    best = (s,p)
    return best[1], components(H0_true,H0_ref,om,shape,flos_model,best[1],setup)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="/mnt/data/G102_fast_run")
    ap.add_argument("--H0-true", type=float, default=73.04)
    ap.add_argument("--H0-ref", type=float, default=67.40)
    ap.add_argument("--omega-m", type=float, default=0.30)
    args = ap.parse_args()
    root = Path(args.outdir); results = root/"results"; plots = root/"plots"
    results.mkdir(parents=True, exist_ok=True); plots.mkdir(parents=True, exist_ok=True)
    H0_true, H0_ref, om = args.H0_true, args.H0_ref, args.omega_m
    setup = {
        "z_sn": np.linspace(0.01,2.0,80), "sigma_mu":0.10,
        "z_bao": np.array([0.38,0.51,0.61,1.48,2.33]),
        "sigma_bao_DM_frac":0.02, "sigma_bao_DH_frac":0.03,
        "zstar":1089.0, "sigma_cmb_D_frac":0.005,
        "z_hz":np.array([0.1,0.2,0.35,0.5,0.7,1.0,1.5,2.0]), "sigma_H_frac":0.07,
    }
    shapes = ["linear","saturating","log"]  # omit stam_su by default; it overgrows at high z in G101
    models = ["constant","one_plus_z","log"]
    fit_sets = {
        "distance_only": ["SN","BAO_DM","CMB_D"],
        "all_probes": ["SN","BAO_DM","CMB_D","BAO_DH","H_z"],
        "SN_BAO": ["SN","BAO_DM"],
        "CMB_only": ["CMB_D"],
    }
    rows = []
    for shape in shapes:
        for model in models:
            for fitset, include in fit_sets.items():
                p, comps = grid_fit(H0_true,H0_ref,om,shape,model,setup,include)
                row = {"shape":shape, "flos_model":model, "fit_set":fitset,
                       "params":";".join(f"{x:.6g}" for x in p)}
                for i,x in enumerate(p): row[f"p{i}"] = float(x)
                for k,(chi2,N) in comps.items():
                    row[f"chi2_{k}"] = chi2
                    row[f"N_{k}"] = N
                    row[f"chi2_per_{k}"] = chi2/N
                rows.append(row)
    csv_path = results/"G102_structure_dependent_flos_grid.csv"
    fields = sorted(set().union(*(r.keys() for r in rows)))
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w=csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    best_dist = min([r for r in rows if r["fit_set"]=="distance_only"], key=lambda r:r["chi2_per_TOTAL"])
    best_all = min([r for r in rows if r["fit_set"]=="all_probes"], key=lambda r:r["chi2_per_TOTAL"])
    print("="*100)
    print("G102 FAST: structure-dependent f_los")
    print("="*100)
    print(f"A0={A0:.12f}, H0_true={H0_true}, H0_ref={H0_ref}, b={b_mpc(H0_true):.6f} Mpc")
    print()
    print(f"{'shape':<12}{'model':<12}{'fit':<14}{'params':<20}{'total':>8}{'SN':>8}{'BAO_DM':>9}{'BAO_DH':>9}{'CMB':>9}{'Hz':>8}")
    print("-"*110)
    for r in rows:
        if r["fit_set"] in ("distance_only","all_probes"):
            print(f"{r['shape']:<12}{r['flos_model']:<12}{r['fit_set']:<14}{r['params']:<20}"
                  f"{r['chi2_per_TOTAL']:>8.2f}{r['chi2_per_SN']:>8.2f}{r['chi2_per_BAO_DM']:>9.2f}"
                  f"{r['chi2_per_BAO_DH']:>9.2f}{r['chi2_per_CMB_D']:>9.2f}{r['chi2_per_H_z']:>8.2f}")
    print()
    print("Best distance-only:", best_dist)
    print("Best all-probes:", best_all)
    # plots
    dist_rows=[r for r in rows if r["fit_set"]=="distance_only"]
    labels=[f"{r['shape']}\n{r['flos_model']}" for r in dist_rows]
    vals=[r["chi2_per_TOTAL"] for r in dist_rows]
    plt.figure(figsize=(12,6)); plt.bar(range(len(vals)), vals)
    plt.xticks(range(len(vals)), labels, rotation=65, ha="right", fontsize=8)
    plt.yscale("log"); plt.ylabel("chi2/N total"); plt.title("G102 distance-only fits")
    plt.grid(axis="y", alpha=0.3); plt.tight_layout()
    p1=plots/"G102_chi2_comparison.png"; plt.savefig(p1,dpi=200); plt.close()
    z=np.linspace(0,5,400)
    plt.figure(figsize=(10,6))
    for name,r in [("best distance-only",best_dist),("best all-probes",best_all)]:
        p=np.array([float(x) for x in r["params"].split(";")])
        plt.plot(z, flos_z(z,r["flos_model"],p), label=f"{name}: {r['shape']}/{r['flos_model']}")
    plt.axhline(0,color="black",lw=0.8); plt.xlabel("z"); plt.ylabel("f_los(z)")
    plt.title("G102 best structure-dependent f_los(z)"); plt.grid(alpha=0.3); plt.legend()
    p2=plots/"G102_best_flos_of_z.png"; plt.tight_layout(); plt.savefig(p2,dpi=200); plt.close()
    zres=np.linspace(0.01,2.0,300)
    ref=mu(DL(zres, DC(zres,H0_ref,om)))
    plt.figure(figsize=(10,6))
    for name,r in [("best distance-only",best_dist),("best all-probes",best_all)]:
        p=np.array([float(x) for x in r["params"].split(";")])
        mod=mu(DL(zres, biased_DC(zres,H0_true,om,r["shape"],r["flos_model"],p)))
        plt.plot(zres, mod-ref, label=f"{name}: {r['shape']}/{r['flos_model']}")
    plt.axhline(0,color="black",lw=0.8); plt.xlabel("z")
    plt.ylabel("mu_STAM_biased - mu_ref_lowerH0")
    plt.title("G102 residuals vs lower-H0 reference"); plt.grid(alpha=0.3); plt.legend(fontsize=8)
    p3=plots/"G102_best_residuals.png"; plt.tight_layout(); plt.savefig(p3,dpi=200); plt.close()
    # summary
    summary=results/"G102_structure_dependent_flos_summary.md"
    md=f"""# G102 — structure-dependent f_los test

## Purpose

G102 replaces the single global `f_los` from G101 with simple redshift-dependent forms:

```text
constant:      f_los(z)=f0
one_plus_z:    f_los(z)=f0+f1 z/(1+z)
log:           f_los(z)=f0+f1 ln(1+z)
```

## Locked convention

```text
D_C_obs(z) = D_C_intrinsic(z; H0_true) + b*f_los(z)*shape(z)
D_L_obs(z) = (1+z)*D_C_obs(z)
Direct H(z) probes are not distance-biased.
```

## Constants

- `A0 = {A0:.12f}`
- `H0_true = {H0_true:.4f}`
- `H0_ref = {H0_ref:.4f}`
- `b = {b_mpc(H0_true):.6f}` Mpc = `{b_mpc(H0_true)*MPC_TO_MLY:.6f}` Mly

## Best distance-only fit

- shape: `{best_dist['shape']}`
- f_los model: `{best_dist['flos_model']}`
- params: `{best_dist['params']}`
- chi2/N total: `{best_dist['chi2_per_TOTAL']:.3f}`
- SN: `{best_dist['chi2_per_SN']:.3f}`
- BAO_DM: `{best_dist['chi2_per_BAO_DM']:.3f}`
- CMB_D: `{best_dist['chi2_per_CMB_D']:.3f}`

## Best all-probes fit

- shape: `{best_all['shape']}`
- f_los model: `{best_all['flos_model']}`
- params: `{best_all['params']}`
- chi2/N total: `{best_all['chi2_per_TOTAL']:.3f}`
- SN: `{best_all['chi2_per_SN']:.3f}`
- BAO_DM: `{best_all['chi2_per_BAO_DM']:.3f}`
- BAO_DH: `{best_all['chi2_per_BAO_DH']:.3f}`
- CMB_D: `{best_all['chi2_per_CMB_D']:.3f}`
- H_z: `{best_all['chi2_per_H_z']:.3f}`

## Interpretation

If a structure-dependent `f_los(z)` improves SN/BAO_DM/CMB_D but BAO_DH and H_z remain strained, this supports the interpretation that the mechanism is a distance-bias layer, not an expansion-rate fix.

This remains a toy compressed likelihood. A real test must use actual SN covariance, BAO likelihoods, CMB acoustic angle, chronometers, structure growth, and locked priors.

## Files

- CSV: `{csv_path}`
- Plot: `{p1}`
- Plot: `{p2}`
- Plot: `{p3}`
"""
    summary.write_text(md, encoding="utf-8")
    print(f"CSV written: {csv_path}")
    print(f"Summary written: {summary}")
    print(f"Plot written: {p1}")
    print(f"Plot written: {p2}")
    print(f"Plot written: {p3}")

if __name__ == "__main__":
    main()

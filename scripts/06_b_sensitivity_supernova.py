#!/usr/bin/env python3
"""Supernova catalog discrepancy and b-sensitivity audit for STAM Model-A.

This script answers two separate questions:

1. Are Union3/Pantheon and DES discrepant in their z -> Hubble-ordinate relation?
2. Can a modest change in STAM's b parameter bring Union3, Pantheon, and DES closer?

Important:
- The default locked candidate is b=354.95.
- Pantheon/Union3 define the reference offset.
- DES is evaluated as a holdout; it is not allowed to define its own b.
- A fixed-offset diagnostic is also reported, but should not be treated as a valid
  model fit because it can make DES look better by worsening Pantheon/Union3.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None


C_MLY_PER_MPC = 3.261563776
H_STAM = 0.000243635
L_STAM = C_MLY_PER_MPC / H_STAM
B_LOCKED_DEFAULT = 354.95


def d_adj_mly(z: np.ndarray, b: float) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    return L_STAM * z * (1.0 + 0.5 * z) + b * z


def mu_stam(z: np.ndarray, b: float) -> np.ndarray:
    d_mpc = d_adj_mly(z, b) / C_MLY_PER_MPC
    return 5.0 * np.log10(d_mpc) + 25.0


def median_abs_dev(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x)
    return float(np.nanmedian(np.abs(x - med)))


def robust_summary(x: np.ndarray) -> dict[str, float | int]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return {
        "n": int(len(x)),
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "std": float(np.std(x, ddof=1)) if len(x) > 1 else 0.0,
        "rmse": float(np.sqrt(np.mean(x * x))),
        "mad": median_abs_dev(x),
        "q05": float(np.quantile(x, 0.05)),
        "q16": float(np.quantile(x, 0.16)),
        "q84": float(np.quantile(x, 0.84)),
        "q95": float(np.quantile(x, 0.95)),
    }


def load_catalogs(union3_path: Path, pantheon_path: Path, des_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    union_raw = pd.read_csv(union3_path)
    pant_raw = pd.read_csv(pantheon_path)
    des_raw = pd.read_csv(des_path)

    union = pd.DataFrame(
        {
            "catalog": "Union3",
            "id": union_raw["bin"].astype(str) if "bin" in union_raw else np.arange(len(union_raw)).astype(str),
            "z": pd.to_numeric(union_raw["z"], errors="coerce"),
            "y": pd.to_numeric(union_raw["mb"], errors="coerce"),
            "err": np.nan,
            "ordinate": "mb",
        }
    )

    pantheon = pd.DataFrame(
        {
            "catalog": "Pantheon",
            "id": pant_raw["CID"].astype(str) if "CID" in pant_raw else np.arange(len(pant_raw)).astype(str),
            "z": pd.to_numeric(pant_raw["zHD"], errors="coerce"),
            "y": pd.to_numeric(pant_raw["MU_SH0ES"], errors="coerce"),
            "err": pd.to_numeric(pant_raw["MU_SH0ES_ERR_DIAG"], errors="coerce") if "MU_SH0ES_ERR_DIAG" in pant_raw else np.nan,
            "ordinate": "MU_SH0ES",
        }
    )

    des = pd.DataFrame(
        {
            "catalog": "DES",
            "id": des_raw["CID"].astype(str) if "CID" in des_raw else np.arange(len(des_raw)).astype(str),
            "z": pd.to_numeric(des_raw["zHD"], errors="coerce"),
            "y": pd.to_numeric(des_raw["MU"], errors="coerce"),
            "err": pd.to_numeric(des_raw["MUERR"], errors="coerce") if "MUERR" in des_raw else np.nan,
            "ordinate": "MU",
        }
    )

    def clean(df: pd.DataFrame) -> pd.DataFrame:
        return df[np.isfinite(df["z"]) & np.isfinite(df["y"]) & (df["z"] > 0)].copy()

    return clean(union), clean(pantheon), clean(des)


def cat_resid(df: pd.DataFrame, b: float, offset: float = 0.0) -> np.ndarray:
    return df["y"].to_numpy(dtype=float) - (mu_stam(df["z"].to_numpy(dtype=float), b) + offset)


def reference_offset(b: float, pantheon: pd.DataFrame, union3: pd.DataFrame) -> float:
    """Catalog-balanced reference intercept from Pantheon and Union3 only."""
    med_p = np.nanmedian(cat_resid(pantheon, b, 0.0))
    med_u = np.nanmedian(cat_resid(union3, b, 0.0))
    return float(np.nanmedian([med_p, med_u]))


def balanced_rmse(resids: dict[str, np.ndarray], names: tuple[str, ...]) -> float:
    return float(np.mean([np.sqrt(np.mean(resids[name] ** 2)) for name in names]))


def balanced_median_abs(resids: dict[str, np.ndarray], names: tuple[str, ...]) -> float:
    return float(np.mean([np.median(np.abs(resids[name])) for name in names]))


def make_binned_empirical(union3: pd.DataFrame, pantheon: pd.DataFrame, des: pd.DataFrame, width: float, zmin: float, zmax: float) -> pd.DataFrame:
    edges = np.arange(zmin, zmax + width * 0.5, width)
    rows = []

    for lo, hi in zip(edges[:-1], edges[1:]):
        row: dict[str, float | int] = {
            "bin_lo": float(lo),
            "bin_hi": float(hi),
            "z_center": float((lo + hi) / 2.0),
        }

        for name, df in [("Pantheon", pantheon), ("Union3", union3), ("DES", des)]:
            d = df[(df["z"] >= lo) & (df["z"] < hi)]
            row[f"{name}_n"] = int(len(d))
            row[f"{name}_median"] = float(np.median(d["y"])) if len(d) else np.nan
            row[f"{name}_mean"] = float(np.mean(d["y"])) if len(d) else np.nan

        ref_vals = []
        for name in ("Pantheon", "Union3"):
            val = row[f"{name}_median"]
            if np.isfinite(val):
                ref_vals.append(val)

        row["Ref_balanced_median"] = float(np.median(ref_vals)) if ref_vals else np.nan
        row["DES_minus_Ref"] = (
            row["DES_median"] - row["Ref_balanced_median"]
            if np.isfinite(row["DES_median"]) and np.isfinite(row["Ref_balanced_median"])
            else np.nan
        )
        row["Union3_minus_Pantheon"] = (
            row["Union3_median"] - row["Pantheon_median"]
            if np.isfinite(row["Union3_median"]) and np.isfinite(row["Pantheon_median"])
            else np.nan
        )
        rows.append(row)

    return pd.DataFrame(rows)


def interpolate_reference_curve(binned: pd.DataFrame, z: np.ndarray) -> np.ndarray:
    ref = binned.dropna(subset=["Ref_balanced_median"]).copy()
    x = ref["z_center"].to_numpy(dtype=float)
    y = ref["Ref_balanced_median"].to_numpy(dtype=float)
    order = np.argsort(x)
    return np.interp(z, x[order], y[order], left=np.nan, right=np.nan)


def write_plots(out: Path, fine: pd.DataFrame, candidate: pd.DataFrame, des_eval: pd.DataFrame, ranked: pd.DataFrame, b_locked: float) -> None:
    if plt is None:
        return

    fig = plt.figure(figsize=(8, 5))
    plt.scatter(des_eval["z"], des_eval["des_minus_emp_ref"], s=8, alpha=0.45, label="DES - empirical Pantheon/Union3 ref")
    plt.scatter(des_eval["z"], des_eval["des_minus_stam_locked_b"], s=8, alpha=0.45, label=f"DES - STAM b={b_locked:g}")
    plt.axhline(0.0, linewidth=1)
    plt.xlabel("zHD")
    plt.ylabel("residual [mag]")
    plt.title("DES empirical discrepancy vs STAM locked-b residual")
    plt.legend()
    fig.savefig(out / "des_empirical_discrepancy_vs_stam_residual.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig = plt.figure(figsize=(8, 5))
    plt.plot(fine["b"], fine["balanced_rmse_all3"], label="balanced RMSE all 3")
    plt.plot(fine["b"], fine["des_rmse"], label="DES RMSE")
    plt.axvline(b_locked, linewidth=1, linestyle="--", label=f"b={b_locked:g}")
    plt.xlabel("b")
    plt.ylabel("RMSE [mag]")
    plt.title("b sensitivity with Pantheon/Union3 offset recalibrated")
    plt.legend()
    fig.savefig(out / "b_sensitivity_rmse_reference_offset.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig = plt.figure(figsize=(8, 5))
    plt.plot(fine["b"], fine["des_median"])
    plt.axhline(0.0, linewidth=1)
    plt.axvline(b_locked, linewidth=1, linestyle="--")
    plt.xlabel("b")
    plt.ylabel("DES median residual [mag]")
    plt.title("Changing b does not remove DES offset under reference calibration")
    fig.savefig(out / "des_median_residual_vs_b_reference_offset.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig = plt.figure(figsize=(8, 5))
    for cat in ["Union3", "Pantheon", "DES"]:
        plt.plot(candidate["b"], candidate[f"{cat}_median"], marker="o", label=cat)
    plt.axhline(0.0, linewidth=1)
    plt.axvline(b_locked, linewidth=1, linestyle="--")
    plt.xlabel("b")
    plt.ylabel("median residual [mag]")
    plt.title("Catalog median residuals at candidate b values")
    plt.legend()
    fig.savefig(out / "candidate_b_catalog_median_residuals.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    ranked = ranked.copy()
    ranked["cum_power_frac"] = np.cumsum(ranked["resid_centered_by_des_median"] ** 2) / np.sum(ranked["resid_centered_by_des_median"] ** 2)
    fig = plt.figure(figsize=(8, 5))
    plt.plot(np.arange(1, len(ranked) + 1), ranked["cum_power_frac"])
    plt.axvline(math.ceil(0.05 * len(ranked)), linewidth=1, linestyle="--", label="top 5%")
    plt.xlabel("number of largest centered DES residuals")
    plt.ylabel("cumulative residual-power fraction")
    plt.title("DES outlier tail after removing DES median offset")
    plt.legend()
    fig.savefig(out / "des_outlier_cumulative_power.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--union3", type=Path, default=Path("data/union3_bins.csv"))
    parser.add_argument("--pantheon", type=Path, default=Path("data/pantheon.csv"))
    parser.add_argument("--des", type=Path, default=Path("data/des.csv"))
    parser.add_argument("--b", type=float, default=B_LOCKED_DEFAULT)
    parser.add_argument("--output-dir", type=Path, default=Path("results/b_sensitivity_supernova"))
    parser.add_argument("--grid-min", type=float, default=0.0)
    parser.add_argument("--grid-max", type=float, default=2000.0)
    parser.add_argument("--grid-step", type=float, default=0.1)
    parser.add_argument("--bin-width", type=float, default=0.05)
    parser.add_argument("--clean-output", action="store_true")
    args = parser.parse_args()

    if args.clean_output and args.output_dir.exists():
        shutil.rmtree(args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    union3, pantheon, des = load_catalogs(args.union3, args.pantheon, args.des)

    zmin = max(float(union3["z"].min()), float(pantheon["z"].min()), float(des["z"].min()))
    zmax = min(float(union3["z"].max()), float(pantheon["z"].max()), float(des["z"].max()))

    union_c = union3[(union3["z"] >= zmin) & (union3["z"] <= zmax)].copy()
    pant_c = pantheon[(pantheon["z"] >= zmin) & (pantheon["z"] <= zmax)].copy()
    des_c = des[(des["z"] >= zmin) & (des["z"] <= zmax)].copy()

    offset_locked = reference_offset(args.b, pant_c, union_c)
    res_locked = {
        "Union3": cat_resid(union_c, args.b, offset_locked),
        "Pantheon": cat_resid(pant_c, args.b, offset_locked),
        "DES": cat_resid(des_c, args.b, offset_locked),
    }

    all_resid = []
    for name, df in [("Union3", union_c), ("Pantheon", pant_c), ("DES", des_c)]:
        d = df.copy()
        d["model_mu_locked_b_ref_offset"] = mu_stam(d["z"].to_numpy(dtype=float), args.b) + offset_locked
        d["residual_locked_b_ref_offset"] = d["y"] - d["model_mu_locked_b_ref_offset"]
        all_resid.append(d)
    all_resid_df = pd.concat(all_resid, ignore_index=True)
    all_resid_df.to_csv(args.output_dir / "catalog_residuals_locked_b.csv", index=False)

    binned = make_binned_empirical(
        union_c,
        pant_c,
        des_c,
        width=args.bin_width,
        zmin=max(0.025, zmin - args.bin_width / 2.0),
        zmax=zmax + args.bin_width,
    )
    binned.to_csv(args.output_dir / "binned_empirical_discrepancy.csv", index=False)

    des_eval = des_c.copy()
    des_eval["emp_ref_mu"] = interpolate_reference_curve(binned, des_eval["z"].to_numpy(dtype=float))
    des_eval = des_eval[np.isfinite(des_eval["emp_ref_mu"])].copy()
    des_eval["des_minus_emp_ref"] = des_eval["y"] - des_eval["emp_ref_mu"]
    des_eval["stam_mu_locked_b"] = mu_stam(des_eval["z"].to_numpy(dtype=float), args.b) + offset_locked
    des_eval["des_minus_stam_locked_b"] = des_eval["y"] - des_eval["stam_mu_locked_b"]
    des_eval["stam_resid_minus_empirical_discrepancy"] = des_eval["des_minus_stam_locked_b"] - des_eval["des_minus_emp_ref"]
    des_eval.to_csv(args.output_dir / "des_pointwise_empirical_vs_stam_residuals.csv", index=False)

    b_values = np.arange(args.grid_min, args.grid_max + args.grid_step / 2.0, args.grid_step)
    grid_rows = []
    fixed_rows = []

    for b in b_values:
        off = reference_offset(float(b), pant_c, union_c)
        r_u = cat_resid(union_c, float(b), off)
        r_p = cat_resid(pant_c, float(b), off)
        r_d = cat_resid(des_c, float(b), off)
        resids = {"Union3": r_u, "Pantheon": r_p, "DES": r_d}
        grid_rows.append(
            {
                "b": float(b),
                "offset_reference": off,
                "balanced_rmse_all3": balanced_rmse(resids, ("Union3", "Pantheon", "DES")),
                "balanced_medabs_all3": balanced_median_abs(resids, ("Union3", "Pantheon", "DES")),
                "des_rmse": float(np.sqrt(np.mean(r_d ** 2))),
                "des_median": float(np.median(r_d)),
                "des_mean": float(np.mean(r_d)),
                "pantheon_median": float(np.median(r_p)),
                "union3_median": float(np.median(r_u)),
                "max_abs_catalog_median": float(max(abs(np.median(r_u)), abs(np.median(r_p)), abs(np.median(r_d)))),
            }
        )

        # Diagnostic only: hold the locked-b offset fixed while varying b.
        r_u_f = cat_resid(union_c, float(b), offset_locked)
        r_p_f = cat_resid(pant_c, float(b), offset_locked)
        r_d_f = cat_resid(des_c, float(b), offset_locked)
        fixed_rows.append(
            {
                "b": float(b),
                "fixed_offset": offset_locked,
                "balanced_rmse": float(
                    np.mean([
                        np.sqrt(np.mean(r_u_f ** 2)),
                        np.sqrt(np.mean(r_p_f ** 2)),
                        np.sqrt(np.mean(r_d_f ** 2)),
                    ])
                ),
                "des_rmse": float(np.sqrt(np.mean(r_d_f ** 2))),
                "des_median": float(np.median(r_d_f)),
                "pantheon_median": float(np.median(r_p_f)),
                "union3_median": float(np.median(r_u_f)),
            }
        )

    fine = pd.DataFrame(grid_rows)
    fixed = pd.DataFrame(fixed_rows)
    fine.to_csv(args.output_dir / "b_sensitivity_reference_offset_fine_grid.csv", index=False)
    fixed.to_csv(args.output_dir / "b_sensitivity_fixed_offset_diagnostic_fine_grid.csv", index=False)

    candidate_bs = [0.0, 250.0, 300.0, args.b, 400.0, 461.3626922, 500.0, 650.0, 800.0, 1000.0, 1335.412792]
    cand_rows = []
    for b in candidate_bs:
        off = reference_offset(float(b), pant_c, union_c)
        resids = {
            "Union3": cat_resid(union_c, float(b), off),
            "Pantheon": cat_resid(pant_c, float(b), off),
            "DES": cat_resid(des_c, float(b), off),
        }
        row = {"b": float(b), "reference_offset": off}
        for name in ["Union3", "Pantheon", "DES"]:
            row[f"{name}_median"] = float(np.median(resids[name]))
            row[f"{name}_mean"] = float(np.mean(resids[name]))
            row[f"{name}_rmse"] = float(np.sqrt(np.mean(resids[name] ** 2)))
            row[f"{name}_median_abs"] = float(np.median(np.abs(resids[name])))
        row["balanced_rmse_all3"] = balanced_rmse(resids, ("Union3", "Pantheon", "DES"))
        row["balanced_medabs_all3"] = balanced_median_abs(resids, ("Union3", "Pantheon", "DES"))
        cand_rows.append(row)
    candidate = pd.DataFrame(cand_rows)
    candidate.to_csv(args.output_dir / "candidate_b_table_reference_offset.csv", index=False)

    des_r = des_c.copy()
    des_r["resid_locked_b"] = res_locked["DES"]
    des_median = float(np.median(des_r["resid_locked_b"]))
    des_r["resid_centered_by_des_median"] = des_r["resid_locked_b"] - des_median
    des_r["abs_centered_resid"] = np.abs(des_r["resid_centered_by_des_median"])
    des_r["abs_resid_locked_b"] = np.abs(des_r["resid_locked_b"])
    ranked = des_r.sort_values("abs_centered_resid", ascending=False).reset_index(drop=True)
    ranked_abs = des_r.sort_values("abs_resid_locked_b", ascending=False).reset_index(drop=True)
    ranked.head(200).to_csv(args.output_dir / "des_top200_centered_outliers_locked_b.csv", index=False)
    ranked_abs.head(200).to_csv(args.output_dir / "des_top200_absolute_residuals_locked_b.csv", index=False)

    n_des = len(des_r)
    power_centered = float(np.sum(des_r["resid_centered_by_des_median"] ** 2))
    power_abs = float(np.sum(des_r["resid_locked_b"] ** 2))
    top5pct = math.ceil(0.05 * n_des)

    valid_des_bins = binned.dropna(subset=["DES_minus_Ref"])
    valid_up_bins = binned.dropna(subset=["Union3_minus_Pantheon"])

    summary = {
        "inputs": {
            "union3": {"file": str(args.union3), "ordinate_used": "mb", "n_total": int(len(union3)), "n_common_overlap": int(len(union_c))},
            "pantheon": {"file": str(args.pantheon), "ordinate_used": "MU_SH0ES", "z_used": "zHD", "n_total": int(len(pantheon)), "n_common_overlap": int(len(pant_c))},
            "des": {"file": str(args.des), "ordinate_used": "MU", "z_used": "zHD", "n_total": int(len(des)), "n_common_overlap": int(len(des_c))},
            "common_overlap_z_min": zmin,
            "common_overlap_z_max": zmax,
        },
        "locked_model": {
            "b_locked": args.b,
            "reference_offset_mag_balanced_median_Pantheon_Union3": offset_locked,
            "L_Mly": L_STAM,
            "C_Mly_per_Mpc": C_MLY_PER_MPC,
        },
        "catalog_residuals_at_locked_b_reference_offset": {
            "Union3": robust_summary(res_locked["Union3"]),
            "Pantheon": robust_summary(res_locked["Pantheon"]),
            "DES": robust_summary(res_locked["DES"]),
        },
        "empirical_discrepancy": {
            "binned_DES_minus_balanced_Pantheon_Union3_reference": robust_summary(valid_des_bins["DES_minus_Ref"].to_numpy(dtype=float)),
            "binned_DES_minus_balanced_reference_z_ge_0p10": robust_summary(valid_des_bins[valid_des_bins["z_center"] >= 0.10]["DES_minus_Ref"].to_numpy(dtype=float)),
            "binned_Union3_minus_Pantheon": robust_summary(valid_up_bins["Union3_minus_Pantheon"].to_numpy(dtype=float)),
            "binned_Union3_minus_Pantheon_z_ge_0p10": robust_summary(valid_up_bins[valid_up_bins["z_center"] >= 0.10]["Union3_minus_Pantheon"].to_numpy(dtype=float)),
            "pointwise_DES_minus_empirical_reference": robust_summary(des_eval["des_minus_emp_ref"].to_numpy(dtype=float)),
            "pointwise_DES_minus_STAM_locked_b": robust_summary(des_eval["des_minus_stam_locked_b"].to_numpy(dtype=float)),
            "pointwise_STAM_residual_minus_empirical_discrepancy": robust_summary(des_eval["stam_resid_minus_empirical_discrepancy"].to_numpy(dtype=float)),
            "correlation_empirical_DES_discrepancy_vs_STAM_residual": float(
                np.corrcoef(des_eval["des_minus_emp_ref"], des_eval["des_minus_stam_locked_b"])[0, 1]
            ),
        },
        "b_sensitivity_reference_recalibrated_each_b": {
            "best_balanced_rmse_all3": fine.loc[fine["balanced_rmse_all3"].idxmin()].to_dict(),
            "best_DES_rmse": fine.loc[fine["des_rmse"].idxmin()].to_dict(),
            "best_balanced_median_abs_all3": fine.loc[fine["balanced_medabs_all3"].idxmin()].to_dict(),
            "b_locked_row": fine.iloc[(fine["b"] - args.b).abs().argmin()].to_dict(),
            "des_median_residual_range_grid": {"min": float(fine["des_median"].min()), "max": float(fine["des_median"].max())},
            "des_rmse_range_grid": {"min": float(fine["des_rmse"].min()), "max": float(fine["des_rmse"].max())},
        },
        "b_sensitivity_fixed_offset_diagnostic_not_recommended_for_fitting": {
            "offset_fixed_from_locked_b": offset_locked,
            "b_that_zeroes_DES_median_with_fixed_offset": fixed.loc[fixed["des_median"].abs().idxmin()].to_dict(),
            "best_balanced_rmse_with_fixed_offset": fixed.loc[fixed["balanced_rmse"].idxmin()].to_dict(),
        },
        "des_outlier_tail": {
            "top5pct_n": top5pct,
            "absolute_residual_power_top5pct_fraction": float(np.sum(ranked_abs.head(top5pct)["resid_locked_b"] ** 2) / power_abs),
            "centered_residual_power_top5pct_fraction": float(np.sum(ranked.head(top5pct)["resid_centered_by_des_median"] ** 2) / power_centered),
        },
        "interpretation": {
            "short": "Changing b modestly does not bring DES into alignment with Pantheon/Union3 under a fair reference-calibrated test. The DES offset is catalog/intercept/outlier-structure dominated, not solved by a slight b adjustment.",
            "fixed_offset_warning": "Large b can reduce DES if the b=354.95 offset is artificially held fixed, but this worsens Pantheon/Union3 and should not be used as a valid model fit.",
        },
    }

    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    write_plots(args.output_dir, fine, candidate, des_eval, ranked, args.b)

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote outputs to {args.output_dir}")


if __name__ == "__main__":
    main()

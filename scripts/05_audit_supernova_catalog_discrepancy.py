#!/usr/bin/env python3
"""Audit Union3/Pantheon vs DES z-to-magnitude discrepancies for STAM Model-A.

This script is intentionally not a DES fit. It locks STAM Model-A's b value
(default 354.95), learns only a global STAM magnitude offset from Pantheon and
Union3, and then compares DES against:

1. an empirical Pantheon reference curve;
2. the locked-b STAM Model-A curve;
3. Union3-centered binned catalog medians.

Expected input files by default:

    data/union3_bins.csv     columns: bin,z,mb
    data/pantheon.csv        Pantheon-style columns: zHD,MU_SH0ES,mB
    data/des.csv             DES SNANA-style columns: VARNAMES row, SN rows, zHD,MU,mB

Outputs are written under results/supernova_discrepancy by default.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import zipfile
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

C_MLY_PER_MPC = 3.261563776
H_STAM = 0.000243635
L_STAM = C_MLY_PER_MPC / H_STAM
B_LOCKED_DEFAULT = 354.95


def ffloat(value: object) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except Exception:
        return math.nan


def stam_distance_mly(z: float, b: float, L: float = L_STAM) -> float:
    return L * z * (1.0 + 0.5 * z) + b * z


def stam_mu(z: float, b: float, L: float = L_STAM) -> float:
    d_mpc = stam_distance_mly(z, b, L) / C_MLY_PER_MPC
    return 5.0 * math.log10(d_mpc) + 25.0


def infer_b_from_mu(z: float, mu: float, L: float = L_STAM) -> float:
    d_mpc = 10.0 ** ((mu - 25.0) / 5.0)
    d_mly = d_mpc * C_MLY_PER_MPC
    return (d_mly - L * z * (1.0 + 0.5 * z)) / z


def read_union3(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            z = ffloat(row.get("z"))
            y = ffloat(row.get("mb"))
            if math.isfinite(z) and math.isfinite(y) and z > 0:
                rows.append(
                    {
                        "catalog": "Union3",
                        "id": row.get("bin", ""),
                        "z": z,
                        "y": y,
                        "y_col": "mb",
                        "raw": row,
                    }
                )
    return rows


def read_pantheon(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            z = ffloat(row.get("zHD") or row.get("zCMB") or row.get("z"))
            y = ffloat(row.get("MU_SH0ES"))
            err = ffloat(row.get("MU_SH0ES_ERR_DIAG"))
            if math.isfinite(z) and math.isfinite(y) and z > 0 and y > 0:
                rows.append(
                    {
                        "catalog": "Pantheon",
                        "id": row.get("CID", ""),
                        "z": z,
                        "y": y,
                        "err": err if math.isfinite(err) and err > 0 else math.nan,
                        "raw_mB": ffloat(row.get("mB")),
                        "y_col": "MU_SH0ES",
                        "raw": row,
                    }
                )
    return rows


def read_des(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader)
        if header and header[0].strip().upper().startswith("VARNAMES"):
            header = header[1:]
        for row in reader:
            if not row:
                continue
            if row[0].strip().upper().startswith("SN"):
                row = row[1:]
            if len(row) < len(header):
                continue
            raw = dict(zip(header, row))
            z = ffloat(raw.get("zHD") or raw.get("zCMB") or raw.get("zHEL"))
            y = ffloat(raw.get("MU"))
            if math.isfinite(z) and math.isfinite(y) and z > 0 and y > 0:
                rows.append(
                    {
                        "catalog": "DES",
                        "id": raw.get("CID", ""),
                        "z": z,
                        "y": y,
                        "err": ffloat(raw.get("MUERR")),
                        "raw_mB": ffloat(raw.get("mB")),
                        "mures": ffloat(raw.get("MURES")),
                        "y_col": "MU",
                        "raw": raw,
                    }
                )
    return rows


def catalog_summary(name: str, rows: list[dict[str, object]]) -> dict[str, object]:
    z = np.array([float(r["z"]) for r in rows], dtype=float)
    y = np.array([float(r["y"]) for r in rows], dtype=float)
    return {
        "catalog": name,
        "n": int(len(rows)),
        "z_min": float(np.min(z)),
        "z_median": float(np.median(z)),
        "z_max": float(np.max(z)),
        "y_col": rows[0].get("y_col", "") if rows else "",
        "y_min": float(np.min(y)),
        "y_median": float(np.median(y)),
        "y_max": float(np.max(y)),
    }


def robust_empirical_fit(rows: list[dict[str, object]], degree: int = 3) -> dict[str, object]:
    """Fit y = 5log10(z) + polynomial(z), with simple MAD clipping."""
    z = np.array([float(r["z"]) for r in rows], dtype=float)
    y = np.array([float(r["y"]) for r in rows], dtype=float)
    base = 5.0 * np.log10(z)
    x = np.column_stack([np.ones_like(z)] + [z**k for k in range(1, degree + 1)])
    mask = np.ones(len(z), dtype=bool)
    coef = np.zeros(degree + 1)
    for _ in range(6):
        coef = np.linalg.lstsq(x[mask], (y - base)[mask], rcond=None)[0]
        residual = y - (base + x @ coef)
        med = np.median(residual[mask])
        mad = 1.4826 * np.median(np.abs(residual[mask] - med))
        if not np.isfinite(mad) or mad <= 0:
            break
        new_mask = np.abs(residual - med) < 3.5 * mad
        if int(np.sum(new_mask)) == int(np.sum(mask)):
            break
        mask = new_mask
    residual = y - (base + x @ coef)
    return {
        "degree": degree,
        "coef": coef.tolist(),
        "n_input": int(len(z)),
        "n_used_after_clip": int(np.sum(mask)),
        "fit_residual_median_mag": float(np.median(residual[mask])),
        "fit_residual_std_mag": float(np.std(residual[mask], ddof=1)),
    }


def empirical_predict(z_values: Iterable[float], fit: dict[str, object]) -> np.ndarray:
    z = np.array(list(z_values), dtype=float)
    coef = np.array(fit["coef"], dtype=float)
    degree = int(fit["degree"])
    x = np.column_stack([np.ones_like(z)] + [z**k for k in range(1, degree + 1)])
    return 5.0 * np.log10(z) + x @ coef


def stats(values: Iterable[float]) -> dict[str, object]:
    arr = np.array(list(values), dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return {"n": 0, "median": math.nan, "mean": math.nan, "std": math.nan, "rmse": math.nan}
    return {
        "n": int(len(arr)),
        "median": float(np.median(arr)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0,
        "rmse": float(np.sqrt(np.mean(arr**2))),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    keys = list(rows[0].keys())
    for row in rows:
        for key in row.keys():
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def build_union_centered_bins(union: list[dict[str, object]], zmax: float) -> list[tuple[float, float, float]]:
    centers = np.array([float(r["z"]) for r in union if float(r["z"]) <= zmax], dtype=float)
    bins: list[tuple[float, float, float]] = []
    for i, center in enumerate(centers):
        if i == 0:
            left = max(0.0, center - (centers[i + 1] - center) / 2.0)
        else:
            left = (centers[i - 1] + center) / 2.0
        if i == len(centers) - 1:
            right = min(zmax, center + (center - centers[i - 1]) / 2.0)
        else:
            right = (center + centers[i + 1]) / 2.0
        bins.append((float(center), float(left), float(right)))
    return bins


def median_or_nan(values: list[float]) -> float:
    return float(np.median(values)) if values else math.nan


def power_fractions(values: np.ndarray) -> dict[str, float]:
    power = values * values
    total = float(np.sum(power))
    order = np.argsort(power)[::-1]
    out: dict[str, float] = {}
    for k in (1, 5, 10, 20, 50, 100, 200):
        kk = min(k, len(values))
        out[f"top_{k}_power_fraction"] = float(np.sum(power[order[:kk]]) / total) if total > 0 else math.nan
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--union3", type=Path, default=Path("data/union3_bins.csv"))
    parser.add_argument("--pantheon", type=Path, default=Path("data/pantheon.csv"))
    parser.add_argument("--des", type=Path, default=Path("data/des.csv"))
    parser.add_argument("--b", type=float, default=B_LOCKED_DEFAULT)
    parser.add_argument("--degree", type=int, default=3)
    parser.add_argument("--output-dir", type=Path, default=Path("results/supernova_discrepancy"))
    args = parser.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    union = read_union3(args.union3)
    pantheon = read_pantheon(args.pantheon)
    des = read_des(args.des)

    zmin = max(min(float(r["z"]) for r in des), 0.025)
    zmax = min(max(float(r["z"]) for r in des), max(float(r["z"]) for r in pantheon), max(float(r["z"]) for r in union))

    pantheon_overlap = [r for r in pantheon if zmin <= float(r["z"]) <= zmax]
    union_overlap = [r for r in union if zmin <= float(r["z"]) <= zmax]
    des_overlap = [r for r in des if zmin <= float(r["z"]) <= zmax]

    # Empirical reference curve from Pantheon only, used to test Union3 and DES consistency.
    pantheon_fit = robust_empirical_fit(pantheon_overlap, degree=args.degree)
    union_empirical_resid = np.array(
        [float(r["y"]) - empirical_predict([float(r["z"])], pantheon_fit)[0] for r in union_overlap],
        dtype=float,
    )
    des_empirical_resid = np.array(
        [float(r["y"]) - empirical_predict([float(r["z"])], pantheon_fit)[0] for r in des_overlap],
        dtype=float,
    )

    # STAM locked-b offset learned from Pantheon and Union3 only, equal catalog weight.
    pantheon_stam_resid = np.array([float(r["y"]) - stam_mu(float(r["z"]), args.b) for r in pantheon_overlap], dtype=float)
    union_stam_resid = np.array([float(r["y"]) - stam_mu(float(r["z"]), args.b) for r in union_overlap], dtype=float)
    pantheon_offset = float(np.median(pantheon_stam_resid))
    union_offset = float(np.median(union_stam_resid))
    equal_catalog_stam_offset = float((pantheon_offset + union_offset) / 2.0)

    des_stam_resid = np.array(
        [float(r["y"]) - (stam_mu(float(r["z"]), args.b) + equal_catalog_stam_offset) for r in des_overlap],
        dtype=float,
    )
    residual_difference = des_stam_resid - des_empirical_resid

    # Binned checks centered on Union3 bins.
    bin_rows: list[dict[str, object]] = []
    for center, left, right in build_union_centered_bins(union, zmax):
        u = [r for r in union if abs(float(r["z"]) - center) < 1e-12]
        p = [r for r in pantheon if left <= float(r["z"]) < right]
        d = [r for r in des if left <= float(r["z"]) < right]
        y_u = float(u[0]["y"]) if u else math.nan
        y_p = median_or_nan([float(r["y"]) for r in p])
        y_d = median_or_nan([float(r["y"]) for r in d])
        ref_vals = [v for v in (y_p, y_u) if math.isfinite(v)]
        y_ref = median_or_nan(ref_vals)
        d_stam_bin = [float(r["y"]) - (stam_mu(float(r["z"]), args.b) + equal_catalog_stam_offset) for r in d]
        bin_rows.append(
            {
                "bin_center_z": center,
                "z_left": left,
                "z_right": right,
                "n_pantheon": len(p),
                "n_des": len(d),
                "union3_mb": y_u,
                "pantheon_median_mu": y_p,
                "des_median_mu": y_d,
                "reference_equal_catalog_mu": y_ref,
                "union3_minus_pantheon_med": y_u - y_p if math.isfinite(y_u) and math.isfinite(y_p) else math.nan,
                "des_minus_pantheon_med": y_d - y_p if math.isfinite(y_d) and math.isfinite(y_p) else math.nan,
                "des_minus_union3": y_d - y_u if math.isfinite(y_d) and math.isfinite(y_u) else math.nan,
                "des_minus_reference_equal": y_d - y_ref if math.isfinite(y_d) and math.isfinite(y_ref) else math.nan,
                "des_stam_locked_b_residual_median": median_or_nan(d_stam_bin),
            }
        )

    # DES pointwise table and outliers.
    des_point_rows: list[dict[str, object]] = []
    for r, empirical_res, stam_res in zip(des_overlap, des_empirical_resid, des_stam_resid):
        z = float(r["z"])
        y = float(r["y"])
        b_inf = infer_b_from_mu(z, y - equal_catalog_stam_offset)
        raw = r.get("raw", {}) if isinstance(r.get("raw"), dict) else {}
        des_point_rows.append(
            {
                "CID": r.get("id", ""),
                "z": z,
                "MU": y,
                "MUERR": r.get("err", math.nan),
                "raw_mB": r.get("raw_mB", math.nan),
                "DES_minus_Pantheon_empirical_curve_mag": empirical_res,
                "DES_minus_STAM_locked_b_mag": stam_res,
                "STAM_residual_minus_empirical_discrepancy_mag": stam_res - empirical_res,
                "b_inferred_after_reference_offset": b_inf,
                "DES_SNANA_MURES": r.get("mures", math.nan),
                "HOST_LOGMASS": raw.get("HOST_LOGMASS", ""),
                "x1": raw.get("x1", ""),
                "c": raw.get("c", ""),
                "FITPROB": raw.get("FITPROB", ""),
                "CUTMASK": raw.get("CUTMASK", ""),
                "IDSAMPLE": raw.get("IDSAMPLE", ""),
                "IZBIN": raw.get("IZBIN", ""),
            }
        )

    des_median_stam_offset = float(np.median(des_stam_resid))
    centered_stam_resid = des_stam_resid - des_median_stam_offset
    for row, centered in zip(des_point_rows, centered_stam_resid):
        row["DES_minus_STAM_after_DES_median_offset_removed_mag"] = centered
        row["abs_DES_minus_STAM_after_DES_median_offset_removed_mag"] = abs(centered)
        row["abs_DES_minus_STAM_locked_b_mag"] = abs(float(row["DES_minus_STAM_locked_b_mag"]))

    top_centered = sorted(des_point_rows, key=lambda r: float(r["abs_DES_minus_STAM_after_DES_median_offset_removed_mag"]), reverse=True)[:50]
    top_abs = sorted(des_point_rows, key=lambda r: float(r["abs_DES_minus_STAM_locked_b_mag"]), reverse=True)[:50]

    summary = {
        "locked_b": args.b,
        "L": L_STAM,
        "ordinate_used": {
            "Union3": "mb column, treated as distance-modulus-like Hubble-diagram ordinate",
            "Pantheon": "MU_SH0ES",
            "DES": "MU",
        },
        "catalog_summary": [catalog_summary("Union3", union), catalog_summary("Pantheon", pantheon), catalog_summary("DES", des)],
        "common_overlap_range": {"z_min": zmin, "z_max": zmax},
        "pantheon_empirical_fit": pantheon_fit,
        "union3_vs_pantheon_empirical_curve": stats(union_empirical_resid),
        "des_vs_pantheon_empirical_curve": stats(des_empirical_resid),
        "stam_reference_offset_mag": {
            "pantheon_median": pantheon_offset,
            "union3_median": union_offset,
            "equal_catalog_average_used": equal_catalog_stam_offset,
        },
        "des_vs_STAM_locked_b": stats(des_stam_resid),
        "STAM_DES_residual_minus_DES_empirical_discrepancy": stats(residual_difference),
        "correlation_empirical_discrepancy_vs_STAM_residual": float(np.corrcoef(des_empirical_resid, des_stam_resid)[0, 1]),
        "des_b_inferred_after_reference_offset": stats([float(r["b_inferred_after_reference_offset"]) for r in des_point_rows]),
        "outlier_power_absolute_STAM_residual": power_fractions(des_stam_resid),
        "outlier_power_after_DES_median_offset_removed": power_fractions(centered_stam_resid),
    }

    write_csv(out / "catalog_summary.csv", summary["catalog_summary"])  # type: ignore[arg-type]
    write_csv(out / "union_centered_binned_discrepancy.csv", bin_rows)
    write_csv(out / "des_pointwise_empirical_and_stam_residuals.csv", des_point_rows)
    write_csv(out / "des_top50_outliers_after_des_median_offset.csv", top_centered)
    write_csv(out / "des_top50_absolute_locked_b_residuals.csv", top_abs)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    # Plots.
    z_des = np.array([float(r["z"]) for r in des_overlap], dtype=float)
    plt.figure(figsize=(8, 5))
    plt.axhline(0.0, linewidth=1)
    plt.scatter(z_des, des_empirical_resid, s=10, alpha=0.6, label="DES - Pantheon empirical curve")
    plt.scatter(z_des, des_stam_resid, s=10, alpha=0.45, label="DES - STAM locked b")
    plt.xlabel("DES z")
    plt.ylabel("residual [mag]")
    plt.title("DES empirical discrepancy vs STAM locked-b residual")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out / "des_empirical_vs_stam_pointwise.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.axhline(0.0, linewidth=1)
    plt.scatter(z_des, residual_difference, s=10, alpha=0.6)
    plt.xlabel("DES z")
    plt.ylabel("STAM residual - empirical discrepancy [mag]")
    plt.title("How closely locked STAM residual matches DES catalog discrepancy")
    plt.tight_layout()
    plt.savefig(out / "stam_residual_minus_empirical_discrepancy.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.axhline(0.0, linewidth=1)
    plt.scatter(z_des, centered_stam_resid, s=10, alpha=0.6)
    plt.xlabel("DES z")
    plt.ylabel("residual after DES median offset removed [mag]")
    plt.title("DES residual shape/outlier tail after removing global DES offset")
    plt.tight_layout()
    plt.savefig(out / "des_centered_stam_residual_outlier_tail.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8, 5))
    bin_z = [float(r["bin_center_z"]) for r in bin_rows]
    plt.axhline(0.0, linewidth=1)
    plt.plot(bin_z, [float(r["union3_minus_pantheon_med"]) for r in bin_rows], marker="o", label="Union3 - Pantheon bin median")
    plt.plot(bin_z, [float(r["des_minus_pantheon_med"]) for r in bin_rows], marker="o", label="DES - Pantheon bin median")
    plt.xlabel("Union3-centered z bin")
    plt.ylabel("median difference [mag]")
    plt.title("Binned catalog discrepancies")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out / "union_centered_binned_catalog_discrepancies.png", dpi=180)
    plt.close()

    zip_path = out / "supernova_discrepancy_outputs.zip"
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(out.rglob("*")):
            if path.is_file() and path.name != zip_path.name:
                zf.write(path, path.relative_to(out))

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote outputs to {out}")
    print(f"Wrote zip to {zip_path}")


if __name__ == "__main__":
    main()

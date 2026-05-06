#!/usr/bin/env python3
"""STAM Model-A supernova light-curve time-dilation proxy diagnostic.

This script uses available SALT2 x1 columns from Pantheon and DES.

Important:
  x1 is a fitted light-curve shape / standardization parameter.
  It is not raw observed-frame light-curve duration.

Therefore this script is a proxy/catalog diagnostic, not a final test of
observed supernova time dilation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

C_MLY_PER_MPC = 3.261563776
H_STAM = 0.000243635
L_MPC = 1.0 / H_STAM


def D_adj0_mpc(z):
    z = np.asarray(z, dtype=float)
    return L_MPC * z * (1.0 + 0.5 * z)


def mu_from_mpc(d_mpc):
    return 5.0 * np.log10(np.asarray(d_mpc, dtype=float)) + 25.0


def mu_model_adj0(z):
    return mu_from_mpc(D_adj0_mpc(z))


def load_catalog(name: str, path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if name.lower() == "pantheon":
        out = pd.DataFrame(
            {
                "catalog": "Pantheon",
                "object_id": df["CID"].astype(str),
                "z": pd.to_numeric(df["zHD"], errors="coerce"),
                "mu_obs": pd.to_numeric(df["MU_SH0ES"], errors="coerce"),
                "x1": pd.to_numeric(df["x1"], errors="coerce"),
                "x1err": pd.to_numeric(df["x1ERR"], errors="coerce"),
                "c": pd.to_numeric(df["c"], errors="coerce"),
            }
        )
    elif name.lower() == "des":
        out = pd.DataFrame(
            {
                "catalog": "DES",
                "object_id": df["CID"].astype(str),
                "z": pd.to_numeric(df["zHD"], errors="coerce"),
                "mu_obs": pd.to_numeric(df["MU"], errors="coerce"),
                "x1": pd.to_numeric(df["x1"], errors="coerce"),
                "x1err": pd.to_numeric(df["x1ERR"], errors="coerce"),
                "c": pd.to_numeric(df["c"], errors="coerce"),
            }
        )
    else:
        raise ValueError(name)

    out = out[np.isfinite(out["z"]) & (out["z"] > 0) & np.isfinite(out["x1"])].copy()
    out["mu_stam_adj0"] = mu_model_adj0(out["z"])
    out["resid_mu_adj0"] = out["mu_obs"] - out["mu_stam_adj0"]
    out["log1pz"] = np.log1p(out["z"])
    return out


def corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 3:
        return np.nan, np.nan
    pearson = float(np.corrcoef(x, y)[0, 1])
    spearman = float(pd.Series(x).rank().corr(pd.Series(y).rank()))
    return pearson, spearman


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pantheon", type=Path, default=Path("data/pantheon.csv"))
    parser.add_argument("--des", type=Path, default=Path("data/des.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/sn_lightcurve_time_proxy"))
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    frames = []
    if args.pantheon.exists():
        frames.append(load_catalog("Pantheon", args.pantheon))
    if args.des.exists():
        frames.append(load_catalog("DES", args.des))
    all_df = pd.concat(frames, ignore_index=True)

    common = all_df[(all_df["z"] >= 0.01) & (all_df["z"] <= 1.14418)].copy()

    rows = []
    for cat, g in common.groupby("catalog"):
        pz, sz = corr(g["z"], g["x1"])
        plog, slog = corr(g["log1pz"], g["x1"])
        pres, sres = corr(g["x1"], g["resid_mu_adj0"])
        rows.append(
            {
                "catalog": cat,
                "n": int(len(g)),
                "x1_median": float(g["x1"].median()),
                "x1_mean": float(g["x1"].mean()),
                "x1_std": float(g["x1"].std(ddof=1)),
                "corr_x1_z_pearson": pz,
                "corr_x1_z_spearman": sz,
                "corr_x1_log1pz_pearson": plog,
                "corr_x1_log1pz_spearman": slog,
                "corr_x1_resid_adj0_pearson": pres,
                "corr_x1_resid_adj0_spearman": sres,
            }
        )
    summary_df = pd.DataFrame(rows)

    all_df.to_csv(args.output_dir / "sn_lightcurve_proxy_pointwise.csv", index=False)
    summary_df.to_csv(args.output_dir / "sn_lightcurve_proxy_summary.csv", index=False)

    future_schema = pd.DataFrame(
        [
            {"column": "z", "meaning": "redshift"},
            {"column": "width_observed_days", "meaning": "observed-frame light-curve width or duration"},
            {"column": "width_rest_or_template_days", "meaning": "rest-frame/template width"},
            {"column": "duration_ratio", "meaning": "width_observed_days / width_rest_or_template_days"},
            {"column": "compare_to_standard", "meaning": "duration_ratio / (1+z)"},
            {"column": "compare_to_STAM_candidate", "meaning": "candidate accumulation-rate/stretch law once defined"},
        ]
    )
    future_schema.to_csv(args.output_dir / "future_raw_lightcurve_duration_schema.csv", index=False)

    summary = {
        "test": "STAM Model-A supernova light-curve time-dilation proxy diagnostic",
        "scope": "Uses SALT2 x1 only; not a raw observed-duration test.",
        "interpretation": (
            "A direct STAM time-dilation test requires raw observed-frame light-curve widths. "
            "This script only checks the available catalog stretch-like proxy x1."
        ),
        "summary": rows,
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    print("PASS: SN light-curve proxy diagnostic complete.")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote outputs to {args.output_dir}")


if __name__ == "__main__":
    main()

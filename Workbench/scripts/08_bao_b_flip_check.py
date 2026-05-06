#!/usr/bin/env python3
"""Check whether b explains the BAO low/high mapping flip.

This script is deliberately small. It compares:
  - D_geo/(1+z)
  - D_adj/(1+z) with b=0
  - D_adj/(1+z) with b=354.95
  - D_adj/(1+z) with historical bridge values

It fits only one nuisance scale r_d_eff for each mapping, matching the earlier
BAO geometric split test.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

C_MLY_PER_MPC = 3.261563776
H_STAM = 0.000243635
L_STAM = C_MLY_PER_MPC / H_STAM

B_VALUES = {
    "no_b": 0.0,
    "pantheon_union_bridge": 354.95,
    "original_retained_bridge": 461.3626922,
    "des_style_posthoc_bridge": 1335.412792,
}


def d_geo(z: np.ndarray) -> np.ndarray:
    return L_STAM * z * (1.0 + 0.15 * z)


def d_adj(z: np.ndarray, b: float) -> np.ndarray:
    return L_STAM * z * (1.0 + 0.5 * z) + b * z


def load_bao() -> pd.DataFrame:
    files = [
        Path("data/bao/sdss_dr17_bao_only_dm.csv"),
        Path("data/bao/desi_dr1_bao_dm.csv"),
    ]
    frames = []
    for path in files:
        if not path.exists():
            raise FileNotFoundError(f"Missing BAO input file: {path}")
        frames.append(pd.read_csv(path))
    df = pd.concat(frames, ignore_index=True)
    df = df[df["observable"].str.lower().isin(["dm_over_rd", "d_m_over_r_d", "dm/rd", "dm_over_rs"])]
    df = df.dropna(subset=["z", "value", "error"]).copy()
    df = df[(df["z"] > 0) & (df["error"] > 0)].copy()

    # Match earlier script: inverse-variance combine duplicate rounded redshifts within each survey.
    df["z_key"] = df["z"].round(3)
    rows = []
    for (survey, z_key), g in df.groupby(["survey", "z_key"], sort=True):
        w = 1.0 / np.square(g["error"].to_numpy(float))
        rows.append(
            {
                "survey": survey,
                "sample": "+".join(g["sample"].astype(str)),
                "z": float(np.sum(w * g["z"].to_numpy(float)) / np.sum(w)),
                "value": float(np.sum(w * g["value"].to_numpy(float)) / np.sum(w)),
                "error": float(math.sqrt(1.0 / np.sum(w))),
            }
        )
    return pd.DataFrame(rows).sort_values(["survey", "z"]).reset_index(drop=True)


def fit_one_scale(df: pd.DataFrame, x: np.ndarray) -> dict[str, float]:
    y = df["value"].to_numpy(float)
    sigma = df["error"].to_numpy(float)
    w = 1.0 / sigma**2
    alpha = float(np.sum(w * x * y) / np.sum(w * x * x))
    pred = alpha * x
    residual = y - pred
    chi2 = float(np.sum((residual / sigma) ** 2))
    dof = max(1, len(df) - 1)
    return {
        "n": int(len(df)),
        "chi2_dof": chi2 / dof,
        "r_d_eff_mpc": (1.0 / alpha) / C_MLY_PER_MPC,
        "weighted_rmse": float(np.sqrt(np.sum(w * residual**2) / np.sum(w))),
        "median_residual": float(np.median(residual)),
    }


def run() -> pd.DataFrame:
    data = load_bao()
    scenarios = {
        "combined_all_z": data,
        "combined_z_le_1p6": data[data["z"] <= 1.6].copy(),
    }

    rows = []
    for scenario, df in scenarios.items():
        z = df["z"].to_numpy(float)

        result = fit_one_scale(df, d_geo(z) / (1.0 + z))
        rows.append({"scenario": scenario, "mapping": "geo_over_1pz", "b": None, **result})

        for name, b in B_VALUES.items():
            result = fit_one_scale(df, d_adj(z, b=b) / (1.0 + z))
            rows.append({"scenario": scenario, "mapping": f"adj_over_1pz_{name}", "b": b, **result})

    out = pd.DataFrame(rows)
    out_dir = Path("results/bao_b_flip_check")
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / "bao_b_flip_check_summary.csv", index=False)

    print(out.to_string(index=False))
    print(f"\nWrote {out_dir / 'bao_b_flip_check_summary.csv'}")
    return out


if __name__ == "__main__":
    run()

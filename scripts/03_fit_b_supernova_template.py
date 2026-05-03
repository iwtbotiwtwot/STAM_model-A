#!/usr/bin/env python3
"""Template for fitting STAM_model-A b from a supernova/distance catalog.

This script does not include a dataset. Provide a CSV with either:

    z, D_obs

or:

    z, mu

where D_obs is in the same distance scale as L_STAM, and mu is standard distance modulus.

Usage:

    python scripts/03_fit_b_supernova_template.py data/my_catalog.csv --name my_catalog
"""

from __future__ import annotations

from pathlib import Path
import argparse
import json

import numpy as np
import pandas as pd

from stam_model_a.constants import L_STAM
from stam_model_a.cosmology import D_adj, distance_modulus_to_mly, infer_b_from_distance


def fit_b_least_squares(z: np.ndarray, d_obs: np.ndarray, L: float = L_STAM) -> float:
    """Fit b in D_obs ~= Lz(1+0.5z)+bz by least squares."""
    base = L * z * (1.0 + 0.5 * z)
    numerator = np.sum(z * (d_obs - base))
    denominator = np.sum(z**2)
    if denominator <= 0:
        raise ValueError("cannot fit b with zero redshift leverage")
    return float(numerator / denominator)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("--name", default="catalog")
    parser.add_argument("--output-dir", type=Path, default=Path("results/fits"))
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    if "z" not in df.columns:
        raise ValueError("CSV must contain column 'z'")
    z = df["z"].to_numpy(dtype=float)
    if np.any(z <= 0):
        raise ValueError("all z values must be positive for b fitting")

    if "D_obs" in df.columns:
        d_obs = df["D_obs"].to_numpy(dtype=float)
        distance_source = "D_obs"
    elif "mu" in df.columns:
        d_obs = distance_modulus_to_mly(df["mu"].to_numpy(dtype=float))
        distance_source = "mu_converted_to_Mly"
    else:
        raise ValueError("CSV must contain either 'D_obs' or 'mu'")

    b_fit = fit_b_least_squares(z, d_obs, L=L_STAM)
    pred = D_adj(z, b=b_fit, L=L_STAM)
    residual = d_obs - pred
    b_each = infer_b_from_distance(z, d_obs, L=L_STAM)

    result = {
        "catalog_name": args.name,
        "input_file": str(args.csv),
        "distance_source": distance_source,
        "n": int(len(df)),
        "L": L_STAM,
        "b_fit": b_fit,
        "mean_residual": float(np.mean(residual)),
        "mean_abs_residual": float(np.mean(np.abs(residual))),
        "rmse_residual": float(np.sqrt(np.mean(residual**2))),
        "b_inferred_mean": float(np.mean(b_each)),
        "b_inferred_std": float(np.std(b_each, ddof=1)) if len(b_each) > 1 else 0.0,
        "label": "fit_not_prediction",
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.output_dir / f"{args.name}_b_fit.json"
    out_path.write_text(json.dumps(result, indent=2, sort_keys=True))

    table = pd.DataFrame({
        "z": z,
        "D_obs": d_obs,
        "D_pred": pred,
        "residual": residual,
        "b_inferred": b_each,
    })
    table_path = args.output_dir / f"{args.name}_b_fit_rows.csv"
    table.to_csv(table_path, index=False)

    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"Wrote {out_path}")
    print(f"Wrote {table_path}")


if __name__ == "__main__":
    main()

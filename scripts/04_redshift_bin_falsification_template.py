#!/usr/bin/env python3
"""Template for redshift-bin b stability tests.

Provide a CSV with z and D_obs, or z and mu. The script reports inferred b by redshift bin.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import json

import numpy as np
import pandas as pd

from stam_model_a.constants import L_STAM
from stam_model_a.cosmology import distance_modulus_to_mly, infer_b_from_distance


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("--name", default="catalog")
    parser.add_argument("--bins", type=int, default=8)
    parser.add_argument("--output-dir", type=Path, default=Path("results/fits"))
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    if "z" not in df.columns:
        raise ValueError("CSV must contain column 'z'")
    z = df["z"].to_numpy(dtype=float)
    if np.any(z <= 0):
        raise ValueError("all z values must be positive")

    if "D_obs" in df.columns:
        d_obs = df["D_obs"].to_numpy(dtype=float)
    elif "mu" in df.columns:
        d_obs = distance_modulus_to_mly(df["mu"].to_numpy(dtype=float))
    else:
        raise ValueError("CSV must contain either 'D_obs' or 'mu'")

    b_each = infer_b_from_distance(z, d_obs, L=L_STAM)
    work = pd.DataFrame({"z": z, "D_obs": d_obs, "b_inferred": b_each})
    work = work.sort_values("z").reset_index(drop=True)
    work["bin"] = pd.qcut(work["z"], q=min(args.bins, len(work)), duplicates="drop")

    grouped = work.groupby("bin", observed=True).agg(
        n=("z", "size"),
        z_min=("z", "min"),
        z_max=("z", "max"),
        z_mean=("z", "mean"),
        b_mean=("b_inferred", "mean"),
        b_std=("b_inferred", "std"),
        b_median=("b_inferred", "median"),
    ).reset_index()
    grouped["bin"] = grouped["bin"].astype(str)

    # Simple linear trend diagnostic: b_inferred = slope*z + intercept.
    slope, intercept = np.polyfit(work["z"].to_numpy(), work["b_inferred"].to_numpy(), deg=1)

    result = {
        "catalog_name": args.name,
        "input_file": str(args.csv),
        "n": int(len(work)),
        "bins_requested": args.bins,
        "bins_used": int(len(grouped)),
        "b_global_mean": float(work["b_inferred"].mean()),
        "b_global_std": float(work["b_inferred"].std(ddof=1)) if len(work) > 1 else 0.0,
        "b_vs_z_slope": float(slope),
        "b_vs_z_intercept": float(intercept),
        "label": "post_hoc_diagnostic_not_prediction",
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / f"{args.name}_redshift_bin_b_summary.json"
    csv_path = args.output_dir / f"{args.name}_redshift_bin_b_table.csv"
    rows_path = args.output_dir / f"{args.name}_rowwise_b_inferred.csv"

    json_path.write_text(json.dumps(result, indent=2, sort_keys=True))
    grouped.to_csv(csv_path, index=False)
    work.to_csv(rows_path, index=False)

    print(json.dumps(result, indent=2, sort_keys=True))
    print(grouped.to_string(index=False))
    print(f"Wrote {json_path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {rows_path}")


if __name__ == "__main__":
    main()

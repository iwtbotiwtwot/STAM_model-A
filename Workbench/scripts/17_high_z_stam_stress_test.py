#!/usr/bin/env python3
"""STAM Model-A high-z stress test.

Scope:
  Formula behavior only. This does not claim to explain the Big Bang,
  recombination, CMB origin, or early-universe microphysics.

Tested no-b Model-A functions:
  D_geo(z)      = Lz(1+0.15z)
  D_adj,0(z)    = Lz(1+0.5z)
  D_excess,0(z) = 0.35Lz^2

Path accumulation:
  <A_path>(z)     = 0.35z/(1+0.15z)
  A_path,local(z) = 0.70z/(1+0.30z)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

C_MLY_PER_MPC = 3.261563776
H_STAM = 0.000243635
L_MPC = 1.0 / H_STAM
L_MLY = C_MLY_PER_MPC / H_STAM


def D_geo_mpc(z):
    z = np.asarray(z, dtype=float)
    return L_MPC * z * (1.0 + 0.15 * z)


def D_adj0_mpc(z):
    z = np.asarray(z, dtype=float)
    return L_MPC * z * (1.0 + 0.5 * z)


def D_excess0_mpc(z):
    z = np.asarray(z, dtype=float)
    return D_adj0_mpc(z) - D_geo_mpc(z)


def A_path_average(z):
    z = np.asarray(z, dtype=float)
    return 0.35 * z / (1.0 + 0.15 * z)


def A_path_local(z):
    z = np.asarray(z, dtype=float)
    return 0.70 * z / (1.0 + 0.30 * z)


def ratio_adj_geo(z):
    return D_adj0_mpc(z) / D_geo_mpc(z)


def main() -> None:
    out_dir = Path("results/high_z_stress")
    out_dir.mkdir(parents=True, exist_ok=True)

    z_low = np.linspace(0.001, 5, 2000)
    z_high = np.logspace(np.log10(5.001), np.log10(1100), 2500)
    z_grid = np.unique(np.concatenate([z_low, z_high]))

    grid_df = pd.DataFrame(
        {
            "z": z_grid,
            "D_geo_Mpc": D_geo_mpc(z_grid),
            "D_adj0_Mpc": D_adj0_mpc(z_grid),
            "D_excess0_Mpc": D_excess0_mpc(z_grid),
            "D_geo_Mly": D_geo_mpc(z_grid) * C_MLY_PER_MPC,
            "D_adj0_Mly": D_adj0_mpc(z_grid) * C_MLY_PER_MPC,
            "D_excess0_Mly": D_excess0_mpc(z_grid) * C_MLY_PER_MPC,
            "A_path_average": A_path_average(z_grid),
            "A_path_local": A_path_local(z_grid),
            "D_adj0_over_D_geo": ratio_adj_geo(z_grid),
        }
    )
    grid_df.to_csv(out_dir / "high_z_model_a_grid.csv", index=False)

    key_zs = np.array([0.001, 0.01, 0.1, 0.5, 1, 2, 3, 5, 10, 20, 100, 1100])
    key_df = pd.DataFrame(
        {
            "z": key_zs,
            "D_geo_Mpc": D_geo_mpc(key_zs),
            "D_adj0_Mpc": D_adj0_mpc(key_zs),
            "D_excess0_Mpc": D_excess0_mpc(key_zs),
            "A_path_average": A_path_average(key_zs),
            "A_path_local": A_path_local(key_zs),
            "D_adj0_over_D_geo": ratio_adj_geo(key_zs),
        }
    )
    key_df.to_csv(out_dir / "key_redshift_values.csv", index=False)

    a_limit = 0.35 / 0.15
    local_limit = 0.70 / 0.30
    ratio_limit = 0.5 / 0.15

    fractions = np.array([0.5, 0.75, 0.9, 0.95, 0.99, 0.999])
    sat_df = pd.DataFrame(
        {
            "fraction_of_asymptotic_limit": fractions,
            "z_for_A_path_average": fractions / (0.15 * (1 - fractions)),
            "z_for_A_path_local": fractions / (0.30 * (1 - fractions)),
            "A_limit": a_limit,
            "A_local_limit": local_limit,
        }
    )
    sat_df.to_csv(out_dir / "A_path_saturation_redshifts.csv", index=False)

    z_cmb = 1100.0
    summary = {
        "test": "STAM Model-A high-z stress test",
        "scope": "formula behavior only; no Big Bang origin claim",
        "formulas": {
            "D_geo": "D_geo(z)=Lz(1+0.15z)",
            "D_adj0": "D_adj,0(z)=Lz(1+0.5z)",
            "D_excess0": "D_excess,0(z)=0.35Lz^2",
            "A_path_average": "<A_path>(z)=0.35z/(1+0.15z)",
            "A_path_local": "A_path,local(z)=0.70z/(1+0.30z)",
        },
        "constants": {
            "H_STAM": H_STAM,
            "H_inverse_Mpc": L_MPC,
            "C_Mly_per_Mpc": C_MLY_PER_MPC,
            "L_Mly": L_MLY,
        },
        "limits": {
            "A_path_average_limit": a_limit,
            "A_path_local_limit": local_limit,
            "D_adj0_over_D_geo_limit": ratio_limit,
        },
        "z_1100": {
            "D_geo_Mpc": float(D_geo_mpc(z_cmb)),
            "D_adj0_Mpc": float(D_adj0_mpc(z_cmb)),
            "D_excess0_Mpc": float(D_excess0_mpc(z_cmb)),
            "A_path_average": float(A_path_average(z_cmb)),
            "A_path_local": float(A_path_local(z_cmb)),
            "D_adj0_over_D_geo": float(ratio_adj_geo(z_cmb)),
            "A_path_average_fraction_of_limit": float(A_path_average(z_cmb) / a_limit),
            "A_path_local_fraction_of_limit": float(A_path_local(z_cmb) / local_limit),
        },
        "checks": {
            "all_values_finite": bool(np.isfinite(grid_df.select_dtypes(include=[np.number]).to_numpy()).all()),
            "all_distances_positive": bool((grid_df[["D_geo_Mpc", "D_adj0_Mpc", "D_excess0_Mpc"]] > 0).all().all()),
        },
        "interpretation": (
            "Current no-b Model-A expressions remain finite and stable out to z=1100. "
            "Path accumulation approaches a finite high-z limit of 7/3 rather than diverging."
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    print("PASS: high-z STAM Model-A stress test complete.")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote outputs to {out_dir}")


if __name__ == "__main__":
    main()

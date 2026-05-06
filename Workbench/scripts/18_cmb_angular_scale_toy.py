#!/usr/bin/env python3
"""STAM Model-A CMB angular-scale toy test.

Scope:
  Toy angular-scale check only.
  This is not a model of Big Bang origin, recombination, CMB microphysics,
  or the CMB power spectrum.

Reference:
  Planck reports 100 theta_* ≈ 1.0411, so theta_* ≈ 0.010411 rad.

Test:
  Use current no-b Model-A distance layers at z≈1090 and ask whether simple
  angular mappings reproduce theta ≈ r_s / D_M.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

C_MLY_PER_MPC = 3.261563776
H_STAM = 0.000243635
L_MPC = 1.0 / H_STAM

Z_STAR = 1090.0
THETA_STAR_OBS = 1.0411 / 100.0
SOUND_HORIZONS_MPC = {
    "r_s_star_144.4Mpc_toy": 144.4,
    "r_d_drag_147.09Mpc_reference": 147.09,
}


def D_geo_mpc(z):
    z = np.asarray(z, dtype=float)
    return L_MPC * z * (1.0 + 0.15 * z)


def D_adj0_mpc(z):
    z = np.asarray(z, dtype=float)
    return L_MPC * z * (1.0 + 0.5 * z)


def D_excess0_mpc(z):
    return D_adj0_mpc(z) - D_geo_mpc(z)


def candidate_distances(z):
    d_geo = float(D_geo_mpc(z))
    d_adj = float(D_adj0_mpc(z))
    d_exc = float(D_excess0_mpc(z))
    return {
        "D_geo": d_geo,
        "D_adj0": d_adj,
        "D_excess0": d_exc,
        "D_geo_over_1pz": d_geo / (1.0 + z),
        "D_adj0_over_1pz": d_adj / (1.0 + z),
        "D_excess0_over_1pz": d_exc / (1.0 + z),
        "D_geo_over_1pz2": d_geo / (1.0 + z) ** 2,
        "D_adj0_over_1pz2": d_adj / (1.0 + z) ** 2,
    }


def main() -> None:
    out_dir = Path("results/cmb_angular_toy")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for z_label, z in [("z_star_1090", Z_STAR), ("z_1100", 1100.0)]:
        distances = candidate_distances(z)
        for r_name, r_s in SOUND_HORIZONS_MPC.items():
            d_required = r_s / THETA_STAR_OBS
            for mapping, d in distances.items():
                theta = r_s / d
                rows.append(
                    {
                        "z_label": z_label,
                        "z": z,
                        "sound_horizon_name": r_name,
                        "sound_horizon_Mpc": r_s,
                        "mapping": mapping,
                        "D_mapping_Mpc": d,
                        "D_required_Mpc_for_observed_theta": d_required,
                        "theta_pred_rad": theta,
                        "theta_obs_rad": THETA_STAR_OBS,
                        "theta_ratio_pred_over_obs": theta / THETA_STAR_OBS,
                        "D_ratio_mapping_over_required": d / d_required,
                        "abs_log10_theta_ratio": abs(math.log10(theta / THETA_STAR_OBS)),
                    }
                )

    results_df = pd.DataFrame(rows)
    results_df.to_csv(out_dir / "cmb_angular_mapping_results.csv", index=False)

    best_df = (
        results_df.sort_values("abs_log10_theta_ratio")
        .groupby(["z_label", "sound_horizon_name"], as_index=False)
        .first()
    )
    best_df.to_csv(out_dir / "cmb_angular_best_mappings.csv", index=False)

    best = best_df[
        (best_df["z_label"] == "z_star_1090")
        & (best_df["sound_horizon_name"] == "r_s_star_144.4Mpc_toy")
    ].iloc[0].to_dict()

    summary = {
        "test": "STAM Model-A CMB angular-scale toy test",
        "scope": "toy only; no CMB origin or recombination claim",
        "observational_reference": {
            "Planck_2018_100_theta_star": 1.0411,
            "theta_star_obs_rad": THETA_STAR_OBS,
            "l_acoustic_approx_pi_over_theta": math.pi / THETA_STAR_OBS,
        },
        "best_simple_mapping_at_z1090": {
            "mapping": best["mapping"],
            "D_mapping_Mpc": best["D_mapping_Mpc"],
            "D_required_Mpc": best["D_required_Mpc_for_observed_theta"],
            "theta_pred_rad": best["theta_pred_rad"],
            "theta_obs_rad": best["theta_obs_rad"],
            "theta_ratio_pred_over_obs": best["theta_ratio_pred_over_obs"],
            "D_ratio_mapping_over_required": best["D_ratio_mapping_over_required"],
        },
        "interpretation": (
            "No tested simple STAM distance-layer mapping at z~1090 reproduces the observed CMB acoustic angular scale. "
            "A CMB-facing distance/visibility rule is needed."
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    print("PASS: CMB angular-scale toy test complete.")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote outputs to {out_dir}")


if __name__ == "__main__":
    main()

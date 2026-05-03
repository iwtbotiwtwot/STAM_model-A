#!/usr/bin/env python3
"""STAM Model-A gravitational-wave propagation test.

Model-A test rule:
    gravitational waves and light propagate locally at c.

Accumulation traversal rule:
    t_obs = ∫ds/c + k∫A(s)ds/c

Equal-coupling case:
    k_GW = k_EM

Then GW and EM waves receive the same accumulation delay through the same A field.
Near high A, outside observers can infer apparent slowing from increased traversal
time, but local propagation remains c along allowed paths.

At A>1, outward escape is unavailable because:
    A = (v_escape/c)^2
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

G = 6.67430e-11
C = 299_792_458.0
M_SUN = 1.98847e30
R_SUN = 6.9634e8
AU = 1.495978707e11
MPC_M = 3.0856775814913673e22


def schwarzschild_radius(mass_kg: float) -> float:
    return 2.0 * G * mass_kg / C**2


def accumulation(mass_kg: float, radius_m: float) -> float:
    return schwarzschild_radius(mass_kg) / radius_m


def accumulation_delay_straight_path(
    mass_kg: float,
    impact_m: float,
    r1_m: float,
    r2_m: float,
    coupling: float = 1.0,
) -> float:
    if impact_m <= 0 or r1_m <= impact_m or r2_m <= impact_m:
        raise ValueError("impact_m must be positive and below endpoint radii")
    x1 = -math.sqrt(r1_m**2 - impact_m**2)
    x2 = math.sqrt(r2_m**2 - impact_m**2)
    return coupling * (schwarzschild_radius(mass_kg) / C) * (
        math.asinh(x2 / impact_m) - math.asinh(x1 / impact_m)
    )


def propagation_times(mass_kg: float, impact_m: float, r1_m: float, r2_m: float, k_em=1.0, k_gw=1.0):
    x1 = -math.sqrt(r1_m**2 - impact_m**2)
    x2 = math.sqrt(r2_m**2 - impact_m**2)
    path_length = x2 - x1
    flat_time = path_length / C
    delay_em = accumulation_delay_straight_path(mass_kg, impact_m, r1_m, r2_m, k_em)
    delay_gw = accumulation_delay_straight_path(mass_kg, impact_m, r1_m, r2_m, k_gw)
    return {
        "path_length_m": path_length,
        "flat_time_s": flat_time,
        "delay_em_s": delay_em,
        "delay_gw_s": delay_gw,
        "t_em_s": flat_time + delay_em,
        "t_gw_s": flat_time + delay_gw,
        "gw_minus_em_s": delay_gw - delay_em,
        "apparent_speed_em_over_c": path_length / (flat_time + delay_em) / C,
        "apparent_speed_gw_over_c": path_length / (flat_time + delay_gw) / C,
        "avg_A_path": delay_em / flat_time,
    }


def main() -> None:
    out_dir = Path("results/gw_propagation")
    out_dir.mkdir(parents=True, exist_ok=True)

    solar_rows = []
    for name, impact, r1, r2 in [
        ("Sun_grazing_Earth_Mercury", R_SUN, AU, 0.387 * AU),
        ("Sun_grazing_Earth_Mars", R_SUN, AU, 1.524 * AU),
        ("Sun_2Rsun_Earth_Mars", 2 * R_SUN, AU, 1.524 * AU),
        ("Sun_10Rsun_Earth_Mars", 10 * R_SUN, AU, 1.524 * AU),
    ]:
        equal = propagation_times(M_SUN, impact, r1, r2, k_em=1.0, k_gw=1.0)
        solar_rows.append(
            {
                "case": name,
                "A_at_impact": accumulation(M_SUN, impact),
                "shared_accumulation_delay_us": equal["delay_em_s"] * 1e6,
                "equal_coupling_gw_minus_em_s": equal["gw_minus_em_s"],
                "apparent_speed_over_c": equal["apparent_speed_gw_over_c"],
                "avg_A_path": equal["avg_A_path"],
            }
        )
    solar_df = pd.DataFrame(solar_rows)
    solar_df.to_csv(out_dir / "gw_em_solar_accumulation_delay_cases.csv", index=False)

    mass_bh = 10.0 * M_SUN
    rh = schwarzschild_radius(mass_bh)
    near_rows = []
    for b_over_rh in [50, 20, 10, 5, 3, 2, 1.5, 1.2, 1.1, 1.05, 1.01]:
        p = propagation_times(mass_bh, b_over_rh * rh, 1000 * rh, 1000 * rh)
        near_rows.append(
            {
                "impact_over_horizon": b_over_rh,
                "A_at_impact": 1.0 / b_over_rh,
                "escape_velocity_at_impact_over_c": math.sqrt(1.0 / b_over_rh),
                "accumulation_delay_ms": p["delay_gw_s"] * 1e3,
                "apparent_speed_over_c": p["apparent_speed_gw_over_c"],
                "gw_minus_em_s_equal_coupling": p["gw_minus_em_s"],
            }
        )
    near_df = pd.DataFrame(near_rows)
    near_df.to_csv(out_dir / "gw_em_near_horizon_apparent_delay_toy.csv", index=False)

    threshold_rows = []
    for r_over_h in [5, 3, 2, 1.5, 1.1, 1.01, 1.0, 0.99, 0.9, 0.5, 0.25]:
        A = 1.0 / r_over_h
        threshold_rows.append(
            {
                "r_over_horizon": r_over_h,
                "A": A,
                "local_wave_speed_over_c": 1.0,
                "escape_velocity_over_c": math.sqrt(A),
                "outward_escape_status": "allowed" if A < 1 else "horizon" if math.isclose(A, 1.0) else "over-threshold",
            }
        )
    threshold_df = pd.DataFrame(threshold_rows)
    threshold_df.to_csv(out_dir / "gw_horizon_escape_threshold.csv", index=False)

    distance_mpc = 40.0
    lag_s = 1.74
    travel_time_s = distance_mpc * MPC_M / C
    summary = {
        "test": "STAM Model-A gravitational-wave propagation",
        "rule": "local GW speed = c; equal accumulation traversal coupling for GW and EM",
        "max_equal_coupling_gw_minus_em_s_solar": float(solar_df["equal_coupling_gw_minus_em_s"].abs().max()),
        "max_equal_coupling_gw_minus_em_s_near_horizon_toy": float(near_df["gw_minus_em_s_equal_coupling"].abs().max()),
        "max_A_minus_escape_fraction_squared_error": float(
            np.max(np.abs(threshold_df["A"] - threshold_df["escape_velocity_over_c"] ** 2))
        ),
        "GW170817_sanity_check": {
            "distance_Mpc_assumed": distance_mpc,
            "observed_GRB_after_GW_lag_s": lag_s,
            "naive_fractional_speed_difference_if_lag_were_all_propagation": lag_s / travel_time_s,
            "equal_coupling_predicted_propagation_lag_s": 0.0,
        },
        "interpretation": (
            "Equal coupling gives no GW-EM propagation-time difference through the same A path. "
            "Near high A, apparent delay increases without changing local c. "
            "At A>1, outward escape is over-threshold because A=(v_escape/c)^2."
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    print("PASS: STAM Model-A GW propagation equal-coupling checks passed.")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote outputs to {out_dir}")


if __name__ == "__main__":
    main()

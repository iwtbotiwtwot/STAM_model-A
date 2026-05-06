#!/usr/bin/env python3
"""STAM Model-A black-hole decay / horizon recession test.

Core STAM rules:
    A(r,t)=2GM(t)/(c^2 r)
    r_h(t)=2GM(t)/c^2

Mass balance:
    dM/dt = Mdot_in - P_out/c^2

Therefore:
    dr_h/dt = (2G/c^2)dM/dt

If P_out > Mdot_in c^2, then dM/dt < 0 and the A=1 horizon recedes.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

G = 6.67430e-11
C = 299_792_458.0
HBAR = 1.054571817e-34
M_SUN = 1.98847e30
YEAR = 365.25 * 86400.0
DAY = 86400.0
KM = 1000.0


def schwarzschild_radius(mass_kg: float) -> float:
    return 2.0 * G * mass_kg / C**2


def accumulation_at_radius(mass_kg: float, radius_m: float) -> float:
    return schwarzschild_radius(mass_kg) / radius_m


def mass_rate(mdot_in_kg_s: float, p_out_w: float) -> float:
    return mdot_in_kg_s - p_out_w / C**2


def horizon_rate(dmdt_kg_s: float) -> float:
    return 2.0 * G * dmdt_kg_s / C**2


def power_required_for_mass_fraction_loss(
    mass_kg: float,
    fraction_loss: float,
    duration_s: float,
    mdot_in_kg_s: float = 0.0,
) -> float:
    target_dmdt = -(fraction_loss * mass_kg) / duration_s
    return (mdot_in_kg_s - target_dmdt) * C**2


def hawking_power_w(mass_kg: float) -> float:
    """Optional reference: Hawking-radiation power for neutral non-rotating BH."""
    return HBAR * C**6 / (15360.0 * math.pi * G**2 * mass_kg**2)


def main() -> None:
    out_dir = Path("results/black_hole_decay")
    out_dir.mkdir(parents=True, exist_ok=True)

    m0 = 10.0 * M_SUN
    rh0 = schwarzschild_radius(m0)

    mass_fractions = np.array([1.25, 1.10, 1.01, 1.00, 0.99, 0.95, 0.90, 0.75, 0.50, 0.25, 0.10])
    rows = []
    for f in mass_fractions:
        mass = f * m0
        rh = schwarzschild_radius(mass)
        rows.append(
            {
                "mass_fraction_M_over_M0": f,
                "mass_solar": mass / M_SUN,
                "horizon_radius_km": rh / KM,
                "horizon_radius_fraction_rh_over_rh0": rh / rh0,
                "A_at_initial_horizon_radius": accumulation_at_radius(mass, rh0),
                "initial_horizon_radius_status_now": (
                    "inside/over-threshold"
                    if accumulation_at_radius(mass, rh0) > 1
                    else "horizon"
                    if math.isclose(accumulation_at_radius(mass, rh0), 1.0)
                    else "outside/sub-threshold"
                ),
            }
        )
    threshold_df = pd.DataFrame(rows)
    threshold_df.to_csv(out_dir / "mass_fraction_horizon_recession.csv", index=False)

    scenarios = [
        {
            "scenario": "exact_balance_no_horizon_motion",
            "mass_solar": 10.0,
            "mdot_in_kg_s": 1.0e10,
            "p_out_w": 1.0e10 * C**2,
            "duration_s": YEAR,
        },
        {
            "scenario": "net_growth_accretion_dominates",
            "mass_solar": 10.0,
            "mdot_in_kg_s": 1.0e12,
            "p_out_w": 1.0e28,
            "duration_s": YEAR,
        },
        {
            "scenario": "net_decay_power_dominates",
            "mass_solar": 10.0,
            "mdot_in_kg_s": 1.0e12,
            "p_out_w": 1.0e32,
            "duration_s": YEAR,
        },
        {
            "scenario": "toy_5_percent_mass_loss_in_one_year",
            "mass_solar": 10.0,
            "mdot_in_kg_s": 0.0,
            "p_out_w": power_required_for_mass_fraction_loss(m0, 0.05, YEAR),
            "duration_s": YEAR,
        },
    ]

    rate_rows = []
    for s in scenarios:
        mass = s["mass_solar"] * M_SUN
        rh = schwarzschild_radius(mass)
        dmdt = mass_rate(s["mdot_in_kg_s"], s["p_out_w"])
        drdt = horizon_rate(dmdt)
        mass_end = max(0.0, mass + dmdt * s["duration_s"])
        rh_end = schwarzschild_radius(mass_end) if mass_end > 0 else 0.0
        rate_rows.append(
            {
                **s,
                "dMdt_kg_s": dmdt,
                "net_status": "growth" if dmdt > 0 else "decay" if dmdt < 0 else "balanced",
                "dr_h_dt_m_s": drdt,
                "dr_h_dt_m_per_year": drdt * YEAR,
                "initial_horizon_radius_km": rh / KM,
                "final_mass_solar": mass_end / M_SUN,
                "final_horizon_radius_km": rh_end / KM,
                "horizon_change_km_over_duration": (rh_end - rh) / KM,
                "A_at_initial_horizon_after_duration": rh_end / rh if rh > 0 else np.nan,
            }
        )
    rate_df = pd.DataFrame(rate_rows)
    rate_df.to_csv(out_dir / "mass_balance_rate_scenarios.csv", index=False)

    hawking_rows = []
    for name, mass in [
        ("10 solar mass", 10.0 * M_SUN),
        ("1 solar mass", 1.0 * M_SUN),
        ("1e15 kg primordial-scale", 1.0e15),
        ("1e12 kg small BH", 1.0e12),
        ("1e11 kg small BH", 1.0e11),
    ]:
        p = hawking_power_w(mass)
        dmdt = -p / C**2
        hawking_rows.append(
            {
                "object": name,
                "mass_kg": mass,
                "horizon_radius_km": schwarzschild_radius(mass) / KM,
                "hawking_power_W_optional_reference": p,
                "hawking_dMdt_kg_s": dmdt,
                "hawking_drh_dt_m_per_year": horizon_rate(dmdt) * YEAR,
            }
        )
    pd.DataFrame(hawking_rows).to_csv(out_dir / "optional_hawking_reference_examples.csv", index=False)

    checks = {
        "horizon_radius_scales_linearly_with_mass_max_error": float(
            np.max(np.abs(threshold_df["horizon_radius_fraction_rh_over_rh0"] - threshold_df["mass_fraction_M_over_M0"]))
        ),
        "A_at_initial_horizon_equals_mass_fraction_max_error": float(
            np.max(np.abs(threshold_df["A_at_initial_horizon_radius"] - threshold_df["mass_fraction_M_over_M0"]))
        ),
        "decay_scenarios_have_negative_drhdt": bool((rate_df[rate_df["net_status"] == "decay"]["dr_h_dt_m_s"] < 0).all()),
        "growth_scenarios_have_positive_drhdt": bool((rate_df[rate_df["net_status"] == "growth"]["dr_h_dt_m_s"] > 0).all()),
    }

    summary = {
        "test": "STAM Model-A black-hole decay and horizon recession",
        "formulas": {
            "A": "A(r,t)=2GM(t)/(c^2 r)",
            "horizon": "r_h(t)=2GM(t)/c^2",
            "mass_balance": "dM/dt=Mdot_in-P_out/c^2",
            "horizon_rate": "dr_h/dt=(2G/c^2)dM/dt",
            "decay_condition": "P_out > Mdot_in c^2",
        },
        "checks": checks,
        "ten_solar_mass_initial_horizon_radius_km": rh0 / KM,
        "ten_solar_mass_after_5_percent_loss_horizon_radius_km": schwarzschild_radius(0.95 * m0) / KM,
        "A_at_old_horizon_after_5_percent_loss": accumulation_at_radius(0.95 * m0, rh0),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    print("PASS: STAM Model-A black-hole decay / horizon recession checks complete.")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote outputs to {out_dir}")


if __name__ == "__main__":
    main()

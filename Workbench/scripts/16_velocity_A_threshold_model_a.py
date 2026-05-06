#!/usr/bin/env python3
"""STAM Model-A moving mass-energy A=1 threshold.

If total energy sources A, then a moving object has:

    E = gamma m c^2
    A_v(R) = gamma * 2Gm/(c^2 R)

At the object's effective radius R:

    A_v = 1

gives:

    gamma_threshold = R / Rs_rest

where:

    Rs_rest = 2Gm/c^2

This script asks at what mass/radius/velocity the moving object's A reaches 1.
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
M_EARTH = 5.9722e24
R_EARTH = 6_378_137.0
R_SUN = 6.9634e8


def schwarzschild_radius_rest(mass_kg: float) -> float:
    return 2.0 * G * mass_kg / C**2


def gamma_from_beta(beta: float) -> float:
    if beta >= 1:
        return math.inf
    return 1.0 / math.sqrt(1.0 - beta**2)


def beta_from_gamma(gamma: float) -> float:
    if gamma < 1:
        return 0.0
    return math.sqrt(1.0 - 1.0 / gamma**2)


def one_minus_beta_approx(gamma: float) -> float:
    return 1.0 / (2.0 * gamma**2)


def threshold_row(name: str, mass_kg: float, radius_m: float) -> dict[str, float | str]:
    rs = schwarzschild_radius_rest(mass_kg)
    A_rest = rs / radius_m
    gamma_thr = radius_m / rs
    beta_thr = beta_from_gamma(gamma_thr) if gamma_thr < 1e154 else 1.0
    return {
        "object": name,
        "mass_kg": mass_kg,
        "radius_m": radius_m,
        "rest_schwarzschild_radius_m": rs,
        "A_rest_at_radius": A_rest,
        "gamma_threshold_for_A_eq_1": gamma_thr,
        "beta_threshold_v_over_c": beta_thr,
        "one_minus_beta_threshold_approx": one_minus_beta_approx(gamma_thr),
        "check_gamma_times_A_rest": gamma_thr * A_rest,
    }


def main() -> None:
    out_dir = Path("results/velocity_A_threshold")
    out_dir.mkdir(parents=True, exist_ok=True)

    objects = [
        ("1 kg, 0.1 m radius", 1.0, 0.1),
        ("100 kg person-scale, 0.5 m radius", 100.0, 0.5),
        ("1000 kg craft, 2 m radius", 1_000.0, 2.0),
        ("10000 kg craft, 5 m radius", 10_000.0, 5.0),
        ("Earth", M_EARTH, R_EARTH),
        ("Sun", M_SUN, R_SUN),
        ("Neutron star toy: 1.4 Msun, 12 km radius", 1.4*M_SUN, 12_000.0),
    ]

    threshold_df = pd.DataFrame([threshold_row(*obj) for obj in objects])
    threshold_df.to_csv(out_dir / "object_A1_velocity_thresholds.csv", index=False)

    betas = [0.0, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999]
    scan_rows = []
    for name, mass, radius in objects:
        rs = schwarzschild_radius_rest(mass)
        A0 = rs / radius
        for beta in betas:
            gamma = gamma_from_beta(beta)
            scan_rows.append(
                {
                    "object": name,
                    "beta_v_over_c": beta,
                    "gamma": gamma,
                    "A_rest": A0,
                    "A_moving": gamma * A0,
                    "kinetic_A_added": (gamma - 1.0) * A0,
                    "threshold_reached": gamma * A0 >= 1.0,
                }
            )
    scan_df = pd.DataFrame(scan_rows)
    scan_df.to_csv(out_dir / "A_vs_velocity_scan.csv", index=False)

    checks = {
        "max_abs_gamma_threshold_times_A_rest_minus_1": float(
            np.max(np.abs(threshold_df["check_gamma_times_A_rest"] - 1.0))
        ),
        "all_gamma_thresholds_ge_1": bool((threshold_df["gamma_threshold_for_A_eq_1"] >= 1.0).all()),
    }

    summary = {
        "test": "STAM Model-A moving mass-energy A=1 threshold",
        "formulas": {
            "A_rest": "A0(R)=Rs_rest/R=2Gm/(c^2 R)",
            "moving_A": "A_v(R)=gamma*A0(R)",
            "threshold": "A_v=1 => gamma_threshold=R/Rs_rest",
            "velocity": "v_threshold/c=sqrt(1-1/gamma_threshold^2)",
        },
        "checks": checks,
        "interpretation": (
            "If total energy sources A, then motion increases the object's effective accumulation contribution by gamma. "
            "A=1 occurs when gamma Rs_rest/R = 1."
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    print("PASS: moving mass-energy A=1 threshold calculation complete.")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote outputs to {out_dir}")


if __name__ == "__main__":
    main()

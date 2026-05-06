#!/usr/bin/env python3
"""STAM Model-A galactic cumulative A scratch test.

Purpose:
  1. Check whether individually tiny A contributions can sum into a galaxy-scale field.
  2. Compare a linear visible-matter A sum with an optional exploratory collective A-envelope.

Core linear Model-A:
    A_total(x)=Σ 2Gm_i/(c²|x-x_i|)

Claim boundary:
  The linear sum is Model-A-native.
  The collective-envelope term is exploratory and not yet a formal Model-A claim.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import mpmath as mp
    HAS_MPMATH = True
except Exception:
    HAS_MPMATH = False

G = 6.67430e-11
C = 299_792_458.0
M_SUN = 1.98847e30
KPC = 3.0856775814913673e19

M_DISK = 6.0e10 * M_SUN
M_GAS = 1.0e10 * M_SUN
M_BULGE = 1.0e10 * M_SUN
M_TOTAL_TOY = M_DISK + M_GAS + M_BULGE
R_D_DISK = 3.0 * KPC
R_D_GAS = 7.0 * KPC
BULGE_SOFTENING = 0.5 * KPC
DISK_SOFTENING = 0.05 * KPC

V_EDGE_TARGET = 220_000.0
R_REF = 20.0 * KPC
R_C_ENV = 1.0 * KPC
R_A_ENV = 100.0 * KPC


def schwarzschild_radius(mass_kg: float) -> float:
    return 2.0 * G * mass_kg / C**2


def A_point(mass_kg: float, radius_m: float) -> float:
    return schwarzschild_radius(mass_kg) / radius_m


def v_circ_compact_from_A(A: float) -> float:
    return C * math.sqrt(max(A, 0.0) / 2.0)


def exponential_disk_elements(M: float, Rd: float, r_max: float = 40*KPC, n_r: int = 96, n_phi: int = 240):
    edges = np.linspace(0.0, r_max, n_r + 1)
    r_mid = 0.5 * (edges[:-1] + edges[1:])
    dr = edges[1:] - edges[:-1]
    x = r_max / Rd
    norm = Rd**2 * (1.0 - math.exp(-x) * (1.0 + x))
    ring_mass = M * (r_mid * np.exp(-r_mid/Rd) * dr) / norm
    phi = np.linspace(0.0, 2.0*np.pi, n_phi, endpoint=False)
    xs, ys, ms = [], [], []
    for R, m_ring in zip(r_mid, ring_mass):
        xs.append(R * np.cos(phi))
        ys.append(R * np.sin(phi))
        ms.append(np.full(n_phi, m_ring / n_phi))
    return np.concatenate(xs), np.concatenate(ys), np.concatenate(ms)


def disk_accumulation_and_accel_at_R(R, x_src, y_src, m_src, softening):
    dx = R - x_src
    dy = -y_src
    dist2 = dx*dx + dy*dy + softening*softening
    dist = np.sqrt(dist2)
    A = np.sum(2.0 * G * m_src / (C**2 * dist))
    gx = -np.sum(G * m_src * dx / (dist2 * dist))
    return float(A), float(max(-gx, 0.0))


def bulge_accel_plummer(R, M=M_BULGE, a=BULGE_SOFTENING):
    return G * M * R / (R*R + a*a)**1.5


def bulge_A_softened(R, M=M_BULGE, a=BULGE_SOFTENING):
    return 2.0 * G * M / (C**2 * math.sqrt(R*R + a*a))


def A_env(R, eps, rc=R_C_ENV, RA=R_A_ENV):
    if R >= RA:
        return 0.0
    return eps * math.log((RA + rc) / (R + rc))


def g_env(R, eps, rc=R_C_ENV, RA=R_A_ENV):
    if R >= RA:
        return 0.0
    return (C**2 / 2.0) * eps / (R + rc)


def main() -> None:
    out_dir = Path("results/galactic_cumulative_A")
    out_dir.mkdir(parents=True, exist_ok=True)

    high_rows = []
    if HAS_MPMATH:
        mp.mp.dps = 1000
        G_mp = mp.mpf("6.67430e-11")
        C_mp = mp.mpf("299792458")
        M_sun_mp = mp.mpf("1.98847e30")
        KPC_mp = mp.mpf("3.0856775814913673e19")
    for radius_kpc in [1, 5, 10, 20, 50, 100]:
        if HAS_MPMATH:
            r_mp = mp.mpf(str(radius_kpc)) * KPC_mp
            A_one = (2*G_mp*M_sun_mp)/(C_mp**2 * r_mp)
            A_1e11 = (2*G_mp*(mp.mpf("1e11")*M_sun_mp))/(C_mp**2 * r_mp)
            high_rows.append(
                {
                    "radius_kpc": radius_kpc,
                    "A_one_solar_mass_1000dps": mp.nstr(A_one, 80),
                    "A_1e11_solar_masses_1000dps": mp.nstr(A_1e11, 80),
                    "A_one_float": A_point(M_SUN, radius_kpc*KPC),
                    "A_1e11_float": A_point(1e11*M_SUN, radius_kpc*KPC),
                }
            )
        else:
            high_rows.append(
                {
                    "radius_kpc": radius_kpc,
                    "A_one_float": A_point(M_SUN, radius_kpc*KPC),
                    "A_1e11_float": A_point(1e11*M_SUN, radius_kpc*KPC),
                }
            )
    high_df = pd.DataFrame(high_rows)
    high_df.to_csv(out_dir / "high_precision_A_scale_check.csv", index=False)

    x_disk, y_disk, m_disk = exponential_disk_elements(M_DISK, R_D_DISK)
    x_gas, y_gas, m_gas = exponential_disk_elements(M_GAS, R_D_GAS)
    x_src = np.concatenate([x_disk, x_gas])
    y_src = np.concatenate([y_disk, y_gas])
    m_src = np.concatenate([m_disk, m_gas])

    disk_rows = []
    for R_kpc in np.linspace(1.0, 40.0, 160):
        R = R_kpc * KPC
        A_diskgas, g_diskgas = disk_accumulation_and_accel_at_R(R, x_src, y_src, m_src, DISK_SOFTENING)
        A_linear = A_diskgas + bulge_A_softened(R)
        g_linear = g_diskgas + bulge_accel_plummer(R)
        v_linear = math.sqrt(max(R*g_linear, 0.0))
        A_compact = A_point(M_TOTAL_TOY, R)
        disk_rows.append(
            {
                "R_kpc": R_kpc,
                "A_linear_visible_disk_bulge_gas": A_linear,
                "A_compact_total_same_mass": A_compact,
                "v_linear_visible_km_s": v_linear / 1000,
                "v_compact_same_mass_km_s": v_circ_compact_from_A(A_compact) / 1000,
            }
        )
    disk_df = pd.DataFrame(disk_rows)
    disk_df.to_csv(out_dir / "linear_cumulative_visible_disk_toy.csv", index=False)

    idx_ref = (disk_df["R_kpc"] - R_REF/KPC).abs().idxmin()
    R_ref_actual = float(disk_df.loc[idx_ref, "R_kpc"]) * KPC
    v_visible_ref = float(disk_df.loc[idx_ref, "v_linear_visible_km_s"]) * 1000
    v_env_needed_sq = max(0.0, V_EDGE_TARGET**2 - v_visible_ref**2)
    eps_fit = 2.0 * v_env_needed_sq * (R_ref_actual + R_C_ENV) / (C**2 * R_ref_actual)

    coll_rows = []
    for _, row in disk_df.iterrows():
        R = float(row["R_kpc"]) * KPC
        g_lin = (float(row["v_linear_visible_km_s"]) * 1000)**2 / R
        g_extra = g_env(R, eps_fit)
        v_total = math.sqrt(max(R*(g_lin + g_extra), 0.0))
        coll_rows.append(
            {
                "R_kpc": row["R_kpc"],
                "epsilon_collective_fit": eps_fit,
                "A_linear_visible": row["A_linear_visible_disk_bulge_gas"],
                "A_collective_envelope": A_env(R, eps_fit),
                "A_total_linear_plus_collective": row["A_linear_visible_disk_bulge_gas"] + A_env(R, eps_fit),
                "v_linear_visible_km_s": row["v_linear_visible_km_s"],
                "v_total_with_collective_km_s": v_total / 1000,
                "v_target_edge_km_s": V_EDGE_TARGET / 1000,
            }
        )
    coll_df = pd.DataFrame(coll_rows)
    coll_df.to_csv(out_dir / "collective_envelope_toy.csv", index=False)

    outer = disk_df[(disk_df["R_kpc"] >= 15) & (disk_df["R_kpc"] <= 35)]
    outer_coll = coll_df[(coll_df["R_kpc"] >= 15) & (coll_df["R_kpc"] <= 35)]

    def slope(x, y):
        return float(np.polyfit(x, y, 1)[0])

    summary = {
        "test": "STAM Model-A galactic cumulative A scratch test",
        "linear_formula": "A_total(x)=sum_i 2Gm_i/(c^2 |x-x_i|)",
        "A_single_sun_at_10kpc": A_point(M_SUN, 10*KPC),
        "A_1e11_suns_at_10kpc": A_point(1e11*M_SUN, 10*KPC),
        "compact_v_for_1e11_suns_at_10kpc_km_s": v_circ_compact_from_A(A_point(1e11*M_SUN, 10*KPC))/1000,
        "linear_visible_outer_median_v_km_s": float(outer["v_linear_visible_km_s"].median()),
        "linear_visible_outer_slope_km_s_per_kpc": slope(outer["R_kpc"], outer["v_linear_visible_km_s"]),
        "collective_envelope_epsilon_fit": eps_fit,
        "collective_outer_median_v_km_s": float(outer_coll["v_total_with_collective_km_s"].median()),
        "collective_outer_slope_km_s_per_kpc": slope(outer_coll["R_kpc"], outer_coll["v_total_with_collective_km_s"]),
        "claim_boundary": "Linear cumulative A is Model-A-native; collective envelope is exploratory only.",
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    print("PASS: galactic cumulative A scratch test complete.")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote outputs to {out_dir}")


if __name__ == "__main__":
    main()

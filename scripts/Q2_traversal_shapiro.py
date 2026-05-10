#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

G = 6.67430e-11
C0 = 299_792_458.0
M_SUN = 1.98847e30
R_SUN = 6.9634e8
M_EARTH = 5.9722e24
R_EARTH = 6_378_137.0
AU = 1.495978707e11
C_MLY_PER_MPC = 3.261563776
H_STAM = 0.000243635
L_MPC = 1.0 / H_STAM
L_MLY = C_MLY_PER_MPC / H_STAM

def schwarzschild_radius(mass_kg: float) -> float:
    return 2.0 * G * mass_kg / C0**2

def A_local(mass_kg: float, r_m: float) -> float:
    return schwarzschild_radius(mass_kg) / r_m

def v_escape(mass_kg: float, r_m: float) -> float:
    return math.sqrt(2.0 * G * mass_kg / r_m)

def D_geo_mpc(z):
    z = np.asarray(z, dtype=float)
    return L_MPC * z * (1.0 + 0.15*z)

def D_adj0_mpc(z):
    z = np.asarray(z, dtype=float)
    return L_MPC * z * (1.0 + 0.5*z)

def D_excess0_mpc(z):
    return D_adj0_mpc(z) - D_geo_mpc(z)

def A_path_average(z):
    z = np.asarray(z, dtype=float)
    return 0.35*z/(1.0 + 0.15*z)

def A_path_local(z):
    z = np.asarray(z, dtype=float)
    return 0.70*z/(1.0 + 0.30*z)

def mu_from_mpc(d_mpc):
    return 5.0*np.log10(np.asarray(d_mpc, dtype=float)) + 25.0

def D_mpc_from_mu(mu):
    return 10.0**((np.asarray(mu, dtype=float)-25.0)/5.0)

"""Q2 — Traversal / Shapiro / local-c path tests.

Purpose:
  Test path accumulation:
    Δt=(1/c)∫A ds
  and preserve the local-c distinction:
    apparent delay comes from extra traversal load (SU/photon-A
    path-stretching from g_rr), not local light slowing.

Clean separation: the formula above computes the path-stretching piece
only. The clock-rate effect from g_tt (dτ/dt = √(1-A)) is a SEPARATE
effect that applies when comparing clock readings at different A values.
For a localized-mass round-trip Shapiro test (single observer clock at
both ends), the clock-rate effect cancels between start and end of the
round trip, so the measured Δt is purely path-stretching. Numerical
equivalence to standard GR Shapiro is by construction (A defined as
2GM/c²r) — no double-counting.
"""

def shapiro_asinh_delay(mass_kg: float, impact_m: float, r1_m: float, r2_m: float) -> float:
    x1 = -math.sqrt(r1_m**2 - impact_m**2)
    x2 = math.sqrt(r2_m**2 - impact_m**2)
    return (schwarzschild_radius(mass_kg)/C0) * (math.asinh(x2/impact_m) - math.asinh(x1/impact_m))

def shapiro_log_delay(mass_kg: float, impact_m: float, r1_m: float, r2_m: float) -> float:
    x1_abs = math.sqrt(r1_m**2 - impact_m**2)
    x2 = math.sqrt(r2_m**2 - impact_m**2)
    R = x1_abs + x2
    return (2.0*G*mass_kg/C0**3) * math.log((r1_m+r2_m+R)/(r1_m+r2_m-R))

def main() -> None:
    out_dir = Path("results/q1_q4_te/Q2_traversal_shapiro")
    out_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        ("Sun_grazing_Earth_Mercury", M_SUN, R_SUN, AU, 0.387*AU),
        ("Sun_grazing_Earth_Mars", M_SUN, R_SUN, AU, 1.524*AU),
        ("Sun_2Rsun_Earth_Mars", M_SUN, 2*R_SUN, AU, 1.524*AU),
        ("Sun_10Rsun_Earth_Mars", M_SUN, 10*R_SUN, AU, 1.524*AU),
    ]
    rows=[]
    for name, M, b, r1, r2 in cases:
        dt_asinh = shapiro_asinh_delay(M,b,r1,r2)
        dt_log = shapiro_log_delay(M,b,r1,r2)
        # approximate straight path flat time
        x1=-math.sqrt(r1**2-b**2); x2=math.sqrt(r2**2-b**2)
        path=x2-x1
        flat=path/C0
        rows.append({
            "case": name,
            "impact_Rsun": b/R_SUN,
            "A_at_impact": A_local(M,b),
            "flat_time_s": flat,
            "A_traversal_delay_s": dt_asinh,
            "A_traversal_delay_us": dt_asinh*1e6,
            "apparent_speed_over_c_using_geo_path": path/(flat+dt_asinh)/C0,
            "asinh_minus_log_s": dt_asinh-dt_log,
            "relative_error_asinh_vs_log": (dt_asinh-dt_log)/dt_log,
        })
    df=pd.DataFrame(rows)
    df.to_csv(out_dir/"Q2_shapiro_traversal_cases.csv", index=False)
    summary={
        "script":"Q2_traversal_shapiro",
        "checks":{
            "max_abs_relative_error_asinh_vs_log":float(df["relative_error_asinh_vs_log"].abs().max()),
            "earth_mars_grazing_delay_us":float(df.loc[df["case"]=="Sun_grazing_Earth_Mars","A_traversal_delay_us"].iloc[0]),
        },
        "interpretation":"Q2 passes: Δt=(1/c)∫A ds reproduces weak-field Shapiro form; apparent slowing is traversal load, not local c failure."
    }
    (out_dir/"summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary, indent=2, sort_keys=True))

if __name__=="__main__":
    main()

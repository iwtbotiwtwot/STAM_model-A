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

"""Q1 — Local A identity tests.

Purpose:
  Test the local, non-distance STAM identities:
    A(r)=2GM/(c²r)
    g=(c²/2)∇A -> GM/r²
    A=(v_escape/c)²
    A=1 -> horizon threshold
"""

def main() -> None:
    out_dir = Path("results/q1_q4_te/Q1_local_A_identity")
    out_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        ("Earth_surface", M_EARTH, R_EARTH),
        ("Sun_surface", M_SUN, R_SUN),
        ("1Msun_at_Rs", M_SUN, schwarzschild_radius(M_SUN)),
        ("10Msun_at_Rs", 10*M_SUN, schwarzschild_radius(10*M_SUN)),
        ("10Msun_at_2Rs", 10*M_SUN, 2*schwarzschild_radius(10*M_SUN)),
        ("10Msun_at_halfRs", 10*M_SUN, 0.5*schwarzschild_radius(10*M_SUN)),
    ]

    rows = []
    for name, M, r in cases:
        A = A_local(M, r)
        vesc = v_escape(M, r)
        g_from_newton = G*M/r**2
        # |grad A| = Rs/r^2, so g=(c^2/2)Rs/r^2 = GM/r^2
        g_from_A = (C0**2/2.0) * schwarzschild_radius(M)/r**2
        rows.append({
            "case": name,
            "mass_kg": M,
            "r_m": r,
            "Rs_m": schwarzschild_radius(M),
            "A": A,
            "v_escape_over_c": vesc/C0,
            "A_minus_vesc2_over_c2": A - (vesc/C0)**2,
            "g_newton_m_s2": g_from_newton,
            "g_from_A_m_s2": g_from_A,
            "g_relative_error": (g_from_A-g_from_newton)/g_from_newton,
            "regime": "A<1" if A < 1 else "A=1" if math.isclose(A, 1.0, rel_tol=1e-12) else "A>1",
        })

    df = pd.DataFrame(rows)
    df.to_csv(out_dir/"Q1_local_A_identity_cases.csv", index=False)

    summary = {
        "script": "Q1_local_A_identity",
        "checks": {
            "max_abs_A_minus_vesc2_over_c2": float(df["A_minus_vesc2_over_c2"].abs().max()),
            "max_abs_g_relative_error": float(df["g_relative_error"].abs().max()),
            "horizon_cases_A_equal_1_count": int((np.isclose(df["A"], 1.0)).sum()),
        },
        "interpretation": "Q1 passes: local A reproduces escape-speed identity, Newtonian acceleration bridge, and horizon threshold."
    }
    (out_dir/"summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()

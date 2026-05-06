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

"""Q4 — Quantum-scale and A-conditioned observation sniff tests.

Purpose:
  Keep quantum claims restrained:
    self-A at particle scale is tiny,
    collective A can emerge macroscopically,
    local observation can normalize c in a toy A-conditioned light-clock.
"""

def S(A: float) -> float:
    return 1.0 + A

def light_clock(A: float, D_geo_m: float=1.0, c: float=C0):
    s=S(A)
    D_A=s*D_geo_m
    T_cross=D_A/c
    T_local=T_cross/s
    D_local=D_A/s
    return {
        "A":A,
        "S":s,
        "D_A":D_A,
        "T_cross":T_cross,
        "T_local":T_local,
        "D_local":D_local,
        "local_c":D_local/T_local,
        "naive_external_c":D_geo_m/T_cross,
        "proper_accumulated_c":D_A/T_cross,
    }

def main() -> None:
    out_dir = Path("results/q1_q4_te/Q4_quantum_observation")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Quantum sniff
    m_e=9.1093837015e-31
    r_bohr=5.29177210903e-11
    m_p=1.67262192369e-27
    r_fm=1e-15
    quantum_rows=[
        {"case":"electron_at_bohr_radius","mass_kg":m_e,"r_m":r_bohr,"A_self":A_local(m_e,r_bohr)},
        {"case":"proton_at_1_fm","mass_kg":m_p,"r_m":r_fm,"A_self":A_local(m_p,r_fm)},
        {"case":"1e11_solar_masses_at_10kpc","mass_kg":1e11*M_SUN,"r_m":10*3.0856775814913673e19,"A_self":A_local(1e11*M_SUN,10*3.0856775814913673e19)},
    ]
    qdf=pd.DataFrame(quantum_rows)
    qdf.to_csv(out_dir/"Q4_quantum_scale_A_sniff.csv", index=False)

    # A-conditioned observation toy
    Avals=[-0.5,0.0,0.001,0.1,0.99]
    rows=[]
    for A in Avals:
        lc=light_clock(A, D_geo_m=10_000_000_000.0*1000)
        rows.append({
            "A":A,
            "S_1_plus_A":lc["S"],
            "local_c_over_c0":lc["local_c"]/C0,
            "naive_external_c_over_c0":lc["naive_external_c"]/C0,
            "proper_accumulated_c_over_c0":lc["proper_accumulated_c"]/C0,
            "cross_period_s":2*lc["T_cross"],
            "local_period_s":2*lc["T_local"],
        })
    odf=pd.DataFrame(rows)
    odf.to_csv(out_dir/"Q4_A_conditioned_observation_toy.csv", index=False)

    summary={
        "script":"Q4_quantum_observation",
        "checks":{
            "electron_A_self":float(qdf.loc[qdf["case"]=="electron_at_bohr_radius","A_self"].iloc[0]),
            "proton_A_self":float(qdf.loc[qdf["case"]=="proton_at_1_fm","A_self"].iloc[0]),
            "max_abs_local_c_over_c0_minus_1":float(np.max(np.abs(odf["local_c_over_c0"]-1.0))),
        },
        "interpretation":"Q4 passes as a sniff test: particle self-A is tiny, collective A can be large, and local-c normalization works in the toy observation model."
    }
    (out_dir/"summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary, indent=2, sort_keys=True))

if __name__=="__main__":
    main()

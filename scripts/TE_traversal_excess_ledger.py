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

"""TE ledger — Traversal Excess ledger.

Purpose:
  Create the missing TE ledger:
    TE(z)=D_excess,0(z)=D_adj,0-D_geo=0.35Lz²
  plus A_path average/local and high-z limit behavior.
"""

def main() -> None:
    out_dir = Path("results/q1_q4_te/TE_traversal_excess_ledger")
    out_dir.mkdir(parents=True, exist_ok=True)

    zs=np.array([0.001,0.01,0.05,0.1,0.3,0.5,1.0,2.0,3.0,5.0,10.0,20.0,100.0,1100.0])
    df=pd.DataFrame({
        "z":zs,
        "D_geo_Mpc":D_geo_mpc(zs),
        "D_adj0_Mpc":D_adj0_mpc(zs),
        "TE_D_excess0_Mpc":D_excess0_mpc(zs),
        "TE_D_excess0_Mly":D_excess0_mpc(zs)*C_MLY_PER_MPC,
        "TE_fraction_of_D_geo_A_path_avg":A_path_average(zs),
        "TE_local_differential_A_path":A_path_local(zs),
        "D_adj0_over_D_geo":D_adj0_mpc(zs)/D_geo_mpc(zs),
    })
    df["A_path_avg_fraction_of_limit_7over3"]=df["TE_fraction_of_D_geo_A_path_avg"]/(7/3)
    df["A_path_local_fraction_of_limit_7over3"]=df["TE_local_differential_A_path"]/(7/3)
    df.to_csv(out_dir/"TE_ledger_key_redshifts.csv", index=False)

    # Saturation table
    fracs=np.array([0.5,0.75,0.9,0.95,0.99,0.999])
    sat=pd.DataFrame({
        "fraction_of_limit":fracs,
        "z_for_A_path_avg":fracs/(0.15*(1-fracs)),
        "z_for_A_path_local":fracs/(0.30*(1-fracs)),
        "limit":7/3,
    })
    sat.to_csv(out_dir/"TE_A_path_saturation_ledger.csv", index=False)

    # Identity checks
    ztest=np.linspace(0.001,10,1000)
    te_direct=D_adj0_mpc(ztest)-D_geo_mpc(ztest)
    te_formula=0.35*L_MPC*ztest**2
    summary={
        "script":"TE_traversal_excess_ledger",
        "definition":"TE(z)=D_excess,0(z)=D_adj,0-D_geo=0.35Lz^2",
        "checks":{
            "max_abs_TE_identity_error_Mpc":float(np.max(np.abs(te_direct-te_formula))),
            "A_path_average_limit":7/3,
            "A_path_local_limit":7/3,
            "D_adj0_over_D_geo_limit":10/3,
        },
        "interpretation":"TE ledger built: traversal excess is the no-b distance excess layer used for A_path."
    }
    (out_dir/"summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary, indent=2, sort_keys=True))

if __name__=="__main__":
    main()

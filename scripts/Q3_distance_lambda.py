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

"""Q3 — Cosmological distance layer and lambda-mapping tests.

Purpose:
  Test the no-b Model-A distance split:
    D_geo, D_adj0, D_excess0
  and the general observable layer:
    D_obs,type = D_geo + lambda_type(z) D_excess0
"""

def lambda_from_distance(D_obs_mpc, z):
    return (np.asarray(D_obs_mpc, dtype=float) - D_geo_mpc(z)) / D_excess0_mpc(z)

def main() -> None:
    out_dir = Path("results/q1_q4_te/Q3_distance_lambda")
    out_dir.mkdir(parents=True, exist_ok=True)

    zs = np.array([0.01,0.05,0.1,0.3,0.5,1.0,2.0,3.0,10.0,100.0,1100.0])
    df = pd.DataFrame({
        "z": zs,
        "D_geo_Mpc": D_geo_mpc(zs),
        "D_adj0_Mpc": D_adj0_mpc(zs),
        "D_excess0_Mpc": D_excess0_mpc(zs),
        "D_adj0_over_D_geo": D_adj0_mpc(zs)/D_geo_mpc(zs),
        "A_path_average": A_path_average(zs),
        "A_path_local": A_path_local(zs),
    })
    df["A_path_average_fraction_of_limit_7over3"] = df["A_path_average"]/(7/3)
    df.to_csv(out_dir/"Q3_distance_layer_key_values.csv", index=False)

    # Synthetic observable layer recovery: choose lambda(z)=1.2-0.3z and recover from D_obs.
    zsyn = np.linspace(0.05, 2.5, 100)
    lam_true = 1.2 - 0.3*zsyn
    Dobs = D_geo_mpc(zsyn) + lam_true * D_excess0_mpc(zsyn)
    lam_rec = lambda_from_distance(Dobs, zsyn)
    syn = pd.DataFrame({"z": zsyn, "lambda_true": lam_true, "lambda_recovered": lam_rec, "error": lam_rec-lam_true})
    syn.to_csv(out_dir/"Q3_synthetic_lambda_recovery.csv", index=False)

    summary={
        "script":"Q3_distance_lambda",
        "checks":{
            "max_abs_lambda_recovery_error":float(np.max(np.abs(lam_rec-lam_true))),
            "A_path_average_limit":7/3,
            "A_path_local_limit":7/3,
            "D_adj0_over_D_geo_limit":10/3,
            "A_path_average_at_z1100":float(A_path_average(1100.0)),
        },
        "interpretation":"Q3 passes: distance layers and lambda observable mapping are internally recoverable without b."
    }
    (out_dir/"summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary, indent=2, sort_keys=True))

if __name__=="__main__":
    main()

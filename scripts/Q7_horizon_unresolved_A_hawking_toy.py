#!/usr/bin/env python3
"""
Q7_horizon_unresolved_A_hawking_toy.py

Accumulation Theory / STAM Model-A
Q7 — horizon-unresolved-A / Hawking-like decay toy

Question:
    Can we make a toy ledger where horizon-locked information becomes unavailable
    to ordinary relational exchange, producing an unresolved-A boundary channel
    that behaves like Hawking-like radiation and causes black-hole mass decay?

Important:
    This is NOT a derivation of Hawking radiation.
    This is NOT quantum field theory in curved spacetime.
    This is a conceptual ledger toy.

Core STAM ingredients:
    A = Rs/r = 2GM/(c^2 r)
    Horizon threshold: A = 1
    r_h = 2GM/c^2

Mass ledger:
    dM/dt = Mdot_in - P_out/c^2

Toy unresolved-horizon channel:
    Let U_h be "unresolved boundary availability":
        U_h = 0 means horizon information is fully active/resolved in the universe ledger.
        U_h = 1 means horizon information is preserved but unavailable to ordinary relational exchange.

    Let P_out scale with:
        P_out ∝ U_h * horizon_factor(M)

Known Hawking-like scaling:
    P_Hawking ∝ 1/M^2

So the toy uses:
    P_out = U_h * K / M^2

Then:
    dM/dt = Mdot_in - U_h*K/M^2

The test:
    Compare mass/horizon evolution for:
        U_h = 0      no unresolved boundary radiation
        U_h = 0.25   weak unresolved channel
        U_h = 0.75   strong unresolved channel
        U_h = 1.0    fully active unresolved boundary channel

What would support the concept internally:
    U_h=0 produces no decay unless inflow/outflow terms exist.
    Higher U_h produces stronger mass loss.
    As M decreases, P_out rises ~1/M^2, making late decay accelerate.
    Horizon radius r_h shrinks with M.

What this proves:
    Only that the proposed ledger is internally coherent as a toy.
"""

from pathlib import Path
import argparse, json, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Use normalized units:
# c=1, G=1/2 gives r_h = M for convenience if desired.
# But we'll keep r_h normalized proportional to M.
def simulate(M0=10.0, U_h=1.0, K=0.05, Mdot_in=0.0, tmax=8000.0, dt=1.0, M_stop=0.5):
    rows = []
    M = float(M0)
    t = 0.0
    step = 0
    while t <= tmax and M > M_stop:
        P_out = U_h * K / (M*M)
        dMdt = Mdot_in - P_out
        r_h = M  # normalized horizon radius
        A_at_fixed_r = M / M0  # A at initial horizon radius r=M0, normalized
        # boundary unresolved "pressure" diagnostic, not physical pressure
        unresolved_flux = P_out

        rows.append({
            "step": step,
            "t": t,
            "M": M,
            "U_h": U_h,
            "P_out": P_out,
            "Mdot_in": Mdot_in,
            "dMdt": dMdt,
            "r_h_norm": r_h,
            "A_at_initial_horizon_radius": A_at_fixed_r,
            "unresolved_flux": unresolved_flux,
        })

        M = max(M + dMdt * dt, 0.0)
        t += dt
        step += 1

    return pd.DataFrame(rows)

def run(outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    scenarios = [
        {"name": "resolved_horizon_no_radiation", "U_h": 0.0, "Mdot_in": 0.0},
        {"name": "weak_unresolved_boundary", "U_h": 0.25, "Mdot_in": 0.0},
        {"name": "strong_unresolved_boundary", "U_h": 0.75, "Mdot_in": 0.0},
        {"name": "fully_unresolved_boundary", "U_h": 1.0, "Mdot_in": 0.0},
        {"name": "fully_unresolved_with_inflow", "U_h": 1.0, "Mdot_in": 0.0002},
    ]

    summary_rows = []
    plt.figure(figsize=(10, 6))

    for sc in scenarios:
        df = simulate(U_h=sc["U_h"], Mdot_in=sc["Mdot_in"])
        df.to_csv(outdir / f"Q7_{sc['name']}.csv", index=False)

        initial = df.iloc[0]
        final = df.iloc[-1]
        half_mass = df[df["M"] <= initial["M"]/2]
        t_half = float(half_mass["t"].iloc[0]) if len(half_mass) else np.nan

        summary_rows.append({
            "scenario": sc["name"],
            "U_h": sc["U_h"],
            "Mdot_in": sc["Mdot_in"],
            "initial_M": float(initial["M"]),
            "final_M": float(final["M"]),
            "mass_lost": float(initial["M"] - final["M"]),
            "initial_P_out": float(initial["P_out"]),
            "final_P_out": float(final["P_out"]),
            "radiation_power_increase_factor": float(final["P_out"]/initial["P_out"]) if initial["P_out"] > 0 else np.nan,
            "final_r_h_norm": float(final["r_h_norm"]),
            "t_half_mass": t_half,
            "ended_at_t": float(final["t"]),
        })

        plt.plot(df["t"], df["M"], label=sc["name"])

        # Individual mass/power plot
        fig = plt.figure(figsize=(10,5))
        ax1 = fig.add_subplot(111)
        ax1.plot(df["t"], df["M"], label="M")
        ax1.plot(df["t"], df["r_h_norm"], "--", label="r_h normalized")
        ax1.set_xlabel("toy time")
        ax1.set_ylabel("M / r_h")
        ax2 = ax1.twinx()
        ax2.plot(df["t"], df["P_out"], label="P_out", alpha=0.6)
        ax2.set_ylabel("P_out")
        ax1.set_title(f"Q7: {sc['name']}")
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1+lines2, labels1+labels2, loc="best")
        fig.tight_layout()
        fig.savefig(outdir / f"Q7_{sc['name']}_mass_power.png", dpi=180)
        plt.close(fig)

    plt.title("Q7 horizon-unresolved-A toy: stronger unresolved boundary gives stronger decay")
    plt.xlabel("toy time")
    plt.ylabel("black-hole mass M")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "Q7_mass_evolution_overlay.png", dpi=180)
    plt.close()

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(outdir / "Q7_summary_metrics.csv", index=False)

    # Scaling check: P_out vs M should follow slope -2 in log-log for fixed U_h.
    df_full = pd.read_csv(outdir / "Q7_fully_unresolved_boundary.csv")
    valid = df_full[(df_full["P_out"] > 0) & (df_full["M"] > 0)]
    slope, intercept = np.polyfit(np.log(valid["M"]), np.log(valid["P_out"]), 1)

    plt.figure(figsize=(7,5))
    plt.scatter(np.log(valid["M"]), np.log(valid["P_out"]), s=4, alpha=0.5)
    plt.plot(np.log(valid["M"]), slope*np.log(valid["M"])+intercept, label=f"slope={slope:.3f}")
    plt.title("Q7 scaling check: P_out vs M")
    plt.xlabel("log(M)")
    plt.ylabel("log(P_out)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "Q7_power_scaling_loglog.png", dpi=180)
    plt.close()

    verdict = {
        "script": "Q7_horizon_unresolved_A_hawking_toy.py",
        "claim_tested": "A horizon-unresolved-A boundary channel can be represented as Hawking-like mass loss in a toy mass ledger.",
        "core_rule": "P_out = U_h*K/M^2",
        "mass_ledger": "dM/dt = Mdot_in - P_out",
        "scaling_check_loglog_slope": float(slope),
        "main_result": [
            "U_h=0 produces no radiation/decay in this toy.",
            "Higher U_h produces stronger mass loss.",
            "P_out increases as M decreases, matching the intended 1/M^2 Hawking-like scaling.",
            "The horizon radius shrinks with mass."
        ],
        "interpretation": "If horizon information is preserved but unavailable to ordinary relational exchange, the boundary can be modeled as an unresolved-A interface with an outgoing radiation channel.",
        "caution": [
            "This is not a Hawking derivation.",
            "This does not prove information becomes universe-unavailable.",
            "This only tests whether the proposed Accumulation Theory ledger is internally coherent as a toy."
        ],
        "barstool_read": "The black hole does not delete the receipt. It locks the receipt at the boundary, and the weird locked boundary pays out tiny radiation while the mass ledger shrinks."
    }
    with open(outdir / "Q7_verdict.json", "w") as f:
        json.dump(verdict, f, indent=2)

    print("Q7 horizon-unresolved-A / Hawking-like decay toy complete.")
    print(f"Output directory: {outdir}")
    print()
    print(summary.to_string(index=False))
    print()
    print(f"Scaling check log(P_out) vs log(M): slope ≈ {slope:.3f} (target -2)")
    print()
    print("Main result:")
    print("  U_h=0 gives no decay in this toy.")
    print("  Higher unresolved boundary availability gives stronger mass loss.")
    print("  P_out rises as M falls, and r_h shrinks with M.")
    print()
    print("Limit:")
    print("  This is a conceptual ledger toy, not a derivation of Hawking radiation.")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="results/Q7_horizon_unresolved_A_hawking_toy")
    args = p.parse_args()
    run(args.outdir)

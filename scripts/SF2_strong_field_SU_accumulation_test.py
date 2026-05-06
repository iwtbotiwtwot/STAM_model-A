#!/usr/bin/env python3
"""
SF2_strong_field_SU_accumulation_test.py

Accumulation Theory / STAM Model-A
SF2 — Strong-field SU accumulation test

Purpose:
    Treat SU as a native unit of accumulation, not as distance itself.

Weak-field STAM path rule:
    SU_weak ∝ ∫ A ds

Strong-field A-metric radial traversal:
    ds^2 = -(1-A)c^2dt^2 + dr^2/(1-A) + r^2dΩ^2

For radial null travel:
    dr/dt = c(1-A)

Flat-space radial time:
    dt_flat = dr/c

A-metric radial coordinate time:
    dt_metric = dr/[c(1-A)]

Extra traversal load:
    dt_metric - dt_flat = [1/(1-A) - 1] dr/c
                          = [A/(1-A)] dr/c

Therefore define a strong-field accumulation load:
    A_load,strong = A/(1-A)

Weak limit:
    A << 1:
        A/(1-A) ≈ A

Horizon behavior:
    A -> 1:
        A/(1-A) -> infinity

In normalized units Rs=1:
    A = 1/r

    SU_weak(r1 -> r2) = ∫ A dr/Rs
                      = ln(r2/r1)

    SU_strong(r1 -> r2) = ∫ A/(1-A) dr/Rs
                        = ln((r2-Rs)/(r1-Rs))

This script compares weak SU and strong SU across field strength.
"""

from pathlib import Path
import argparse, json, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RS = 1.0

def A_of_r(r):
    return RS / np.asarray(r, dtype=float)

def A_load_weak(A):
    return np.asarray(A, dtype=float)

def A_load_strong(A):
    A = np.asarray(A, dtype=float)
    return A / (1.0 - A)

def SU_weak_radial(r1, r2):
    return math.log(r2/r1)

def SU_strong_radial(r1, r2):
    return math.log((r2-RS)/(r1-RS))

def run(outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # A-grid up to near horizon
    A_grid = np.linspace(0.001, 0.999, 2000)
    r_grid = RS / A_grid

    grid = pd.DataFrame({
        "A": A_grid,
        "r_over_Rs": r_grid,
        "A_load_weak": A_load_weak(A_grid),
        "A_load_strong": A_load_strong(A_grid),
        "strong_over_weak": A_load_strong(A_grid)/A_load_weak(A_grid),
    })
    grid.to_csv(outdir / "SF2_A_load_grid.csv", index=False)

    # Radial intervals from r1 to r2=100Rs
    r2 = 100.0
    r1_values = [50.0, 20.0, 10.0, 5.0, 3.0, 2.0, 1.5, 1.2, 1.1, 1.05, 1.01, 1.001]
    rows = []
    for r1 in r1_values:
        suw = SU_weak_radial(r1, r2)
        sus = SU_strong_radial(r1, r2)
        A1 = 1.0/r1
        rows.append({
            "r1_over_Rs": r1,
            "r2_over_Rs": r2,
            "A_at_r1": A1,
            "SU_weak_ln_r2_over_r1": suw,
            "SU_strong_ln_r2minus1_over_r1minus1": sus,
            "SU_strong_minus_weak": sus - suw,
            "SU_strong_over_weak": sus/suw if suw != 0 else np.nan,
            "A_load_at_r1_weak": A1,
            "A_load_at_r1_strong": A1/(1-A1),
            "load_strong_over_weak_at_r1": 1/(1-A1),
        })
    intervals = pd.DataFrame(rows)
    intervals.to_csv(outdir / "SF2_radial_SU_intervals.csv", index=False)

    # Fixed outer radius, sweep inner radius continuously.
    r1_sweep = np.concatenate([
        np.linspace(2.0, 50.0, 400),
        np.linspace(1.05, 2.0, 400),
        np.linspace(1.001, 1.05, 400),
    ])
    r1_sweep = np.unique(np.sort(r1_sweep))
    sweep_rows = []
    for r1 in r1_sweep:
        suw = SU_weak_radial(r1, r2)
        sus = SU_strong_radial(r1, r2)
        sweep_rows.append({
            "r1_over_Rs": r1,
            "A_at_r1": 1/r1,
            "SU_weak": suw,
            "SU_strong": sus,
            "SU_excess_strong_minus_weak": sus-suw,
            "strong_over_weak": sus/suw if suw != 0 else np.nan,
        })
    sweep = pd.DataFrame(sweep_rows)
    sweep.to_csv(outdir / "SF2_radial_SU_sweep.csv", index=False)

    # Plots
    plt.figure(figsize=(8,5))
    plt.plot(grid["A"], grid["A_load_weak"], label="weak load = A")
    plt.plot(grid["A"], np.minimum(grid["A_load_strong"], 20), label="strong load = A/(1-A), clipped at 20")
    plt.axvline(1/3, linestyle=":", label="A=1/3 ISCO")
    plt.axvline(2/3, linestyle=":", label="A=2/3 photon sphere")
    plt.axvline(1, linestyle=":", label="A=1 horizon")
    plt.title("SF2: SU accumulation load per path element")
    plt.xlabel("A")
    plt.ylabel("accumulation load")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "SF2_A_load_weak_vs_strong.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8,5))
    plt.plot(sweep["A_at_r1"], sweep["SU_weak"], label="SU_weak = ∫A dr/Rs")
    plt.plot(sweep["A_at_r1"], sweep["SU_strong"], label="SU_strong = ∫A/(1-A) dr/Rs")
    plt.axvline(1/3, linestyle=":", label="A=1/3")
    plt.axvline(2/3, linestyle=":", label="A=2/3")
    plt.axvline(1, linestyle=":", label="A=1")
    plt.title("SF2: radial SU from r1 to 100Rs")
    plt.xlabel("A at inner radius r1")
    plt.ylabel("SU count, normalized")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "SF2_radial_SU_weak_vs_strong.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8,5))
    plt.plot(sweep["A_at_r1"], sweep["strong_over_weak"])
    plt.axvline(1/3, linestyle=":", label="A=1/3")
    plt.axvline(2/3, linestyle=":", label="A=2/3")
    plt.axvline(1, linestyle=":", label="A=1")
    plt.title("SF2: strong SU / weak SU")
    plt.xlabel("A at inner radius r1")
    plt.ylabel("SU_strong / SU_weak")
    plt.ylim(0, min(8, np.nanmax(sweep["strong_over_weak"])))
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "SF2_strong_over_weak_SU_ratio.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8,5))
    plt.semilogy(sweep["A_at_r1"], sweep["SU_strong"])
    plt.axvline(1/3, linestyle=":", label="A=1/3")
    plt.axvline(2/3, linestyle=":", label="A=2/3")
    plt.axvline(1, linestyle=":", label="A=1")
    plt.title("SF2: strong SU divergence near horizon")
    plt.xlabel("A at inner radius r1")
    plt.ylabel("SU_strong, log scale")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "SF2_SU_strong_divergence_log.png", dpi=180)
    plt.close()

    # Check weak-limit agreement at small A.
    weak_checks = []
    for A in [0.001, 0.01, 0.05, 0.1, 1/3, 2/3, 0.9, 0.99]:
        weak = A
        strong = A/(1-A)
        weak_checks.append({
            "A": A,
            "A_load_weak": weak,
            "A_load_strong": strong,
            "strong_minus_weak": strong-weak,
            "strong_over_weak": strong/weak if weak else np.nan,
        })
    weak_df = pd.DataFrame(weak_checks)
    weak_df.to_csv(outdir / "SF2_load_check_values.csv", index=False)

    verdict = {
        "script": "SF2_strong_field_SU_accumulation_test.py",
        "purpose": "Treat SU as a native unit of accumulation and extend it into strong-field A-metric traversal.",
        "weak_SU": "SU_weak ∝ ∫A ds",
        "strong_SU": "SU_strong ∝ ∫A/(1-A) ds",
        "normalized_radial_results": {
            "SU_weak": "ln(r2/r1)",
            "SU_strong": "ln((r2-Rs)/(r1-Rs))"
        },
        "main_result": [
            "Strong SU reduces to weak SU when A is small.",
            "Strong SU grows faster than weak SU as A rises.",
            "Strong SU diverges as A approaches 1.",
            "This gives SU a strong-field meaning: count of accumulated traversal load, not distance itself."
        ],
        "interpretation": "In weak/mid-field cosmology, spreadsheet SU can be treated as accumulated A-units. Near horizons, the A-metric suggests replacing A with A/(1-A) for traversal load.",
        "caution": "This is a strong-field extension of the SU concept using the Schwarzschild A-metric baseline. It is not yet a new metric modification."
    }
    with open(outdir / "SF2_verdict.json", "w") as f:
        json.dump(verdict, f, indent=2)

    print("SF2 strong-field SU accumulation test complete.")
    print(f"Output directory: {outdir}")
    print()
    print("Load check values:")
    print(weak_df.to_string(index=False))
    print()
    print("Radial SU intervals:")
    print(intervals.to_string(index=False))
    print()
    print("Main result:")
    print("  Weak SU:   ∫A dr/Rs = ln(r2/r1)")
    print("  Strong SU: ∫A/(1-A) dr/Rs = ln((r2-Rs)/(r1-Rs))")
    print("  Strong SU agrees in weak field and diverges near A=1.")
    print("  This supports SU as accumulated traversal load rather than distance itself.")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="results/SF2_strong_field_SU_accumulation")
    args = p.parse_args()
    run(args.outdir)

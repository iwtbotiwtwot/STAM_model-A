#!/usr/bin/env python3
"""
G96_BNS_engine_time_scaling.py

STAM / Model-A BNS engine-time scaling test.

Purpose
=======
This script tests the STAM substance-velocity-cap / high-A engine-time
prediction against binary-neutron-star events with electromagnetic counterparts.

Interpretation
==============
The observed GW-to-EM delay is NOT treated as a long-distance propagation-speed
difference between gravitational waves and light.

STAM reading:
    GW and light propagate essentially together after release.
    The observed delay is local engine time near the merger, caused by
    high-A substance dynamics near the critical doughnut / shell threshold.

Locked model inputs
===================
A_0 = 1 / (12 pi)

Threshold equation:
    A_0 = x(1 - x) / (1 + x)

where:
    x = R_s / (2 r_threshold)

Then:
    r_threshold / R_s = 1 / (2x)

Engine time:
    tau_engine = (5/8) * (r_threshold/R_s)^4 * (R_s/c)

Since R_s/c = 2GM/c^3, the prediction is linear in total mass:
    tau_engine = slope * (M_total / M_sun)

Outputs
=======
- Prints the derived constants and event comparison table.
- Writes:
    results/G96_BNS_engine_time_scaling_summary.md
    results/G96_BNS_engine_time_scaling_events.csv
    plots/G96_BNS_engine_time_scaling.png

Usage
=====
Run default built-in GW170817 test:

    python G96_BNS_engine_time_scaling.py

Optionally supply a CSV:

    python G96_BNS_engine_time_scaling.py --events path/to/events.csv

CSV columns:
    event,total_mass_msun,observed_delay_s,delay_uncertainty_s,notes

Only event,total_mass_msun,observed_delay_s are required.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import List, Dict, Any

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------

G = 6.67430e-11                  # m^3 kg^-1 s^-2
c = 299_792_458.0                # m/s
M_sun = 1.98847e30               # kg

A0 = 1.0 / (12.0 * math.pi)


def solve_threshold_x(A0_value: float = A0):
    """Solve A0 = x(1-x)/(1+x), return physical root(s)."""
    # A0(1+x) = x - x^2
    # x^2 + (A0 - 1)x + A0 = 0
    a = 1.0
    b = A0_value - 1.0
    ccoef = A0_value
    disc = b*b - 4*a*ccoef
    if disc < 0:
        raise ValueError("No real threshold roots.")
    roots = [(-b - math.sqrt(disc)) / (2*a),
             (-b + math.sqrt(disc)) / (2*a)]
    return roots


def choose_outer_threshold_root(roots):
    """Pick the smaller x root, corresponding to larger r/Rs."""
    roots_pos = [r for r in roots if r > 0]
    if not roots_pos:
        raise ValueError("No positive root.")
    return min(roots_pos)


def engine_slope_seconds_per_msun(A0_value: float = A0):
    roots = solve_threshold_x(A0_value)
    x = choose_outer_threshold_root(roots)
    r_over_Rs = 1.0 / (2.0 * x)
    Rs_over_c_per_msun = 2.0 * G * M_sun / c**3
    slope = (5.0 / 8.0) * (r_over_Rs ** 4) * Rs_over_c_per_msun
    return slope, x, r_over_Rs, roots


def predict_delay(total_mass_msun: float, slope: float) -> float:
    return slope * total_mass_msun


DEFAULT_EVENTS = [
    {
        "event": "GW170817 / GRB 170817A",
        "total_mass_msun": 2.70,
        "observed_delay_s": 1.74,
        "delay_uncertainty_s": "",
        "notes": "Canonical BNS + gamma-ray counterpart; observed delay often quoted as ~1.74 s. Total mass rounded to 2.70 Msun for STAM locked test.",
    }
]


def load_events(path: Path | None):
    if path is None:
        return DEFAULT_EVENTS

    rows = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"event", "total_mass_msun", "observed_delay_s"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"CSV missing required columns: {sorted(missing)}")
        for row in reader:
            rows.append(row)
    return rows


def normalize_event(row: Dict[str, Any]):
    return {
        "event": str(row.get("event", "")).strip(),
        "total_mass_msun": float(row.get("total_mass_msun")),
        "observed_delay_s": float(row.get("observed_delay_s")),
        "delay_uncertainty_s": row.get("delay_uncertainty_s", ""),
        "notes": str(row.get("notes", "")).strip(),
    }


def main():
    parser = argparse.ArgumentParser(description="G96 BNS engine-time scaling test")
    parser.add_argument("--events", type=str, default=None,
                        help="Optional CSV of BNS+EM events.")
    parser.add_argument("--outdir", type=str, default=".",
                        help="Output root directory. Defaults to current directory.")
    args = parser.parse_args()

    root = Path(args.outdir).resolve()
    results_dir = root / "results"
    plots_dir = root / "plots"
    results_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    slope, x, r_over_Rs, roots = engine_slope_seconds_per_msun(A0)

    print("=" * 88)
    print("G96: STAM BNS engine-time scaling test")
    print("=" * 88)
    print()
    print("Interpretation:")
    print("  GW and light propagate essentially together after release.")
    print("  Delay is local engine time near the merger, not cosmic travel-speed slippage.")
    print()
    print(f"A0 = 1/(12π) = {A0:.12f}")
    print("Threshold equation: A0 = x(1-x)/(1+x)")
    print(f"Roots x = {roots[0]:.12f}, {roots[1]:.12f}")
    print(f"Chosen outer-threshold root x = {x:.12f}")
    print(f"r_threshold/Rs = {r_over_Rs:.6f}")
    print(f"Engine-time slope = {slope:.6f} s / M_sun")
    print()

    raw_events = load_events(Path(args.events) if args.events else None)
    events = [normalize_event(r) for r in raw_events]

    rows = []
    for ev in events:
        pred = predict_delay(ev["total_mass_msun"], slope)
        resid = ev["observed_delay_s"] - pred
        pct_err = 100.0 * abs(resid) / ev["observed_delay_s"] if ev["observed_delay_s"] != 0 else float("nan")
        rows.append({
            **ev,
            "predicted_delay_s": pred,
            "residual_obs_minus_pred_s": resid,
            "absolute_percent_error": pct_err,
        })

    print(f"{'Event':<28}{'Mtot(Msun)':>12}{'Obs(s)':>12}{'Pred(s)':>12}{'Obs-Pred':>12}{'Abs % err':>12}")
    print("-" * 88)
    for r in rows:
        print(f"{r['event']:<28}{r['total_mass_msun']:>12.4f}{r['observed_delay_s']:>12.4f}"
              f"{r['predicted_delay_s']:>12.4f}{r['residual_obs_minus_pred_s']:>12.4f}"
              f"{r['absolute_percent_error']:>12.3f}")
    print()

    # CSV output
    csv_path = results_dir / "G96_BNS_engine_time_scaling_events.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "event", "total_mass_msun", "observed_delay_s", "predicted_delay_s",
            "residual_obs_minus_pred_s", "absolute_percent_error",
            "delay_uncertainty_s", "notes"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in fieldnames})

    # Plot
    masses = np.linspace(1.8, 4.0, 300)
    delays = slope * masses

    plt.figure(figsize=(9, 6))
    plt.plot(masses, delays, label=f"STAM prediction: τ = {slope:.3f} s/Msun × M")
    for r in rows:
        plt.scatter([r["total_mass_msun"]], [r["observed_delay_s"]], s=70)
        plt.annotate(r["event"], (r["total_mass_msun"], r["observed_delay_s"]),
                     xytext=(8, 8), textcoords="offset points", fontsize=8)
        # vertical residual
        plt.plot([r["total_mass_msun"], r["total_mass_msun"]],
                 [r["predicted_delay_s"], r["observed_delay_s"]],
                 linestyle="--", linewidth=1)
    plt.xlabel("Total binary mass (M_sun)")
    plt.ylabel("GW-to-EM delay / engine time (seconds)")
    plt.title("G96 STAM BNS engine-time scaling")
    plt.grid(alpha=0.3)
    plt.legend()
    plot_path = plots_dir / "G96_BNS_engine_time_scaling.png"
    plt.tight_layout()
    plt.savefig(plot_path, dpi=200)
    plt.close()

    # Markdown summary
    md_path = results_dir / "G96_BNS_engine_time_scaling_summary.md"
    md = []
    md.append("# G96 — BNS engine-time scaling test\n")
    md.append("## Interpretation\n")
    md.append("The observed GW-to-EM delay is treated as **local engine time near the merger**, not as a long-distance propagation-speed difference. GW and light propagate essentially together after release.\n")
    md.append("## Locked derivation\n")
    md.append("```text\n")
    md.append("A0 = 1/(12π)\n")
    md.append("A0 = x(1-x)/(1+x), where x = Rs/(2r_threshold)\n")
    md.append("r_threshold/Rs = 1/(2x)\n")
    md.append("tau_engine = (5/8)(r_threshold/Rs)^4 (Rs/c)\n")
    md.append("tau_engine = slope × (M_total/M_sun)\n")
    md.append("```\n")
    md.append(f"- A0 = `{A0:.12f}`\n")
    md.append(f"- x roots = `{roots[0]:.12f}`, `{roots[1]:.12f}`\n")
    md.append(f"- chosen root = `{x:.12f}`\n")
    md.append(f"- r_threshold/Rs = `{r_over_Rs:.6f}`\n")
    md.append(f"- slope = `{slope:.6f}` s/M_sun\n")
    md.append("## Event comparison\n")
    md.append("| Event | M_total (M_sun) | Observed delay (s) | STAM predicted (s) | Obs - Pred (s) | Abs % error |\n")
    md.append("|---|---:|---:|---:|---:|---:|\n")
    for r in rows:
        md.append(f"| {r['event']} | {r['total_mass_msun']:.4f} | {r['observed_delay_s']:.4f} | "
                  f"{r['predicted_delay_s']:.4f} | {r['residual_obs_minus_pred_s']:.4f} | "
                  f"{r['absolute_percent_error']:.3f}% |\n")
    md.append("\n## Falsification handle\n")
    md.append("Future BNS + EM counterpart events should fall approximately on the linear mass-scaling relation if the STAM doughnut / substance-engine-time mechanism is correct. If events with reliable EM launch-time interpretation do not follow this scaling, this mechanism is wrong or incomplete.\n")
    md.append("## Files\n")
    md.append(f"- CSV: `{csv_path}`\n")
    md.append(f"- Plot: `{plot_path}`\n")
    md_path.write_text("".join(md), encoding="utf-8")

    print(f"CSV written: {csv_path}")
    print(f"Summary written: {md_path}")
    print(f"Plot written: {plot_path}")


if __name__ == "__main__":
    main()

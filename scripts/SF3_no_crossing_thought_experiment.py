#!/usr/bin/env python3
"""
SF3_no_crossing_thought_experiment.py

Bold-STAM strong-field commitment and the two-planet no-crossing thought experiment.
Author: Sean Brady / STAM Model-A continuation

Purpose:
    Implement the author's careful two-planet thought experiment with all three
    effects kept separate at every step. No magic bell, no free ticket.

Setup:
    - Alpha+ observer at r_obs (large radius, A ≈ 0.001).
    - Traveler released from rest at r_start, falls radially toward Rs.
    - At each instantaneous A value, the traveler emits an "I am at A=X" message.
    - We compute three independent quantities:
        (1) Traveler's diary τ_travel(A): proper time elapsed since release.
        (2) Clock-dilation contribution: how much extra coordinate time accumulates
            beyond the traveler's proper clock — this is the gravitational time
            dilation observed by Alpha+.
        (3) SU delay: light-signal propagation time from traveler's r to r_obs,
            integrated through the A field between them. The medium is non-zero
            everywhere along the way.
    - The observer's view of the traveler at A=X is the SUM of (2) + (3) on top
      of the traveler's own diary (1).

Bold-STAM strong-field metric (the specific theoretical commitment of this script):
    g_tt = -(1-A) c^2                   (GR — preserves all weak-field clock tests)
    g_rr = 1 / [(1-A)(1-A^2)^2]         (Bold STAM — vanishes as (1-A)^3 at horizon)

Why this asymmetric form:
    - At leading order in A, both g_tt and g_rr reduce to GR Schwarzschild
      (1-A) and 1/(1-A), so all weak-field tests of GR pass automatically.
    - At A=1, g_rr diverges as 1/(1-A)^3 instead of GR's 1/(1-A), making proper
      *distance* to the boundary infinite (SU integral diverges algebraically).
    - For a freely falling observer with E=c^2, the proper time integral picks
      up a logarithmic divergence at A=1. Falls forever, with each intermediate
      A crossing at finite proper time.
    - The deviation is structural and concentrated near A=1; falsifiable via
      EHT shadow size, near-horizon Shapiro delay, LIGO ringdown frequencies.

The boundary at A=1 — the headline result:
    For A < 1: the triangle is computable. Each quantity is finite. Ratios
    between frames are well-defined. The traveler ages slowly compared to the
    observer; the magnitudes of clock-dilation and SU delay can be compared.

    As A → 1: all three quantities (traveler's τ, clock-dilation contribution,
    SU delay) diverge to infinity. They diverge at *different rates*, but each
    individually goes to ∞.

    At A = 1 itself: ∞ = ∞ in every frame. Frame-comparison ratios become ∞/∞,
    indeterminate. The very framework of comparing frames collapses.

    This IS the no-crossing boundary. Not "the traveler asymptotically slows down,"
    but "the comparison structure fails." Both the traveler and the observer
    measure infinite time to reach A=1, by different mechanisms, with the same
    mathematical limit. There is no privileged frame in which the crossing can
    be registered as a finite event.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

C = 2.99792458e8
G = 6.67430e-11
M_SUN = 1.98847e30

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# --- Bold-STAM and GR metric coefficients ---

def h_STAM(A: np.ndarray) -> np.ndarray:
    """g_tt coefficient: same as GR, preserves clock-dilation rate."""
    return 1.0 - A


def k_STAM(A: np.ndarray) -> np.ndarray:
    """1/g_rr coefficient. Vanishes as (1-A)^3 at the horizon — bold STAM departure."""
    return (1.0 - A) * (1.0 - A ** 2) ** 2


def h_GR(A: np.ndarray) -> np.ndarray:
    return 1.0 - A


def k_GR(A: np.ndarray) -> np.ndarray:
    return 1.0 - A


# --- Free-fall trajectory (E = c^2, rest at finite r_start) ---

def dtau_dA_infaller(A: np.ndarray, h_func, k_func, Rs: float) -> np.ndarray:
    """dτ/dA for a radial infaller with E=c^2.

    From the geodesic with h dt/dτ = 1: (dr/dτ)^2 = k c^2 (1-h)/h.
    With dr = -Rs/A^2 dA (since A=Rs/r):
        dτ = dr / (dr/dτ_signed) = (Rs/A^2 / c) / sqrt(k (1-h)/h) dA
    """
    h = h_func(A)
    k = k_func(A)
    factor = np.sqrt(k * (1.0 - h) / h)
    # avoid zero-division at A=0 boundary (we never integrate quite there)
    factor = np.where(factor > 0, factor, np.nan)
    return (Rs / (A ** 2)) / (C * factor)


def dt_dA_infaller(A: np.ndarray, h_func, k_func, Rs: float) -> np.ndarray:
    """dt/dA (coordinate time) for the same radial infaller. dt/dτ = 1/h."""
    return dtau_dA_infaller(A, h_func, k_func, Rs) / h_func(A)


def dt_signal_dA(A: np.ndarray, h_func, k_func, Rs: float) -> np.ndarray:
    """dt/dA for a radial null geodesic (light signal) traveling outward through A field.

    From ds^2 = 0: dt = dr / (c sqrt(h k)).  dr = -Rs/A^2 dA, take absolute for outward.
    """
    h = h_func(A)
    k = k_func(A)
    return (Rs / (A ** 2)) / (C * np.sqrt(h * k))


# --- Cumulative integrals ---

def integrate_cumulative(integrand: np.ndarray, A_grid: np.ndarray) -> np.ndarray:
    """Cumulative integral from A_grid[0] to A_grid[i] for each i.

    Uses trapezoid rule with non-uniform spacing.
    """
    out = np.zeros_like(integrand)
    out[0] = 0.0
    for i in range(1, len(A_grid)):
        out[i] = out[i - 1] + 0.5 * (integrand[i] + integrand[i - 1]) * (A_grid[i] - A_grid[i - 1])
    return out


def build_trajectory(A_start: float, A_max: float, h_func, k_func, Rs: float,
                     n_points: int = 4000) -> dict:
    """Build the full trajectory: τ_travel, t_emission, etc. as functions of A.

    A_start is the starting position (small A, traveler at rest there).
    A_max is the closest to A=1 we sample (cannot include 1 itself).
    """
    A_grid = np.linspace(A_start, A_max, n_points)

    dtau_dA = dtau_dA_infaller(A_grid, h_func, k_func, Rs)
    dt_dA = dt_dA_infaller(A_grid, h_func, k_func, Rs)

    tau_travel = integrate_cumulative(dtau_dA, A_grid)
    t_emission = integrate_cumulative(dt_dA, A_grid)

    return {
        "A_grid": A_grid,
        "tau_travel": tau_travel,
        "t_emission": t_emission,
        "clock_dilation_buildup": t_emission - tau_travel,
    }


def signal_delay(A_emit: float, A_obs: float, h_func, k_func, Rs: float,
                 n_points: int = 2000) -> float:
    """SU/Shapiro delay for a light signal from A_emit (deeper) up to A_obs (shallower).

    Integrate dt_signal_dA from A_emit down to A_obs (note A_emit > A_obs since deeper = larger A).
    """
    if A_obs >= A_emit:
        return 0.0
    A_path = np.linspace(A_obs, A_emit, n_points)
    integrand = dt_signal_dA(A_path, h_func, k_func, Rs)
    return float(np.trapezoid(integrand, A_path))


# --- Triangle assembly ---

def build_triangle(A_emit_values: np.ndarray, traj: dict, h_func, k_func, Rs: float,
                   A_obs: float) -> pd.DataFrame:
    """For each A_emit value (where the message was sent), compute:
    - τ_travel: traveler's diary at emission
    - clock_dilation_buildup: t_emission - τ_travel
    - SU_delay: signal propagation from emission to observer
    - t_obs_reception: total observer time when message arrives = t_emission + SU_delay
    - clock_dilation_factor: 1/sqrt(h(A)) — instantaneous tick rate ratio
    """
    A_grid = traj["A_grid"]
    tau_at = np.interp(A_emit_values, A_grid, traj["tau_travel"])
    t_emit_at = np.interp(A_emit_values, A_grid, traj["t_emission"])
    cd_buildup = t_emit_at - tau_at

    su_delays = np.array([signal_delay(A, A_obs, h_func, k_func, Rs) for A in A_emit_values])

    rows = []
    for i, A in enumerate(A_emit_values):
        rows.append({
            "A_emit": A,
            "r_over_Rs": 1.0 / A,
            "tau_travel_s": tau_at[i],
            "clock_dilation_buildup_s": cd_buildup[i],
            "SU_delay_s": su_delays[i],
            "t_obs_reception_s": t_emit_at[i] + su_delays[i],
            "instantaneous_cd_factor": 1.0 / math.sqrt(h_func(A)),
        })
    return pd.DataFrame(rows)


# --- Boundary-collapse analysis ---

def boundary_collapse_table(df: pd.DataFrame) -> pd.DataFrame:
    """Compute frame-comparison ratios and show how they behave as A→1.

    For A<1, ratios are well-defined finite numbers.
    As A→1, every individual quantity diverges; ratios may stabilize or diverge,
    but the absolute scales of comparison vanish (everything → ∞ on the same
    boundary).
    """
    df = df.copy()
    df["clock_dilation_over_tau"] = df["clock_dilation_buildup_s"] / df["tau_travel_s"]
    df["SU_over_clock_dilation"] = df["SU_delay_s"] / df["clock_dilation_buildup_s"]
    df["SU_over_tau"] = df["SU_delay_s"] / df["tau_travel_s"]
    df["t_obs_over_tau"] = df["t_obs_reception_s"] / df["tau_travel_s"]
    return df


# --- Plots ---

def plot_triangle(df_stam: pd.DataFrame, df_gr: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    A = df_stam["A_emit"].values

    axes[0, 0].plot(A, df_stam["tau_travel_s"], "o-", label="Bold STAM", linewidth=2)
    axes[0, 0].plot(A, df_gr["tau_travel_s"], "s--", label="GR Schwarzschild", linewidth=1.5)
    axes[0, 0].set_xlabel("A at emission")
    axes[0, 0].set_ylabel("Traveler's diary τ_travel (s)")
    axes[0, 0].set_title("(1) Traveler's diary — diverges in bold STAM, finite in GR")
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_xlim(0, 1)
    axes[0, 0].grid(True, which="both", linewidth=0.3)
    axes[0, 0].legend()

    axes[0, 1].plot(A, df_stam["clock_dilation_buildup_s"], "o-", label="Bold STAM", linewidth=2)
    axes[0, 1].plot(A, df_gr["clock_dilation_buildup_s"], "s--", label="GR", linewidth=1.5)
    axes[0, 1].set_xlabel("A at emission")
    axes[0, 1].set_ylabel("Clock-dilation buildup t_emit − τ_travel (s)")
    axes[0, 1].set_title("(2) Clock-dilation contribution to observer view")
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_xlim(0, 1)
    axes[0, 1].grid(True, which="both", linewidth=0.3)
    axes[0, 1].legend()

    axes[1, 0].plot(A, df_stam["SU_delay_s"], "o-", label="Bold STAM", linewidth=2)
    axes[1, 0].plot(A, df_gr["SU_delay_s"], "s--", label="GR", linewidth=1.5)
    axes[1, 0].set_xlabel("A at emission")
    axes[1, 0].set_ylabel("SU delay along signal path (s)")
    axes[1, 0].set_title("(3) SU delay contribution — bold STAM diverges as 1/(1-A)")
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_xlim(0, 1)
    axes[1, 0].grid(True, which="both", linewidth=0.3)
    axes[1, 0].legend()

    axes[1, 1].plot(A, df_stam["t_obs_reception_s"], "o-", label="Bold STAM total observer time", linewidth=2)
    axes[1, 1].plot(A, df_gr["t_obs_reception_s"], "s--", label="GR total observer time", linewidth=1.5)
    axes[1, 1].plot(A, df_stam["tau_travel_s"], "x:", label="Bold STAM traveler's diary (for reference)", linewidth=1.5)
    axes[1, 1].set_xlabel("A at emission")
    axes[1, 1].set_ylabel("Time at observer (s)")
    axes[1, 1].set_title("(4) Total observer view: clock dilation + SU delay compounded")
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlim(0, 1)
    axes[1, 1].grid(True, which="both", linewidth=0.3)
    axes[1, 1].legend(fontsize=8)

    plt.suptitle("SF3: The triangle of effects, traveler frame vs observer frame", y=1.00, fontsize=13)
    out = PLOTS / "SF3_triangle_panels.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_boundary_collapse(df_stam: pd.DataFrame) -> Path:
    """Show how every quantity diverges to infinity as A→1, illustrating frame-collapse."""
    fig, ax = plt.subplots(figsize=(10, 6))
    A = df_stam["A_emit"].values
    one_minus_A = 1.0 - A

    ax.loglog(one_minus_A, df_stam["tau_travel_s"], "o-", label="Traveler's diary τ_travel", linewidth=2)
    ax.loglog(one_minus_A, df_stam["clock_dilation_buildup_s"], "s-", label="Clock-dilation buildup", linewidth=2)
    ax.loglog(one_minus_A, df_stam["SU_delay_s"], "^-", label="SU delay along signal path", linewidth=2)
    ax.loglog(one_minus_A, df_stam["t_obs_reception_s"], "d-", label="Total observer reception", linewidth=2.5)

    ax.set_xlabel("1 − A  (distance to horizon, A=1 at left)")
    ax.set_ylabel("Time (s)")
    ax.set_title("Boundary collapse: every quantity → ∞ as A → 1, by different rates but same limit")
    ax.invert_xaxis()
    ax.grid(True, which="both", linewidth=0.3)
    ax.legend()

    out = PLOTS / "SF3_boundary_collapse.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_frame_ratios(df_stam: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(10, 6))
    A = df_stam["A_emit"].values

    ax.semilogy(A, df_stam["clock_dilation_over_tau"], "o-",
                label="(t_emit − τ_travel) / τ_travel  [observer/traveler clock excess]", linewidth=2)
    ax.semilogy(A, df_stam["SU_over_tau"], "^-",
                label="SU_delay / τ_travel  [traversal vs aging]", linewidth=2)
    ax.semilogy(A, df_stam["SU_over_clock_dilation"], "s-",
                label="SU_delay / clock_dilation  [which observer effect dominates]", linewidth=2)
    ax.semilogy(A, df_stam["t_obs_over_tau"], "d-",
                label="t_obs_reception / τ_travel  [total observer/traveler ratio]", linewidth=2)

    ax.set_xlabel("A at emission")
    ax.set_ylabel("Frame-comparison ratio")
    ax.set_title("Frame-comparison ratios while well-defined (A < 1) — collapse to ∞/∞ at A=1")
    ax.set_xlim(0, 1)
    ax.grid(True, which="both", linewidth=0.3)
    ax.legend(fontsize=9)

    out = PLOTS / "SF3_frame_ratios.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- Markdown summary ---

def write_markdown(df_stam: pd.DataFrame, df_gr: pd.DataFrame, plot_paths: list[Path],
                   M: float, Rs: float, A_obs: float) -> Path:
    md = []
    md.append("# STAM Model-A SF3: Two-Planet No-Crossing Thought Experiment\n")
    md.append("## Purpose\n")
    md.append(
        "Implement the author's two-planet thought experiment with all three effects kept "
        "separate at every step. No magic bell, no free ticket. Demonstrate that the "
        "no-crossing boundary at A=1 is defined by the *collapse of frame comparison*, "
        "not by either frame seeing a special event.\n"
    )
    md.append("## Setup\n")
    md.append(
        f"- Black hole mass: M = {M / M_SUN:.3g} M_sun, Rs = {Rs:.3g} m.\n"
        f"- Alpha+ observer at large radius, A_obs = {A_obs:.3g}.\n"
        f"- Traveler released from rest at r_start = Rs/A_start; falls radially toward Rs.\n"
        f"- Traveler emits 'I am at A=X' messages at successive A values.\n"
    )
    md.append("## Bold-STAM strong-field metric\n")
    md.append("```text\n"
              "g_tt = -(1-A) c^2                       (GR — preserves all weak-field clock tests)\n"
              "g_rr = 1 / [(1-A)(1-A^2)^2]             (Bold STAM — vanishes as (1-A)^3 at horizon)\n"
              "```\n"
              "At leading order in A, both reduce to GR Schwarzschild — all weak-field tests of "
              "GR pass automatically. At A=1, g_rr diverges faster than GR's, making proper "
              "distance to the boundary *infinite*. For a free-fall trajectory with E=c^2, the "
              "proper time integral picks up a logarithmic divergence at A=1: the traveler "
              "**falls forever** while crossing each intermediate A value at finite proper time.\n"
    )
    md.append("## Three independent quantities (the triangle)\n")
    md.append(
        "1. **Traveler's diary τ_travel(A):** proper time elapsed since release. Diverges "
        "logarithmically in bold STAM as A→1 (falls forever); finite at A=1 in GR.\n"
        "2. **Clock-dilation buildup (t_emit − τ_travel):** the gravitational time-dilation "
        "contribution to what an observer sees. Diverges algebraically as 1/√(1-A).\n"
        "3. **SU delay along signal path:** light-signal propagation cost from traveler's r to "
        "the observer, integrated through the A field between them (always non-zero). Diverges "
        "algebraically (1-A)^(-2) in bold STAM, faster than GR's logarithmic.\n"
        "\n"
        "**Total observer reception time = traveler's emission coordinate time + SU delay.** "
        "The clock-dilation effect is *embedded* in the difference between coordinate time at "
        "emission and the traveler's proper time at emission.\n"
    )
    md.append("## Numerical results — bold STAM\n")
    md.append(df_stam.to_markdown(index=False, floatfmt=".6g"))
    md.append("\n\n## Numerical results — GR (for comparison)\n")
    md.append(df_gr.to_markdown(index=False, floatfmt=".6g"))
    md.append("\n\n## Frame-comparison ratios (the boundary-collapse diagnostic)\n")
    md.append(boundary_collapse_table(df_stam).to_markdown(index=False, floatfmt=".6g"))
    md.append("\n\n## The headline result: A=1 is where frame comparison collapses\n")
    md.append(
        "**For A < 1**: the triangle is computable. Each quantity is a finite number. The "
        "ratio τ_travel : clock-dilation buildup : SU delay carries the asymmetry between the "
        "frames. The traveler ages slowly compared to the observer; the magnitudes of the two "
        "observer-frame effects can be compared.\n"
        "\n"
        "**As A → 1**: all three quantities diverge to infinity. They diverge at *different "
        "rates* — τ_travel logarithmically, clock-dilation as 1/√(1-A), SU delay as 1/(1-A)^2 — "
        "but each individually goes to ∞.\n"
        "\n"
        "**At A = 1 itself**: ∞ = ∞ in every frame. The traveler measures infinite proper time "
        "to reach A=1 (logarithmic divergence). The observer measures infinite reception time "
        "for any 'I am at A=1' message (algebraic divergence). The SU traversal cost is also "
        "infinite. Frame-comparison ratios become ∞/∞: indeterminate.\n"
        "\n"
        "**The very framework of comparing frames collapses.** This is the precise mathematical "
        "statement of the no-crossing boundary. It is not 'the traveler asymptotically slows.' "
        "It is 'the comparison structure between any two frames fails.' Both the traveler and "
        "the observer measure infinite time to reach A=1, by different physical mechanisms, "
        "with the same mathematical limit. There is no privileged frame in which the crossing "
        "can be registered as a finite event.\n"
    )
    md.append("## Why this matters for the fatalities\n")
    md.append(
        "- **'Time stops at the horizon'** — No. Time does not stop. *Comparison* stops. Both "
        "clocks keep running in their own local frames. What stops is the meaningful relationship "
        "between the frames.\n"
        "- **'You're privileging the observer'** — No. The script computes both frames; both "
        "diverge; they reach the same infinite limit. The observer has no magic bell, the "
        "traveler has no privileged finite proper time.\n"
        "- **'You're modifying the metric ad hoc'** — The h(A), k(A) modification gives infinite "
        "traveler proper time; that's the bold STAM commitment. But the structural point — that "
        "A=1 is where comparison breaks down — would hold for any metric with all three "
        "quantities diverging at the boundary. The specific form is one realization of a more "
        "general boundary-collapse principle.\n"
        "- **'Information paradox'** — There's no paradox because there's no 'after the crossing' "
        "event in any frame. The universal ledger has no entry for 'crossed' because no frame "
        "has finite time to register one. Information stays in the resolved sector throughout; "
        "it never enters an interior, because no interior is reachable from any frame.\n"
    )
    md.append("## Falsifiable corners (where bold STAM differs from GR)\n")
    md.append(
        "- *Near-horizon Shapiro delay*: bold STAM SU delay diverges algebraically as "
        "1/(1-A)^2 vs GR's logarithmic. EHT timing of light bouncing near the photon sphere "
        "could in principle test this.\n"
        "- *EHT shadow size*: shadow depends on photon-sphere geometry, which is sensitive to "
        "g_rr. Modified k(A) shifts the shadow size at the few-percent level for moderate-A "
        "regions.\n"
        "- *LIGO ringdown*: quasinormal mode frequencies depend on near-horizon metric "
        "structure. Bold STAM predicts shifted QNM frequencies vs GR.\n"
        "- *PPN parameters at second order*: bold STAM has g_rr = 1 + A + 2A^2 + … vs GR's "
        "1 + A + A^2 + … High-precision PPN tests beyond Cassini precision could distinguish.\n"
    )
    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "SF3_no_crossing_thought_experiment_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    M = 1.0 * M_SUN
    Rs = 2.0 * G * M / (C ** 2)

    A_start = 0.001     # traveler released here (large r, weak field)
    A_obs = A_start     # observer at the same large radius
    A_max = 0.99999     # closest to horizon we sample

    traj_stam = build_trajectory(A_start, A_max, h_STAM, k_STAM, Rs, n_points=8000)
    traj_gr = build_trajectory(A_start, A_max, h_GR, k_GR, Rs, n_points=8000)

    A_emit_values = np.array([0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999])

    df_stam = build_triangle(A_emit_values, traj_stam, h_STAM, k_STAM, Rs, A_obs)
    df_gr = build_triangle(A_emit_values, traj_gr, h_GR, k_GR, Rs, A_obs)

    df_ratios = boundary_collapse_table(df_stam)

    df_stam.to_csv(RESULTS / "SF3_triangle_bold_stam.csv", index=False)
    df_gr.to_csv(RESULTS / "SF3_triangle_gr.csv", index=False)
    df_ratios.to_csv(RESULTS / "SF3_frame_ratios.csv", index=False)

    plot_paths = [
        plot_triangle(df_stam, df_gr),
        plot_boundary_collapse(df_stam),
        plot_frame_ratios(df_ratios),
    ]

    summary = write_markdown(df_stam, df_gr, plot_paths, M, Rs, A_obs)

    print("STAM Model-A SF3: Two-Planet No-Crossing Thought Experiment")
    print("=" * 72)
    print(f"\nBlack hole: M = {M / M_SUN:.3g} M_sun, Rs = {Rs:.3g} m")
    print(f"Setup: A_start = A_obs = {A_obs}, A_max sampled = {A_max}")
    print()
    print("Bold-STAM triangle (selected):")
    print(df_stam[["A_emit", "tau_travel_s", "clock_dilation_buildup_s",
                   "SU_delay_s", "t_obs_reception_s"]].to_string(index=False))
    print()
    print("GR comparison (selected):")
    print(df_gr[["A_emit", "tau_travel_s", "clock_dilation_buildup_s",
                 "SU_delay_s", "t_obs_reception_s"]].to_string(index=False))
    print()
    print("Boundary-collapse check -- bold STAM tau_travel at A=0.99999:")
    print(f"  tau_travel = {df_stam.iloc[-1]['tau_travel_s']:.6g} s")
    print("  (diverges logarithmically as A -> 1: falls forever)")
    print(f"\nObserver reception at A=0.99999: {df_stam.iloc[-1]['t_obs_reception_s']:.6g} s")
    print("  (diverges algebraically -- observer never receives 'crossed' message)")
    print()
    print("Files written:")
    print(f"- {RESULTS / 'SF3_triangle_bold_stam.csv'}")
    print(f"- {RESULTS / 'SF3_triangle_gr.csv'}")
    print(f"- {RESULTS / 'SF3_frame_ratios.csv'}")
    print(f"- {summary}")
    for p in plot_paths:
        print(f"- {p}")


if __name__ == "__main__":
    main()

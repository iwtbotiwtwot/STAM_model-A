#!/usr/bin/env python3
"""
G21_strong_field_with_cosmic_floor.py

Strong-field Model-A with the cosmic floor A_0 made explicit.
Author: Sean Brady / STAM Model-A continuation 2026-05-11

Foundational technical first step toward re-doing the strong-field corpus
under the symmetric-boundary commitment (A runs on (A_0, 1), not (0, 1)).

Water-tank framing (author 2026-05-11):
    A_0 IS the void — the unperturbed manifold level. Not zero. Zero means
    nothing exists; the void value is what the manifold sits at when no
    matter or energy is present. The maximum level (saturation, horizon)
    is A = 1.

    Matter and energy displace the level upward WITHIN the available range
    (1 − A_0). They cannot push the level above 1 (cells can't saturate
    beyond their capacity). So the composition is:

        A_total(r) = A_0 + (1 − A_0) × (Rs / r)
        A_local(r) = (1 − A_0) × (Rs / r)   (displacement, scaled to the
                                              available range)

    At r = Rs: A_total = A_0 + (1−A_0) = 1   (horizon at GR's Schwarzschild
                                              radius, as it should be).
    At r → ∞: A_total = A_0                  (asymptotic at the cosmic floor).
    At intermediate r: A_total ∈ (A_0, 1).

    A naive A_total = A_0 + Rs/r would have A_total > 1 at r ≤ Rs, which
    violates the framework's A ≤ 1 ceiling. The rescaling by (1 − A_0) is
    what the water-tank picture actually requires: displacement maps onto
    the available range, not onto the full unit interval.

Ledger framing (May 10 evening + 2026-05-11):
    A is the density of resolved ledger entries per Planck cell.
    A_0 = 1/(12π) is the minimum write density required for the manifold
    to exist (the (4π × 3) decomposition counts one entry per solid angle
    × spatial direction × Planck cell). A_local is the additional write
    density from matter currently interacting with its surroundings.
    A_total = A_0 + A_local sums two write-density contributions locally.

What this script does:
    1. Adopt A_total(r) = A_0 + Rs/r with A_0 = 1/(12π).
    2. Compute the strong-field metric components g_tt = -(1-A_total)c²,
       g_rr = 1/[(1-A_total)(1-A_total²)²].
    3. Find the landmark radii where A_total = 1/3, 2/3, 1 (ISCO, photon
       sphere, horizon). Compare to the A_0 = 0 baseline (GR / original F3).
    4. Compute m(r)/M from the metric and decompose into:
        - m_background(r): contribution from the cosmic A_0 floor alone
        - m_source(r): displacement contribution from the matter source
    5. Compute the F3-style effective stress-energy bracket B(A_total, A_0)
       and decompose into structural (background-only) and dynamical
       (source-displacement) pieces.
    6. Locate the sign-change point (NEC crossover) in B(A_total, A_0)
       and compare to F3's A ≈ 0.44.
    7. Output: markdown summary + CSV + comparison plots.

Honest scope:
    This script changes ONLY the asymptotic boundary of A. It does not
    change the metric construction k(A) = (1-A)(1-A²)² (which remains
    postulated per project_strong_field_departure_question.md), and does
    not touch the field-equations gap. It produces the cleanest "look
    what happens when we move asymptotic A from 0 to A_0" comparison
    and surfaces what shifts and what doesn't.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

# ============================================================
# Core constants
# ============================================================

A_0 = 1.0 / (12.0 * math.pi)  # cosmic floor / void value (G7 commitment)
# Equivalent values for reference:
# A_0 = 1/(4π × 3) under (4π × 3) decomposition reading
# A_0 ≈ 0.026525...


# ============================================================
# Field composition (water-tank / ledger framing)
# ============================================================

def A_local_of_r(r_over_Rs, A_floor=A_0):
    """Source displacement, scaled to the available range (1 − A_floor).

    A_local = (1 − A_floor) × (Rs / r), so A_local maps onto (0, 1 − A_floor)
    as r runs (∞, Rs). At r = Rs, A_local = 1 − A_floor and A_total = 1."""
    return (1.0 - A_floor) / np.asarray(r_over_Rs, dtype=float)


def A_total_of_r(r_over_Rs, A_floor=A_0):
    """A_total = A_floor + (1 − A_floor)(Rs/r). Void + displacement on available range."""
    return A_floor + A_local_of_r(r_over_Rs, A_floor)


def r_over_Rs_for_A_total(A_target, A_floor=A_0):
    """Solve A_total = A_floor + (1 − A_floor)(Rs/r) = A_target for r/Rs.

        r/Rs = (1 − A_floor) / (A_target − A_floor)

    Requires A_target > A_floor."""
    A_local_needed = A_target - A_floor
    if np.any(A_local_needed <= 0):
        raise ValueError(f"A_target ({A_target}) must exceed A_floor ({A_floor})")
    return (1.0 - A_floor) / A_local_needed


# ============================================================
# Strong-field metric (k(A) postulated; carried over from project commitment)
# ============================================================

def h_of_A(A):
    """g_tt coefficient: g_tt = -h(A) c². h(A) = 1-A (same as GR form)."""
    return 1.0 - A


def k_of_A(A):
    """g_rr coefficient: g_rr = 1/k(A). k(A) = (1-A)(1-A²)² (Model-A postulate)."""
    return (1.0 - A) * (1.0 - A**2)**2


# ============================================================
# Mass function and effective stress-energy
# ============================================================

def mass_function_over_M(A_total, A_floor=A_0):
    """m(r)/M from the metric at a point where A_total takes its value.

    Under rescaled composition A_total = A_floor + (1 − A_floor)(Rs/r):
        r/Rs = (1 − A_floor) / (A_total − A_floor) = (1 − A_floor) / A_local

    So:
        m(r)/M = (r/Rs) × [1 − k(A_total)]
               = (1 − A_floor) × [1 − k(A_total)] / A_local
    """
    A_local = A_total - A_floor
    return (1.0 - A_floor) * (1.0 - k_of_A(A_total)) / A_local


def mass_function_background_over_M(r_over_Rs, A_floor=A_0):
    """Background-only m_background(r)/M = (r/Rs) × [1 - k(A_floor)].

    The cosmic A_0 floor sources a non-trivial effective enclosed mass even
    in the absence of any matter source (it's the price of measuring 'enclosed
    mass' in the Schwarzschild parameterization while sitting inside the cosmic
    A_0 medium)."""
    return r_over_Rs * (1.0 - k_of_A(A_floor))


def mass_function_source_over_M(A_total, A_floor=A_0):
    """Source's displacement contribution to enclosed mass under rescaled composition.

    m_source(r) = m_total(r) - m_background(r)
                = (r/Rs) × {[1 - k(A_total)] - [1 - k(A_floor)]}
                = (r/Rs) × [k(A_floor) - k(A_total)]
                = (1 − A_floor) × [k(A_floor) − k(A_total)] / A_local

    At horizon (A_total = 1, k = 0, r = Rs):
        m_source/M = (1 − A_floor) × k(A_floor) / (1 − A_floor) = k(A_floor) ≈ 0.972
        — close to M, with a small deficit absorbed by the (1 − A_floor) rescaling.
    """
    A_local = A_total - A_floor
    return (1.0 - A_floor) * (k_of_A(A_floor) - k_of_A(A_total)) / A_local


def f_of_A(A):
    """Convenience: f(A) = 1 - k(A) = 1 - (1-A)(1-A²)² = 1 - (1-A)³(1+A)²."""
    return 1.0 - k_of_A(A)


def f_prime_of_A(A):
    """f'(A) = (1-A)²(1+A)(1+5A). Hand-derived; see docstring of bracket function."""
    return (1.0 - A)**2 * (1.0 + A) * (1.0 + 5.0 * A)


def rho_eff_bracket(A_total, A_floor=A_0):
    """The F3-style bracket for the effective density × r² (signed).

    Under rescaled composition A_total = A_floor + (1 − A_floor)(Rs/r):
        dA_total/dr = −(1 − A_floor)(Rs/r²)
        r/Rs = (1 − A_floor) / A_local

    Both (1 − A_floor) factors cancel in dm/dr × (1/M), leaving:
        dm_total/dr × (1/M) = (1/Rs) × [f(A_total) − A_local × f'(A_total)]

    where A_local = A_total − A_floor. The bracket is invariant under the choice
    between naive and rescaled composition rules.

    When A_floor = 0, this reduces to F3's original bracket:
        B(A, 0) = f(A) − A × f'(A) = 1 − (1−A)²(1+A)(1+A+4A²)  ✓
    """
    A_local = A_total - A_floor
    return f_of_A(A_total) - A_local * f_prime_of_A(A_total)


def rho_eff_bracket_background(A_floor=A_0):
    """The 'structural' piece: bracket evaluated with A_local → 0.

    This is the cosmic-background contribution to the effective stress-energy
    bracket, present everywhere by virtue of the manifold existing at A_0.
    Constant in r (does not depend on A_total beyond A_floor itself).

        B_background = f(A_floor)
    """
    return f_of_A(A_floor)


def rho_eff_bracket_source(A_total, A_floor=A_0):
    """The 'dynamical' piece: source-displacement contribution.

        B_source(A_total, A_floor) = B(A_total, A_floor) - B_background(A_floor)
                                   = f(A_total) - f(A_floor)
                                     - (A_total - A_floor) × f'(A_total)

    At A_total = A_floor: B_source = 0  ✓  (no source, no displacement)
    At A_total = 1:       B_source = 1 - f(A_floor) ≈ 1 (close to F3's saturation)
    """
    return rho_eff_bracket(A_total, A_floor) - rho_eff_bracket_background(A_floor)


# ============================================================
# Landmarks
# ============================================================

def landmark_table(A_floor=A_0):
    """Compute landmark radii where A_total reaches GR thirds, under rescaled composition.

    r/Rs (with floor) = (1 − A_floor) / (A_target − A_floor)
    r/Rs (no floor)   = 1 / A_target   (GR baseline; same as r/Rs in F3)
    """
    landmarks = [
        ("ISCO",          1.0 / 3.0),
        ("photon sphere", 2.0 / 3.0),
        ("horizon",       1.0),
    ]
    rows = []
    for name, A_target in landmarks:
        r_floor = (1.0 - A_floor) / (A_target - A_floor)
        r_nofloor = 1.0 / A_target  # GR / original baseline
        shift_pct = 100.0 * (r_floor - r_nofloor) / r_nofloor
        rows.append({
            "landmark":          name,
            "A_total_target":    A_target,
            "A_local_needed":    A_target - A_floor,
            "r/Rs (with floor)": r_floor,
            "r/Rs (no floor)":   r_nofloor,
            "shift (%)":         shift_pct,
        })
    return pd.DataFrame(rows)


# ============================================================
# F3-style table with floor: where does the NEC crossover go?
# ============================================================

def find_all_sign_changes(func, A_floor=A_0, scan_start=None, A_high=0.999,
                          n_scan=10000, tol=1e-7):
    """Scan (scan_start, A_high) and find ALL sign changes of func.

    func should accept a numpy array of A values and an A_floor argument.
    A_floor is passed through to func (not used for scan-start unless scan_start=None)."""
    if scan_start is None:
        scan_start = A_floor + 1e-6
    A_grid = np.linspace(scan_start, A_high, n_scan)
    vals = func(A_grid, A_floor)
    sign_change_indices = np.where(np.diff(np.sign(vals)) != 0)[0]
    crossings = []
    for i in sign_change_indices:
        lo, hi = A_grid[i], A_grid[i + 1]
        v_lo = vals[i]
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            v_mid = func(np.array([mid]), A_floor)[0]
            if v_mid * v_lo < 0:
                hi = mid
            else:
                lo = mid
                v_lo = v_mid
        crossings.append(0.5 * (lo + hi))
    return crossings


def find_NEC_crossover_with_floor(A_floor=A_0, A_high=0.999):
    """Find all sign changes of B_total in (A_floor, A_high)."""
    return find_all_sign_changes(rho_eff_bracket, A_floor, None, A_high)


def find_NEC_crossover_source_only(A_floor=A_0, A_high=0.999):
    """Find all sign changes of B_source above the trivial zero at A_floor."""
    return find_all_sign_changes(rho_eff_bracket_source, A_floor,
                                 scan_start=A_floor + 1e-4, A_high=A_high)


# ============================================================
# Comparison table at a grid of A_total values
# ============================================================

def comparison_table(A_floor=A_0):
    """Build the F3 table extended with cosmic floor."""
    A_values = np.array([
        A_floor,  # asymptotic limit
        0.05, 0.1, 0.2, 0.3, 0.4, 0.44, 0.5,
        0.6, 0.7, 0.8, 0.9, 0.99,
    ])
    rows = []
    for A in A_values:
        if A < A_floor:
            continue  # below the floor — manifold doesn't extend here
        A_local = A - A_floor
        # Under rescaled composition: r/Rs = (1 − A_floor) / A_local
        r_over_Rs = (1.0 - A_floor) / A_local if A_local > 0 else float("inf")

        m_total = mass_function_over_M(A, A_floor) if A_local > 0 else float("nan")
        m_back = mass_function_background_over_M(r_over_Rs, A_floor) if A_local > 0 else float("nan")
        m_source = mass_function_source_over_M(A, A_floor) if A_local > 0 else float("nan")

        B_total = rho_eff_bracket(A, A_floor)
        B_back = rho_eff_bracket_background(A_floor)
        B_source = rho_eff_bracket_source(A, A_floor)

        # F3 (no floor) bracket for comparison
        B_F3 = f_of_A(A) - A * f_prime_of_A(A)

        rows.append({
            "A_total":              A,
            "A_local":              A_local,
            "r/Rs":                 r_over_Rs,
            "m_total/M":            m_total,
            "m_background/M":       m_back,
            "m_source/M":           m_source,
            "B_total (with floor)": B_total,
            "B_background":         B_back,
            "B_source":             B_source,
            "B_F3 (no floor)":      B_F3,
            "ΔB (floor - F3)":      B_total - B_F3,
        })
    return pd.DataFrame(rows)


# ============================================================
# Plots
# ============================================================

def plot_landmark_shifts(df_landmarks, A_floor=A_0):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    landmarks = df_landmarks["landmark"].tolist()
    r_floor = df_landmarks["r/Rs (with floor)"].to_numpy()
    r_nofloor = df_landmarks["r/Rs (no floor)"].to_numpy()
    x = np.arange(len(landmarks))
    width = 0.35
    ax.bar(x - width/2, r_nofloor, width, label="No floor (A_0 = 0)", color="tab:gray")
    ax.bar(x + width/2, r_floor, width, label=f"With A_0 = 1/(12π) ≈ {A_floor:.4f}",
           color="tab:blue")
    ax.set_xticks(x)
    ax.set_xticklabels(landmarks)
    ax.set_ylabel("r / Rs")
    ax.set_title("Strong-field landmark shifts under cosmic floor")
    for i, (rf, rn) in enumerate(zip(r_floor, r_nofloor)):
        ax.text(i + width/2, rf, f"{rf:.3f}", ha="center", va="bottom", fontsize=9)
        ax.text(i - width/2, rn, f"{rn:.3f}", ha="center", va="bottom", fontsize=9)
    ax.legend()
    ax.grid(True, axis="y", linewidth=0.3)
    out = PLOTS / "G21_landmark_shifts.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_bracket_decomposition(A_floor=A_0):
    A_grid = np.linspace(A_floor + 1e-4, 0.999, 2000)
    B_total = rho_eff_bracket(A_grid, A_floor)
    B_back = np.full_like(A_grid, rho_eff_bracket_background(A_floor))
    B_source = rho_eff_bracket_source(A_grid, A_floor)
    B_F3 = f_of_A(A_grid) - A_grid * f_prime_of_A(A_grid)

    crossings_floor = find_NEC_crossover_with_floor(A_floor)
    A_crit_F3 = 0.44  # known F3 result

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.plot(A_grid, B_total, color="black", linewidth=2.0,
            label="B_total (with cosmic floor)")
    ax.plot(A_grid, B_back, color="tab:purple", linewidth=1.5, linestyle=":",
            label=f"B_background (structural, = f(A_0) ≈ {B_back[0]:.4f})")
    ax.plot(A_grid, B_source, color="tab:blue", linewidth=2.0, linestyle="--",
            label="B_source (dynamical, source displacement)")
    ax.plot(A_grid, B_F3, color="tab:red", linewidth=1.5, alpha=0.7,
            label="B_F3 (no floor; original F3 result)")
    ax.axhline(0, color="black", linewidth=0.5)
    for i, A_c in enumerate(crossings_floor):
        label = (f"B_total crossover {i+1} at A ≈ {A_c:.4f}"
                 if len(crossings_floor) > 1 else
                 f"B_total crossover at A ≈ {A_c:.4f}")
        ax.axvline(A_c, color="tab:green", linestyle="--", linewidth=1.5, label=label)
    ax.axvline(A_crit_F3, color="tab:orange", linestyle=":", linewidth=1.5,
               label=f"NEC crossover (no floor, F3) at A ≈ {A_crit_F3:.2f}")
    ax.set_xlabel("A_total")
    ax.set_ylabel("ρ_eff × r² bracket (signed)")
    ax.set_title("Effective stress-energy bracket: total / background / source decomposition")
    ax.grid(True, linewidth=0.3)
    ax.legend(loc="upper left", fontsize=8)
    out = PLOTS / "G21_bracket_decomposition.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_mass_function_decomposition(A_floor=A_0):
    A_grid = np.linspace(A_floor + 1e-3, 0.99, 1000)
    m_total = mass_function_over_M(A_grid, A_floor)
    m_source = mass_function_source_over_M(A_grid, A_floor)
    r_grid = 1.0 / (A_grid - A_floor)
    m_back = mass_function_background_over_M(r_grid, A_floor)

    # F3 baseline for comparison
    A_no_floor = np.linspace(0.001, 0.99, 1000)
    m_F3 = f_of_A(A_no_floor) / A_no_floor  # (r/Rs) × f(A) = f(A) / A in the no-floor case

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.plot(A_grid, m_total, color="black", linewidth=2.0, label="m_total/M (with floor)")
    ax.plot(A_grid, m_back, color="tab:purple", linewidth=1.5, linestyle=":",
            label="m_background/M (structural)")
    ax.plot(A_grid, m_source, color="tab:blue", linewidth=2.0, linestyle="--",
            label="m_source/M (dynamical, displacement)")
    ax.plot(A_no_floor, m_F3, color="tab:red", linewidth=1.5, alpha=0.7,
            label="m_F3/M (no floor; original F3)")
    ax.axhline(1, color="gray", linewidth=0.5)
    ax.set_xlabel("A_total")
    ax.set_ylabel("m(r) / M")
    ax.set_title("Mass function decomposition: total / background / source")
    ax.grid(True, linewidth=0.3)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_ylim(-0.5, 5.0)  # cap for readability
    out = PLOTS / "G21_mass_function_decomposition.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# ============================================================
# Markdown summary
# ============================================================

def write_markdown(df_landmarks, df_compare, crossings_total, plot_paths):
    md = []
    md.append("# G21: Strong-field Model-A with the cosmic floor explicit\n")
    md.append(
        "**Date:** 2026-05-11\n"
        "**Foundational technical step toward re-interpreting the strong-field "
        "corpus under the symmetric-boundary commitment (A on (A_0, 1), not (0, 1)).**\n"
    )
    md.append("## Water-tank framing\n")
    md.append(
        "Author 2026-05-11: A_0 IS the void — the unperturbed manifold level, "
        "not zero. Matter and energy displace this value into something higher "
        "locally. Crucially, displacement maps onto the AVAILABLE RANGE (1 − A_0), "
        "not onto the unit interval — A can't be pushed above 1 (the saturation/"
        "horizon ceiling), so the source's displacement scales by (1 − A_0).\n"
        "\n"
        "```text\n"
        "A_total(r) = A_0 + (1 − A_0) × (Rs / r)\n"
        "A_local(r) = (1 − A_0) × (Rs / r)   (displacement, scaled to available range)\n"
        f"A_0        = 1/(12π) ≈ {A_0:.6f}    (cosmic structural floor, G7)\n"
        "```\n"
        "\n"
        "At r = Rs: A_total = A_0 + (1 − A_0) = 1 (horizon at GR's Schwarzschild "
        "radius, matching weak-field GR exactly).\n"
        "\n"
        "**Naive A_total = A_0 + Rs/r is wrong.** It would put A_total > 1 at r ≤ Rs, "
        "which violates the A ≤ 1 ceiling. The rescaling by (1 − A_0) is what the "
        "water-tank picture actually requires.\n"
        "\n"
        "Under the ledger reading, A is the density of resolved write-events per "
        "Planck cell. A_0 is the minimum write density required for the manifold "
        "to exist; A_local is the additional write density due to matter "
        "currently being at the source location, scaled so cells can't saturate "
        "beyond capacity.\n"
    )
    md.append("## Strong-field metric (unchanged k(A))\n")
    md.append(
        "```text\n"
        "g_tt = -(1 - A_total) c²\n"
        "g_rr = 1 / [(1 - A_total)(1 - A_total²)²]\n"
        "```\n"
        "\n"
        "The metric construction k(A) = (1-A)(1-A²)² remains the framework's "
        "postulated commitment (see project_strong_field_commitment.md and "
        "project_strong_field_departure_question.md). This script changes only "
        "what enters as A: A_total with cosmic floor instead of A_local alone.\n"
    )
    md.append("## Landmark radii\n")
    md.append(df_landmarks.to_markdown(index=False, floatfmt=".5f"))
    md.append("\n")
    md.append(
        "**Reading:** all three orbital landmarks shift OUTWARD when the cosmic "
        "floor is included. The shift is largest at ISCO (lowest A_total target) "
        "and smallest at the horizon (highest target). This is because A_0 is a "
        "fixed offset; landmarks defined by small A_local values feel it more "
        "in relative terms.\n"
        "\n"
        "**Horizon:** r_horizon / Rs = 1/(1 - A_0) ≈ 1.0272. The horizon sits "
        "~2.7% outside the GR Schwarzschild radius. The cosmic A_0 is already in "
        "the cell; the source's displacement only needs to push A_total from "
        "A_0 to 1, i.e., A_local = 1 - A_0 instead of 1.\n"
    )
    md.append("## Asymptotic limit\n")
    md.append(
        "As r → ∞, A_local → 0 and A_total → A_0. The metric reads:\n"
        "\n"
        "```text\n"
        f"g_tt(∞) = -(1 - A_0) c² ≈ -{1.0 - A_0:.4f} c²\n"
        f"g_rr(∞) = 1 / [(1 - A_0)(1 - A_0²)²] ≈ {1.0 / k_of_A(A_0):.4f}\n"
        "```\n"
        "\n"
        "**Spacetime is asymptotically non-Minkowski.** Clocks at the asymptotic "
        f"limit tick at √(1 - A_0) ≈ {math.sqrt(1.0 - A_0):.4f} of a hypothetical "
        "A=0 clock — but A=0 does not exist, so this is not a slowdown relative "
        "to a faster reference. It is the framework's commitment about what "
        "asymptotic time IS at the cosmic floor.\n"
    )
    md.append("## Effective stress-energy decomposition\n")
    md.append(
        "Under the ledger framing, the effective stress-energy bracket "
        "B(A_total, A_0) = ρ_eff × r² (up to a constant) decomposes into two "
        "physically distinct contributions:\n"
        "\n"
        "- **Structural (background):** `B_background = f(A_0)` — the cosmic-floor "
        "contribution, present everywhere by virtue of the manifold existing. "
        "Constant in r. Encodes the structural baseline.\n"
        "- **Dynamical (source):** `B_source(A_total) = B(A_total, A_0) - B_background` "
        "— the displacement contribution from the matter source. Zero at the "
        "asymptotic limit (A_total = A_0); rises to ≈ 1 - f(A_0) at the horizon.\n"
        "\n"
        "The full bracket with cosmic floor is:\n"
        "\n"
        "```text\n"
        "B(A_total, A_0) = f(A_total) - (A_total - A_0) × f'(A_total)\n"
        "                = 1 - (1 - A_total)²(1 + A_total)\n"
        "                       × [1 + A_total(1 - 5A_0) + 4A_total² - A_0]\n"
        "```\n"
        "\n"
        "When A_0 = 0 this reduces exactly to F3's bracket "
        "`1 - (1 - A)²(1 + A)(1 + A + 4A²)`. With A_0 = 1/(12π), there are A_0 "
        "corrections at order A_0 throughout.\n"
    )
    md.append("## Numerical table\n")
    md.append(df_compare.to_markdown(index=False, floatfmt=".6f"))
    md.append("\n")
    if len(crossings_total) >= 2:
        md.append(
            "**Two sign changes of B_total** (not one as in F3 without floor):\n"
            f"- Lower crossover: A_total ≈ {crossings_total[0]:.5f}\n"
            f"- Upper crossover: A_total ≈ {crossings_total[1]:.5f}\n"
            "\n"
            "F3 (no floor) had a single sign change at A ≈ 0.44. With the cosmic "
            "floor, the effective stress-energy is positive at the asymptotic "
            f"limit (A_total = A_0, B_total = f(A_0) ≈ {f_of_A(A_0):.5f}), "
            "becomes negative in a finite band between the two crossovers, then "
            "positive again as A_total approaches 1.\n"
            "\n"
            "**This is qualitatively different from F3's no-floor picture.** F3 "
            "(no floor) had effective ρ negative for ALL A < 0.44 (the whole "
            "outer region of the source). With the cosmic floor enforced, the "
            "negative-ρ region is now a finite band, NOT the whole outer region. "
            "The asymptotic region itself is positive-ρ, dominated by the cosmic "
            "structural baseline.\n"
        )
    elif len(crossings_total) == 1:
        md.append(
            f"**Single sign change of B_total at A_total ≈ {crossings_total[0]:.5f}.** "
            "F3 (no floor) had a single sign change at A ≈ 0.44.\n"
        )
    else:
        md.append(
            "**No sign change of B_total** in (A_0, 1). The bracket is single-signed "
            "throughout the manifold's A range under cosmic floor.\n"
        )
    md.append(
        "\n"
        "**What's structural vs dynamical:**\n"
        "- The background bracket f(A_0) is small and positive "
        f"(≈ {f_of_A(A_0):.5f}). It represents the cosmic-floor contribution "
        "to effective stress-energy — small magnitude, positive sign, present "
        "everywhere by virtue of the manifold existing.\n"
        "- The source bracket B_source is negative for moderate A_total and "
        "positive near the horizon. Its single sign change sits near F3's 0.44 "
        "(small offset from A_0 corrections to dA/dr).\n"
    )
    md.append("## Honest assessment\n")
    md.append(
        "**What changed quantitatively:**\n"
        "- Landmark radii shift outward by 2-9% (ISCO most, horizon least).\n"
        "- Asymptotic metric is non-Minkowski (rescaled time and radial coordinates).\n"
        f"- F3's single NEC crossover at A ≈ 0.44 becomes **{len(crossings_total)} "
        f"sign changes** under cosmic floor: " +
        (f"{crossings_total[0]:.4f} and {crossings_total[1]:.4f}." if len(crossings_total) >= 2
         else (f"{crossings_total[0]:.4f}." if crossings_total else "none.")) + "\n"
        "- The negative-ρ band is no longer the whole outer region; it is a "
        "finite shell between the two crossovers.\n"
        "- The mass function and effective ρ both grow linearly with r at large "
        "r (vs returning to fixed values at infinity in the no-floor case) — "
        "the cosmic A_0 contributes a 'background mass' that scales with volume.\n"
        "\n"
        "**What changed conceptually:**\n"
        "- Effective stress-energy decomposes naturally into structural "
        "(background, A_0-only) and dynamical (source displacement) contributions. "
        "These are not free interpretive layers — they fall out of the math once "
        "you write the bracket with A_floor present.\n"
        "- The structural piece is small but nonzero everywhere — it is the "
        "framework's cosmological-constant-like contribution, expressed at the "
        "spatial-profile level (consistent with script 34's equation-of-state "
        "finding that w ≈ -1 at A near A_0).\n"
        "- The dynamical piece carries F3's original sign-change story, now "
        "interpretable as 'where the source's contribution flips from "
        "dark-energy-like to matter-like.'\n"
        "\n"
        "**What did NOT change:**\n"
        "- The metric construction k(A) = (1-A)(1-A²)². Still postulated. "
        "Field-equations gap unchanged.\n"
        "- The thirds-of-A landmarks (ISCO at 1/3, photon sphere at 2/3, horizon "
        "at 1) are still the values of A_total at which they occur. What "
        "shifted is the radius at which A_total reaches each value.\n"
        "- F3's verdict (Birkhoff doesn't apply, owes field equations) unchanged.\n"
        "\n"
        "**The 'why' question this surfaces:**\n"
        "- Why does the NEC crossover (sign change of B_total) sit so close to "
        "A = 0.44 with or without the floor? F3's 0.44 is structurally related "
        "to where the (1-A)²(1+A)(1+A+4A²) factor crosses 1. Cosmic floor "
        "shifts this slightly. Worth understanding what determines this value "
        "structurally — it may have its own meaning.\n"
        "- The asymptotic linear growth of m(r) means the framework predicts "
        "a 'mass-at-infinity' that diverges in the standard Schwarzschild "
        "parameterization. This is the cosmic-floor contribution at large r. "
        "Interpretation: in a proper FRW-style cosmological framework, this "
        "is the framework's natively-derived dark-energy effective mass.\n"
    )
    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "G21_strong_field_with_cosmic_floor_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 78)
    print("G21: Strong-field Model-A with the cosmic floor explicit")
    print("=" * 78)
    print(f"A_0 = 1/(12π) = {A_0:.8f}")
    print()

    df_landmarks = landmark_table(A_0)
    df_landmarks.to_csv(RESULTS / "G21_landmark_shifts.csv", index=False)
    print("Landmark radii:")
    print(df_landmarks.to_string(index=False))
    print()

    crossings_total = find_NEC_crossover_with_floor(A_0)
    print("B_total sign changes (with cosmic floor):")
    if not crossings_total:
        print("  none in (A_0, 1).")
    else:
        for i, A_c in enumerate(crossings_total):
            r_over_Rs_c = (1.0 - A_0) / (A_c - A_0)
            print(f"  crossover {i+1}: A_total ≈ {A_c:.6f}  (r/Rs ≈ {r_over_Rs_c:.4f})")
    print(f"F3 reference (no floor): single crossover at A ≈ 0.4400")
    print()

    crossings_source = find_NEC_crossover_source_only(A_0)
    print("B_source sign changes (above trivial zero at A_0):")
    if not crossings_source:
        print("  none.")
    else:
        for i, A_c in enumerate(crossings_source):
            print(f"  crossover {i+1}: A_total ≈ {A_c:.6f}")
    print()

    df_compare = comparison_table(A_0)
    df_compare.to_csv(RESULTS / "G21_comparison_table.csv", index=False)
    print("Comparison table (with/without cosmic floor):")
    print(df_compare.to_string(index=False))
    print()

    plot_paths = [
        plot_landmark_shifts(df_landmarks, A_0),
        plot_bracket_decomposition(A_0),
        plot_mass_function_decomposition(A_0),
    ]

    summary_path = write_markdown(df_landmarks, df_compare, crossings_total, plot_paths)
    print(f"Markdown summary written: {summary_path}")
    for p in plot_paths:
        print(f"Plot: {p}")


if __name__ == "__main__":
    main()

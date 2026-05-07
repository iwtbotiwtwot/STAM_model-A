#!/usr/bin/env python3
"""
Q9_resolution_temperature_derivation.py

STAM Model-A: Path 1 derivation of boundary temperature from
resolution-event statistics.
Author: Sean Brady / STAM Model-A continuation

Purpose:
    The Q8 script postulated k_B T = (1 / (4 pi)) * hbar * c * |grad A|.
    This Q9 script asks how much of that follows from STAM-native premises
    alone, without importing the standard surface-gravity-to-temperature
    relation.

Premises used (each labeled STAM-native or basic-quantum):
    P1 (basic-quantum): vacuum is full of unresolved A-fluctuations of
        energy E with lifetime tau(E) = hbar / E. (Time-energy uncertainty.)
    P2 (basic-quantum): the bulk density of fluctuation states per unit
        volume per unit energy is rho(E) = E^2 / (2 pi^2 hbar^3 c^3).
        (Standard 3D massless scalar density of states. We use it because
        STAM has not committed to a different density of states; an
        alternative STAM-specific density would change the prefactor but
        not the scaling.)
    P3 (STAM-native): the resolution rate per fluctuation in a region
        with gradient |grad A| is r_res = c |grad A|. This is the only
        rate scale that can be built from c, |grad A|, and no extra
        dimensional inputs.
    P4 (STAM-native): at an A=1 boundary, only outward-oriented
        resolutions complete. Inward resolutions have no committed
        substrate to anchor to.

What this script tries to extract:
    (a) The temperature SCALE that comes out of premises P1-P4 with no
        further input.
    (b) The SHAPE of the resolved-flux spectrum predicted by P1-P4, and
        how it compares to a Planckian at the Q8 temperature.
    (c) Whether the 1/(4 pi) prefactor falls out cleanly, or whether an
        additional physical principle is required.

Honest reporting:
    The script reports what it actually finds. If the prefactor does not
    come out as 1/(4 pi), the residual mismatch is reported as the
    O(1) gap that Path 2 (SU imaginary-time argument) would need to close.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

C = 2.99792458e8
HBAR = 1.054571817e-34
KB = 1.380649e-23
G = 6.67430e-11
M_SUN = 1.98847e30

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# --- Premises as functions ---

def fluctuation_lifetime(E: np.ndarray) -> np.ndarray:
    """P1: tau(E) = hbar / E."""
    return HBAR / E


def density_of_states(E: np.ndarray) -> np.ndarray:
    """P2: standard 3D massless density of states per unit volume per unit energy."""
    return (E ** 2) / (2.0 * (math.pi ** 2) * (HBAR ** 3) * (C ** 3))


def resolution_rate(grad_A: float) -> float:
    """P3: only rate scale that can be built from c and |grad A|."""
    return C * grad_A


def resolution_probability_per_lifetime(E: np.ndarray, grad_A: float) -> np.ndarray:
    """Probability that a fluctuation of energy E commits during its lifetime tau(E).

    P_res = r_res * tau = (c |grad A|) * (hbar / E) = hbar c |grad A| / E.

    Saturates at 1: a probability cannot exceed unity. When P_res > 1, the
    fluctuation almost certainly commits during its lifetime.
    """
    P = HBAR * C * grad_A / E
    return np.minimum(P, 1.0)


def commit_rate_per_unit_volume_per_unit_energy(E: np.ndarray, grad_A: float) -> np.ndarray:
    """Number of resolved fluctuations per unit volume per unit time per unit energy.

    = (density of states at E) * (probability each commits per lifetime)
        / (lifetime)
    = rho(E) * P_res(E) / tau(E)
    = rho(E) * (P_res(E) * E / hbar)
    """
    rho = density_of_states(E)
    P = resolution_probability_per_lifetime(E, grad_A)
    return rho * P * E / HBAR


# --- Derived results ---

def characteristic_resolved_energy(grad_A: float) -> float:
    """The crossover energy where P_res = 1: E_cross = hbar c |grad A|."""
    return HBAR * C * grad_A


def predicted_spectrum_normalised(E_grid: np.ndarray, grad_A: float) -> np.ndarray:
    """Outgoing flux spectrum, normalized by its peak for shape comparison."""
    s = commit_rate_per_unit_volume_per_unit_energy(E_grid, grad_A)
    return s / np.max(s)


def planckian_normalised(E_grid: np.ndarray, T: float) -> np.ndarray:
    """Planck distribution dN/dE proportional to E^2 / (exp(E/kT) - 1), normalized by peak."""
    x = E_grid / (KB * T)
    expm1 = np.expm1(x)
    safe = np.where(expm1 > 0, expm1, np.nan)
    p = (E_grid ** 2) / safe
    return p / np.nanmax(p)


def find_knee_energy(E_grid: np.ndarray, grad_A: float) -> float:
    """Locate the knee in S(E) — the energy where P_res transitions from saturated to falling.

    By construction this is at E = hbar c |grad A|. Verifying it numerically as a sanity check.
    """
    P = HBAR * C * grad_A / E_grid
    idx = int(np.argmin(np.abs(P - 1.0)))
    return float(E_grid[idx])


def integrate_total_emission_rate(E_grid: np.ndarray, spectrum: np.ndarray) -> float:
    """Integrate spectrum over E_grid. If the integral depends strongly on the upper bound,
    the spectrum has no UV cutoff in its native form."""
    return float(np.trapezoid(spectrum, E_grid))


# --- Comparison cases ---

def case_solar_mass_bh() -> dict:
    M = 1.0 * M_SUN
    Rs = 2.0 * G * M / (C ** 2)
    grad_A = 1.0 / Rs
    T_q8 = HBAR * C * grad_A / (4.0 * math.pi * KB)  # the Q8 postulate value
    return {
        "label": "Solar-mass black hole",
        "grad_A": grad_A,
        "T_Q8": T_q8,
        "Rs_m": Rs,
    }


def case_micro_bh() -> dict:
    M = 1e15  # kg, asteroid-mass primordial BH
    Rs = 2.0 * G * M / (C ** 2)
    grad_A = 1.0 / Rs
    T_q8 = HBAR * C * grad_A / (4.0 * math.pi * KB)
    return {
        "label": "10^15 kg primordial BH",
        "grad_A": grad_A,
        "T_Q8": T_q8,
        "Rs_m": Rs,
    }


def analyse_case(case: dict) -> dict:
    grad_A = case["grad_A"]
    T_Q8 = case["T_Q8"]

    E_cross = characteristic_resolved_energy(grad_A)
    E_grid_short = np.logspace(np.log10(E_cross * 1e-3), np.log10(E_cross * 1e3), 4000)
    E_grid_long = np.logspace(np.log10(E_cross * 1e-3), np.log10(E_cross * 1e6), 4000)

    s_short = commit_rate_per_unit_volume_per_unit_energy(E_grid_short, grad_A)
    s_long = commit_rate_per_unit_volume_per_unit_energy(E_grid_long, grad_A)

    knee_E = find_knee_energy(E_grid_short, grad_A)
    integral_short = integrate_total_emission_rate(E_grid_short, s_short)
    integral_long = integrate_total_emission_rate(E_grid_long, s_long)

    return {
        **case,
        "E_cross_J": E_cross,
        "knee_E_J": knee_E,
        "knee_over_E_cross": knee_E / E_cross,
        "knee_over_kT_Q8": knee_E / (KB * T_Q8),
        "T_Q8_K": T_Q8,
        "integral_grid_to_1e3_E_cross": integral_short,
        "integral_grid_to_1e6_E_cross": integral_long,
        "integral_ratio_long_over_short": integral_long / integral_short,
        "_E_grid": E_grid_short,
        "_s_pred": s_short,
    }


# --- Plots ---

def plot_spectrum_comparison(result: dict) -> Path:
    grad_A = result["grad_A"]
    T_Q8 = result["T_Q8"]
    E_grid = result["_E_grid"]
    s_pred = result["_s_pred"]

    s_pred_n = s_pred / np.max(s_pred)
    p_q8 = planckian_normalised(E_grid, T_Q8)

    plt.figure(figsize=(8, 5))
    plt.plot(E_grid, s_pred_n, label="Path 1 prediction (resolution statistics)", linewidth=2.2)
    plt.plot(E_grid, p_q8, "--", label=f"Planckian at T_Q8 = {T_Q8:.3g} K", linewidth=1.8)

    plt.axvline(KB * T_Q8, color="gray", linestyle=":", alpha=0.7,
                label=f"k_B T_Q8 = {KB * T_Q8:.3g} J")
    plt.axvline(HBAR * C * grad_A, color="purple", linestyle=":", alpha=0.7,
                label=f"hbar c |grad A| = {HBAR * C * grad_A:.3g} J")

    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(1e-6, 2)
    plt.xlabel("Fluctuation energy E (J)")
    plt.ylabel("Spectrum (normalized to peak)")
    plt.title(f"Path 1 spectrum vs Planckian — {result['label']}")
    plt.legend(fontsize=8)
    plt.grid(True, which="both", linewidth=0.3)

    safe = result["label"].replace(" ", "_").replace("^", "").replace("/", "_")
    out = PLOTS / f"Q9_spectrum_{safe}.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_knee_scaling_check() -> Path:
    """Confirm that the Path 1 knee tracks |grad A|, even though the spectrum is not Planckian."""
    grad_A_values = np.logspace(-30, 5, 40)
    knee_energies = HBAR * C * grad_A_values  # by construction
    kT_Q8_values = HBAR * C * grad_A_values / (4.0 * math.pi)

    plt.figure(figsize=(8, 5))
    plt.plot(grad_A_values, knee_energies, "o", label="Path 1 knee energy = hbar c |grad A|", markersize=5)
    plt.plot(grad_A_values, kT_Q8_values, "--",
             label="Q8 thermal scale k_B T_Q8 = hbar c |grad A| / (4 pi)", linewidth=1.8)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("|grad A| (1/m)")
    plt.ylabel("Energy (J)")
    plt.title("Path 1 knee and Q8 thermal scale both track |grad A|, separated by factor 4 pi")
    plt.legend()
    plt.grid(True, which="both", linewidth=0.3)

    out = PLOTS / "Q9_knee_vs_kT_scaling.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- Markdown summary ---

def write_markdown(results: list[dict], plot_paths: list[Path]) -> Path:
    md = []
    md.append("# STAM Model-A Q9: Path 1 Derivation of Boundary Temperature\n")
    md.append("## Purpose\n")
    md.append(
        "The Q8 script postulated `k_B T = hbar c |grad A| / (4 pi)` and showed it reproduces "
        "Schwarzschild, Unruh, and de Sitter temperatures. Q9 asks how much of that follows "
        "from STAM-native premises *alone*, without importing the standard "
        "surface-gravity-to-temperature relation.\n"
    )
    md.append("## Premises\n")
    md.append(
        "- **P1 (basic-quantum):** Vacuum holds unresolved A-fluctuations of energy E with "
        "lifetime tau = hbar / E.\n"
        "- **P2 (basic-quantum):** Bulk density of fluctuation states is rho(E) = "
        "E^2 / (2 pi^2 hbar^3 c^3) per unit volume per unit energy.\n"
        "- **P3 (STAM-native):** Resolution rate per fluctuation = c |grad A|. This is the "
        "*only* rate scale that can be built from c, |grad A|, and no extra dimensional inputs.\n"
        "- **P4 (STAM-native):** At an A=1 boundary, only outward-oriented resolutions complete "
        "(no committed substrate inward to anchor to).\n"
    )
    md.append("## Predicted outgoing-flux spectrum (Path 1)\n")
    md.append(
        "Combining P1-P4: rate of resolved fluctuations per unit volume per unit time per "
        "unit energy is\n"
        "```text\n"
        "S(E) = rho(E) * min(1, hbar c |grad A| / E) * E / hbar\n"
        "     ~ E^3       for E < hbar c |grad A|   (P_res saturated)\n"
        "     ~ E^2 |grad A|   for E > hbar c |grad A|   (P_res falling as 1/E)\n"
        "```\n"
        "Both branches are POWER LAWS that grow with E. There is no exponential cutoff. "
        "The Planckian distribution at T_Q8 = hbar c |grad A| / (4 pi k_B) has a peak near "
        "k_B T_Q8 and falls exponentially above; the Path 1 spectrum has neither feature.\n"
    )
    md.append("## Case-by-case results\n")
    df_rows = []
    for r in results:
        df_rows.append({
            "case": r["label"],
            "grad_A_per_m": r["grad_A"],
            "knee_E_J": r["knee_E_J"],
            "knee / (hbar c |grad A|)": r["knee_over_E_cross"],
            "knee / (k_B T_Q8)": r["knee_over_kT_Q8"],
            "integral up to 1e3 E_cross": r["integral_grid_to_1e3_E_cross"],
            "integral up to 1e6 E_cross": r["integral_grid_to_1e6_E_cross"],
            "ratio (1e6 / 1e3)": r["integral_ratio_long_over_short"],
        })
    df = pd.DataFrame(df_rows)
    md.append(df.to_markdown(index=False, floatfmt=".6g"))
    md.append(
        "\n\nThe `ratio (1e6 / 1e3)` column shows that extending the integration upper limit "
        "by a factor of 1000 changes the integrated emission by orders of magnitude. This "
        "is the diagnostic that the spectrum has no UV cutoff in Path 1's premises — the "
        "total emission rate is not finite without additional input.\n"
    )
    md.append("## What fell out, what didn't — honest report\n")
    md.append(
        "- **Characteristic energy scale:** YES. The Path 1 spectrum has a structural feature "
        "(the knee where probability saturates) at E = hbar c |grad A|. This matches the "
        "Q8 thermal scale up to a factor of 4 pi; both quantities track |grad A| with the "
        "same dependence.\n"
        "- **Thermal Planckian shape:** NO. The Path 1 spectrum is two power laws (E^3 below "
        "the knee, E^2 above), with no exponential suppression. It is not a thermal "
        "distribution.\n"
        "- **A meaningful temperature:** NOT FROM PATH 1 ALONE. Without an exponential UV "
        "cutoff, no finite mean energy and no well-defined temperature can be extracted. The "
        "first-moment estimator depends entirely on where the integration grid is truncated, "
        "which means it has no physical meaning.\n"
        "- **Prefactor 1/(4 pi):** NO. Path 1 cannot fix the prefactor because it cannot fix "
        "the temperature.\n"
    )
    md.append("## What Path 1 actually accomplished\n")
    md.append(
        "Path 1 successfully established two pieces of structure from STAM-native premises:\n"
        "1. A characteristic energy scale `hbar c |grad A|` exists at any A=1 boundary, with "
        "no help from QFT-on-curved-spacetime.\n"
        "2. The spectrum has a knee at that scale, where the resolution probability saturates.\n"
    )
    md.append(
        "What Path 1 *fails* to establish is the exponential UV cutoff that would make the "
        "spectrum thermal. This is exactly the structural feature that Path 2 would need "
        "to provide. The natural candidate is the strong-field redshift / clock-rate relation: "
        "fluctuations with very high local energy have very small asymptotic energy, and a "
        "proper accounting of the redshift between local and asymptotic frequencies should "
        "produce the missing exponential cutoff via the SU logarithmic divergence near the "
        "boundary.\n"
    )
    md.append("## Plain-language summary\n")
    md.append(
        "STAM-native ingredients (uncertainty + resolution rate + asymmetric boundary) are "
        "enough to identify the right *energy scale* for thermal emission at an A=1 boundary, "
        "but they are NOT enough to produce a thermal spectrum or to fix the 1/(4 pi) "
        "prefactor. The missing ingredient is whatever produces the exponential cutoff at "
        "high energies — most likely the strong-field redshift relation, which would "
        "naturally connect to Path 2 via the SU integral.\n"
    )
    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "Q9_resolution_temperature_derivation_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    cases = [case_solar_mass_bh(), case_micro_bh()]
    results = [analyse_case(c) for c in cases]

    plot_paths = [plot_spectrum_comparison(r) for r in results]
    plot_paths.append(plot_knee_scaling_check())

    summary = write_markdown(results, plot_paths)

    print("STAM Model-A Q9: Path 1 Derivation of Boundary Temperature")
    print("=" * 64)
    for r in results:
        print(f"\n{r['label']}:")
        print(f"  |grad A|              = {r['grad_A']:.6g} 1/m")
        print(f"  knee energy           = {r['knee_E_J']:.6g} J")
        print(f"  knee / (hbar c |grad A|) = {r['knee_over_E_cross']:.6g}  (should be ~1)")
        print(f"  knee / (k_B T_Q8)     = {r['knee_over_kT_Q8']:.6g}  (should be ~4 pi = {4*math.pi:.4g})")
        print(f"  integral [E_c*1e-3, E_c*1e3] = {r['integral_grid_to_1e3_E_cross']:.6g}")
        print(f"  integral [E_c*1e-3, E_c*1e6] = {r['integral_grid_to_1e6_E_cross']:.6g}")
        print(f"  ratio (extending UV by 1000x) = {r['integral_ratio_long_over_short']:.6g}")
    print("\nHonest finding: Path 1 spectrum has no UV cutoff, so no meaningful T")
    print("can be extracted from it alone. The knee is at hbar c |grad A| (matches")
    print("Q8 scale up to 4 pi), but the spectrum is power-law, not Planckian.")
    print("\nFiles written:")
    print(f"- {summary}")
    for p in plot_paths:
        print(f"- {p}")


if __name__ == "__main__":
    main()

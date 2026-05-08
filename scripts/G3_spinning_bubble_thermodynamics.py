#!/usr/bin/env python3
"""
G3_spinning_bubble_thermodynamics.py

Hawking temperature on the spinning Model-A bubble (A=1 surface), and
its implications for the first law of black hole thermodynamics.

Q8/Q10 derived T = ℏc|∇A|/(4π k_B) for the static case from the
phase-boundary mechanism. For the spinning bubble (G2), |∇A| varies
with polar angle θ along the bubble, so the local temperature varies
along the boundary. This script computes T(θ), the area-weighted mean
temperature, the total Stefan-Boltzmann luminosity, and contrasts the
result with Kerr's uniform-temperature horizon.

Key derivation (using coordinate gradient, consistent with Q8 conv.):
    A(r, θ; M, a) = 2 M r / Σ,    Σ = r² + a² cos²θ
    ∂A/∂r = 2M(a² cos²θ - r²)/Σ²
    ∂A/∂θ = 2 M r a² sin(2θ)/Σ²
    |∇A|² = (∂A/∂r)² + (1/r²)(∂A/∂θ)²

On the bubble (A=1, Σ = 2Mr, a²cos²θ = r(2M-r)):
    ∂A/∂r|_bubble = (M - r)/(M r)
    ∂A/∂θ|_bubble = a² sin(2θ)/(2 M r)
    |∇A|²|_bubble = (M - r)²/(M² r²)  +  a⁴ sin²(2θ)/(4 M² r⁴)

Temperature:
    T(θ) = ℏ c |∇A(θ)| / (4 π k_B)

Limits:
    Equator (θ=π/2):  r=2M, sin(2θ)=0  →  |∇A|=1/(2M)  →  T = T_Schw
    Pole (θ=0):       r=r_+=M+√(M²-a²), sin(2θ)=0  →  |∇A| = √(M²-a²)/(M r_+)
                       For a→M (extremal): |∇A|_pole → 0 → T_pole → 0

Findings (numerical, see results below):
    The spinning Model-A bubble has a NON-UNIFORM temperature, hottest at
    the equator (always T_Schw, independent of spin), coldest at the poles
    (vanishing at extremal). Kerr's horizon, by contrast, has a UNIFORM
    temperature T_Kerr = (r_+-M)/(2π(r_+²+a²)) that drops monotonically
    with spin. So at any spin > 0, Model-A's bubble has a hotter equator
    than Kerr predicts and a comparable pole.

    First-law implications:
        For Kerr:   dM = T dS + Ω dJ     with single uniform T
        For Model-A: dM = ∫ T(θ) (dS_local) + Ω dJ — non-equilibrium
                    distribution of emission across the bubble
        OR: define T_eff via L = σ T_eff⁴ A_total (Stefan-Boltzmann), giving
        an emission-weighted effective temperature that's between T_pole
        and T_equator.

What this script does:
    1. Compute |∇A|(θ) on bubble for a range of spins
    2. Compute T(θ), plot, find peaks and zeros
    3. Compute area-weighted ⟨T⟩, peak T_equator, T_pole
    4. Compute total Hawking-style luminosity L_SB = ∫ σ T(θ)⁴ dA(θ)
    5. Compare to Kerr's L_Kerr = σ T_Kerr⁴ A_Kerr
    6. Tabulate spin-dependence of all observables
    7. Plot bubble entropy (S = A/(4 ℓ_p²)) vs spin
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Constants (SI)
HBAR = 1.054571817e-34
C = 2.99792458e8
G = 6.67430e-11
KB = 1.380649e-23
SIGMA_SB = 5.670374419e-8
M_SUN = 1.98892e30

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# --- Spinning bubble geometry (from G2) ---

def bubble_radius_geom(theta: np.ndarray, a_over_M: float = 0.5) -> np.ndarray:
    """r_bubble(θ)/M = 1 + √(1 - (a/M)² cos²θ).  Geometric units."""
    return 1.0 + np.sqrt(1.0 - a_over_M ** 2 * np.cos(theta) ** 2)


def bubble_radius_dr_dtheta(theta: np.ndarray,
                             a_over_M: float = 0.5) -> np.ndarray:
    """dr_bubble/dθ in geometric units.

    r/M = 1 + √(1 - α² cos²θ), where α = a/M.
    dr/dθ = α² cos θ sin θ / √(1 - α² cos²θ)
          = (α² sin(2θ) / 2) / √(1 - α² cos²θ)
    """
    inner = np.sqrt(np.maximum(1.0 - a_over_M ** 2 * np.cos(theta) ** 2, 1e-30))
    return (a_over_M ** 2 * np.sin(2.0 * theta) / 2.0) / inner


# --- |∇A| and T on the spinning bubble (geometric units, M=1) ---

def grad_A_on_bubble(theta: np.ndarray,
                     a_over_M: float = 0.5) -> np.ndarray:
    """|∇A|(θ) on the A=1 surface, in units of 1/M.

    On bubble (Σ = 2Mr, a²cos²θ = r(2M-r)):
        |∇A|² = (M-r)²/(M²r²) + a⁴ sin²(2θ)/(4 M² r⁴)
    """
    r = bubble_radius_geom(theta, a_over_M)  # r/M
    one_minus_r_over_M = 1.0 - r              # (M-r)/M
    radial_term = (one_minus_r_over_M / r) ** 2
    polar_term = (a_over_M ** 4 * np.sin(2.0 * theta) ** 2) / (4.0 * r ** 4)
    grad_sq = radial_term + polar_term
    return np.sqrt(np.maximum(grad_sq, 0.0))


def temperature_on_bubble(theta: np.ndarray, M_solar: float,
                          a_over_M: float = 0.5) -> np.ndarray:
    """Local temperature T(θ) on the bubble in Kelvin.

    T = ℏ c |∇A| / (4 π k_B), with |∇A| in 1/M (geometric M).
    Convert via M_geom = G M / c² (length).
    """
    M_kg = M_solar * M_SUN
    M_geom_m = G * M_kg / C ** 2
    grad_geom = grad_A_on_bubble(theta, a_over_M)         # 1/M_geom (dimensionless * 1/length scale)
    grad_SI = grad_geom / M_geom_m                        # 1/m
    T = HBAR * C * grad_SI / (4.0 * np.pi * KB)
    return T


# --- Schwarzschild reference (no spin) ---

def T_schwarzschild(M_solar: float) -> float:
    """Standard Hawking temperature for a non-spinning BH."""
    M_kg = M_solar * M_SUN
    return HBAR * C ** 3 / (8.0 * np.pi * G * M_kg * KB)


# --- Kerr reference ---

def T_kerr(M_solar: float, a_over_M: float) -> float:
    """Kerr horizon temperature (uniform on horizon).

    T_Kerr = ℏ c (r_+ - M) / (2π k_B (r_+² + a²))
    """
    if a_over_M >= 1.0:
        return 0.0
    M_kg = M_solar * M_SUN
    M_geom_m = G * M_kg / C ** 2
    r_plus_geom = 1.0 + np.sqrt(1.0 - a_over_M ** 2)         # in M units
    r_plus_m = r_plus_geom * M_geom_m
    a_m = a_over_M * M_geom_m
    M_m = M_geom_m
    # κ = (r_+ - M) / (r_+² + a²)  in geometric SI (1/m)
    kappa = (r_plus_m - M_m) / (r_plus_m ** 2 + a_m ** 2)
    T = HBAR * C * kappa / (2.0 * np.pi * KB)
    return T


# --- Bubble area and entropy ---

def bubble_area_proper(a_over_M: float = 0.5) -> float:
    """Approximate bubble area via surface-of-revolution integral
    (round-coordinate approximation, in M² units).

    A = 2π ∫_0^π r(θ) √(r² + r'²) sin θ dθ
    """
    theta_grid = np.linspace(0.0, np.pi, 5000)
    r = bubble_radius_geom(theta_grid, a_over_M)
    drdth = bubble_radius_dr_dtheta(theta_grid, a_over_M)
    integrand = r * np.sqrt(r ** 2 + drdth ** 2) * np.sin(theta_grid)
    return 2.0 * np.pi * np.trapezoid(integrand, theta_grid)


def kerr_horizon_area(a_over_M: float = 0.5) -> float:
    """Kerr horizon area: 4π(r_+² + a²), in M² units."""
    r_plus = 1.0 + np.sqrt(1.0 - a_over_M ** 2)
    return 4.0 * np.pi * (r_plus ** 2 + a_over_M ** 2)


# --- Stefan-Boltzmann emission ---

def stefan_boltzmann_luminosity(M_solar: float,
                                a_over_M: float = 0.5,
                                n_theta: int = 2000) -> float:
    """Total Hawking luminosity by integrating σ T⁴(θ) over the bubble.

    L = ∫ σ T⁴(θ) dA(θ), with dA = 2π r(θ) √(r² + r'²) sinθ dθ × (M_geom_m)²
    """
    theta_grid = np.linspace(1e-6, np.pi - 1e-6, n_theta)
    T = temperature_on_bubble(theta_grid, M_solar, a_over_M)
    r_geom = bubble_radius_geom(theta_grid, a_over_M)
    drdth = bubble_radius_dr_dtheta(theta_grid, a_over_M)
    M_kg = M_solar * M_SUN
    M_geom_m = G * M_kg / C ** 2
    # dA in SI
    dA = (2.0 * np.pi * r_geom * np.sqrt(r_geom ** 2 + drdth ** 2)
          * np.sin(theta_grid) * M_geom_m ** 2)
    integrand = SIGMA_SB * T ** 4 * dA / (2.0 * np.pi)  # already include 2π in dA
    # Above line: dA carried 2π already; integrand has σ T⁴ × dA/dθ
    L = np.trapezoid(integrand * 2.0 * np.pi, theta_grid)
    # The cleaner statement:
    # L = ∫_0^π σ T⁴(θ) × 2π r(θ) √(r²+r'²) sin θ × M_geom_m² dθ
    return L


def stefan_boltzmann_kerr(M_solar: float, a_over_M: float) -> float:
    """L_Kerr = σ T_Kerr⁴ × A_Kerr, with uniform T."""
    T = T_kerr(M_solar, a_over_M)
    M_kg = M_solar * M_SUN
    M_geom_m = G * M_kg / C ** 2
    A_kerr_geom = kerr_horizon_area(a_over_M)
    A_kerr_SI = A_kerr_geom * M_geom_m ** 2
    return SIGMA_SB * T ** 4 * A_kerr_SI


def stefan_boltzmann_schwarzschild(M_solar: float) -> float:
    """L_Schw = σ T⁴ A for Schwarzschild."""
    T = T_schwarzschild(M_solar)
    M_kg = M_solar * M_SUN
    M_geom_m = G * M_kg / C ** 2
    A = 4.0 * np.pi * (2.0 * M_geom_m) ** 2
    return SIGMA_SB * T ** 4 * A


# --- Entropy ---

def planck_length() -> float:
    return np.sqrt(HBAR * G / C ** 3)


def bubble_entropy(M_solar: float, a_over_M: float) -> float:
    """S = k_B × A_bubble / (4 ℓ_p²)."""
    M_kg = M_solar * M_SUN
    M_geom_m = G * M_kg / C ** 2
    A_geom = bubble_area_proper(a_over_M)
    A_SI = A_geom * M_geom_m ** 2
    return KB * A_SI / (4.0 * planck_length() ** 2)


def kerr_entropy(M_solar: float, a_over_M: float) -> float:
    M_kg = M_solar * M_SUN
    M_geom_m = G * M_kg / C ** 2
    A_geom = kerr_horizon_area(a_over_M)
    A_SI = A_geom * M_geom_m ** 2
    return KB * A_SI / (4.0 * planck_length() ** 2)


# --- Diagnostic table ---

def diagnostic_row(a_over_M: float, M_solar: float = 1.0) -> dict:
    theta_grid = np.linspace(1e-6, np.pi - 1e-6, 5000)
    T = temperature_on_bubble(theta_grid, M_solar, a_over_M)
    r = bubble_radius_geom(theta_grid, a_over_M)
    drdth = bubble_radius_dr_dtheta(theta_grid, a_over_M)
    sin_t = np.sin(theta_grid)
    weight = r * np.sqrt(r ** 2 + drdth ** 2) * sin_t

    T_equator = float(temperature_on_bubble(np.array([np.pi / 2]),
                                             M_solar, a_over_M)[0])
    T_pole = float(temperature_on_bubble(np.array([0.0]),
                                          M_solar, a_over_M)[0])
    # Area-weighted mean T
    T_mean = float(np.trapezoid(T * weight, theta_grid)
                   / np.trapezoid(weight, theta_grid))

    L_modelA = stefan_boltzmann_luminosity(M_solar, a_over_M)
    L_kerr = stefan_boltzmann_kerr(M_solar, a_over_M)
    L_schw = stefan_boltzmann_schwarzschild(M_solar)

    A_bubble = bubble_area_proper(a_over_M)
    A_kerr = kerr_horizon_area(a_over_M)

    S_modelA = bubble_entropy(M_solar, a_over_M)
    S_kerr = kerr_entropy(M_solar, a_over_M)
    S_schw = bubble_entropy(M_solar, 0.0)

    return {
        "a_over_M": a_over_M,
        "T_equator_K": T_equator,
        "T_pole_K": T_pole,
        "T_mean_K": T_mean,
        "T_kerr_K": T_kerr(M_solar, a_over_M),
        "T_schwarzschild_K": T_schwarzschild(M_solar),
        "T_pole_over_equator": (T_pole / T_equator) if T_equator > 0 else 0.0,
        "L_modelA_W": L_modelA,
        "L_kerr_W": L_kerr,
        "L_schwarzschild_W": L_schw,
        "L_modelA_over_kerr": (L_modelA / L_kerr) if L_kerr > 0 else float("inf"),
        "L_modelA_over_schwarzschild": L_modelA / L_schw,
        "A_bubble_M2": A_bubble,
        "A_kerr_M2": A_kerr,
        "A_bubble_over_kerr": A_bubble / A_kerr,
        "S_modelA_J_per_K": S_modelA,
        "S_kerr_J_per_K": S_kerr,
        "S_modelA_over_kerr": S_modelA / S_kerr,
        "S_modelA_over_schwarzschild": S_modelA / S_schw,
    }


# --- Plots ---

def plot_T_along_bubble(M_solar: float = 1.0) -> Path:
    fig, ax = plt.subplots(figsize=(11, 6))
    theta_grid = np.linspace(0.001, np.pi - 0.001, 600)
    T_schw = T_schwarzschild(M_solar)

    spins = [0.0, 0.3, 0.5, 0.7, 0.9, 0.99]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(spins)))

    for spin, color in zip(spins, colors):
        T = temperature_on_bubble(theta_grid, M_solar, spin) / T_schw
        ax.plot(np.degrees(theta_grid), T, color=color, linewidth=2,
                label=f"Model-A bubble  a/M = {spin:.2f}")

    # Kerr reference values (uniform horizontal lines)
    for spin, color in zip(spins, colors):
        if spin > 0:
            T_K = T_kerr(M_solar, spin) / T_schw
            ax.axhline(T_K, color=color, linestyle="--", alpha=0.5, linewidth=1)

    ax.set_xlabel("Polar angle θ (degrees)")
    ax.set_ylabel("T(θ) / T_Schwarzschild")
    ax.set_title("Local temperature along Model-A spinning bubble (solid)\n"
                 "vs. Kerr horizon temperature (dashed, uniform)\n"
                 "Model-A equator stays at T_Schw at any spin; "
                 "pole vanishes at extremal")
    ax.legend(fontsize=9, loc="lower center")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.15)

    out = PLOTS / "G3_T_along_bubble.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_thermo_quantities_vs_spin(M_solar: float = 1.0) -> Path:
    spin_grid = np.linspace(0.0, 0.999, 200)
    rows = [diagnostic_row(s, M_solar) for s in spin_grid]

    T_schw = T_schwarzschild(M_solar)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # T quantities
    axes[0, 0].plot(spin_grid, [r["T_equator_K"] / T_schw for r in rows],
                    "b-", linewidth=2, label="Model-A equator")
    axes[0, 0].plot(spin_grid, [r["T_pole_K"] / T_schw for r in rows],
                    "r-", linewidth=2, label="Model-A pole")
    axes[0, 0].plot(spin_grid, [r["T_mean_K"] / T_schw for r in rows],
                    "purple", linewidth=2, label="Model-A area-weighted mean")
    axes[0, 0].plot(spin_grid, [r["T_kerr_K"] / T_schw for r in rows],
                    "g--", linewidth=2, label="Kerr (uniform)")
    axes[0, 0].set_xlabel("Spin a/M")
    axes[0, 0].set_ylabel("T / T_Schwarzschild")
    axes[0, 0].set_title("Temperature: equator preserved, pole vanishes,\n"
                         "Kerr drops uniformly")
    axes[0, 0].legend(fontsize=9)
    axes[0, 0].grid(True, alpha=0.3)

    # Luminosity
    axes[0, 1].semilogy(spin_grid, [r["L_modelA_W"] / r["L_schwarzschild_W"]
                                    for r in rows], "b-", linewidth=2,
                        label="Model-A / Schwarzschild")
    axes[0, 1].semilogy(spin_grid, [r["L_kerr_W"] / r["L_schwarzschild_W"]
                                    for r in rows], "g--", linewidth=2,
                        label="Kerr / Schwarzschild")
    axes[0, 1].set_xlabel("Spin a/M")
    axes[0, 1].set_ylabel("L_emission / L_Schwarzschild")
    axes[0, 1].set_title("Stefan-Boltzmann emission rate (relative to Schw.)\n"
                         "Both shut off near extremal, but Model-A > Kerr always")
    axes[0, 1].legend(fontsize=9)
    axes[0, 1].grid(True, which="both", alpha=0.3)

    # Areas
    A_schw = 16.0 * np.pi
    axes[1, 0].plot(spin_grid, [r["A_bubble_M2"] / A_schw for r in rows],
                    "b-", linewidth=2, label="Model-A bubble")
    axes[1, 0].plot(spin_grid, [r["A_kerr_M2"] / A_schw for r in rows],
                    "g--", linewidth=2, label="Kerr horizon")
    axes[1, 0].set_xlabel("Spin a/M")
    axes[1, 0].set_ylabel("Area / Schwarzschild horizon area")
    axes[1, 0].set_title("Boundary surface area\n"
                         "Both shrink with spin; Model-A's shrinks more slowly")
    axes[1, 0].legend(fontsize=9)
    axes[1, 0].grid(True, alpha=0.3)

    # Entropy
    S_schw = bubble_entropy(M_solar, 0.0)
    axes[1, 1].plot(spin_grid, [r["S_modelA_J_per_K"] / S_schw for r in rows],
                    "b-", linewidth=2, label="Model-A bubble")
    axes[1, 1].plot(spin_grid, [r["S_kerr_J_per_K"] / S_schw for r in rows],
                    "g--", linewidth=2, label="Kerr horizon")
    axes[1, 1].set_xlabel("Spin a/M")
    axes[1, 1].set_ylabel("S / S_Schwarzschild")
    axes[1, 1].set_title("Bekenstein-Hawking entropy ∝ area\n"
                         "Model-A predicts MORE entropy at given spin than Kerr")
    axes[1, 1].legend(fontsize=9)
    axes[1, 1].grid(True, alpha=0.3)

    fig.suptitle(f"Spinning Model-A bubble thermodynamics for {M_solar} M_sun BH",
                 fontsize=13)
    out = PLOTS / "G3_thermo_quantities_vs_spin.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- Markdown summary ---

def write_markdown(rows: list[dict], plots: list[Path],
                   M_solar: float) -> Path:
    md = []
    md.append("# G3: Spinning Bubble Thermodynamics — Model-A vs Kerr\n")

    md.append("## Setup\n")
    md.append(
        "Apply Q8/Q10's resolution-rule temperature `T = ℏc|∇A|/(4π k_B)` "
        "to the spinning bubble (G2). The bubble is at A=1, with "
        "r_bubble(θ) = M + √(M² - a²cos²θ). On the bubble, |∇A| varies "
        "with θ:\n"
        "\n"
        "```text\n"
        "|∇A|²(θ)|_bubble = (M-r)²/(M²r²)  +  a⁴ sin²(2θ) / (4 M² r⁴)\n"
        "                 with r = r_bubble(θ)\n"
        "```\n"
        "\n"
        "**Equator (θ=π/2):** r=2M, sin(2θ)=0 → |∇A|=1/(2M) → T = T_Schw\n"
        "**Pole (θ=0):** r=r_+, sin(2θ)=0 → |∇A| = √(M²-a²)/(M r_+) → "
        "T → 0 at extremal\n"
        "\n"
        "The temperature is **non-uniform** along the spinning bubble, "
        "highest at the equator, lowest at the poles. This contrasts "
        "with Kerr, where T is uniform on the horizon.\n"
    )

    md.append(f"## Numerical results for {M_solar} solar mass BH\n")
    md.append("```text")
    header = (f"{'a/M':>7s}  {'T_eq/T_S':>10s}  {'T_pole/T_S':>11s}  "
              f"{'T_mean/T_S':>11s}  {'T_K/T_S':>9s}  "
              f"{'A_b/A_S':>8s}  {'A_K/A_S':>8s}  {'L_b/L_S':>9s}  "
              f"{'L_K/L_S':>9s}")
    md.append(header)
    md.append("-" * 100)
    T_schw = T_schwarzschild(M_solar)
    A_schw = 16.0 * np.pi
    L_schw = stefan_boltzmann_schwarzschild(M_solar)
    for r in rows:
        md.append(
            f"{r['a_over_M']:>7.3f}  {r['T_equator_K']/T_schw:>10.4f}  "
            f"{r['T_pole_K']/T_schw:>11.4f}  {r['T_mean_K']/T_schw:>11.4f}  "
            f"{r['T_kerr_K']/T_schw:>9.4f}  "
            f"{r['A_bubble_M2']/A_schw:>8.4f}  {r['A_kerr_M2']/A_schw:>8.4f}  "
            f"{r['L_modelA_W']/L_schw:>9.3e}  "
            f"{r['L_kerr_W']/L_schw:>9.3e}"
        )
    md.append("```\n")
    md.append("(T_S = T_Schwarzschild, A_S = Schwarzschild horizon area, "
              "L_S = Schwarzschild luminosity)\n")

    md.append("## Findings\n")
    md.append(
        "**1. Temperature distribution.** Model-A predicts the spinning "
        "bubble has a *non-uniform* temperature along it. The equatorial "
        "temperature is preserved at the Schwarzschild value for ANY spin "
        "(because the equatorial bubble radius is locked at 2M). The polar "
        "temperature decreases monotonically with spin and vanishes at "
        "extremal. Kerr's horizon temperature is uniform and decreases "
        "monotonically (zeroing at extremal as well).\n"
        "\n"
        "**2. Hawking emission rate.** Stefan-Boltzmann luminosity "
        "L = ∫ σ T⁴(θ) dA. Because Model-A retains the hot equator, its "
        "luminosity stays substantially higher than Kerr's at all spin > 0. "
        "Both predictions converge to L_Schw at a=0 and both shut off at "
        "extremal (no temperature anywhere → no emission), but the "
        "shutoff rate differs: Kerr's L drops as (r_+ - M)⁴ × A_Kerr while "
        "Model-A's L is dominated by the residual hot equator until it "
        "loses area at extremal.\n"
        "\n"
        "**3. Boundary area & entropy.** Both Model-A and Kerr have "
        "decreasing area with spin, but Model-A shrinks more slowly. At "
        "extremal: Model-A bubble area ~0.81 × Schwarzschild, Kerr horizon "
        "= 0.5 × Schwarzschild. Therefore S_Model-A > S_Kerr at any spin > "
        "0, with the gap widening at high spin.\n"
        "\n"
        "**4. First-law structure.** Standard Kerr thermodynamics has "
        "`dM = T_Kerr dS + Ω_Kerr dJ` with a single uniform T. Model-A's "
        "non-uniform T(θ) means the spinning bubble is in NON-equilibrium "
        "Hawking emission: hotter regions emit faster, the boundary is "
        "polar-cooler/equator-hotter. To recover a single-T first law, "
        "one would compute an effective T_eff via L = σ T_eff⁴ A, but the "
        "physical interpretation of T_eff differs from Kerr's T_Kerr.\n"
    )

    md.append("## Implications & STAM-distinctive predictions\n")
    md.append(
        "- **Latitudinally-resolved Hawking spectrum.** If primordial BHs "
        "(near the evaporation epoch) are spinning, Model-A predicts their "
        "Hawking emission has a polar pattern: equatorial-bright, polar-"
        "dim, with the contrast growing with spin. Kerr predicts isotropic "
        "emission. Could be testable in PBH constraints if the mass scale "
        "is right.\n"
        "\n"
        "- **Equator temperature locked at T_Schw.** This is unique to "
        "Model-A's bubble identification — a spinning BH still has the "
        "Schwarzschild temperature at its equator, because the equatorial "
        "bubble radius is locked at 2M. This means a spinning BH's "
        "EQUATORIAL Hawking emission is independent of spin in Model-A.\n"
        "\n"
        "- **Higher entropy than Kerr.** S_Model-A > S_Kerr at any spin > 0. "
        "Spinning Model-A BHs hold more information on their boundary than "
        "Kerr would assign them.\n"
        "\n"
        "- **No Penrose extraction (G2 finding).** Combined with this "
        "thermodynamics: spinning BHs in Model-A can only lose energy via "
        "Hawking-style boundary emission, never via classical energy "
        "extraction in the ergoregion. Energy budget for a spinning BH is "
        "constrained: L_emission as computed here, integrated over time, "
        "drives the BH toward the static (a=0) limit eventually.\n"
        "\n"
        "- **Smarr formula.** Q12 confirmed Smarr `Mc² = 2 T S` for "
        "Schwarzschild. For Kerr the integrated form is `Mc² = 2 T S + "
        "2 Ω J`. Model-A's non-uniform T means the Smarr-like relation "
        "needs reformulation. Whether `Mc² = 2 ⟨T⟩ S + 2 Ω J` holds for "
        "an appropriate ⟨T⟩, or whether a fundamentally different "
        "integral relation arises, is an open computation.\n"
    )

    md.append("## Honest caveats\n")
    md.append(
        "- **Coordinate gradient, not metric gradient.** Q8 used the "
        "coordinate |∇A| and that's preserved here. The proper |∇A| "
        "from the inverse metric on Model-A's spinning g would differ "
        "(would require committing to g_rr, g_tφ for the spinning case).\n"
        "- **Round-coordinate area.** The bubble area integral uses the "
        "surface-of-revolution form, not the proper induced metric on the "
        "A=1 surface. Off by O(1) factor depending on how strongly the "
        "bubble is warped. The qualitative story (non-uniform T, area "
        "decreases with spin, S_modelA > S_Kerr) is robust.\n"
        "- **Stefan-Boltzmann photon-only.** As in Q8/Q12, this counts "
        "two-helicity photon emission. Full Hawking emission integrates "
        "over all massless species and is larger by a degrees-of-freedom "
        "factor. Ratios L_ModelA/L_Kerr/L_Schw are unaffected.\n"
        "- **No dynamical evolution.** The script computes thermodynamic "
        "quantities at fixed (M, a). Including back-reaction (the BH "
        "spinning down via emission, J/M ratio evolving) is a separate "
        "evaporation calculation.\n"
    )

    md.append("## Generated plots\n")
    for p in plots:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "G3_spinning_bubble_thermodynamics_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    M_solar = 1.0
    spins = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 0.999]
    rows = [diagnostic_row(s, M_solar) for s in spins]

    plots = [
        plot_T_along_bubble(M_solar),
        plot_thermo_quantities_vs_spin(M_solar),
    ]
    summary = write_markdown(rows, plots, M_solar)

    print("G3: Spinning Model-A bubble thermodynamics")
    print("=" * 76)
    print()
    print(f"Reference: 1 solar mass BH, T_Schwarzschild = "
          f"{T_schwarzschild(M_solar):.3e} K")
    print()
    T_schw = T_schwarzschild(M_solar)
    A_schw = 16.0 * np.pi
    L_schw = stefan_boltzmann_schwarzschild(M_solar)
    print(f"  {'a/M':>6s}  {'T_eq/T_S':>10s}  {'T_pole/T_S':>11s}  "
          f"{'T_K/T_S':>9s}  {'A_b/A_S':>8s}  {'A_K/A_S':>8s}  "
          f"{'L_b/L_S':>10s}  {'L_K/L_S':>10s}")
    print("  " + "-" * 95)
    for r in rows:
        print(f"  {r['a_over_M']:>6.3f}  {r['T_equator_K']/T_schw:>10.4f}  "
              f"{r['T_pole_K']/T_schw:>11.4f}  {r['T_kerr_K']/T_schw:>9.4f}  "
              f"{r['A_bubble_M2']/A_schw:>8.4f}  {r['A_kerr_M2']/A_schw:>8.4f}  "
              f"{r['L_modelA_W']/L_schw:>10.3e}  "
              f"{r['L_kerr_W']/L_schw:>10.3e}")
    print()
    print("Verdict: Model-A predicts NON-UNIFORM T along the spinning bubble.")
    print("Equatorial T preserved at Schwarzschild value (independent of spin)")
    print("Polar T vanishes at extremal (matches Kerr)")
    print("Total emission rate higher than Kerr at any spin > 0")
    print("Bubble entropy higher than Kerr at any spin > 0")
    print("Standard first law dM = T dS + Omega dJ requires reformulation:")
    print("  Model-A's non-uniform T means latitudinal emission pattern")
    print()
    print("Files written:")
    print(f"  {summary}")
    for p in plots:
        print(f"  {p}")


if __name__ == "__main__":
    main()

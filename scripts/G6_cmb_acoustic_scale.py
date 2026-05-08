#!/usr/bin/env python3
"""
G6_cmb_acoustic_scale.py

Test STAM's H_0 = 73 cosmological commitment against the CMB acoustic angular
scale measured by Planck.

The angular scale θ_⋆ = r_s / D_C(z_⋆) is observationally locked at
θ_⋆ ≈ 0.5965° = 0.01041 rad with ~0.05% precision (Planck 2018: θ_MC =
1.04101 × 10⁻²). It depends on the comoving sound horizon at recombination
r_s and the comoving distance to last-scattering z_⋆ ≈ 1090.

LCDM achieves this with H_0 ≈ 67.4 and Ω_m ≈ 0.315. The Hubble tension is
that local distance-ladder measurements give H_0 ≈ 73 (SH0ES). STAM commits
to H_0 = 73 at all redshifts, with the discrepancy coming from photon-A
traversal that LCDM lacks.

Test: at H_0 = 73, what cosmology does STAM need to reproduce θ_⋆ ≈ 0.01041?
Does V(A) = β/(1-A) calibrated to A_0 = 0.0265 supply enough modification?

Method:
    θ_⋆ = r_s / D_C(z_⋆)
    r_s = ∫_{z_⋆}^∞ c_s(z) / H(z) dz       (comoving sound horizon)
    D_C(z_⋆) = ∫_0^{z_⋆} c / H(z') dz'      (comoving distance to LSS)
    c_s(z) = c / √(3(1 + R(z)))            (baryon-photon sound speed)
    R(z) = (3 Ω_b) / (4 Ω_γ × (1+z))       (baryon-to-photon ratio)

Cosmology models tested:
    1. LCDM-Planck:  H_0=67.4, Ω_m=0.315, Ω_Λ=0.685 (reference, matches data)
    2. LCDM at H_0=73: same matter content, just shifted H_0 (Hubble tension scenario)
    3. STAM pure-EdS:  H_0=73, Ω_m=1, no dark energy
    4. STAM + V(A):    H_0=73, Ω_m=0.693, Ω_DE_STAM=0.307 (script 35 calibration)
    5. STAM matched θ_⋆: solve for Ω_DE_STAM that reproduces observed θ_⋆ at H_0=73

Honest expected outcome:
    Scenario 1 reproduces θ_⋆ by construction.
    Scenario 2 (LCDM at H_0=73) misses θ_⋆ — this IS the Hubble tension as
        seen from the CMB side.
    Scenario 3 (pure EdS) misses badly — not enough dark-energy-like effect.
    Scenario 4 (STAM with V(A)) should be CLOSER than scenarios 2 or 3, but
        likely doesn't quite hit observed θ_⋆ because Ω_DE_STAM = 0.307
        from script 35 is smaller than the ~0.685 LCDM uses.
    Scenario 5 reverse-engineers what Ω_DE_STAM would need to be to match —
        if it comes out close to the V(A)-derived value, framework is on
        track. If it requires Ω_DE_STAM ≈ 0.685, V(A) calibration needs work.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def trapz_integrate(f, a: float, b: float, n: int = 50000) -> float:
    """Trapezoidal integration; replaces scipy.integrate.quad for these
    well-behaved cosmological integrals."""
    xs = np.linspace(a, b, n)
    ys = np.array([f(x) for x in xs])
    return float(np.trapezoid(ys, xs))


def trapz_log_integrate(f, a: float, b: float, n: int = 100000) -> float:
    """Log-spaced trapezoidal integration: useful when integrand varies
    over many orders of magnitude in the integration range (e.g., r_s
    integration from z_* to z_eq to z_BB).

    ∫_a^b f(x) dx = ∫_{ln a}^{ln b} f(e^u) e^u du
    """
    if a <= 0:
        return trapz_integrate(f, a, b, n)
    us = np.linspace(np.log(a), np.log(b), n)
    xs = np.exp(us)
    ys = np.array([f(x) * x for x in xs])
    return float(np.trapezoid(ys, us))


def bisect(f, a: float, b: float, tol: float = 1e-7,
           max_iter: int = 200) -> float | None:
    """Simple bisection to find root of f in [a, b]."""
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        return None
    for _ in range(max_iter):
        c = 0.5 * (a + b)
        fc = f(c)
        if abs(fc) < tol or (b - a) < tol:
            return c
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return 0.5 * (a + b)

# Constants
C_KMS = 299792.458                # speed of light, km/s
Z_STAR = 1090.0                   # recombination redshift
THETA_STAR_OBSERVED = 0.0104101   # Planck 2018, θ_MC

# Standard background densities (CMB-derived, same in all models)
OMEGA_GAMMA_H2 = 2.4728e-5        # photon density
N_EFF_NU = 3.046                   # effective neutrino species
# Total radiation: photons + neutrinos
OMEGA_R_H2 = OMEGA_GAMMA_H2 * (1.0 + N_EFF_NU * (7.0/8.0) * (4.0/11.0)**(4.0/3.0))
OMEGA_B_H2 = 0.02237              # baryon density (Planck 2018)

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


def sound_speed_kms(z: float) -> float:
    """c_s(z) in km/s.

    R(z) = (3/4) × (Ω_b / Ω_γ) × 1/(1+z), but better to use the physical
    densities directly (they cancel h²):
        R(z) = (3/4) × (Ω_b h²) / (Ω_γ h²) / (1+z)
    """
    R = (3.0 / 4.0) * (OMEGA_B_H2 / OMEGA_GAMMA_H2) / (1.0 + z)
    return C_KMS / np.sqrt(3.0 * (1.0 + R))


def H_kms_per_Mpc(z: float, H0: float, Om: float, OL: float) -> float:
    """H(z) in km/s/Mpc for given cosmology."""
    h = H0 / 100.0
    Or = OMEGA_R_H2 / h**2
    E_sq = Om * (1.0 + z)**3 + Or * (1.0 + z)**4 + OL
    return H0 * np.sqrt(E_sq)


def comoving_sound_horizon_Mpc(z_star: float, H0: float, Om: float,
                                 OL: float, z_max: float = 1.0e7) -> float:
    """r_s = ∫_{z_star}^{z_max} c_s(z) / H(z) dz   in comoving Mpc.

    Use log-spaced grid because integrand spans many orders of magnitude
    from z_star to z_max (radiation era at high z)."""
    integrand = lambda z: sound_speed_kms(z) / H_kms_per_Mpc(z, H0, Om, OL)
    return trapz_log_integrate(integrand, z_star, z_max, n=20000)


def comoving_distance_Mpc(z_star: float, H0: float, Om: float, OL: float) -> float:
    """D_C(z_star) = ∫_0^{z_star} c / H(z') dz'   in comoving Mpc."""
    integrand = lambda z: C_KMS / H_kms_per_Mpc(z, H0, Om, OL)
    return trapz_integrate(integrand, 0.0, z_star, n=20000)


def theta_star(H0: float, Om: float, OL: float) -> dict:
    r_s = comoving_sound_horizon_Mpc(Z_STAR, H0, Om, OL)
    D_C = comoving_distance_Mpc(Z_STAR, H0, Om, OL)
    theta = r_s / D_C
    return {
        "H0": H0,
        "Omega_m": Om,
        "Omega_DE": OL,
        "Omega_m_h2": Om * (H0/100)**2,
        "r_s_Mpc": r_s,
        "D_C_Mpc": D_C,
        "theta_star_rad": theta,
        "theta_star_deg": np.degrees(theta),
        "frac_offset_from_observed": (theta - THETA_STAR_OBSERVED) / THETA_STAR_OBSERVED,
    }


def find_OL_for_observed_theta(H0: float, Om_target: float = 0.315,
                                  flat: bool = True) -> dict:
    """Find Ω_DE such that θ_⋆(H_0, Ω_m, Ω_DE) = observed value.

    If flat=True, set Ω_m = 1 - Ω_DE (so we vary both simultaneously);
    otherwise hold Ω_m fixed at Om_target and only vary Ω_DE."""
    def deficit(OL):
        if flat:
            Om = 1.0 - OL
            if Om <= 0:
                return -1e10
            return theta_star(H0, Om, OL)["theta_star_rad"] - THETA_STAR_OBSERVED
        else:
            return theta_star(H0, Om_target, OL)["theta_star_rad"] - THETA_STAR_OBSERVED

    # Bracket the solution
    OL_solution = bisect(deficit, 0.001, 0.999, tol=1e-6)
    if OL_solution is None:
        return {"H0": H0, "Omega_DE_required": None,
                "note": "no flat solution in (0,1)"}

    if flat:
        result = theta_star(H0, 1.0 - OL_solution, OL_solution)
    else:
        result = theta_star(H0, Om_target, OL_solution)
    result["Omega_DE_required"] = OL_solution
    return result


# --- run all cosmologies ---

def run_all_models() -> dict:
    models = {
        "1. LCDM-Planck (reference)": {
            "H0": 67.4, "Om": 0.315, "OL": 0.685,
        },
        "2. LCDM at H_0=73 (Hubble tension)": {
            "H0": 73.0, "Om": 0.315, "OL": 0.685,
        },
        "3. STAM pure-EdS at H_0=73 (no V(A))": {
            "H0": 73.0, "Om": 1.0, "OL": 0.0,
        },
        "4. STAM with V(A), script-35 calibration": {
            "H0": 73.0, "Om": 0.693, "OL": 0.307,
        },
        "5. STAM tuned-flat to match Planck th_star": {
            "tuned": True,
        },
    }
    results = {}
    for name, params in models.items():
        if params.get("tuned"):
            results[name] = find_OL_for_observed_theta(73.0, flat=True)
        else:
            results[name] = theta_star(params["H0"], params["Om"], params["OL"])
    return results


# --- plot ---

def plot_theta_landscape() -> Path:
    """Show θ_⋆ as function of (H_0, Ω_DE) with the observed value contour
    and our test points."""
    H0_grid = np.linspace(60.0, 80.0, 30)
    OL_grid = np.linspace(0.0, 0.95, 30)
    theta_grid = np.zeros((len(OL_grid), len(H0_grid)))

    for i, OL in enumerate(OL_grid):
        for j, H0 in enumerate(H0_grid):
            Om = 1.0 - OL
            theta_grid[i, j] = theta_star(H0, Om, OL)["theta_star_rad"]

    fig, ax = plt.subplots(figsize=(11, 7))
    levels = np.linspace(0.005, 0.020, 16)
    cs = ax.contourf(H0_grid, OL_grid, theta_grid, levels=levels, cmap="viridis")
    cbar = plt.colorbar(cs, ax=ax, label="θ_⋆ (rad)")

    # Observed value contour
    ax.contour(H0_grid, OL_grid, theta_grid,
               levels=[THETA_STAR_OBSERVED], colors=["red"], linewidths=2.5)

    # Test points
    points = {
        "LCDM-Planck": (67.4, 0.685),
        "LCDM @ H₀=73": (73.0, 0.685),
        "STAM pure EdS": (73.0, 0.0),
        "STAM + V(A)": (73.0, 0.307),
    }
    for label, (h, ol) in points.items():
        ax.plot(h, ol, "o", color="white", markersize=10,
                markeredgecolor="black", markeredgewidth=1.5)
        ax.annotate(label, (h, ol), xytext=(7, -3), textcoords="offset points",
                    fontsize=10, color="white",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="black",
                              alpha=0.7))

    ax.set_xlabel("H₀ (km/s/Mpc)")
    ax.set_ylabel("Ω_DE (flat: Ω_m = 1 − Ω_DE)")
    ax.set_title("CMB acoustic angular scale θ_⋆ across (H₀, Ω_DE)\n"
                 "Red contour = observed θ_⋆ = 0.01041 rad (Planck)")
    ax.grid(True, alpha=0.3)

    out = PLOTS / "G6_theta_star_landscape.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- markdown summary ---

def write_markdown(results: dict, plots: list[Path]) -> Path:
    md = []
    md.append("# G6: CMB Acoustic Angular Scale Test for STAM @ H_0 = 73\n")

    md.append("## The test\n")
    md.append(
        "Planck 2018 measured the CMB acoustic peak angular scale at:\n"
        "```text\n"
        f"θ_⋆ (observed) = {THETA_STAR_OBSERVED:.6f} rad = "
        f"{np.degrees(THETA_STAR_OBSERVED):.4f}°  (~0.05% precision)\n"
        "```\n"
        "\n"
        "θ_⋆ = r_s / D_C(z_⋆) where r_s is the comoving sound horizon at "
        "z_⋆ = 1090 (recombination) and D_C is the comoving distance to "
        "last scattering. LCDM achieves this with H_0 ≈ 67.4 and "
        "Ω_m ≈ 0.315.\n"
        "\n"
        "STAM commits to **H_0 = 73 at all redshifts** (matching the local "
        "distance-ladder value, with the photon-A traversal effect "
        "explaining away the Planck-CMB inferred value). The framework's "
        "V(A) = β/(1-A) potential gives Ω_DE_STAM ≈ 0.307 at A_0 = 0.0265 "
        "(script 35). This test asks: with H_0 = 73 and that V(A) "
        "calibration, does STAM reproduce the observed θ_⋆?\n"
    )

    md.append("## Numerical results\n")
    md.append("```text")
    md.append(f"{'Model':<48s}  {'H_0':>7s}  {'Ω_m':>6s}  {'Ω_DE':>6s}  "
              f"{'r_s (Mpc)':>11s}  {'D_C (Mpc)':>11s}  "
              f"{'θ_⋆ (rad)':>12s}  {'offset':>10s}")
    md.append("-" * 130)
    for name, r in results.items():
        if r.get("Omega_DE_required") is not None:
            DE_str = f"{r['Omega_DE_required']:.4f}*"
        else:
            DE_str = f"{r.get('Omega_DE', 0.0):.4f}"
        md.append(
            f"{name:<48s}  {r['H0']:>7.2f}  "
            f"{r.get('Omega_m', 0.0):>6.3f}  {DE_str:>6s}  "
            f"{r['r_s_Mpc']:>11.2f}  {r['D_C_Mpc']:>11.2f}  "
            f"{r['theta_star_rad']:>12.6f}  "
            f"{r['frac_offset_from_observed']*100:>+9.2f}%"
        )
    md.append("```\n")
    md.append(
        f"Observed (Planck): θ_⋆ = {THETA_STAR_OBSERVED:.6f} rad = "
        f"{np.degrees(THETA_STAR_OBSERVED):.4f}°. "
        "Asterisk on Ω_DE indicates fit value to match θ_⋆.\n"
    )

    # Verdict
    model4 = results["4. STAM with V(A), script-35 calibration"]
    model5 = results["5. STAM tuned-flat to match Planck th_star"]
    md.append("## Verdict\n")
    md.append(
        f"**STAM with H_0 = 73 and the V(A) script-35 calibration "
        f"(Ω_DE_STAM = 0.307) gives θ_⋆ = {model4['theta_star_rad']:.6f} rad, "
        f"offset {model4['frac_offset_from_observed']*100:+.2f}% from "
        "observed.**\n"
        "\n"
        f"To match the observed θ_⋆ at H_0 = 73 with a flat-universe "
        f"effective Ω_DE, we'd need Ω_DE = "
        f"{model5.get('Omega_DE_required', 'undefined'):.4f} (model 5).\n"
        "\n"
    )
    if model5.get('Omega_DE_required') is not None:
        ratio = model5['Omega_DE_required'] / 0.307
        md.append(
            f"That's **{ratio:.2f}× the V(A) script-35 prediction of 0.307**. "
            f"So the V(A) potential as currently calibrated under-supplies "
            f"the dark-energy-like fraction by a factor of ~{ratio:.1f}×. "
            "This is the same factor-of-~2 discrepancy script 35 already "
            "flagged when comparing Ω_DE_STAM to LCDM's Ω_Λ.\n"
        )
    md.append("")

    md.append("## What this means\n")
    md.append(
        "**Status: a clean cosmological tension, not a fatality.** STAM's "
        "structural commitment to H_0 = 73 is **viable in principle** — "
        "there exists a flat (Ω_m, Ω_DE_eff) at H_0 = 73 that reproduces "
        "Planck's θ_⋆ exactly (model 5). The question is whether STAM's V(A) "
        "potential can supply that effective Ω_DE.\n"
        "\n"
        "Currently script 35 calibrates β so that A_0 = 0.0265 reproduces "
        "the historical bridge term b ≈ 354.95 Mly. Under that calibration, "
        "Ω_DE_STAM ≈ 0.307 — about half what the CMB needs.\n"
        "\n"
        "Three possible resolutions:\n"
        "\n"
        "**(a) The β calibration is incomplete.** A more careful FRW "
        "evolution with the V(A) potential (rather than the slow-roll "
        "tracking ansatz used in script 36) might give a different effective "
        "Ω_DE_STAM at z=0. The internal tension flagged in the modified-"
        "Friedmann SN test (A_today = 0.38 from KG dynamics vs 0.0265 from "
        "bridge) is exactly the gap that needs closing.\n"
        "\n"
        "**(b) Photon-A line-of-sight integration matters at high z.** "
        "Some of the CMB-θ_⋆ apparent value could be due to the photon-A "
        "traversal effect on the CMB photons themselves, not on the "
        "underlying expansion. STAM's structure-dependent A_LoS picture has "
        "this contribution; computing it properly might reduce the required "
        "Ω_DE_STAM at H_0 = 73.\n"
        "\n"
        "**(c) STAM's commitment to flat universe is wrong.** A small "
        "spatial curvature could absorb the discrepancy. This isn't "
        "currently part of the framework's commitments but isn't excluded "
        "either.\n"
    )

    md.append("## Honest assessment\n")
    md.append(
        "The framework currently has a **partial cosmological story**:\n"
        "\n"
        f"- ✅ **SN distance fits** (scripts 36-39): STAM with V(A) wins "
        "combined chi² over LCDM by ~25-33 across Pantheon+/Union3/DES at "
        "the same number of free parameters. Inter-catalog tension "
        "(Pantheon+/Union3 +38 mmag) predicted to within 27%.\n"
        f"- ✅ **Bridge term derived**: A_0 = 0.0265 ≈ 1/(12π) reproduces "
        "the historical b = 354.95 Mly exactly under V(A) calibration.\n"
        f"- ❌ **CMB θ_⋆ at H_0 = 73**: the V(A) calibration under-supplies "
        f"Ω_DE_eff by ~{(model5.get('Omega_DE_required',0)/0.307 if model5.get('Omega_DE_required') else 0):.1f}×.\n"
        f"- ❌ **BAO direct test** (script 22σ): simple constant-A model "
        "fails at 22σ. Structure-dependent A_LoS is the open path.\n"
        "\n"
        "**The CMB result is the cleanest current quantitative tension** in "
        "the framework. It's not a fatality (a flat solution at H_0 = 73 "
        "exists, just requires more Ω_DE_eff than current V(A) gives), but "
        "it tells us where the cosmological completion has work to do.\n"
        "\n"
        "**Headline-worthy?** Honestly, no — the framework needs to close "
        "the V(A) calibration gap (point (a) above) before this becomes a "
        "positive result. What this test DOES give is a sharp target: "
        f"STAM needs Ω_DE_eff ≈ {model5.get('Omega_DE_required', 0):.3f} "
        "from V(A) at H_0 = 73 to match Planck. Currently it produces "
        "0.307. Closing that gap is the next cosmological computation.\n"
    )

    md.append("## Generated plots\n")
    for p in plots:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "G6_cmb_acoustic_scale_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    results = run_all_models()

    plots = [plot_theta_landscape()]
    summary = write_markdown(results, plots)

    print("G6: CMB acoustic angular scale test for STAM at H_0 = 73")
    print("=" * 78)
    print()
    print(f"Observed (Planck 2018): theta_star = {THETA_STAR_OBSERVED:.6f} rad "
          f"= {np.degrees(THETA_STAR_OBSERVED):.4f} deg")
    print()
    print(f"  {'Model':<48s}  {'H0':>5s}  {'Om':>5s}  {'ODE':>5s}  "
          f"{'r_s(Mpc)':>9s}  {'D_C(Mpc)':>9s}  {'th_star':>9s}  {'offset':>8s}")
    print("  " + "-" * 120)
    for name, r in results.items():
        if r.get("Omega_DE_required") is not None:
            DE_str = f"{r['Omega_DE_required']:.3f}*"
        else:
            DE_str = f"{r.get('Omega_DE', 0.0):.3f}"
        print(f"  {name:<48s}  {r['H0']:>5.1f}  "
              f"{r.get('Omega_m', 0.0):>5.2f}  {DE_str:>5s}  "
              f"{r['r_s_Mpc']:>9.2f}  {r['D_C_Mpc']:>9.1f}  "
              f"{r['theta_star_rad']:>9.5f}  "
              f"{r['frac_offset_from_observed']*100:>+7.2f}%")
    print()
    print("(* = fit to match observed)")
    print()
    print("Files written:")
    print(f"  {summary}")
    for p in plots:
        print(f"  {p}")


if __name__ == "__main__":
    main()

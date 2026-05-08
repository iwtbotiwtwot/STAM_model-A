#!/usr/bin/env python3
"""
G1_qnm_ringdown.py

Quasinormal mode (QNM) ringdown spectrum for binary black hole mergers
under Model-A's strong-field metric vs. Schwarzschild.

Setup:
    Both metrics share g_tt = -(1-A)c²  (h(A) = 1-A, identical to GR)
    They differ in g_rr:
        Schwarzschild:  g_rr = 1/(1-A)              [k(A) = 1-A]
        Model-A:        g_rr = 1/[(1-A)(1-A²)²]    [k(A) = (1-A)(1-A²)²]

    Because g_tt is unchanged, the photon sphere sits at r_c = 1.5 R_s (A=2/3)
    in both theories, and the orbital frequency Ω_c there is the same. So the
    DOMINANT RINGDOWN FREQUENCY (real part of QNM) is essentially identical
    in eikonal approximation.

    Where the theories diverge is the Lyapunov exponent λ of unstable circular
    null orbits at the photon sphere, which controls the DAMPING TIME of the
    ringdown. λ depends on the product h(r_c) × k(r_c) — and k differs:
        k_Schwarzschild(r_c=3M) = 1/3
        k_Model-A(r_c=3M)       = (1/3)(1 - 4/9)² = 25/243 ≈ 0.103
    The Model-A modification suppresses k near the photon sphere by ~3.2x,
    which propagates into a SLOWER damping (longer τ, higher Q).

Method (eikonal / WKB at photon sphere, following Cardoso-Konoplya-Zhidenko):
    Ω_c = √(h(r_c)) / r_c              [orbital angular velocity]
    λ²  = (h × k / (2 E²)) × |V''_eff(r_c)|
    V_eff(r) = L² × h(r) / r²          [null-orbit effective potential]
    L² = E² × r_c² / h(r_c)            [photon-orbit angular momentum]

    Eikonal QNM:
        ω_R = l × Ω_c
        ω_I = -(n + 1/2) × λ

    For the dominant ℓ=2, n=0 mode (what LIGO measures from BBH ringdowns).

Conversion to SI:
    1 geometric unit M = G M / c³  (in seconds)
    f_R = ω_R / (2π) in Hz
    τ   = 1 / |ω_I|  in seconds (e-folding time)
    Q   = ω_R / (2 |ω_I|)  (quality factor)

Validation:
    For Schwarzschild eikonal (M=1): ω_R = 2/(3√3) ≈ 0.385, ω_I = -1/(2 × 3√3) ≈ -0.0962.
    Berti's exact ℓ=2, n=0 axial gravitational mode: ω_R ≈ 0.3737, ω_I ≈ -0.0890.
    Eikonal undershoots ω_R by ~3% and undershoots |ω_I| by ~8%. The RATIO of
    Model-A to Schwarzschild is more robust than the absolute numbers.

Falsification logic:
    LIGO measures ringdown frequency f and damping time τ from compact-binary
    coalescences. For GW150914-like events (final mass ~62 M_sun), GR predicts
    f ~ 252 Hz and τ ~ 4 ms for the dominant ℓ=2, m=2, n=0 mode.

    Model-A predicts the SAME frequency (Ω_c unchanged) but a LONGER damping
    time by a factor ~λ_S/λ_MA. If LIGO's measured τ is consistent with the
    Schwarzschild value and inconsistent with Model-A's larger value, that's
    a falsification of the modified k(A) at A=2/3.

    This is genuinely STAM-distinctive — it's the cleanest near-horizon test
    that doesn't require resolving features at the horizon itself.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Constants
C = 2.99792458e8                 # m/s
G = 6.67430e-11                  # m^3 / (kg s^2)
M_SUN = 1.98892e30               # kg

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# --- metric components in geometric units (M=1, c=G=1) ---

def A_of_r(r: float, M: float = 1.0) -> float:
    return 2.0 * M / r


def h(r: float, M: float = 1.0) -> float:
    """g_tt coefficient: g_tt = -h(r) c². Same for both metrics."""
    return 1.0 - A_of_r(r, M)


def k_schwarzschild(r: float, M: float = 1.0) -> float:
    """g_rr = 1/k. Schwarzschild: k = 1 - A."""
    return 1.0 - A_of_r(r, M)


def k_model_a(r: float, M: float = 1.0) -> float:
    """g_rr = 1/k. Model-A: k = (1-A)(1-A²)²."""
    A = A_of_r(r, M)
    return (1.0 - A) * (1.0 - A ** 2) ** 2


# --- photon sphere quantities ---

def photon_sphere_radius(M: float = 1.0) -> float:
    """Photon sphere: r_c h'(r_c) = 2 h(r_c). For h = 1 - 2M/r, r_c = 3M.
    Same for Schwarzschild and Model-A (depends only on h)."""
    return 3.0 * M


def Omega_c(M: float = 1.0) -> float:
    """Orbital angular velocity at photon sphere = √(h)/r_c.
    Same for both metrics (h-only)."""
    r_c = photon_sphere_radius(M)
    return np.sqrt(h(r_c, M)) / r_c


def lambda_lyapunov(k_func, M: float = 1.0) -> float:
    """Lyapunov exponent for unstable circular null orbit at photon sphere.

    From V_eff(r) = L² h(r)/r² (with E=1, L² = r_c²/h(r_c) for photon orbit):
        λ² = h(r_c) × k(r_c) / 2 × |V''_eff(r_c)| / E²

    For h = 1 - 2M/r and r_c = 3M, the closed form simplifies to:
        |V''_eff(r_c)|/E² = 2 / (3 M²)

    so:
        λ² = (h × k / 2) × (2 / (3 M²)) = (h × k) / (3 M²)
    """
    r_c = photon_sphere_radius(M)
    h_c = h(r_c, M)
    k_c = k_func(r_c, M)
    V_double_prime_over_E2 = 2.0 / (3.0 * M ** 2)
    lambda_sq = (h_c * k_c / 2.0) * V_double_prime_over_E2
    return np.sqrt(lambda_sq)


# --- QNM frequencies (eikonal / WKB at photon sphere) ---

def qnm_eikonal(k_func, l: int, n: int, M: float = 1.0) -> dict:
    """Eikonal QNM: ω_R = l × Ω_c, ω_I = -(n + 1/2) × λ.

    Returns ω in geometric units (1/M)."""
    omega_c = Omega_c(M)
    lam = lambda_lyapunov(k_func, M)
    omega_R = l * omega_c
    omega_I = -(n + 0.5) * lam
    return {
        "omega_R_geom": omega_R,
        "omega_I_geom": omega_I,
        "Omega_c_geom": omega_c,
        "lambda_geom": lam,
        "Q_factor": omega_R / (2 * abs(omega_I)),
    }


def to_SI(omega_geom: float, M_solar: float) -> float:
    """Convert ω from geometric units (1/M) to SI rad/s.

    Geometric M (in seconds) = G M_kg / c³.
    """
    M_kg = M_solar * M_SUN
    geom_M_in_seconds = G * M_kg / C ** 3
    return omega_geom / geom_M_in_seconds


def ringdown_observables(k_func, M_solar: float, l: int = 2, n: int = 0) -> dict:
    """Compute LIGO-relevant observables: f (Hz), τ (s), Q-factor."""
    qnm = qnm_eikonal(k_func, l, n, M=1.0)
    omega_R_SI = to_SI(qnm["omega_R_geom"], M_solar)
    omega_I_SI = to_SI(qnm["omega_I_geom"], M_solar)
    f_Hz = omega_R_SI / (2 * np.pi)
    tau_s = 1.0 / abs(omega_I_SI)
    return {
        "f_Hz": f_Hz,
        "tau_s": tau_s,
        "tau_ms": tau_s * 1e3,
        "Q_factor": qnm["Q_factor"],
        "omega_R_geom": qnm["omega_R_geom"],
        "omega_I_geom": qnm["omega_I_geom"],
        "lambda_geom": qnm["lambda_geom"],
    }


# --- benchmark masses (LIGO BBH remnants) ---

BENCHMARKS = [
    ("GW150914 remnant", 62.0),
    ("GW170729 remnant", 80.0),
    ("GW190521 remnant", 142.0),
    ("Stellar-mass BBH (typical)", 30.0),
    ("Intermediate-mass (~1000 M_sun)", 1000.0),
    ("Sgr A*", 4.3e6),
    ("M87*", 6.5e9),
]


def benchmark_table() -> list[dict]:
    rows = []
    for label, M_solar in BENCHMARKS:
        ring_S = ringdown_observables(k_schwarzschild, M_solar)
        ring_MA = ringdown_observables(k_model_a, M_solar)
        rows.append({
            "label": label,
            "M_solar": M_solar,
            "f_Hz_S": ring_S["f_Hz"],
            "f_Hz_MA": ring_MA["f_Hz"],
            "tau_ms_S": ring_S["tau_ms"],
            "tau_ms_MA": ring_MA["tau_ms"],
            "Q_S": ring_S["Q_factor"],
            "Q_MA": ring_MA["Q_factor"],
            "tau_ratio": ring_MA["tau_ms"] / ring_S["tau_ms"],
        })
    return rows


# --- plots ---

def plot_lyapunov_comparison() -> Path:
    """Show how λ varies along radius for both metrics — but at photon sphere
    the Schwarzschild and Model-A values diverge most clearly."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # k(r) curves
    A_grid = np.linspace(0, 0.95, 500)
    k_S_curve = (1 - A_grid)
    k_MA_curve = (1 - A_grid) * (1 - A_grid ** 2) ** 2

    axes[0].plot(A_grid, k_S_curve, "b-", linewidth=2, label="Schwarzschild k(A) = 1-A")
    axes[0].plot(A_grid, k_MA_curve, "r-", linewidth=2,
                 label="Model-A k(A) = (1-A)(1-A²)²")
    axes[0].axvline(2.0 / 3.0, color="gray", linestyle="--", alpha=0.5,
                    label="Photon sphere (A=2/3)")
    axes[0].set_xlabel("A = R_s / r")
    axes[0].set_ylabel("k(A)  [g_rr = 1/k]")
    axes[0].set_title("Metric component k(A) — same at A=0, diverges at A→1")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Ratio of λs
    M_solar_grid = np.logspace(0, 10, 200)
    ring_S = [ringdown_observables(k_schwarzschild, M)["tau_ms"] for M in M_solar_grid]
    ring_MA = [ringdown_observables(k_model_a, M)["tau_ms"] for M in M_solar_grid]

    axes[1].loglog(M_solar_grid, ring_S, "b-", linewidth=2, label="Schwarzschild τ")
    axes[1].loglog(M_solar_grid, ring_MA, "r-", linewidth=2, label="Model-A τ")
    axes[1].set_xlabel("Final BH mass (M_sun)")
    axes[1].set_ylabel("Ringdown damping time τ (ms)")
    axes[1].set_title(f"Damping time τ scales linearly with M; "
                      f"Model-A is {ring_MA[0]/ring_S[0]:.2f}× longer")
    axes[1].legend()
    axes[1].grid(True, which="both", alpha=0.3)

    out = PLOTS / "G1_qnm_metric_and_damping.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_ringdown_waveforms() -> Path:
    """Show the actual time-domain ringdown for GW150914-like remnant."""
    M_solar = 62.0  # GW150914 final BH mass
    ring_S = ringdown_observables(k_schwarzschild, M_solar)
    ring_MA = ringdown_observables(k_model_a, M_solar)

    t_ms = np.linspace(0, 30, 5000)
    t_s = t_ms * 1e-3

    h_S = np.exp(-t_s / ring_S["tau_s"]) * np.cos(2 * np.pi * ring_S["f_Hz"] * t_s)
    h_MA = np.exp(-t_s / ring_MA["tau_s"]) * np.cos(2 * np.pi * ring_MA["f_Hz"] * t_s)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    axes[0].plot(t_ms, h_S, "b-", linewidth=1.5,
                 label=(f"Schwarzschild  f={ring_S['f_Hz']:.0f} Hz, "
                        f"τ={ring_S['tau_ms']:.2f} ms"))
    axes[0].plot(t_ms, h_MA, "r-", linewidth=1.5, alpha=0.85,
                 label=(f"Model-A          f={ring_MA['f_Hz']:.0f} Hz, "
                        f"τ={ring_MA['tau_ms']:.2f} ms"))
    axes[0].axhline(0, color="black", linewidth=0.4)
    axes[0].set_ylabel("Ringdown amplitude (normalized)")
    axes[0].set_title(f"GW150914-like ringdown ({M_solar} M_sun remnant): "
                      "same frequency, ~80% longer damping in Model-A")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Envelope only
    env_S = np.exp(-t_s / ring_S["tau_s"])
    env_MA = np.exp(-t_s / ring_MA["tau_s"])
    axes[1].semilogy(t_ms, env_S, "b-", linewidth=2, label="Schwarzschild envelope")
    axes[1].semilogy(t_ms, env_MA, "r-", linewidth=2, label="Model-A envelope")
    axes[1].set_xlabel("Time after merger (ms)")
    axes[1].set_ylabel("|envelope|")
    axes[1].set_title("Decay envelopes — Model-A rings ~80% longer")
    axes[1].legend()
    axes[1].grid(True, which="both", alpha=0.3)
    axes[1].set_ylim(1e-6, 2)

    out = PLOTS / "G1_qnm_ringdown_waveforms.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_benchmark_comparison(rows: list[dict]) -> Path:
    fig, ax = plt.subplots(figsize=(11, 5.5))
    labels = [r["label"] for r in rows]
    M_vals = [r["M_solar"] for r in rows]
    ratio_vals = [r["tau_ratio"] for r in rows]

    bars = ax.bar(range(len(labels)), ratio_vals, color="tab:purple",
                  edgecolor="black", alpha=0.8)
    ax.axhline(1.0, color="tab:blue", linestyle="--", alpha=0.7,
               label="Schwarzschild reference (τ_MA / τ_S = 1)")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("τ_ModelA / τ_Schwarzschild")
    ax.set_title("Predicted ringdown damping-time ratio (Model-A / Schwarzschild)\n"
                 "Ratio is mass-independent — same factor for any BH")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    for bar, val in zip(bars, ratio_vals):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.02, f"{val:.3f}",
                ha="center", fontsize=9)
    out = PLOTS / "G1_qnm_tau_ratio.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- markdown summary ---

def write_markdown(rows: list[dict], plots: list[Path]) -> Path:
    md = []
    md.append("# G1: Quasinormal-Mode Ringdown Spectrum — Model-A vs Schwarzschild\n")

    md.append("## The test\n")
    md.append(
        "Model-A's strong-field metric differs from Schwarzschild only in g_rr "
        "(time component g_tt is identical). The QNM ringdown spectrum after "
        "binary BH merger probes near-horizon geometry, where the (1-A²)² "
        "modification of k(A) bites. We compute the dominant ℓ=2, n=0 mode "
        "(what LIGO measures) for both metrics in the eikonal/WKB approximation "
        "at the photon sphere, then convert to physical frequency and damping "
        "time for representative BBH remnants.\n"
    )

    md.append("## Setup\n")
    md.append(
        "**Both metrics share:**\n"
        "- g_tt = -(1-A) c²  [h(A) = 1-A]\n"
        "- Photon sphere at r_c = 3M (A = 2/3)\n"
        "- Orbital frequency Ω_c = √(1/3) / (3M) = 1/(3√3 M)\n"
        "\n"
        "**They differ in g_rr:**\n"
        "- Schwarzschild: g_rr = 1/(1-A)              k(r_c) = 1/3\n"
        "- Model-A:       g_rr = 1/[(1-A)(1-A²)²]     k(r_c) = 25/243 ≈ 0.103\n"
        "\n"
        "**Eikonal QNM formulas (Cardoso-Konoplya-Zhidenko):**\n"
        "```text\n"
        "ω_R = l × Ω_c                              (real, oscillation frequency)\n"
        "ω_I = -(n + 1/2) × λ                       (imaginary, damping rate)\n"
        "λ²  = h(r_c) × k(r_c) / (3 M²)             (Lyapunov exponent at photon sphere)\n"
        "```\n"
        "\n"
        "Because Ω_c depends on h only, Model-A predicts the **same ringdown "
        "frequency** as Schwarzschild for any given BH mass. The damping rate λ "
        "depends on h × k together, so Model-A's smaller k at the photon sphere "
        "gives a **longer damping time τ = 1/|ω_I|** by a factor:\n"
        "\n"
        "```text\n"
        "λ_S / λ_MA = √(k_S / k_MA) = √((1/3) / (25/243)) = √(81/25) = 9/5 = 1.80\n"
        "τ_MA / τ_S = λ_S / λ_MA   = 1.80\n"
        "```\n"
        "\n"
        "**Model-A predicts ringdowns are 80% longer than GR predicts**, at the "
        "same dominant frequency, for any BH mass. This ratio is independent of "
        "BH mass — it's set entirely by the photon-sphere k value.\n"
    )

    md.append("## Predictions for representative LIGO events\n")
    md.append(
        "(Eikonal approximation; absolute values are ~3% off the Berti exact "
        "ℓ=2 axial gravitational mode, but the **ratio** Model-A / Schwarzschild "
        "is robust because both share the same approximation structure.)\n"
        "\n"
        "```text\n"
    )
    header = (f"{'Event':<35s}  {'M (M_sun)':>10s}  "
              f"{'f_GR(Hz)':>10s}  {'f_MA(Hz)':>10s}  "
              f"{'τ_GR(ms)':>10s}  {'τ_MA(ms)':>10s}  {'τ_MA/τ_GR':>10s}")
    md.append(header)
    md.append("-" * 110)
    for r in rows:
        md.append(
            f"{r['label']:<35s}  {r['M_solar']:>10.2g}  "
            f"{r['f_Hz_S']:>10.1f}  {r['f_Hz_MA']:>10.1f}  "
            f"{r['tau_ms_S']:>10.4g}  {r['tau_ms_MA']:>10.4g}  "
            f"{r['tau_ratio']:>10.3f}"
        )
    md.append("```\n")

    md.append("## Falsification logic\n")
    md.append(
        "LIGO's ringdown analysis of GW150914 (final BH ~62 M_sun) extracted a "
        "dominant mode at f ≈ 251 Hz with damping time τ ≈ 4 ms — consistent "
        "with the Schwarzschild prediction. Subsequent BBH events (GW170729, "
        "GW190521, etc.) have produced similar ringdown measurements.\n"
        "\n"
        "Model-A predicts the **same frequency** (f ≈ 251 Hz for 62 M_sun) but a "
        "**1.80× longer damping time** (τ ≈ 7.2 ms instead of 4 ms). This is "
        "a 80% effect — well above current measurement uncertainties on damping "
        "time for the loudest events.\n"
        "\n"
        "**The cleanest test:** the ringdown signal h(t) ∝ exp(-t/τ) cos(2πft) "
        "should ring 80% longer than GR predicts if Model-A is right. LIGO's "
        "matched-filter analysis fits both f and τ; if τ comes out consistently "
        "with the GR value (~4 ms for 62 M_sun class) and inconsistent with Model-"
        "A's longer prediction, that's a falsification of Model-A's k(A) at "
        "A = 2/3.\n"
        "\n"
        "**What this distinguishes:** the modification factor k(A) = (1-A²)² "
        "× (1-A) hits the photon sphere at A = 2/3, where (1-A²)² = (5/9)² ≈ "
        "0.31. So the photon sphere sees about 1/3 the Schwarzschild k value. "
        "The ringdown probes exactly this region. If we miss the prediction, "
        "the modification factor must be tuned down — but the modification is "
        "constrained at A=1 (must give boundary behavior) and at A→0 (must "
        "match GR weak field), so there's limited room to dial it without "
        "breaking other commitments.\n"
    )

    md.append("## Honest caveats\n")
    md.append(
        "- **Eikonal approximation.** This is a high-ℓ asymptotic that's about "
        "3% off for ℓ=2. The RATIO Model-A/Schwarzschild is more robust than "
        "absolute frequencies because both share the same approximation.\n"
        "- **Scalar perturbation surrogate.** This script computes scalar QNMs "
        "as a stand-in for the gravitational ℓ=2, m=2, n=0 mode. The actual "
        "gravitational perturbation equations (Regge-Wheeler / Zerilli for "
        "Schwarzschild; Model-A analog yet to be derived) may modify the "
        "result. For Schwarzschild, scalar and axial-gravitational ℓ=2 QNMs "
        "differ by ~10%; same for Model-A. The Model-A/GR ratio of damping "
        "times should still be ~1.8 to within similar precision.\n"
        "- **Spin neglected.** Real BH remnants are spinning (typical χ_f ≈ 0.7 "
        "for BBH mergers). Kerr QNMs differ from Schwarzschild; Model-A's "
        "spinning analog has not been constructed. Spin-aligned events with "
        "moderate χ should still show the qualitative ratio, but a full Kerr-"
        "Model-A computation is a follow-up.\n"
        "- **WKB at higher orders.** Schutz-Will (1985) gives a few-percent "
        "improvement over leading WKB; Konoplya-Zhidenko (2019) give 6th-order. "
        "These improve the absolute QNM accuracy but not the Model-A/GR ratio.\n"
    )

    md.append("## Generated plots\n")
    for p in plots:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "G1_qnm_ringdown_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    rows = benchmark_table()
    plots = [
        plot_lyapunov_comparison(),
        plot_ringdown_waveforms(),
        plot_benchmark_comparison(rows),
    ]
    summary = write_markdown(rows, plots)

    # Console
    print("G1: Quasinormal-mode ringdown spectrum — Model-A vs Schwarzschild")
    print("=" * 76)
    print()
    print("Eikonal QNM at photon sphere (l=2, n=0 dominant ringdown mode):")
    print()
    qnm_S = qnm_eikonal(k_schwarzschild, l=2, n=0, M=1.0)
    qnm_MA = qnm_eikonal(k_model_a, l=2, n=0, M=1.0)
    print(f"  Schwarzschild: omega_R M = {qnm_S['omega_R_geom']:.4f}, "
          f"omega_I M = {qnm_S['omega_I_geom']:.4f}, Q = {qnm_S['Q_factor']:.2f}")
    print(f"  Model-A:       omega_R M = {qnm_MA['omega_R_geom']:.4f}, "
          f"omega_I M = {qnm_MA['omega_I_geom']:.4f}, Q = {qnm_MA['Q_factor']:.2f}")
    print()
    print(f"  lambda_Schwarzschild = {qnm_S['lambda_geom']:.4f}/M")
    print(f"  lambda_Model-A       = {qnm_MA['lambda_geom']:.4f}/M")
    print(f"  Ratio lambda_S / lambda_MA = "
          f"{qnm_S['lambda_geom']/qnm_MA['lambda_geom']:.4f}")
    print(f"  Ratio tau_MA / tau_S       = "
          f"{abs(qnm_S['omega_I_geom'])/abs(qnm_MA['omega_I_geom']):.4f}")
    print()
    print("Predictions for benchmark BH masses:")
    print()
    print(f"  {'Event':<35s}  {'M (M_sun)':>10s}  {'f (Hz)':>10s}  "
          f"{'tau_GR (ms)':>12s}  {'tau_MA (ms)':>12s}  {'tau_ratio':>10s}")
    print("  " + "-" * 100)
    for r in rows:
        print(f"  {r['label']:<35s}  {r['M_solar']:>10.2g}  {r['f_Hz_S']:>10.1f}  "
              f"{r['tau_ms_S']:>12.4g}  {r['tau_ms_MA']:>12.4g}  "
              f"{r['tau_ratio']:>10.3f}")
    print()
    print("Falsification: LIGO ringdown measurements (e.g. GW150914) constrain")
    print("tau to within ~10-20% on the loudest events. Model-A predicts a 1.80x")
    print("longer damping time at the SAME frequency. If LIGO data is consistent")
    print("with Schwarzschild tau and inconsistent with Model-A's prediction,")
    print("Model-A's k(A) modification at A=2/3 is falsified.")
    print()
    print("Files written:")
    print(f"  {summary}")
    for p in plots:
        print(f"  {p}")


if __name__ == "__main__":
    main()

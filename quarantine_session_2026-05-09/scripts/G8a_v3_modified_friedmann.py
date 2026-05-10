"""
G8a - Model-A V_3 modified Friedmann cosmology

V_3 selected by G8 under the symmetric-boundary commitment:
    V_3(A) = alpha/A + beta/(1-A)

with
    alpha/beta = [A_0 / (1-A_0)]^2     (structural minimum at A_0)
    A_0 = 1/(12 pi)                    (G7 commitment)

This is the V_3 retrofit of script 36 (which used V_1 = beta/(1-A)).
V_1 is rejected by G8 because it does not diverge at A=0 and therefore
does not satisfy the symmetric-boundary commitment.

Calibration of beta_tilde:
    Under V_1, beta was set by V'(A_0) = kappa rho_m_today, giving
    beta_tilde = Omega_m * (1-A_0)^2. Under V_3, V'(A_0) = 0 by construction
    (A_0 is the minimum). So we calibrate beta_tilde by demanding that
    Model-A's modified Friedmann closes with h^2(today) = 1, given
    Omega_m_total = 0.315 (Model-A internal closure, not anchored to any
    other framework).

Tracking ansatz:
    Slow-roll equilibrium: Y'(A_eq) = Omega_m * a^-3
    Tracking velocity:     omega = -3 h Y'(A_eq) / Y''(A_eq)

The script:
    1. Solves for beta_tilde such that h^2(today) = 1 with V_3 + tracking.
    2. Computes the equilibrium A_eq(a) numerically across cosmic history.
    3. Tabulates Omega_m, Omega_DE_pot, Omega_DE_kin, Omega_DE_total today.
    4. Computes effective dark-energy w_eff(z).
    5. Verifies tracking by integrating the full Klein-Gordon + Friedmann
       system with RK4.
    6. Plots H(z), density components, w_eff, and tracking-vs-numerical.

Author: STAM Model-A research line, 2026-05-09. Following the two-tables
framing: this script reports Model-A's predictions on its own table.
"""

from __future__ import annotations
import math
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import brentq

# ===================================================================
# Model-A V_3 cosmology inputs (no tuning)
# ===================================================================
A0 = 1.0 / (12.0 * math.pi)        # G7 commitment
OMEGA_M = 0.315                     # total matter content (Model-A internal)
OMEGA_R = 9.2e-5                    # radiation today
H0_KMSMPC = 73.04                   # Model-A's H_0 (plot label only)

# Structural ratio fixed by V_3 minimum at A_0
ALPHA_OVER_BETA = (A0 / (1.0 - A0))**2

# ===================================================================
# V_3 potential and derivatives
# ===================================================================
def Y_pot(A, beta_tilde):
    """V_3(A) / (kappa rho_crit) = alpha_tilde/A + beta_tilde/(1-A)."""
    alpha_tilde = beta_tilde * ALPHA_OVER_BETA
    return alpha_tilde / A + beta_tilde / (1.0 - A)


def Yp(A, beta_tilde):
    """Derivative dY/dA."""
    alpha_tilde = beta_tilde * ALPHA_OVER_BETA
    return -alpha_tilde / A**2 + beta_tilde / (1.0 - A)**2


def Ypp(A, beta_tilde):
    """Second derivative d^2Y/dA^2."""
    alpha_tilde = beta_tilde * ALPHA_OVER_BETA
    return 2.0 * alpha_tilde / A**3 + 2.0 * beta_tilde / (1.0 - A)**3


# ===================================================================
# Tracking equilibrium A_eq(a)
# ===================================================================
def A_eq_at(a, beta_tilde):
    """Solve Y'(A_eq) = Omega_m * a^-3 for A_eq.

    Y'(A_0) = 0, Y'(A=1) -> +inf. Root in (A_0, 1) for any rho_m > 0.
    """
    rho_m = OMEGA_M * a**-3
    f = lambda A: Yp(A, beta_tilde) - rho_m
    return brentq(f, A0 + 1e-12, 1.0 - 1e-9)


def h2_tracking(a, beta_tilde):
    """h^2(a) under V_3 modified Friedmann with tracking ansatz.

    Tracking: omega = -3 h Y'(A_eq) / Y''(A_eq) = -3 h rho_m / Y''(A_eq)
    Friedmann: h^2 = rho_m + rho_r + omega^2/6 + Y(A_eq)
    Solve algebraically for h^2:
        h^2 = (rho_m + rho_r + Y(A_eq)) / (1 - K(a))
        K(a) = (3/2) rho_m^2 / Y''(A_eq)^2
    """
    A_eq = A_eq_at(a, beta_tilde)
    rho_m = OMEGA_M * a**-3
    rho_r = OMEGA_R * a**-4
    Y_val = Y_pot(A_eq, beta_tilde)
    Ypp_val = Ypp(A_eq, beta_tilde)
    K = (3.0 / 2.0) * rho_m**2 / Ypp_val**2
    return (rho_m + rho_r + Y_val) / (1.0 - K)


def h2_components(a, beta_tilde):
    """Decompose h^2(a) into matter, radiation, kinetic, potential."""
    h2 = h2_tracking(a, beta_tilde)
    A_eq = A_eq_at(a, beta_tilde)
    rho_m = OMEGA_M * a**-3
    rho_r = OMEGA_R * a**-4
    Y_val = Y_pot(A_eq, beta_tilde)
    Ypp_val = Ypp(A_eq, beta_tilde)
    omega = -3.0 * np.sqrt(h2) * rho_m / Ypp_val
    Om_kin = omega**2 / 6.0
    return {
        'h2': h2,
        'A_eq': A_eq,
        'omega': omega,
        'Omega_m': rho_m,
        'Omega_r': rho_r,
        'Omega_DE_pot': Y_val,
        'Omega_DE_kin': Om_kin,
        'Omega_DE_total': Y_val + Om_kin,
    }


# ===================================================================
# Calibrate beta_tilde by h^2(today) = 1
# ===================================================================
def closure_residual(beta_tilde):
    return h2_tracking(1.0, beta_tilde) - 1.0


print("=" * 72)
print("Model-A V_3 modified Friedmann (G8a)")
print("=" * 72)
print(f"A_0 (G7 structural)              = 1/(12 pi) = {A0:.8f}")
print(f"alpha/beta (V_3 minimum at A_0)  = {ALPHA_OVER_BETA:.8f}")
print(f"Omega_m (today, Model-A total)   = {OMEGA_M}")
print(f"Omega_r (today)                  = {OMEGA_R:.2e}")

# Find beta_tilde via h(today)=1 closure
# Bracket: at small beta the field has little potential energy, h^2 < 1.
# At large beta, h^2 > 1. Look for sign change.
betas_probe = np.logspace(-3, 1.5, 50)
residuals = []
for b in betas_probe:
    try:
        residuals.append(closure_residual(b))
    except Exception:
        residuals.append(np.nan)
residuals = np.array(residuals)

# Find adjacent indices that bracket zero
sign_changes = np.where(np.sign(residuals[:-1]) != np.sign(residuals[1:]))[0]
if len(sign_changes) == 0:
    raise RuntimeError("h^2(today) = 1 closure bracket not found in beta scan.")
i = sign_changes[0]
beta_low, beta_high = betas_probe[i], betas_probe[i + 1]
BETA_TILDE = brentq(closure_residual, beta_low, beta_high, xtol=1e-10)
ALPHA_TILDE = BETA_TILDE * ALPHA_OVER_BETA

print(f"\nbeta_tilde calibrated by h^2(today) = 1 closure:")
print(f"  beta_tilde  = beta /(kappa rho_crit)  = {BETA_TILDE:.6f}")
print(f"  alpha_tilde = alpha/(kappa rho_crit)  = {ALPHA_TILDE:.6e}")
print(f"  closure residual h^2(1) - 1           = {closure_residual(BETA_TILDE):+.2e}")

# ===================================================================
# Today's Omega budget
# ===================================================================
print("\n--- Closure check at z = 0 (a = 1) ---")
comp_today = h2_components(1.0, BETA_TILDE)
print(f"  Omega_matter             = {comp_today['Omega_m']:.6f}")
print(f"  Omega_radiation          = {comp_today['Omega_r']:.6e}")
print(f"  A_eq (today)             = {comp_today['A_eq']:.8f}")
print(f"  A_eq - A_0               = {comp_today['A_eq'] - A0:+.6e}")
print(f"  omega (today)            = {comp_today['omega']:+.6f}")
print(f"  Omega_DE_potential       = {comp_today['Omega_DE_pot']:.6f}")
print(f"  Omega_DE_kinetic         = {comp_today['Omega_DE_kin']:.6f}")
print(f"  Omega_DE_total           = {comp_today['Omega_DE_total']:.6f}")
print(f"  Omega_total (h^2 today)  = {comp_today['h2']:.6f}")

# ===================================================================
# H(z) tabulation
# ===================================================================
print("\n--- H(z) tabulation ---")
z_grid = np.linspace(0, 5, 600)
a_grid = 1.0 / (1.0 + z_grid)
h_grid = np.array([np.sqrt(h2_tracking(a, BETA_TILDE)) for a in a_grid])

z_check = [0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0]
print(f"  {'z':>6} {'h(z)':>10} {'H(z) [km/s/Mpc]':>20}")
for z in z_check:
    a = 1.0 / (1.0 + z)
    h = np.sqrt(h2_tracking(a, BETA_TILDE))
    print(f"  {z:>6.2f} {h:>10.4f} {h * H0_KMSMPC:>20.2f}")

# ===================================================================
# Effective dark-energy w_eff(z)
# rho_DE = rho_DE_pot + rho_DE_kin
# w_eff = -1 - (1/3) d ln rho_DE / d ln a
# ===================================================================
a_we = np.logspace(-2.5, -0.001, 400)
z_we = 1.0 / a_we - 1.0
rho_DE = np.array([h2_components(a, BETA_TILDE)['Omega_DE_total'] for a in a_we])
ln_a = np.log(a_we)
ln_rho = np.log(rho_DE)
w_eff = -1.0 - (1.0 / 3.0) * np.gradient(ln_rho, ln_a)
w_today = float(np.interp(0.0, z_we[::-1], w_eff[::-1]))
w_z1 = float(np.interp(1.0, z_we[::-1], w_eff[::-1]))
w_z3 = float(np.interp(3.0, z_we[::-1], w_eff[::-1]))
print("\n--- Effective dark-energy w_eff(z) ---")
print(f"  w_eff(z=0)   = {w_today:+.4f}")
print(f"  w_eff(z=1)   = {w_z1:+.4f}")
print(f"  w_eff(z=3)   = {w_z3:+.4f}")

# ===================================================================
# Numerical verification: full coupled KG + Friedmann
# ===================================================================
print("\n--- Numerical KG + Friedmann verification ---")

def rhs_v3(y):
    A, omega, lna = y
    a_ = math.exp(lna)
    Y_v = Y_pot(A, BETA_TILDE)
    Yp_v = Yp(A, BETA_TILDE)
    rho_m = OMEGA_M * a_**-3
    rho_r = OMEGA_R * a_**-4
    h2_ = rho_m + rho_r + omega * omega / 6.0 + Y_v
    h_ = math.sqrt(max(h2_, 1e-30))
    dA = omega
    domega = -3.0 * h_ * omega - 3.0 * (Yp_v - rho_m)
    dlna = h_
    return np.array([dA, domega, dlna])


# Initial condition: tracking from a = 1e-3
lna_init = math.log(1e-3)
a_init = math.exp(lna_init)
A_init = A_eq_at(a_init, BETA_TILDE)
rho_m_init = OMEGA_M * a_init**-3
h_init = math.sqrt(rho_m_init + OMEGA_R * a_init**-4)
omega_init = -3.0 * h_init * Yp(A_init, BETA_TILDE) / Ypp(A_init, BETA_TILDE)

y = np.array([A_init, omega_init, lna_init])
traj_lna = [lna_init]
traj_A = [A_init]
traj_omega = [omega_init]
traj_h = [h_init]

target_dlna = 0.005
max_steps = 200_000
y_today = None
for step in range(max_steps):
    h_now = math.sqrt(
        OMEGA_M * math.exp(-3 * y[2]) + OMEGA_R * math.exp(-4 * y[2])
        + y[1]**2 / 6.0 + Y_pot(y[0], BETA_TILDE)
    )
    dtau = target_dlna / max(h_now, 1e-10)
    k1 = rhs_v3(y)
    k2 = rhs_v3(y + 0.5 * dtau * k1)
    k3 = rhs_v3(y + 0.5 * dtau * k2)
    k4 = rhs_v3(y + dtau * k3)
    y_new = y + (dtau / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    if y_new[2] >= 0.0:
        frac = -y[2] / (y_new[2] - y[2])
        y_today = y + frac * (y_new - y)
        h_today_num = math.sqrt(
            OMEGA_M + OMEGA_R + y_today[1]**2 / 6.0
            + Y_pot(y_today[0], BETA_TILDE)
        )
        traj_lna.append(0.0)
        traj_A.append(y_today[0])
        traj_omega.append(y_today[1])
        traj_h.append(h_today_num)
        break
    y = y_new
    traj_lna.append(y[2])
    traj_A.append(y[0])
    traj_omega.append(y[1])
    traj_h.append(h_now)

traj_lna = np.array(traj_lna)
traj_A = np.array(traj_A)
traj_omega = np.array(traj_omega)
traj_h = np.array(traj_h)
a_num = np.exp(traj_lna)

A_today_tracking = A_eq_at(1.0, BETA_TILDE)
A_today_numerical = traj_A[-1]
print(f"  A(today) numerical    = {A_today_numerical:.8f}")
print(f"  A(today) tracking     = {A_today_tracking:.8f}")
print(f"  abs diff              = {abs(A_today_numerical - A_today_tracking):.3e}")
print(f"  h(today) numerical    = {traj_h[-1]:.6f}")
print(f"  h(today) tracking     = {math.sqrt(h2_tracking(1.0, BETA_TILDE)):.6f}")
print(f"  steps integrated      = {len(traj_lna)}")

# Numerical Omega budget today
A_n = traj_A[-1]
omega_n = traj_omega[-1]
Y_n = Y_pot(A_n, BETA_TILDE)
Om_kin_n = omega_n**2 / 6.0
h2_n = OMEGA_M + OMEGA_R + Om_kin_n + Y_n
print("\n--- Full numerical Omega budget at z=0 ---")
print(f"  A(today)              = {A_n:.6f}")
print(f"  omega(today)          = {omega_n:+.6f}")
print(f"  Omega_DE_potential    = {Y_n:.6f}")
print(f"  Omega_DE_kinetic      = {Om_kin_n:.6f}")
print(f"  Omega_DE_total (num)  = {Y_n + Om_kin_n:.6f}")
print(f"  h^2(today, num)       = {h2_n:.6f}")

# ===================================================================
# Plotting
# ===================================================================
fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# Panel A: H(z)
ax = axes[0, 0]
ax.plot(z_grid, h_grid * H0_KMSMPC, lw=2.2, color='C0',
        label='Model-A V_3 (tracking)')
ax.set_xlabel('Redshift z')
ax.set_ylabel('H(z) [km/s/Mpc]')
ax.set_title(f'Hubble parameter (Model-A, H_0 = {H0_KMSMPC})')
ax.legend(loc='upper left')
ax.grid(alpha=0.3)

# Panel B: A_eq(a) trajectory
ax = axes[0, 1]
ax.plot(a_grid, [A_eq_at(a, BETA_TILDE) for a in a_grid], lw=2.2, color='C2',
        label='A_eq(a) tracking')
ax.axhline(A0, ls=':', color='gray', alpha=0.7,
           label=f'A_0 = 1/(12π) = {A0:.6f}')
ax.set_xscale('log')
ax.set_xlabel('Scale factor a')
ax.set_ylabel('A_eq')
ax.set_title('Tracking equilibrium A_eq(a) — approaches A_0 as a → ∞')
ax.legend()
ax.grid(alpha=0.3)

# Panel C: density fractions over z
ax = axes[1, 0]
Om_frac = np.array([h2_components(a, BETA_TILDE)['Omega_m'] / h2_components(a, BETA_TILDE)['h2']
                    for a in a_grid])
pot_frac = np.array([h2_components(a, BETA_TILDE)['Omega_DE_pot'] / h2_components(a, BETA_TILDE)['h2']
                     for a in a_grid])
kin_frac = np.array([h2_components(a, BETA_TILDE)['Omega_DE_kin'] / h2_components(a, BETA_TILDE)['h2']
                     for a in a_grid])
ax.fill_between(z_grid, 0, Om_frac, label='Ω_m(z)', alpha=0.65, color='C3')
ax.fill_between(z_grid, Om_frac, Om_frac + pot_frac,
                label='Ω_DE_pot(z)', alpha=0.65, color='C0')
ax.fill_between(z_grid, Om_frac + pot_frac, Om_frac + pot_frac + kin_frac,
                label='Ω_DE_kin(z)', alpha=0.65, color='C2')
ax.set_xlabel('Redshift z')
ax.set_ylabel('Fractional density')
ax.set_title('Model-A V_3 density components')
ax.legend(loc='upper right')
ax.grid(alpha=0.3)
ax.set_ylim(0, 1.05)
ax.set_xlim(0, 3)

# Panel D: w_eff(z)
ax = axes[1, 1]
ax.plot(z_we, w_eff, lw=2.2, color='C4', label='Model-A V_3 w_eff(z)')
ax.axhline(-1, ls=':', color='k', alpha=0.5, label='w = -1 reference')
ax.set_xlabel('Redshift z')
ax.set_ylabel('w_eff(z)')
ax.set_title('Effective dark-energy equation of state')
ax.legend(loc='lower right')
ax.grid(alpha=0.3)
ax.set_xlim(0, 3)

plt.suptitle(
    f'Model-A V_3 modified Friedmann | A_0 = 1/(12π), '
    f'β̃ = {BETA_TILDE:.4f} (calibrated by h²(today) = 1) | '
    f'Ω_DE_total = {comp_today["Omega_DE_total"]:.4f}',
    fontsize=12
)
plt.tight_layout(rect=[0, 0, 1, 0.96])

outdir = Path('reports/G8a')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / 'v3_modified_friedmann_summary.png', dpi=130,
            bbox_inches='tight')
plt.savefig(outdir / 'v3_modified_friedmann_summary.pdf', bbox_inches='tight')

# Tracking-vs-numerical check
fig2, ax2 = plt.subplots(figsize=(9, 6))
ax2.plot(np.log10(a_num), traj_A, lw=2.2, color='C0',
         label='Full numerical KG (RK4)')
A_eq_grid = np.array([A_eq_at(a, BETA_TILDE) for a in a_num])
ax2.plot(np.log10(a_num), A_eq_grid, lw=2.2, ls='--', alpha=0.85, color='C1',
         label='Tracking ansatz A_eq(a)')
ax2.axhline(A0, ls=':', color='gray', alpha=0.6,
            label=f'A_0 = 1/(12π) = {A0:.6f}')
ax2.set_xlabel('log10(a)')
ax2.set_ylabel('A')
ax2.set_title('V_3 numerical KG vs tracking ansatz')
ax2.legend(loc='best')
ax2.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(outdir / 'v3_tracking_attractor_check.png', dpi=130,
            bbox_inches='tight')

# CSV summary
import csv
with open(outdir / 'v3_today_budget.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['quantity', 'value', 'units', 'source'])
    w.writerow(['A_0', f'{A0:.10f}', 'dimensionless', 'G7 commitment 1/(12 pi)'])
    w.writerow(['alpha/beta', f'{ALPHA_OVER_BETA:.10e}', 'dimensionless',
                'V_3 structural minimum at A_0'])
    w.writerow(['beta_tilde', f'{BETA_TILDE:.10f}', 'dimensionless',
                'calibrated: h^2(today) = 1 closure'])
    w.writerow(['alpha_tilde', f'{ALPHA_TILDE:.10e}', 'dimensionless',
                'derived: alpha_tilde = beta_tilde * (alpha/beta)'])
    w.writerow(['Omega_m', f'{OMEGA_M}', 'dimensionless', 'input'])
    w.writerow(['Omega_r', f'{OMEGA_R}', 'dimensionless', 'input'])
    w.writerow(['A_eq(today)', f'{comp_today["A_eq"]:.10f}', 'dimensionless',
                'tracking equilibrium today'])
    w.writerow(['omega(today)', f'{comp_today["omega"]:+.10f}', 'dimensionless',
                'tracking velocity today'])
    w.writerow(['Omega_DE_pot(today)', f'{comp_today["Omega_DE_pot"]:.6f}',
                'dimensionless', 'V_3(A_eq) at a=1'])
    w.writerow(['Omega_DE_kin(today)', f'{comp_today["Omega_DE_kin"]:.6f}',
                'dimensionless', 'omega^2/6 at a=1'])
    w.writerow(['Omega_DE_total(today)', f'{comp_today["Omega_DE_total"]:.6f}',
                'dimensionless', 'pot + kin at a=1'])
    w.writerow(['h^2(today)', f'{comp_today["h2"]:.6f}', 'dimensionless',
                'closure: should be 1.000000'])
    w.writerow(['w_eff(z=0)', f'{w_today:+.6f}', 'dimensionless', 'effective EoS'])
    w.writerow(['w_eff(z=1)', f'{w_z1:+.6f}', 'dimensionless', 'effective EoS'])
    w.writerow(['w_eff(z=3)', f'{w_z3:+.6f}', 'dimensionless', 'effective EoS'])

print(f"\nPlots and CSV saved to {outdir.resolve()}")
print("=" * 72)

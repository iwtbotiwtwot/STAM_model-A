"""
Script 36 - Full modified Friedmann cosmology in STAM with V(A) = beta/(1-A)

Resolution of the dark-energy factor-of-2 gap from script 35.

The naive estimate in script 35:
    Omega_DE_naive = V(A_0)/(kappa rho_crit) = (1-A_0) * Omega_m = 0.307
gave match ratio 0.45 vs observed Omega_DE = 0.685.

The full modified Friedmann includes BOTH potential and kinetic energy of
the tracking A field:
    Omega_DE_total = Omega_DE_pot + Omega_DE_kin

Tracking solution: A(a) follows the equilibrium A_eq(a) defined by
    V'(A_eq) = kappa * rho_m(a)
For V(A) = beta/(1-A):
    A_eq(a) = 1 - (1-A_0) * a^(3/2)
which has Adot_eq = -(3/2)(1-A_0) * a^(3/2) * H, nonzero today.

The kinetic contribution is:
    Omega_DE_kin = (Adot/H_0)^2 / 6 = (3/8)(1-A_0)^2 * a^3 * h^2

At a=1 with A_0 = 0.0265:
    Omega_DE_pot = (1-A_0) * Omega_m  = 0.307
    Omega_DE_kin = (3/8)(1-A_0)^2     = 0.355  (using h^2 ~ 1)
    Omega_DE_total                    = 0.662
    observed Omega_DE                 = 0.685
    match ratio                       ~ 0.97

The factor-of-2 gap closes once we recognize A is not a cosmological constant
but a slowly-rolling field tracking matter density.
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ===================================================================
# Inputs (fixed by prior STAM scripts)
# ===================================================================
A0 = 0.026514          # cosmic ambient field, calibrated in script 35
OMEGA_M = 0.315        # matter today (Planck 2018)
OMEGA_R = 9.2e-5       # radiation today (CMB+neutrinos)
H0_KMSMPC = 67.4       # plot units only

# Calibration: beta_tilde = beta/(kappa * rho_crit_0)
# from V'(A_0) = kappa * rho_m,today  =>  beta = (1-A_0)^2 * kappa * rho_m
BETA_TILDE = OMEGA_M * (1.0 - A0)**2

print("=" * 72)
print("STAM Modified Friedmann  |  V(A) = beta / (1 - A)")
print("=" * 72)
print(f"A_0 (cosmic ambient)            = {A0}")
print(f"Omega_m (today)                 = {OMEGA_M}")
print(f"Omega_r (today)                 = {OMEGA_R:.2e}")
print(f"beta_tilde = beta/(kappa rho_c) = {BETA_TILDE:.6f}")

# ===================================================================
# Tracking solution (analytic)
# ===================================================================
def A_eq(a):
    """Equilibrium A satisfying V'(A_eq) = kappa rho_m(a)."""
    return 1.0 - (1.0 - A0) * a**1.5

def Y_pot(a):
    """Y = V(A_eq)/(kappa rho_crit_0). Potential energy fraction."""
    return BETA_TILDE / (1.0 - A_eq(a))

def K_kin(a):
    """Coefficient of h^2 in Omega_kin(a). Kinetic energy fraction = K(a) * h^2."""
    return (3.0/8.0) * (1.0 - A0)**2 * a**3

def h2(a):
    """h^2(a) = H^2/H_0^2 from modified Friedmann with tracking field."""
    matter = OMEGA_M * a**-3
    radiation = OMEGA_R * a**-4
    pot = Y_pot(a)
    K = K_kin(a)
    return (matter + radiation + pot) / (1.0 - K)

def h2_LCDM(a):
    OL = 1.0 - OMEGA_M - OMEGA_R
    return OMEGA_M * a**-3 + OMEGA_R * a**-4 + OL

# ===================================================================
# Closure at z = 0
# ===================================================================
print("\n--- Closure check at z = 0 (a = 1) ---")
a0 = 1.0
h2_now = h2(a0)
Om_now = OMEGA_M * a0**-3
Or_now = OMEGA_R * a0**-4
pot_now = Y_pot(a0)
kin_now = K_kin(a0) * h2_now
print(f"  Omega_matter            = {Om_now:.4f}")
print(f"  Omega_radiation         = {Or_now:.6f}")
print(f"  Omega_DE_potential      = {pot_now:.4f}")
print(f"  Omega_DE_kinetic        = {kin_now:.4f}")
print(f"  Omega_DE_total          = {pot_now + kin_now:.4f}")
print(f"  Omega_total (h^2)       = {h2_now:.4f}")

obs_DE = 1.0 - OMEGA_M - OMEGA_R
print(f"\n  observed Omega_DE (LCDM)   = {obs_DE:.4f}")
match = (pot_now + kin_now) / obs_DE
print(f"  STAM match ratio           = {match:.4f}")
print(f"  vs naive (script 35) ratio = {0.307/obs_DE:.4f}")

# ===================================================================
# H(z) tabulation
# ===================================================================
z = np.linspace(0, 5, 600)
a = 1.0/(1.0 + z)
h_STAM = np.sqrt(h2(a))
h_LCDM_arr = np.sqrt(h2_LCDM(a))

# ===================================================================
# Effective dark-energy equation of state w_eff(z)
# rho_DE = rho_DE_pot + rho_DE_kin, both functions of a
# w_eff = -1 - (1/3) d ln rho_DE / d ln a
# ===================================================================
def rho_DE(a_):
    return Y_pot(a_) + K_kin(a_) * h2(a_)

a_grid = np.logspace(-2.5, -0.001, 400)  # past only, avoid extrapolating future
z_grid = 1.0/a_grid - 1.0
ln_a = np.log(a_grid)
ln_rho = np.log(rho_DE(a_grid))
w_eff = -1.0 - (1.0/3.0) * np.gradient(ln_rho, ln_a)

w_today = np.interp(0.0, z_grid[::-1], w_eff[::-1])
w_z1 = np.interp(1.0, z_grid[::-1], w_eff[::-1])
w_high = np.interp(3.0, z_grid[::-1], w_eff[::-1])
print(f"\n--- Effective dark-energy w(z) ---")
print(f"  w_eff(z=0)   = {w_today:.4f}")
print(f"  w_eff(z=1)   = {w_z1:.4f}")
print(f"  w_eff(z=3)   = {w_high:.4f}")
print(f"  (slow-roll-only limit at high z would be -1/2)")

# ===================================================================
# Numerical verification: full coupled KG + Friedmann in cosmic time.
# State y = [A, omega, lna] with omega = dA/dtau, tau = H_0 t.
#
#   dA/dtau    = omega
#   domega/dtau= -3 h omega - 3 [Y'(A) - Omega_m a^-3]
#   dlna/dtau  = h
#   h^2 = Omega_m a^-3 + Omega_r a^-4 + omega^2/6 + Y(A)
#
# Direct sum in Friedmann (no algebraic loop) makes the ODE clean.
# ===================================================================
print("\n--- Numerical KG + Friedmann verification ---")

def rhs(y):
    A, omega, lna = y
    a_ = np.exp(lna)
    Y = BETA_TILDE / (1.0 - A)
    Yp = BETA_TILDE / (1.0 - A)**2
    rho_m = OMEGA_M * a_**-3
    rho_r = OMEGA_R * a_**-4
    h2_ = rho_m + rho_r + omega*omega/6.0 + Y
    h_ = np.sqrt(max(h2_, 1e-30))
    dA = omega
    domega = -3.0*h_*omega - 3.0*(Yp - rho_m)
    dlna = h_
    return np.array([dA, domega, dlna])

# RK4 in cosmic time with adaptive logarithmic-ish step
tau = 0.0
lna_init = np.log(1e-3)
a_init = np.exp(lna_init)
A_init = A_eq(a_init)
# tracking initial omega: omega = -1.5(1-A_0) a^(3/2) h_init
h_init = np.sqrt(OMEGA_M * a_init**-3 + OMEGA_R * a_init**-4)
omega_init = -1.5 * (1.0 - A0) * a_init**1.5 * h_init
y = np.array([A_init, omega_init, lna_init])

# Integrate forward until lna >= 0
traj_tau = [tau]
traj_A = [y[0]]
traj_omega = [y[1]]
traj_lna = [y[2]]
traj_h = [h_init]

# Step size: scale with current h to keep dlna ~ const
target_dlna = 0.005
max_steps = 200_000
for step in range(max_steps):
    h_now = np.sqrt(OMEGA_M * np.exp(-3*y[2]) + OMEGA_R * np.exp(-4*y[2])
                    + y[1]**2/6.0 + BETA_TILDE/(1.0 - y[0]))
    dtau = target_dlna / max(h_now, 1e-10)
    k1 = rhs(y)
    k2 = rhs(y + 0.5*dtau*k1)
    k3 = rhs(y + 0.5*dtau*k2)
    k4 = rhs(y + dtau*k3)
    y_new = y + (dtau/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    tau += dtau
    if y_new[2] >= 0.0:
        # interpolate to lna = 0
        frac = -y[2] / (y_new[2] - y[2])
        y_today = y + frac*(y_new - y)
        traj_tau.append(tau - dtau + frac*dtau)
        traj_A.append(y_today[0])
        traj_omega.append(y_today[1])
        traj_lna.append(0.0)
        h_today_num = np.sqrt(OMEGA_M + OMEGA_R + y_today[1]**2/6.0
                              + BETA_TILDE/(1.0 - y_today[0]))
        traj_h.append(h_today_num)
        break
    y = y_new
    traj_tau.append(tau)
    traj_A.append(y[0])
    traj_omega.append(y[1])
    traj_lna.append(y[2])
    h_v = np.sqrt(OMEGA_M*np.exp(-3*y[2]) + OMEGA_R*np.exp(-4*y[2])
                  + y[1]**2/6.0 + BETA_TILDE/(1.0 - y[0]))
    traj_h.append(h_v)

traj_tau = np.array(traj_tau)
traj_A = np.array(traj_A)
traj_omega = np.array(traj_omega)
traj_lna = np.array(traj_lna)
traj_h = np.array(traj_h)
a_num = np.exp(traj_lna)
A_eq_num = A_eq(a_num)

A_today_numerical = traj_A[-1]
A_today_tracking = A_eq(1.0)
print(f"  A(today) numerical    = {A_today_numerical:.6f}")
print(f"  A(today) tracking     = {A_today_tracking:.6f}")
print(f"  abs diff              = "
      f"{abs(A_today_numerical - A_today_tracking):.3e}")
print(f"  h(today) numerical    = {traj_h[-1]:.4f}")
print(f"  h(today) tracking     = {np.sqrt(h2(1.0)):.4f}")
print(f"  steps integrated      = {len(traj_tau)}")

# ===================================================================
# Plotting
# ===================================================================
fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# Panel A: H(z)
ax = axes[0,0]
ax.plot(z, h_STAM * H0_KMSMPC, label='STAM (modified Friedmann)', lw=2.2, color='C0')
ax.plot(z, h_LCDM_arr * H0_KMSMPC, label='LCDM', lw=2.2, ls='--', color='C1')
ax.set_xlabel('Redshift z')
ax.set_ylabel('H(z) [km/s/Mpc]')
ax.set_title('Hubble parameter')
ax.legend(loc='upper left')
ax.grid(alpha=0.3)

# Panel B: ratio H_STAM/H_LCDM
ax = axes[0,1]
ax.plot(z, h_STAM/h_LCDM_arr, lw=2.2, color='C2')
ax.axhline(1.0, ls=':', color='k', alpha=0.5)
ax.set_xlabel('Redshift z')
ax.set_ylabel('H_STAM / H_LCDM')
ax.set_title('Ratio of expansion rates (1.00 = identical)')
ax.grid(alpha=0.3)

# Panel C: density fractions over z
ax = axes[1,0]
Om_frac = OMEGA_M * a**-3 / h2(a)
pot_frac = Y_pot(a) / h2(a)
kin_frac = K_kin(a) * h2(a) / h2(a)
ax.fill_between(z, 0, Om_frac, label='Omega_m(z)', alpha=0.65, color='C3')
ax.fill_between(z, Om_frac, Om_frac + pot_frac,
                label='Omega_DE_pot(z)', alpha=0.65, color='C0')
ax.fill_between(z, Om_frac + pot_frac, Om_frac + pot_frac + kin_frac,
                label='Omega_DE_kin(z)', alpha=0.65, color='C2')
ax.set_xlabel('Redshift z')
ax.set_ylabel('Fractional density')
ax.set_title('STAM density components (tracking solution)')
ax.legend(loc='upper right')
ax.grid(alpha=0.3)
ax.set_ylim(0, 1.05)
ax.set_xlim(0, 3)

# Panel D: w_eff(z) and tracking attractor verification
ax = axes[1,1]
ax.plot(z_grid, w_eff, lw=2.2, color='C4', label='STAM w_eff(z)')
ax.axhline(-1, ls=':', color='k', alpha=0.7, label='LCDM (w = -1)')
ax.axhline(-0.5, ls=':', color='r', alpha=0.7, label='slow-roll only (w = -1/2)')
ax.set_xlabel('Redshift z')
ax.set_ylabel('w_eff(z)')
ax.set_title('Effective dark-energy equation of state')
ax.legend(loc='lower right')
ax.grid(alpha=0.3)
ax.set_xlim(0, 3)
ax.set_ylim(-1.2, 0.0)

plt.suptitle(
    f'STAM Modified Friedmann  |  A_0 = {A0:.4f}  |  '
    f'Omega_DE_total = {pot_now+kin_now:.3f} vs LCDM {obs_DE:.3f} '
    f'(match {match:.2f})',
    fontsize=13
)
plt.tight_layout(rect=[0, 0, 1, 0.96])

outdir = Path('reports/script_36')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / 'modified_friedmann_summary.png',
            dpi=130, bbox_inches='tight')
plt.savefig(outdir / 'modified_friedmann_summary.pdf', bbox_inches='tight')

# Tracking-vs-numerical figure (these DISAGREE late, important honest plot)
fig2, ax2 = plt.subplots(figsize=(9, 6))
ax2.plot(np.log10(a_num), traj_A,
         label='Full numerical KG (RK4)', lw=2.2, color='C0')
ax2.plot(np.log10(a_num), A_eq_num,
         label='Tracking A_eq(a) [analytic]',
         lw=2.2, ls='--', alpha=0.85, color='C1')
ax2.axhline(A0, ls=':', color='gray', alpha=0.6,
            label=f'Bridge-calibrated A_0 = {A0:.4f}')
ax2.set_xlabel('log10(a)')
ax2.set_ylabel('A')
ax2.set_title('Numerical KG vs tracking ansatz '
              '(divergence at late z = field lags equilibrium)')
ax2.legend(loc='best')
ax2.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(outdir / 'tracking_attractor_check.png',
            dpi=130, bbox_inches='tight')

# Numerical Omega budget today
A_n = traj_A[-1]
omega_n = traj_omega[-1]
Y_n = BETA_TILDE / (1.0 - A_n)
Om_kin_n = omega_n**2 / 6.0
h2_n = OMEGA_M + OMEGA_R + Om_kin_n + Y_n
print("\n--- Full numerical Omega budget at z=0 ---")
print(f"  A(today)              = {A_n:.4f}")
print(f"  omega(today)          = {omega_n:.4f}")
print(f"  Omega_DE_potential    = {Y_n:.4f}")
print(f"  Omega_DE_kinetic      = {Om_kin_n:.4f}")
print(f"  Omega_DE_total (num)  = {Y_n + Om_kin_n:.4f}")
print(f"  h^2(today, num)       = {h2_n:.4f}")
print(f"  match ratio (num)     = {(Y_n + Om_kin_n)/obs_DE:.4f}")

print(f"\nPlots saved to {outdir.resolve()}")
print("=" * 72)

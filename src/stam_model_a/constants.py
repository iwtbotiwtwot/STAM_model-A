"""Constants and historical candidate parameters for STAM_model-A."""

from __future__ import annotations

# SI constants
G_SI: float = 6.67430e-11  # m^3 kg^-1 s^-2
C_SI: float = 299_792_458.0  # m s^-1

# STAM distance scale constants carried from current notes.
C_MLY_PER_MPC: float = 3.261563776
H_STAM: float = 0.000243635
L_STAM: float = C_MLY_PER_MPC / H_STAM

# Historical candidate catalog/path calibration values.
B_PANTHEON_UNION_STYLE: float = 354.95
B_ORIGINAL_RETAINED: float = 461.3626922
B_DES_STYLE: float = 1335.412792

B_CANDIDATES: dict[str, float] = {
    "pantheon_union_style": B_PANTHEON_UNION_STYLE,
    "original_retained": B_ORIGINAL_RETAINED,
    "des_style": B_DES_STYLE,
}

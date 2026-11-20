"""Introductory example of thermodynamic functions for a selected material.

This script shows the two basic thermodynamic interfaces in STCHLab:
- `funR` for thermal reduction
- `funOX` for re-oxidation

The goal is not to cover every available function, but to show a few common
tasks:
- solve an equilibrium state,
- check the inverse relationship between variables,
- inspect enthalpy,
- and see how an applied voltage shifts the equilibrium.
"""

from stchlab import STCHLab, CeO2

material = CeO2
stch = STCHLab(material)

# Useful reference: `delta_max` is the maximum extent of reduction for the material.
print(f"Thermodynamic functions with {material.name}")
print(f"delta_max = {material.delta_max:.3f}")
print()

# ---------------------------------------------------------------------------
# 1) Thermal reduction
# ---------------------------------------------------------------------------
# Given delta, the extent of reduction, and pO2, solve for the equilibrium
# temperature and then recover the same state through the inverse functions.
delta = 0.01
pO2 = 100.0

T_R = stch.funR.T(delta, pO2)
G_R = stch.funR.G_R(delta, T_R, pO2)
H_R = stch.funR.H_R(delta, T_R)
delta_R = stch.funR.delta(T_R, pO2)
pO2_R = stch.funR.pO2(delta, T_R)

V_R = 0.50
T_R_V = stch.funR.T(delta, pO2, V_R)

# ---------------------------------------------------------------------------
# 2) Re-oxidation
# ---------------------------------------------------------------------------
# Given delta, the extent of reduction, and product molar fraction x, solve for
# the equilibrium temperature and verify the inverse functions.
x = 0.01

T_OX = stch.funOX.T(delta, x)
G_OX = stch.funOX.G_OX(delta, T_OX, x)
H_OX = stch.funOX.H_OX(delta, T_OX)
delta_OX = stch.funOX.delta(T_OX, x)
x_OX = stch.funOX.x(delta, T_OX)

# Thermodynamic lower bound for the product molar fraction at this temperature.
x_lim_900C = stch.funOX.x_lim(900 + 273.15)
V_OX = -0.05
T_OX_V = stch.funOX.T(delta, x, V_OX)

print("Reduction")
print(f"T = {T_R:.3f} K")
print(f"G_R = {G_R:.6e} J/mol")
print(f"H_R = {H_R:.3f} J/mol")
print(f"delta(T, pO2) = {delta_R:.6f}")
print(f"pO2(delta, T) = {pO2_R:.6f} Pa")
print(f"T with V = {V_R:.2f} V -> {T_R_V:.3f} K")
print()
print("Re-oxidation")
print(f"T = {T_OX:.3f} K")
print(f"G_OX = {G_OX:.6e} J/mol")
print(f"H_OX = {H_OX:.3f} J/mol")
print(f"delta(T, x) = {delta_OX:.6f}")
print(f"x(delta, T) = {x_OX:.6f}")
print(f"x_lim at 900 C = {x_lim_900C:.6f}")
print(f"T with V = {V_OX:.2f} V -> {T_OX_V:.3f} K")

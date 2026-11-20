"""Minimal example of counter-current models for a selected material.

The script is split into two parts:
- reactor-level calculations with `funCCR` and `funCCOX`
- full system closure and efficiency with `funCCSYS`

It is meant to show the progression from local reactor quantities such as
`T_R_max` or `T_OX_min` to complete cycle metrics such as `Q_total` and `STH`.
"""

from stchlab import STCHLab, CeO2

material = CeO2
stch = STCHLab(material)

# Minimum complete specification for the counter-current system. These six
# inputs are enough to solve the coupled reduction/oxidation cycle.
data = {
    "delta_OX": 0.01,
    "delta_R": 0.04,
    "pO2_min": 10.0,
    "pO2_max": 1000.0,
    "x_min": 0.001,
    "x_max": 0.1,
}

T_R_max = stch.funCCR.T_R_max(
    data["delta_OX"], data["delta_R"], data["pO2_min"], data["pO2_max"]
)
lmbda = stch.funCCR.lmbda(data["pO2_min"], data["pO2_max"])

T_OX_min = stch.funCCOX.T_OX_min(
    data["delta_OX"], data["delta_R"], data["x_min"], data["x_max"]
)
theta = stch.funCCOX.theta(data["x_min"], data["x_max"])

# Solve the complete counter-current cycle and then evaluate its energy balance.
state = stch.funCCSYS.system(data.copy(), use_db=False)
result = stch.funCCSYS.STH(state, type="Classic")

print(f"Counter-current models with {material.name}")
print()
print("Reactor level")
print("Reduction reactor")
print("Inputs")
print(f"  delta_OX = {data['delta_OX']:.6f}")
print(f"  delta_R = {data['delta_R']:.6f}")
print(f"  pO2_min = {data['pO2_min']:.6f} Pa")
print(f"  pO2_max = {data['pO2_max']:.6f} Pa")
print("Outputs")
print(f"  T_R_max = {T_R_max:.3f} K")
print(f"  lmbda = {lmbda:.6f}")
print()
print("Re-oxidation reactor")
print("Inputs")
print(f"  delta_OX = {data['delta_OX']:.6f}")
print(f"  delta_R = {data['delta_R']:.6f}")
print(f"  x_min = {data['x_min']:.6f}")
print(f"  x_max = {data['x_max']:.6f}")
print("Outputs")
print(f"  T_OX_min = {T_OX_min:.3f} K")
print(f"  theta = {theta:.6f}")
print()
print("System level")
print("Minimum complete system inputs (6 parameters)")
print(f"  delta_OX = {data['delta_OX']:.6f}")
print(f"  delta_R = {data['delta_R']:.6f}")
print(f"  pO2_min = {data['pO2_min']:.6f} Pa")
print(f"  pO2_max = {data['pO2_max']:.6f} Pa")
print(f"  x_min = {data['x_min']:.6f}")
print(f"  x_max = {data['x_max']:.6f}")
print("Outputs")
print(f"  T_R_max = {state['T_R_max']:.3f} K")
print(f"  T_OX_min = {state['T_OX_min']:.3f} K")
print(f"  lmbda = {state['lmbda']:.6f}")
print(f"  theta = {state['theta']:.6f}")
print(f"  Q_total = {result['Q_total'] / 1e3:.3f} kJ/mol H2")
print(f"  STH = {100 * result['STH']:.3f} %")

# Display the energy-balance pie chart for the counter-current system.
stch.funCCSYS.pie(result)

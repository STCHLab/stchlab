"""Minimal example of fixed-bed models for a selected material.

The script first runs one reduction step and one oxidation step at reactor level.
It then closes the full fixed-bed cycle with `funFBSYS.system(...)` and evaluates
the system performance with `funFBSYS.STH(...)`.
"""

from stchlab import STCHLab, CeO2

material = CeO2
stch = STCHLab(material)

# Reduction-step inputs: gas flow, oxide inventory, kinetic cap, and duration.
reduction_data = {
    "pO2_min": 0.2,
    "T_R": 900 + 273.15,
    "V_R": 0.75,
    "ndot_inert": 7.44e-5,
    "n_MOx": 1.0 / 172.12,
    "k_lim_R": 5e-5,
    "time_R": 600,
    "delta_OX": 0.001,
}

reduction_result = stch.funFBR.step(reduction_data)

# Oxidation-step inputs use the reduced state returned by the previous step.
oxidation_data = {
    "x_min": 1e-6,
    "T_OX": 900 + 273.15,
    "ndot_tot": 7.44e-5,
    "n_MOx": 1.0 / 172.12,
    "k_lim_OX": 1e-4,
    "time_OX": 300,
    "delta_R": reduction_result["delta_R"],
}

oxidation_result = stch.funFBOX.step(oxidation_data)

# This dictionary contains both reactor conditions and the information needed to
# close the full fixed-bed cycle consistently.
system_data = {
    "pO2_min": 0.2,
    "T_R": 900 + 273.15,
    "V_R": 0.75,
    "ndot_inert": 7.44e-5,
    "n_MOx": 1.0 / 172.12,
    "x_min": 1e-6,
    "T_OX": 900 + 273.15,
    "ndot_tot": 7.44e-5,
    "k_lim_R": 5e-5,
    "k_lim_OX": 1e-4,
    "time_R": 600,
    "time_OX": 300,
}

state = stch.funFBSYS.system(system_data.copy(), use_db=False)
result = stch.funFBSYS.STH(state)
nH2 = state["n_MOx"] * (state["delta_R"] - state["delta_OX"])
Q_total_kJ_molH2 = result["Q_total"] / nH2 / 1e3

print(f"Fixed-bed models with {material.name}")
print()
print("Reactor level")
print("Reduction reactor")
print("Inputs")
print(f"  pO2_min = {reduction_data['pO2_min']:.6f} Pa")
print(f"  T_R = {reduction_data['T_R']:.3f} K")
print(f"  V_R = {reduction_data['V_R']:.3f} V")
print(f"  ndot_inert = {reduction_data['ndot_inert']:.6e} mol/s")
print(f"  n_MOx = {reduction_data['n_MOx']:.6e} mol")
print(f"  k_lim_R = {reduction_data['k_lim_R']:.6e}")
print(f"  time_R = {reduction_data['time_R']:.3f} s")
print(f"  delta_OX = {reduction_data['delta_OX']:.6f}")
print("Outputs")
print(f"  delta_R = {reduction_result['delta_R']:.6f}")
print(f"  t_R = {reduction_result['t_R']:.3f} s")
print()
print("Re-oxidation reactor")
print("Inputs")
print(f"  x_min = {oxidation_data['x_min']:.6f}")
print(f"  T_OX = {oxidation_data['T_OX']:.3f} K")
print(f"  ndot_tot = {oxidation_data['ndot_tot']:.6e} mol/s")
print(f"  n_MOx = {oxidation_data['n_MOx']:.6e} mol")
print(f"  k_lim_OX = {oxidation_data['k_lim_OX']:.6e}")
print(f"  time_OX = {oxidation_data['time_OX']:.3f} s")
print(f"  delta_R = {oxidation_data['delta_R']:.6f}")
print("Outputs")
print(f"  delta_OX = {oxidation_result['delta_OX']:.6f}")
print(f"  t_OX = {oxidation_result['t_OX']:.3f} s")
print()
print("System level")
print("Minimum complete system inputs")
print(f"  pO2_min = {system_data['pO2_min']:.6f} Pa")
print(f"  T_R = {system_data['T_R']:.3f} K")
print(f"  V_R = {system_data['V_R']:.3f} V")
print(f"  ndot_inert = {system_data['ndot_inert']:.6e} mol/s")
print(f"  n_MOx = {system_data['n_MOx']:.6e} mol")
print(f"  x_min = {system_data['x_min']:.6f}")
print(f"  T_OX = {system_data['T_OX']:.3f} K")
print(f"  ndot_tot = {system_data['ndot_tot']:.6e} mol/s")
print(f"  k_lim_R = {system_data['k_lim_R']:.6e}")
print(f"  k_lim_OX = {system_data['k_lim_OX']:.6e}")
print(f"  time_R = {system_data['time_R']:.3f} s")
print(f"  time_OX = {system_data['time_OX']:.3f} s")
print("Outputs")
print(f"  delta_R = {state['delta_R']:.6f}")
print(f"  delta_OX = {state['delta_OX']:.6f}")
print(f"  accumulated Theta = {state['Theta']:.6f}")
print(f"  accumulated Lmbda = {state['Lmbda']:.6f}")
print(f"  Q_total = {Q_total_kJ_molH2:.3f} kJ/mol H2")
print(f"  STH = {100 * result['STH']:.3f} %")

# Display the energy-balance pie chart for the fixed-bed system.
stch.funFBSYS.pie(result)

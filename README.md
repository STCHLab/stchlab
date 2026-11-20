# STCHLab

STCHLab (Solar Thermochemical Hydrogen in silico Laboratory) is a Python package for
the analysis of two-step non-stoichiometric solar thermochemical redox cycles for
hydrogen production.

The library provides a unified, material-agnostic framework that spans equilibrium
thermodynamics, reactor-level models, and system-level energy balances for both
counter-current and fixed-bed operation. It also supports thermo-electrochemical
operation through electrical bias during reduction and/or re-oxidation.

Materials are represented as self-contained modules with a minimal interface, so new
redox-active materials can be incorporated without modifying the core library.
Reusable code lives in the `stchlab/` package, while runnable examples stay in
`examples/`.

## Nomenclature

STCHLab uses variable names close to scientific notation, with underscores for
subscripts and qualifiers when plain text needs them.

**Variable Names**

| Name | Meaning | Example |
|---|---|---|
| `T` | temperature | `T_ref` |
| `G` | Gibbs energy | `G_R`, `G_OX`, `G_s` |
| `H` | enthalpy | `H_R`, `H_OX`, `H_s` |
| `Q` | heat or energy duty | `Q_total`, `Q_R`, `Q_sep_H2` |
| `W` | useful work or fuel heating value | `W_H2` |
| `V` | voltage or electrical bias | `V`, `V_R`, `V_OX` |
| `eta` | efficiency | `eta` |
| `nu` | recovery or subsystem efficiency | `nu_MOx`, `nu_sep`, `nu_opt` |
| `tau` | transmissivity | `tau` |
| `epsilon` | emissivity | `epsilon` |
| `alpha` | linear heat-loss coefficient | `alpha` |
| `C` | solar concentration ratio | `C` |
| `DNI` | direct normal irradiance | `DNI` |
| `mu` | chemical potential | `mu_solid` |
| `dmu` | derivative of chemical potential with respect to `delta` | `dmu_solid` |
| `p` | pressure | `psys_R` |
| `psys` | system pressure | `psys_R` |
| `pO2` | oxygen partial pressure | `pO2_min` |
| `pstar` | dimensionless oxygen pressure ratio | `pstar_min` |
| `x` | molar fraction | `x_inert` |
| `N` | number of discretization points | `N_R` |
| `MOx` | metal oxide material | `MOx` |
| `delta` | extent of reduction | `delta_min` |
| `lmbda` | sweep ratio | `lmbda` |
| `theta` | product yield | `theta` |
| `Lmbda` | fixed-bed accumulated/effective sweep metric | `Lmbda` |
| `Theta` | fixed-bed accumulated/effective product metric | `Theta` |
| `Ddelta` | reduction swing, `delta_R - delta_OX` | `Ddelta` |
| `ndot` | molar flow rate | `ndot_inert` |
| `k_lim` | kinetic rate limit | `k_lim_R` |
| `frac_Ddelta` | fraction of available fixed-bed delta swing | `frac_Ddelta_R` |

**Subscripts And Qualifiers**

| Subscript | Meaning | Example |
|---|---|---|
| `R` | thermal reduction step | `psys_R`, `N_R` |
| `OX` | re-oxidation step | `psys_OX`, `N_OX` |
| `s` | splitting reactant gas species | `G_s`, `H_s` , `mu_s`|
| `ref` | reference condition | `T_ref` |
| `inert` | inert gas | `x_inert` |
| `sys` | system-level quantity | `psys_R` |
| `amb` | ambient condition | `T_amb` |
| `sep` | separator condition or separation duty | `T_sep_H2`, `Q_sep_H2` |
| `in` | inlet value | `T_R_in` |
| `out` | outlet value | `T_OX_out` |
| `O2` | oxygen species | `pO2_min`, `xO2_max` |
| `H2O` | water or steam species | `ndot_H2O` |
| `H2` | hydrogen product species | `W_H2` |
| `H` | hot heat-exchanger stream | `T_H_in` |
| `C` | cold heat-exchanger stream | `T_C_out` |
| `DNI` | direct normal irradiance | `DNI` |
| `sol` | solar input | `Q_sol` |
| `rec` | receiver input | `Q_rec` |
| `aux` | auxiliary heat | `Q_aux` |
| `rej` | rejected heat | `Q_rej` |
| `ch` | chemical reaction duty | `Q_ch_R` |
| `ideal` | reversible or ideal value | `Q_R_ideal` |
| `used` | energy consumed by downstream uses | `Q_OX_used` |
| `0` | initial guess | `T_0`, `delta_0`, `pO2_0`, `x_0` |
| `lim` | thermodynamic limit | `x_lim` |
| `min` | minimum value | `delta_min` |
| `max` | maximum value | `delta_max` |
| `eff` | effective value used internally | `_delta_min_eff` |
| `MOx` | metal oxide material or phase | `n_MOx` |

Rules used throughout the library:

- `T` always means temperature.
- `N` means a count, for example number of discretization points.
- `pO2` is kept as `pO2`, not `p_O2`.
- Subscripts and qualifiers are written with `_`: `psys_R`, `x_inert`, `T_ref`.
- Function names follow the same subscript convention as variables: `G_R`, `H_OX`, `G_s`, `x_lim`.
- In `funCCR`, `pstar(pO2)` and `pO2(pstar)` are conversion functions between oxygen partial pressure and dimensionless oxygen pressure ratio.
- Chemical notation is preserved when useful: `MOx`, `H2O`, `CO2`, `OX`.
- A term can appear as a variable name and as a subscript if it has both roles, for example `MOx`.
- In `funCCSYS.STH`, `tau` is used for transmissivity so that `T` remains reserved for temperature.

**Counter-Current STH Keys**

The counter-current system energy balance returned by `funCCSYS.STH(...)` uses:

| Category | Keys |
|---|---|
| model and temperatures | `type`, `T_H_in`, `T_H_out`, `T_C_in`, `T_C_out`, `i_DT` |
| useful output and efficiencies | `W_H2`, `eta`, `STH` |
| total inputs and losses | `Q_total`, `Q_sol`, `Q_rec`, `Q_aux`, `Q_rej` |
| reduction and oxidation duties | `Q_R`, `Q_R_ideal`, `Q_OX`, `Q_OX_used` |
| reaction and material heating | `Q_ch_R`, `Q_ch_OX`, `Q_MOx`, `Q_MOx_R`, `Q_MOx_OX`, `Q_inert` |
| separations and steam handling | `Q_sep_H2_ideal`, `Q_sep_H2`, `Q_sep_inert_ideal`, `Q_sep_inert`, `Q_H2O_re`, `Q_H2O_new`, `Q_H2O_lat` |

## Installation

Open a terminal in the repository root before running any install command.
If you only want to run examples or write scripts that use STCHLab, follow
`For Users`. If you plan to edit the library code, follow `For Developers`.

### For Developers

1. Get the repository from GitHub.

If you use git:

- `git clone https://github.com/adlcalle/stchlab.git`

This creates a folder called `stchlab` in your current directory.

If you do not use git:

- download the repository as a ZIP file from GitHub,
- extract it,
- open the extracted folder, usually called `stchlab-main`.

2. Open a terminal inside that folder.

If you used `git clone`, this is usually:

- `cd stchlab`

If you downloaded the ZIP, go to the extracted folder name instead.

You must run the install command from the project root, that is, the folder that
contains `pyproject.toml`, `README.md`, `examples/`, and the inner `stchlab/`
package folder.

Do not go inside the inner `stchlab/` package folder before running the install
command.

3. If you want to use a virtual environment, create and activate it now using the
optional section below.

4. Use an editable install:

- `python -m pip install -e .`

This is the recommended option if you are developing `stchlab`, changing files inside
`stchlab/`, or creating new examples.

With this mode:

- you install the package once in your active environment,
- changes to the source code are picked up automatically,
- you do not need to reinstall after normal code edits.

Use the same Python environment in your IDE and in the terminal.

### For Users

1. Get the repository from GitHub.

If you use git:

- `git clone https://github.com/adlcalle/stchlab.git`

This creates a folder called `stchlab` in your current directory.

If you do not use git:

- download the repository as a ZIP file from GitHub,
- extract it,
- open the extracted folder, usually called `stchlab-main`.

2. Open a terminal inside that folder.

If you used `git clone`, this is usually:

- `cd stchlab`

If you downloaded the ZIP, go to the extracted folder name instead.

You must run the install command from the project root, that is, the folder that
contains `pyproject.toml`, `README.md`, `examples/`, and the inner `stchlab/`
package folder.

Do not go inside the inner `stchlab/` package folder before running the install
command.

3. If you want to use a virtual environment, create and activate it now using the
optional section below.

4. Use a regular install:

- `python -m pip install .`

This is the recommended option if you want to use the package, run the example
scripts, and build your own models without modifying the library code.

With this mode:

- the package is installed in the active environment,
- examples can import `stchlab` normally,
- if the library source code changes later, reinstall the package.

After installing the package, you can create your own working folder anywhere on your
computer and write your own models there; your scripts will be able to import
`stchlab` as long as you run them with the same Python environment where the package
was installed, whether that environment is virtual or not.

Quick sanity check after installation:

- `python -c "from stchlab import STCHLab, CeO2; print(STCHLab(CeO2).funR.T(0.01, 100.0))"`

This should print a finite equilibrium temperature.


### Optional: Create A Virtual Environment

This step is optional. You can skip it and install the package in your current Python
environment if you prefer.

Advantages of using a virtual environment:

- it keeps the dependencies for this project separate from other Python projects,
- it avoids conflicts between package versions,
- it makes it easier to keep a clean and reproducible setup,
- if something breaks, the environment can be deleted and created again easily.

Possible disadvantages:

- it adds one extra setup step before installing the package,
- you need to remember to activate the environment in each new terminal,
- your IDE should be configured to use the same environment.

Create the environment:

- `python -m venv .venv`

Activate it on Windows PowerShell:

- `.\.venv\Scripts\Activate.ps1`

Activate it on macOS/Linux:

- `source .venv/bin/activate`

Optional but recommended:

- `python -m pip install --upgrade pip`

## Examples

All runnable scripts live in `examples/`. After installing the package in the active
Python environment, they can be launched directly from the repository root.

Recommended order if you are new to the package:

1. Run `example_thermodynamic_functions.py` first. It introduces the basic
   thermodynamic functions.
2. Run `example_counter_current_models.py` if your reactor model is
   counter-current.
3. Run `example_fixed_bed_models.py` if your reactor model is fixed-bed.

From the repository root, run an example with:

- `python examples/example_thermodynamic_functions.py`
- `python examples/example_counter_current_models.py`
- `python examples/example_fixed_bed_models.py`

The counter-current and fixed-bed examples also display an energy-balance pie chart.
If a figure window opens, close it to let the script finish.

You do not need to keep your own scripts inside this repository. After installation,
you can create a working folder anywhere on your computer and write a script such as:

```python
from stchlab import STCHLab, CeO2

material = CeO2
stch = STCHLab(material)
```

This creates an STCHLab object for `CeO2`. Most examples then call methods on that
object, such as `stch.funR.T(...)` or `stch.funCCSYS.system(...)`.
By default, the inert sweep gas used in `STH` energy balances is CoolProp `N2`.
Use `inert` to choose another CoolProp fluid, for example:

```python
stch = STCHLab(CeO2, inert="Ar")
```

`STCHLab` supports `oxidant="H2O"` and `oxidant="CO2"`:

```python
co2_stch = STCHLab(CeO2, oxidant="CO2")
```

The system pressures `psys_R` and `psys_OX` must be finite positive values in Pa.
The inert fraction `x_inert` must satisfy `0 <= x_inert < 1`; zero is allowed,
but one is not.

The public API lives under `stchlab/`. Keep new runnable scripts in `examples/` if
they are meant to be shared with the project.

## Thermodynamic Functions

Start by creating an STCHLab object for the material of interest. In the examples
below, `CeO2` is used throughout.

```python
from stchlab import STCHLab, CeO2

material = CeO2
stch = STCHLab(material)
```

Use `funR` for thermal reduction and `funOX` for re-oxidation.
The snippets below build on one another, so run them in order in one script or
Python session.

Common variable names:

| Name | Meaning | Units |
| --- | --- | --- |
| `T` | temperature | K |
| `pO2` | oxygen partial pressure | Pa |
| `x` | product molar fraction in re-oxidation, for example H2 or CO | dimensionless |
| `delta` | extent of reduction, also sometimes called extent of non-stoichiometry | dimensionless |
| `V` | applied electrical bias | V |

The value `material.delta_max` is the largest extent of reduction allowed by the
material model:

```python
print(f"delta_max = {material.delta_max:.3f}")
```

### Reduction

For thermal reduction, choose `delta` and `pO2`, then solve for the equilibrium
temperature:

```python
delta = 0.01
pO2 = 100.0

T_R = stch.funR.T(delta, pO2)
print(f"T_R = {T_R:.3f} K")
```

At equilibrium, the Gibbs free energy expression should be close to zero:

```python
G_R = stch.funR.G_R(delta, T_R, pO2)
print(f"G_R = {G_R:.6e} J/mol")
```

You can also use inverse functions to recover the same state from different inputs:

```python
delta_check = stch.funR.delta(T_R, pO2)
pO2_check = stch.funR.pO2(delta, T_R)

print(f"delta(T_R, pO2) = {delta_check:.6f}")
print(f"pO2(delta, T_R) = {pO2_check:.6f} Pa")
```

If you include an electrical bias, pass it as the third argument:

```python
V_R = 0.50
T_R_V = stch.funR.T(delta, pO2, V_R)

print(f"T_R with V = {V_R:.2f} V -> {T_R_V:.3f} K")
```

### Re-Oxidation

For re-oxidation, use `x`, the product molar fraction, instead of `pO2`:

```python
x = 0.01

T_OX = stch.funOX.T(delta, x)
print(f"T_OX = {T_OX:.3f} K")
```

As above, you can evaluate the free energy and inverse functions:

```python
G_OX = stch.funOX.G_OX(delta, T_OX, x)
delta_check = stch.funOX.delta(T_OX, x)
x_check = stch.funOX.x(delta, T_OX)

print(f"G_OX = {G_OX:.6e} J/mol")
print(f"delta(T_OX, x) = {delta_check:.6f}")
print(f"x(delta, T_OX) = {x_check:.6f}")
```

To check the thermodynamic equilibrium lower bound for `x` at 900 C:

```python
x_lim_900C = stch.funOX.x_lim(900 + 273.15)
print(f"x_lim at 900 C = {x_lim_900C:.6f}")
```

For water splitting, `H2O -> H2 + 0.5 O2`, this is the equilibrium H2 molar
fraction at that temperature. For CO2 splitting, it is the equilibrium CO molar
fraction. In the presence of the reactant gas, product molar fractions below
`x_lim` are not thermodynamically reachable at that temperature.

For thermo-electrochemical re-oxidation, pass the voltage as the third argument:

```python
V_OX = -0.05
T_OX_V = stch.funOX.T(delta, x, V_OX)

print(f"T_OX with V = {V_OX:.2f} V -> {T_OX_V:.3f} K")
```

For a complete runnable version of these snippets, see
`examples/example_thermodynamic_functions.py`.

## Counter-Current Models

Using the same `stch = STCHLab(CeO2)` object, the counter-current helpers can
be used first at reactor level and then at full-system level.

The counter-current reduction and re-oxidation models implemented in STCHLab are
based on de la Calle et al. (2022) for thermal reduction and de la Calle et al.
(2024) for re-oxidation.

The minimum counter-current input dictionary has six values:

| Name | Meaning |
| --- | --- |
| `delta_OX` | extent of reduction after re-oxidation |
| `delta_R` | extent of reduction after reduction |
| `pO2_min` | lowest oxygen partial pressure in the reduction reactor |
| `pO2_max` | highest oxygen partial pressure in the reduction reactor |
| `x_min` | lower product molar fraction limit for re-oxidation |
| `x_max` | upper product molar fraction limit for re-oxidation |

```python
cc_data = {
    "delta_OX": 0.01,
    "delta_R": 0.04,
    "pO2_min": 10.0,
    "pO2_max": 1000.0,
    "x_min": 0.001,
    "x_max": 0.1,
}
```

At reactor level, `funCCR` handles the reduction side:

```python
T_R_max = stch.funCCR.T_R_max(
    cc_data["delta_OX"], cc_data["delta_R"], cc_data["pO2_min"], cc_data["pO2_max"]
)
lmbda = stch.funCCR.lmbda(cc_data["pO2_min"], cc_data["pO2_max"])

print(f"T_R_max = {T_R_max:.3f} K")
print(f"lmbda = {lmbda:.6f}")
```

The re-oxidation side uses `funCCOX`:

```python
T_OX_min = stch.funCCOX.T_OX_min(
    cc_data["delta_OX"], cc_data["delta_R"], cc_data["x_min"], cc_data["x_max"]
)
theta = stch.funCCOX.theta(cc_data["x_min"], cc_data["x_max"])

print(f"T_OX_min = {T_OX_min:.3f} K")
print(f"theta = {theta:.6f}")
```

In counter-current models, `lmbda` is the sweep ratio and `theta` is the
product yield defined by the inlet/outlet composition limits of the reduction
and re-oxidation reactors.
They are not aliases of the fixed-bed `Lmbda` and `Theta` values below.
`funCCSYS.system(...)` reports `ndot_H2O` as the inlet steam flow required by
the re-oxidation mass balance. When `x_inert > 0`, the total re-oxidation gas
flow also includes the inert fraction, and `funCCSYS.STH(...)` folds the inert
recycle sensible heat into `Q_H2O_re`.

For the full counter-current cycle, use `funCCSYS.system(...)` and then compute the
solar-to-hydrogen efficiency with `funCCSYS.STH(...)`:

```python
cc_state = stch.funCCSYS.system(cc_data.copy(), use_db=False)
cc_result = stch.funCCSYS.STH(cc_state, type="Classic")

print(f"Q_total = {cc_result['Q_total'] / 1e3:.3f} kJ/mol H2")
print(f"STH = {100 * cc_result['STH']:.3f} %")
```

To display the energy-balance pie chart:

```python
stch.funCCSYS.pie(cc_result)
```

The public example `example_counter_current_models.py` prints the same inputs and
outputs in a more readable block format.

## Fixed-Bed Models

The fixed-bed helpers follow the same logic: reactor-level calculations first, and then
cycle closure and efficiency at system level.

Fixed-bed inputs look longer than counter-current inputs because this model includes
gas flow, material inventory, kinetic limits, and step duration.
The snippets below build on one another.

For the reduction step, create an input dictionary and call `funFBR.step(...)`:

```python
fb_red = {
    "pO2_min": 0.2,
    "T_R": 900 + 273.15,
    "V_R": 0.75,
    "ndot_inert": 7.44e-5,
    "n_MOx": 1.0 / 172.12,
    "k_lim_R": 5e-5,
    "time_R": 600,
    "delta_OX": 0.001,
}

resR = stch.funFBR.step(fb_red)

print(f"delta_R = {resR['delta_R']:.6f}")
print(f"t_R = {resR['t_R']:.3f} s")
```

The oxidation step uses the reduced state from the previous step:

```python
fb_ox = {
    "x_min": 1e-6,
    "T_OX": 900 + 273.15,
    "ndot_tot": 7.44e-5,
    "n_MOx": 1.0 / 172.12,
    "k_lim_OX": 1e-4,
    "time_OX": 300,
    "delta_R": resR["delta_R"],
}

resOX = stch.funFBOX.step(fb_ox)

print(f"delta_OX = {resOX['delta_OX']:.6f}")
print(f"t_OX = {resOX['t_OX']:.3f} s")
```

For the full fixed-bed system, put the inputs for both steps in one dictionary:

```python
fb_data = {
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
```

Then close the cycle and compute performance:

```python
fb_state = stch.funFBSYS.system(fb_data.copy(), use_db=False)
fb_result = stch.funFBSYS.STH(fb_state)

fb_nH2 = fb_state["n_MOx"] * (fb_state["delta_R"] - fb_state["delta_OX"])
fb_Q_total = fb_result["Q_total"] / fb_nH2 / 1e3

print(f"delta_R = {fb_state['delta_R']:.6f}")
print(f"delta_OX = {fb_state['delta_OX']:.6f}")
print(f"Theta = {fb_state['Theta']:.6f}")
print(f"Lmbda = {fb_state['Lmbda']:.6f}")
print(f"Q_total = {fb_Q_total:.3f} kJ/mol H2")
print(f"STH = {100 * fb_result['STH']:.3f} %")
```

In fixed-bed models, `Lmbda` and `Theta` are accumulated/effective values from
the time-resolved reduction and re-oxidation histories, so they are intentionally
reported with separate names from counter-current `lmbda` and `theta`.

The fixed-bed system energy balance returned by `funFBSYS.STH(...)` uses:

| Category | Keys |
|---|---|
| useful output and efficiencies | `W_H2`, `STH`, `STH_available` |
| total inputs and solar duties | `Q_total`, `Q_total_1`, `Q_total_2`, `Q_total_3`, `Q_solar`, `Q_solar_R`, `Q_solar_W`, `Q_solar_C`, `Q_solar_available` |
| projected/receiver duties | `Q_proj`, `Q_proj_R`, `Q_proj_W`, `Q_proj_C`, `Q_rec_R` |
| reaction and material heating | `Q_ch_R`, `Q_ch_OX`, `Q_MOx`, `Q_inert`, `Q_inert_R`, `Q_inert_W`, `Q_inert_C` |
| separations and steam handling | `Q_sep_H2`, `Q_sep_inert`, `Q_H2O_re`, `Q_H2O_new`, `Q_H2O_lat` |
| losses, auxiliary heat, and diagnostics | `Q_aux`, `Q_rej_C`, `Q_rej_T`, `Q_loss`, `Q_loss_R`, `Q_loss_W`, `Q_loss_C`, `Q_loss_OX`, `Q_loss_2nd_law`, `Q_loss_circ_1`, `Q_loss_circ_2`, `Q_coll`, `Q_miss_C`, `Q_proc`, `Q_avail` |

To display the energy-balance pie chart:

```python
stch.funFBSYS.pie(fb_result)
```

The public example `example_fixed_bed_models.py` prints the reactor-level and
system-level inputs and outputs in the same style as the counter-current example.

The runnable scripts in `examples/` follow the same workflow; the counter-current
and fixed-bed examples also display the energy-balance pie chart, following the
visualization approach of Cheng and de la Calle et al. (2021).

## Database Cache

STCHLab can use SQLite cache files to reuse solved `system()` simulations.
This is useful because some cycle closures are expensive to solve. The cache stores
the solved reactor/system state returned by `funCCSYS.system(...)` or
`funFBSYS.system(...)`; `STH(...)` results are not cached because they are computed
from a solved state and the energy-balance parameters passed to `STH()`.

### What Gets Cached

The system-level helpers use the cache by default:

```python
cc_state = stch.funCCSYS.system(cc_data.copy())
fb_state = stch.funFBSYS.system(fb_data.copy())
```

There are two independent SQLite databases:

- `simulationsCC.db` for counter-current systems,
- `simulationsFB.db` for fixed-bed systems.

Each cache entry is keyed by the `system()` inputs plus the STCHLab context:
material, oxidant, `x_inert`, `psys_R`, `psys_OX`, `N_R`, and `N_OX`. Changing any
of those values creates a different cache key. The global `inert` setting is not
part of the `system()` cache key because it only affects `STH()` energy balances,
not the solved `system()` state.

### Default Cache Files

By default, new cache entries are stored in the user cache folder:

- `~/.stchlab/cache/simulationsCC.db`
- `~/.stchlab/cache/simulationsFB.db`

This keeps normal runs from modifying files inside the installed `stchlab`
package. The database files shipped inside `stchlab/` are reference data used by
the project tests and should not be treated as the live cache during normal use.

### Running Without Cache

To run without the database cache, pass `use_db=False`:

```python
cc_state = stch.funCCSYS.system(cc_data.copy(), use_db=False)
fb_state = stch.funFBSYS.system(fb_data.copy(), use_db=False)
```

This is often the best choice for quick debugging, benchmarking, or parallel runs
where you do not need to reuse results.

### Using A Custom Cache File

To use a specific cache file, pass `db_name`:

```python
cc_state = stch.funCCSYS.system(cc_data.copy(), db_name="my_counter_current_cache.db")
fb_state = stch.funFBSYS.system(fb_data.copy(), db_name="my_fixed_bed_cache.db")
```

If the folder for `db_name` does not exist, STCHLab creates it automatically.

### Changing The Default Cache Folder

To change the default cache folder for a whole run, set the environment variable
`STCHLAB_CACHE_DIR` before starting Python:

Windows PowerShell:

- `$env:STCHLAB_CACHE_DIR = "C:\Users\you\stchlab-cache"`

macOS/Linux:

- `export STCHLAB_CACHE_DIR="$HOME/stchlab-cache"`

### Parallel Sweeps

Large sweeps over many independent operating points can be parallelized outside
STCHLab. For maximum throughput, use `use_db=False` or give each worker a separate
`db_name`, then merge or post-process the results after the sweep.

If avoiding duplicate work matters more than raw write throughput, multiple workers
can share one SQLite cache file with `use_db=True`. STCHLab coordinates shared-cache
work with a small `simulation_claims` table: one worker claims a hash before
computing it, other workers wait for the result, and stale or failed claims are
handled explicitly.

The shared-cache controls are:

```python
fb_state = stch.funFBSYS.system(
    fb_data.copy(),
    db_name="shared_fixed_bed_cache.db",
    stale_seconds=3600,
    retry_failed=False,
    poll_seconds=5.0,
)
```

- `stale_seconds` should be longer than the maximum expected time for one
  simulation. If a worker dies while holding a claim, another worker can reclaim it
  after this time.
- `retry_failed=False` makes repeated failed simulations raise a clear error instead
  of silently retrying forever. Set it to `True` when you intentionally want to
  rerun failed points.
- `poll_seconds` controls how often a waiting worker checks whether the claimed
  simulation has finished.

SQLite is still a file-based database, not a client-server database. For extremely
large runs with heavy write traffic, one cache file per worker remains the most
robust option.

## Creating A New Material With CEF

The simplest way to add a new redox material is to create your own Python file and
define a class that inherits from `CEF`. This lets you generate the thermodynamic
functions required by `STCHLab` from a set of CEF coefficients.

The `CEF` helper in STCHLab is based on the two-parameter compound energy
formalism described by Miller et al. (2022), adapted here as a compact material
model that generates `Cp`, `G_R`, and `H_R` from a fitted set of coefficients.

The material definition needs:

- `name`: material name used by the library and database cache,
- `coef`: a NumPy array with exactly 9 CEF coefficients,
- `sigma_model`: choose `SIGMA_MO2` for fluorite-type `MO2` materials or
  `SIGMA_PEROVSKITE_B` for perovskites,
- `y`: required for the perovskite model; this is the material-composition
  parameter in `AA'B(1-y)B'yO3`, not the re-oxidation product molar fraction `x`,
- `delta_max`: optional maximum extent of reduction if you want to override the
  default value,
- `cp`: optional constant heat capacity in J/mol/K; if omitted, `CEF` uses `80`.

For example, create a file called `my_material.py` and start with the imports:

```python
import numpy as np

from stchlab.materials.CEF import CEF, SIGMA_PEROVSKITE_B
```

Then define the material class. Replace `name`, `coef`, and `y` with the values for
your material:

```python
class MyMaterial(CEF):
    name = "MyMaterial"
    sigma_model = SIGMA_PEROVSKITE_B
    coef = np.array([
        -13.8791047644328,
        778.671753977449,
        -7818.79533842228,
        13.3236010119209,
        1142.68242700902,
        -11852.2821288997,
        0,
        0,
        0,
    ])
    y = 6.7 / 8
```

Finally, expose the module-level names expected by `STCHLab(...)`:

```python
name = MyMaterial.name
coef = MyMaterial.coef
y = MyMaterial.y
delta_max = MyMaterial.delta_max
Cp = MyMaterial.Cp
G_R = MyMaterial.G_R
H_R = MyMaterial.H_R
```

This last block uses the same module-level interface as the built-in material files
such as `CeO2.py` and `BCZM.py`. That is what makes your material directly
compatible with `STCHLab(...)`.

For a fluorite-type `MO2` material, use `SIGMA_MO2` instead. In that case `y` is not
needed, and the maximum extent of reduction is `0.25` by default unless you
override `delta_max` explicitly.

Once your file is defined, you can use it from any working folder:

```python
import my_material

from stchlab import STCHLab

stch = STCHLab(my_material)
```

If you want a template, start from `stchlab/materials/CeO2.py` for an `MO2`
material or `stchlab/materials/BCZM.py` for a perovskite material and replace the
coefficients and metadata with your own.

## Tests

Regression tests live in `tests/` and can be run from the repository root.

- `python tests/test_fast.py`
- `python tests/test_full.py`
- `python tests/test_performance.py`

`test_fast.py` is the deterministic regression test. It checks the main
thermodynamic, reactor-level, system-level, and database-cache workflows with fixed
reference cases.

`test_full.py` is the randomized regression test. It samples operating conditions
inside safe ranges and repeats the same checks with random cases.

`test_performance.py` is a watchdog for large runtime regressions or hangs in
representative workflows; it is not intended as a precise benchmark.

Useful options for the randomized test:

- `python tests/test_full.py --cases 50 --seed 42`
- `python tests/test_full.py --cases 100 --seed 42`

Use `--seed` if you want the same random sequence again.

Useful options for the performance watchdog:

- `python tests/test_performance.py --cases 5 --seed 42 --timeout 20`

## References

- de la Calle, A., Ermanoski, I., and Stechel, E. B. (2022). Towards chemical
  equilibrium in thermochemical water splitting. Part 1: Thermal reduction.
  International Journal of Hydrogen Energy, 47, 10474-10482.
  https://doi.org/10.1016/j.ijhydene.2021.07.167

- de la Calle, A., Ermanoski, I., Miller, J. E., and Stechel, E. B. (2024).
  Towards chemical equilibrium in thermochemical water splitting. Part 2:
  Re-oxidation. International Journal of Hydrogen Energy, 72, 1159-1168.
  https://doi.org/10.1016/j.ijhydene.2024.05.298

- Cheng, W.-H., de la Calle, A., Atwater, H. A., Stechel, E. B., and Xiang, C.
  (2021). Hydrogen from Sunlight and Water: A Side-by-Side Comparison between
  Photoelectrochemical and Solar Thermochemical Water-Splitting. ACS Energy
  Letters, 6, 3096-3113. https://doi.org/10.1021/acsenergylett.1c00758

- Miller, J. E., Babiniec, S. M., Coker, E. N., Loutzenhiser, P. G., Stechel, E. B.,
  and Ambrosini, A. (2022). Modified Calcium Manganites for Thermochemical
  Energy Storage Applications. Frontiers in Energy Research, 10, 774099.
  https://doi.org/10.3389/fenrg.2022.774099

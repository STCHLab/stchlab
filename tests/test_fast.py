"""Fast regression test for the main STCHLab workflows.

Run from the repository root with:

    python tests/test_fast.py
"""

from pathlib import Path
import shutil
import sqlite3
import sys

import matplotlib

matplotlib.use("Agg")

import stchlab.core as core_module
from stchlab import STCHLab, CeO2
from stchlab.core import (
    DEFAULT_CACHE_DIR,
    DEFAULT_CC_DB,
    DEFAULT_FB_DB,
    PACKAGE_CC_DB as CORE_PACKAGE_CC_DB,
    PACKAGE_FB_DB as CORE_PACKAGE_FB_DB,
)


material = CeO2
stch = STCHLab(material)

EQ_TOL = 1e-6
REACTOR_TOL = 1e-3
CACHE_TOL = 1e-10
DB_REPRO_TOL = 1e-5

TMP_CC_DB = Path(__file__).resolve().parent / "__tmp_cc_test.db"
TMP_FB_DB = Path(__file__).resolve().parent / "__tmp_fb_test.db"
TMP_CACHE_DIR = Path(__file__).resolve().parent / "__tmp_cache_dir"
PACKAGE_CC_DB = Path(CORE_PACKAGE_CC_DB)
PACKAGE_FB_DB = Path(CORE_PACKAGE_FB_DB)

PASSED = 0
FAILED = []


def require(condition, message):
    """Raise an assertion with a readable message when a check fails."""
    if not condition:
        raise AssertionError(message)


def approx_equal(value, expected, tol, label):
    """Check two floats within tolerance."""
    require(abs(value - expected) <= tol, f"{label}: |{value} - {expected}| > {tol}")


def compare_dict_values(first, second, keys, tol):
    """Compare selected numeric keys from two dictionaries."""
    for key in keys:
        approx_equal(float(first[key]), float(second[key]), tol, key)


def count_rows(db_path):
    """Count rows in the SQLite simulations table."""
    connection = sqlite3.connect(db_path)
    try:
        cursor = connection.cursor()
        cursor.execute("SELECT COUNT(*) FROM simulations")
        return int(cursor.fetchone()[0])
    finally:
        connection.close()


def table_columns(db_path):
    """Return SQLite simulations table column names."""
    connection = sqlite3.connect(db_path)
    try:
        cursor = connection.cursor()
        cursor.execute("PRAGMA table_info(simulations)")
        return [row[1] for row in cursor.fetchall()]
    finally:
        connection.close()


def fetch_one_row(db_path, query, params=()):
    """Fetch a single row from a SQLite database."""
    connection = sqlite3.connect(db_path)
    try:
        cursor = connection.cursor()
        cursor.execute(query, params)
        row = cursor.fetchone()
        require(row is not None, f"no row returned from {db_path.name}")
        return row
    finally:
        connection.close()


def cleanup_temp_db(db_path):
    """Remove a temporary SQLite file if it exists."""
    try:
        db_path.unlink()
    except FileNotFoundError:
        pass


def cleanup_temp_dir(dir_path):
    """Remove a temporary directory tree if it exists."""
    shutil.rmtree(dir_path, ignore_errors=True)


def run_test(index, title, callback):
    """Run one numbered test and track pass/fail status."""
    global PASSED

    print(f"Test {index}: {title}")
    try:
        detail = callback()
        if detail:
            print(detail)
        print("PASSED")
        print()
        PASSED += 1
    except Exception as exc:
        print(f"FAILED: {exc}")
        print()
        FAILED.append((index, title, str(exc)))


def test_1():
    delta = 0.01
    pO2 = 100.0

    temperature = stch.funR.T(delta, pO2)
    gibbs = stch.funR.G_R(delta, temperature, pO2)
    delta_check = stch.funR.delta(temperature, pO2)
    pO2_check = stch.funR.pO2(delta, temperature)

    approx_equal(gibbs, 0.0, EQ_TOL, "G_R")
    approx_equal(delta_check, delta, EQ_TOL, "delta")
    approx_equal(pO2_check, pO2, EQ_TOL, "pO2")

    return f"  delta={delta}, pO2={pO2} Pa -> T={temperature:.6f} K"


def test_2():
    delta = 0.01
    x = 0.01

    temperature = stch.funOX.T(delta, x)
    gibbs = stch.funOX.G_OX(delta, temperature, x)
    delta_check = stch.funOX.delta(temperature, x)
    x_check = stch.funOX.x(delta, temperature)

    approx_equal(gibbs, 0.0, EQ_TOL, "G_OX")
    approx_equal(delta_check, delta, EQ_TOL, "delta")
    approx_equal(x_check, x, EQ_TOL, "x")

    return f"  delta={delta}, x={x} -> T={temperature:.6f} K"


def test_3():
    delta = 0.01
    pO2 = 100.0
    voltage = 0.5

    temperature = stch.funR.T(delta, pO2, voltage)
    gibbs = stch.funR.G_R(delta, temperature, pO2, voltage)
    delta_check = stch.funR.delta(temperature, pO2, voltage)
    pO2_check = stch.funR.pO2(delta, temperature, voltage)

    approx_equal(gibbs, 0.0, EQ_TOL, "G_R")
    approx_equal(delta_check, delta, EQ_TOL, "delta")
    approx_equal(pO2_check, pO2, EQ_TOL, "pO2")

    return f"  delta={delta}, pO2={pO2} Pa, V={voltage} V -> T={temperature:.6f} K"


def test_4():
    delta = 0.01
    x = 0.01
    voltage = -0.05

    temperature = stch.funOX.T(delta, x, voltage)
    gibbs = stch.funOX.G_OX(delta, temperature, x, voltage)
    delta_check = stch.funOX.delta(temperature, x, voltage)
    x_check = stch.funOX.x(delta, temperature, voltage)

    approx_equal(gibbs, 0.0, EQ_TOL, "G_OX")
    approx_equal(delta_check, delta, EQ_TOL, "delta")
    approx_equal(x_check, x, EQ_TOL, "x")

    return f"  delta={delta}, x={x}, V={voltage} V -> T={temperature:.6f} K"


def test_5():
    delta_ox = 0.01
    delta_r = 0.04
    pO2_min = 10.0
    pO2_max = 1000.0

    T_R_max = stch.funCCR.T_R_max(delta_ox, delta_r, pO2_min, pO2_max)
    lmbda = stch.funCCR.lmbda(pO2_min, pO2_max)
    delta_r_check = stch.funCCR.delta_R(delta_ox, pO2_min, lmbda, T_R_max, delta_r)

    approx_equal(delta_r_check, delta_r, REACTOR_TOL, "delta_R")

    return (
        f"  delta_OX={delta_ox}, delta_R={delta_r}, pO2_min={pO2_min} Pa, "
        f"pO2_max={pO2_max} Pa -> T_R_max={T_R_max:.6f} K"
    )


def test_6():
    delta_ox = 0.01
    delta_r = 0.04
    x_min = 0.001
    x_max = 0.1

    T_OX_min = stch.funCCOX.T_OX_min(delta_ox, delta_r, x_min, x_max)
    theta = stch.funCCOX.theta(x_min, x_max)
    delta_ox_check = stch.funCCOX.delta_OX(delta_r, x_min, theta, T_OX_min, delta_ox)

    approx_equal(delta_ox_check, delta_ox, REACTOR_TOL, "delta_OX")

    return f"  delta_OX={delta_ox}, delta_R={delta_r}, x_min={x_min}, x_max={x_max} -> T_OX_min={T_OX_min:.6f} K"


def test_7():
    data_time = {
        "pO2_min": 0.1,
        "T_R": 850 + 273.15,
        "V_R": 1.0,
        "ndot_inert": 0.01,
        "n_MOx": 1.0,
        "k_lim_R": 5e-4,
        "time_R": 20 * 60,
        "delta_OX": 0.001,
    }

    result_time = stch.funFBR.step(data_time)
    data_frac = {
        "pO2_min": 0.1,
        "T_R": 850 + 273.15,
        "V_R": 1.0,
        "ndot_inert": 0.01,
        "n_MOx": 1.0,
        "k_lim_R": 5e-4,
        "frac_Ddelta_R": result_time["frac_Ddelta_R"],
        "delta_OX": 0.001,
    }
    result_frac = stch.funFBR.step(data_frac)

    approx_equal(result_frac["delta_R"], result_time["delta_R"], REACTOR_TOL, "delta_R")
    for key in ["delta_R", "t_R", "delta_R_max", "frac_Ddelta_R", "Lmbda", "t", "delta", "Ddelta"]:
        require(key in result_time, f"funFBR.step result must include {key}")
    for key in ["deltaR", "tR", "deltaRmax", "tporR", "t_por_R", "lmbda"]:
        require(key not in result_time, f"funFBR.step result must not include old key {key}")

    return f"  T_R={data_time['T_R']:.3f} K, V_R={data_time['V_R']:.3f} V -> delta_R={result_time['delta_R']:.6f}"


def test_8():
    data_time = {
        "x_min": 0.001,
        "T_OX": 850 + 273.15,
        "ndot_inert": 0.01,
        "n_MOx": 1.0,
        "ndot_tot": 0.001,
        "time_OX": 5 * 60,
        "delta_R": 0.1,
    }

    result_time = stch.funFBOX.step(data_time)
    data_frac = {
        "x_min": 0.001,
        "T_OX": 850 + 273.15,
        "ndot_inert": 0.01,
        "n_MOx": 1.0,
        "ndot_tot": 0.001,
        "frac_Ddelta_OX": result_time["frac_Ddelta_OX"],
        "delta_R": 0.1,
    }
    result_frac = stch.funFBOX.step(data_frac)

    approx_equal(result_frac["delta_OX"], result_time["delta_OX"], REACTOR_TOL, "delta_OX")
    for key in ["delta_OX", "t_OX", "delta_OX_min", "frac_Ddelta_OX", "Theta", "t", "delta", "Ddelta"]:
        require(key in result_time, f"funFBOX.step result must include {key}")
    for key in ["deltaOX", "tOX", "deltaOXmin", "tporOX", "t_por_OX", "theta"]:
        require(key not in result_time, f"funFBOX.step result must not include old key {key}")

    return f"  T_OX={data_time['T_OX']:.3f} K, delta_R={data_time['delta_R']:.6f} -> delta_OX={result_time['delta_OX']:.6f}"


def test_9():
    data = {
        "delta_OX": 0.01,
        "delta_R": 0.04,
        "pO2_min": 10.0,
        "pO2_max": 1000.0,
        "x_min": 0.001,
        "x_max": 0.1,
    }

    system_input = data.copy()
    system_input_before = system_input.copy()
    state = stch.funCCSYS.system(system_input, use_db=False)
    result = stch.funCCSYS.STH(state, type="Classic")
    vt_result = stch.funCCSYS.STH(state, type="VT")
    ar_stch = STCHLab(CeO2, inert="Ar")
    ar_state = ar_stch.funCCSYS.system(data.copy(), use_db=False)
    ar_result = ar_stch.funCCSYS.STH(ar_state, type="Classic")
    ox_inert_stch = STCHLab(CeO2, x_inert=0.2)
    ox_inert_state = ox_inert_stch.funCCSYS.system(data.copy(), use_db=False)
    ox_inert_result = ox_inert_stch.funCCSYS.STH(ox_inert_state, type="Classic")
    par = {"C": 2500, "nu_sep": 0.25}
    par_before = par.copy()
    stch.funCCSYS.STH(state, par=par, type="Classic")
    new_result_keys = [
        "type", "W_H2", "eta", "Q_total", "Q_sol", "Q_rec", "Q_aux", "Q_rej",
        "Q_R", "Q_R_ideal", "Q_OX", "Q_OX_used", "Q_ch_R", "Q_ch_OX",
        "Q_MOx", "Q_MOx_R", "Q_MOx_OX", "Q_inert", "Q_sep_H2_ideal",
        "Q_sep_H2", "Q_sep_inert_ideal", "Q_sep_inert", "Q_H2O_re",
        "Q_H2O_new", "Q_H2O_lat", "T_H_in", "T_H_out", "T_C_in", "T_C_out",
    ]
    old_result_keys = [
        "Type", "WH2", "eff", "Qtotal", "Qsol", "Qrec", "Qaux", "Qrej",
        "QR", "QRIdeal", "QOX", "QOXused", "QchR", "QchOX", "QMO",
        "QMO_R", "QMO_OX", "Qinert", "QsH2ideal", "QsH2",
        "Qsinertideal", "Qsinert", "QH2Ore", "QH2Onew", "QH2Olat",
        "THin", "THout", "TCin", "TCout",
    ]

    require(state["delta_R"] > state["delta_OX"], "counter-current system must satisfy delta_R > delta_OX")
    require(state["T_R_max"] > 0.0, "T_R_max must be positive")
    require(state["T_OX_min"] > 0.0, "T_OX_min must be positive")
    require(state["lmbda"] > 0.0, "lmbda must be positive")
    require(state["theta"] > 0.0, "theta must be positive")
    require(result["Q_total"] > 0.0, "Q_total must be positive")
    require(0.0 < result["STH"] < 1.0, "STH must be between 0 and 1")
    require(vt_result["Q_total"] > 0.0, "VT Q_total must be positive")
    require(ar_stch.inert == "Ar", "inert gas should be stored on STCHLab")
    require(ar_result["Q_total"] > 0.0, "Ar counter-current Q_total must be positive")
    expected_ndot_H2O = state["Ddelta"] * (1.0 - stch.x_inert - state["x_min"]) / (state["x_max"] - state["x_min"])
    approx_equal(state["ndot_H2O"], expected_ndot_H2O, CACHE_TOL, "counter-current OX steam flow")
    expected_inert_ndot_H2O = ox_inert_state["Ddelta"] * (1.0 - ox_inert_stch.x_inert - ox_inert_state["x_min"]) / (ox_inert_state["x_max"] - ox_inert_state["x_min"])
    approx_equal(ox_inert_state["ndot_H2O"], expected_inert_ndot_H2O, CACHE_TOL, "counter-current OX steam flow with inert")
    require(ox_inert_result["Q_H2O_re"] > 0.0, "counter-current OX recycle duty with inert must be positive")
    require(ox_inert_result["Q_total"] > 0.0, "counter-current OX inert Q_total must be positive")
    require("i_DT" in vt_result, "VT result must include i_DT")
    require("iDT" not in vt_result, "VT result must not include old key iDT")
    for key in new_result_keys:
        require(key in result, f"funCCSYS.STH result must include {key}")
    for key in old_result_keys:
        require(key not in result, f"funCCSYS.STH result must not include old key {key}")
    require(system_input == system_input_before, "funCCSYS.system must not mutate input data")
    require(par == par_before, "funCCSYS.STH must not mutate par overrides")

    return f"  T_R_max={state['T_R_max']:.6f} K, T_OX_min={state['T_OX_min']:.6f} K, STH={100 * result['STH']:.3f} %"


def test_10():
    data = {
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

    state = stch.funFBSYS.system(data.copy(), use_db=False)
    result = stch.funFBSYS.STH(state)
    ar_stch = STCHLab(CeO2, inert="Ar")
    ar_state = ar_stch.funFBSYS.system(data.copy(), use_db=False)
    ar_result = ar_stch.funFBSYS.STH(ar_state)
    par = {"C": 2500, "dT_dt": 100}
    par_before = par.copy()
    stch.funFBSYS.STH(state, par=par)
    new_result_keys = [
        "STH", "STH_available", "Q_solar_R", "Q_proj_R", "Q_rec_R",
        "Q_inert_R", "Q_inert", "Q_ch_R", "W_el_R", "Q_sep_inert",
        "Q_MOx", "W_H2", "A_ap", "R_ap", "C_R", "C_W_max",
        "Q_solar_W", "Q_ch_OX", "Q_H2O_re", "Q_H2O_new", "Q_H2O_lat",
        "Q_loss_OX", "Q_sep_H2", "Q_proj_C", "Q_solar_C", "Q_rej_C",
        "Q_rej_T", "Q_aux", "Q_solar", "Q_loss", "Q_proj", "Q_coll",
        "Q_miss_C", "Q_total", "Q_total_1", "Q_total_2", "Q_total_3",
        "Q_loss_R", "Q_loss_W", "Q_loss_C", "Q_inert_W", "Q_inert_C",
        "Q_proj_W", "Q_proc", "Q_avail", "Q_solar_available",
        "Q_loss_2nd_law", "Q_loss_circ_1", "Q_loss_circ_2",
    ]
    old_result_keys = [
        "STHT", "QsolarR", "QprojR", "QrecR", "QinertR", "Qinert",
        "QchR", "WelR", "Qsinert", "QMO", "WH2", "CR", "CWmax",
        "QsolarW", "QchOX", "QH2Ore", "QH2Onew", "QH2Olat",
        "QlossOX", "QsH2", "QprojC", "QsolarC", "QrejC", "QrejT",
        "Qaux", "Qsolar", "Qloss", "Qproj", "Qcoll", "QmissC",
        "Qtotal", "Qtotal1", "Qtotal2", "Qtotal3", "QlossR",
        "QlossW", "QlossC", "QinertW", "QinertC", "QprojW",
        "Qproc", "Qavail", "QsolarT", "Qloss2law", "QlossCirc1",
        "QlossCirc2",
    ]
    plot_data = state.copy()
    plot_data_before = plot_data.copy()
    original_show = core_module.plt.show
    try:
        core_module.plt.show = lambda *args, **kwargs: None
        stch.funFBSYS.plot_system(plot_data)
    finally:
        core_module.plt.show = original_show
        core_module.plt.close("all")

    require(state["delta_R"] > state["delta_OX"], "fixed-bed system must satisfy delta_R > delta_OX")
    require(state["Theta"] > 0.0, "Theta must be positive")
    require(state["Lmbda"] > 0.0, "Lmbda must be positive")
    require(result["Q_total"] > 0.0, "Q_total must be positive")
    require(result["STH_available"] > 0.0, "STH_available must be positive")
    require(0.0 < result["STH"] < 1.0, "STH must be between 0 and 1")
    require(ar_result["Q_total"] > 0.0, "Ar fixed-bed Q_total must be positive")
    for key in new_result_keys:
        require(key in result, f"funFBSYS.STH result must include {key}")
    for key in old_result_keys:
        require(key not in result, f"funFBSYS.STH result must not include old key {key}")
    for key in stch.funFBSYS.system_keys:
        require(key in state, f"funFBSYS.system result must include {key}")
    for key in [
        "TR", "TOX", "VR", "VOX", "pO2min", "xmin", "nMOx", "ninert", "ntot",
        "klimR", "klimOX", "timeR", "timeOX", "tporR", "tporOX", "t_por_R", "t_por_OX", "deltaR", "deltaOX",
    ]:
        require(key not in state, f"funFBSYS.system result must not include old key {key}")
    require(par == par_before, "funFBSYS.STH must not mutate par overrides")
    require(plot_data == plot_data_before, "funFBSYS.plot_system must not mutate input data")

    return f"  delta_R={state['delta_R']:.6f}, delta_OX={state['delta_OX']:.6f}, STH={100 * result['STH']:.3f} %"


def test_11():
    cleanup_temp_db(TMP_CC_DB)
    try:
        data = {
            "delta_OX": 0.01,
            "delta_R": 0.04,
            "pO2_min": 10.0,
            "pO2_max": 1000.0,
            "x_min": 0.001,
            "x_max": 0.1,
        }

        first = stch.funCCSYS.system(data.copy(), use_db=True, db_name=str(TMP_CC_DB))
        second = stch.funCCSYS.system(data.copy(), use_db=True, db_name=str(TMP_CC_DB))
        compare_dict_values(first, second, ["delta_OX", "delta_R", "T_R_max", "T_OX_min", "lmbda", "theta"], CACHE_TOL)

        rows = count_rows(TMP_CC_DB)
        require(rows == 1, f"expected 1 cached counter-current row, found {rows}")
        require("inert" not in table_columns(TMP_CC_DB), "inert must not be stored in counter-current system cache")

        return f"  cache file={TMP_CC_DB.name}, rows={rows}"
    finally:
        cleanup_temp_db(TMP_CC_DB)


def test_12():
    cleanup_temp_db(TMP_FB_DB)
    try:
        data = {
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

        first = stch.funFBSYS.system(data.copy(), use_db=True, db_name=str(TMP_FB_DB))
        second = stch.funFBSYS.system(data.copy(), use_db=True, db_name=str(TMP_FB_DB))
        compare_dict_values(first, second, ["delta_R", "delta_OX", "Theta", "Lmbda"], CACHE_TOL)

        rows = count_rows(TMP_FB_DB)
        require(rows == 1, f"expected 1 cached fixed-bed row, found {rows}")
        require("inert" not in table_columns(TMP_FB_DB), "inert must not be stored in fixed-bed system cache")

        return f"  cache file={TMP_FB_DB.name}, rows={rows}"
    finally:
        cleanup_temp_db(TMP_FB_DB)


def test_13():
    row = fetch_one_row(
        PACKAGE_CC_DB,
        """
        SELECT
            pO2_min, pO2_max, x_min, x_max, delta_R, delta_OX,
            T_R_in, T_R_out, T_R_max, T_OX_in, T_OX_out, T_OX_min,
            Ddelta, ndot_inert, ndot_H2O, xO2_min, xO2_max
        FROM simulations
        WHERE MOx = ? AND oxidant = ? AND x_inert = ?
        ORDER BY rowid
        LIMIT 1
        """,
        ("CeO2", "H2O", 0.0),
    )

    (
        pO2_min,
        pO2_max,
        x_min,
        x_max,
        delta_r,
        delta_ox,
        tr_in,
        tr_out,
        tr_max,
        tox_in,
        tox_out,
        tox_min,
        Ddelta,
        ndot_inert,
        ndot_H2O,
        xO2_min,
        xO2_max,
    ) = row

    data = {
        "delta_OX": delta_ox,
        "delta_R": delta_r,
        "pO2_min": pO2_min,
        "pO2_max": pO2_max,
        "x_min": x_min,
        "x_max": x_max,
    }

    state = stch.funCCSYS.system(data.copy(), use_db=False)
    compare_dict_values(
        state,
        {
            "T_R_in": tr_in,
            "T_R_out": tr_out,
            "T_R_max": tr_max,
            "T_OX_in": tox_in,
            "T_OX_out": tox_out,
            "T_OX_min": tox_min,
            "Ddelta": Ddelta,
            "ndot_inert": ndot_inert,
            "ndot_H2O": Ddelta * (1.0 - stch.x_inert - x_min) / (x_max - x_min),
            "xO2_min": xO2_min,
            "xO2_max": xO2_max,
        },
        ["T_R_in", "T_R_out", "T_R_max", "T_OX_in", "T_OX_out", "T_OX_min", "Ddelta", "ndot_inert", "ndot_H2O", "xO2_min", "xO2_max"],
        CACHE_TOL,
    )

    return f"  packaged CC row reproduced: T_R_max={state['T_R_max']:.6f} K, T_OX_min={state['T_OX_min']:.6f} K"


def test_14():
    row = fetch_one_row(
        PACKAGE_FB_DB,
        """
        SELECT
            T_R, T_OX, V_R, V_OX, pO2_min, x_min, n_MOx, ndot_inert, ndot_tot,
            k_lim_R, k_lim_OX, time_R, time_OX, frac_Ddelta_R, frac_Ddelta_OX,
            delta_R, delta_OX, Theta, Lmbda
        FROM simulations
        WHERE MOx = ? AND oxidant = ? AND x_inert = ?
        ORDER BY rowid
        LIMIT 1
        """,
        ("CeO2", "H2O", 0.0),
    )

    (
        T_R_value,
        tox,
        vr,
        vox,
        pO2_min,
        x_min,
        n_MOx,
        ndot_inert,
        ndot_tot,
        klim_r,
        klim_ox,
        time_r,
        time_ox,
        frac_Ddelta_r,
        frac_Ddelta_ox,
        delta_r,
        delta_ox,
        theta,
        Lmbda,
    ) = row

    data = {
        "T_R": T_R_value,
        "T_OX": tox,
        "V_R": vr,
        "V_OX": vox,
        "pO2_min": pO2_min,
        "x_min": x_min,
        "n_MOx": n_MOx,
        "ndot_inert": ndot_inert,
        "ndot_tot": ndot_tot,
        "k_lim_R": klim_r,
        "k_lim_OX": klim_ox,
        "time_R": time_r,
        "time_OX": time_ox,
        "frac_Ddelta_R": frac_Ddelta_r,
        "frac_Ddelta_OX": frac_Ddelta_ox,
    }

    state_recomputed = stch.funFBSYS.system(data.copy(), use_db=False)

    compare_dict_values(
        state_recomputed,
        {"delta_R": delta_r, "delta_OX": delta_ox, "Theta": theta, "Lmbda": Lmbda},
        ["delta_R", "delta_OX", "Theta", "Lmbda"],
        DB_REPRO_TOL,
    )

    return f"  packaged FB row reproduced: delta_R={state_recomputed['delta_R']:.6f}, delta_OX={state_recomputed['delta_OX']:.6f}"


def test_15():
    require(STCHLab(CeO2).inert == "N2", "default inert gas should be N2")
    require(STCHLab(CeO2, inert="Argon").inert == "Argon", "CoolProp inert gas names should be accepted")
    co2_stch = STCHLab(CeO2, oxidant="co2")
    require(co2_stch.oxidant == "CO2", "oxidant should normalize lowercase CO2")

    try:
        STCHLab(CeO2, oxidant="H20")
    except ValueError as exc:
        require("oxidant" in str(exc), "invalid oxidant error should be explicit")
    else:
        raise AssertionError("invalid oxidant should raise ValueError")

    try:
        STCHLab(CeO2, inert="not-a-fluid")
    except ValueError as exc:
        require("inert" in str(exc), "invalid inert error should be explicit")
    else:
        raise AssertionError("invalid inert should raise ValueError")

    return "  lowercase CO2 normalized; inert gas validated"


def test_16():
    impossible_reduction = {
        "pO2_min": 0.12923484877986363,
        "T_R": 1188.98258902546,
        "V_R": 0.5135601695713496,
        "ndot_inert": 0.011870647010743744,
        "n_MOx": 1.0277160178533251,
        "k_lim_R": 0.00010221064050146226,
        "time_R": 998.7020424640032,
        "delta_OX": 0.009961790369998526,
    }

    red_time = stch.funFBR.step(impossible_reduction)
    red_frac_data = impossible_reduction.copy()
    red_frac_data.pop("time_R")
    red_frac_data["frac_Ddelta_R"] = 0.5
    red_frac = stch.funFBR.step(red_frac_data)

    require(red_time.get("status") == "infeasible", "impossible reduction should be flagged")
    require(red_time["delta_R_max"] <= impossible_reduction["delta_OX"], "test case must be reduction-infeasible")
    approx_equal(red_time["delta_R"], impossible_reduction["delta_OX"], CACHE_TOL, "infeasible reduction delta_R")
    approx_equal(red_frac["delta_R"], impossible_reduction["delta_OX"], CACHE_TOL, "infeasible frac_Ddelta_R delta_R")
    approx_equal(red_time["frac_Ddelta_R"], 0.0, CACHE_TOL, "infeasible reduction frac_Ddelta_R")
    approx_equal(red_frac["frac_Ddelta_R"], 0.0, CACHE_TOL, "infeasible frac_Ddelta_R")

    impossible_oxidation = {
        "x_min": 0.1,
        "T_OX": 1400.0,
        "ndot_tot": 1e-3,
        "n_MOx": 1.0,
        "k_lim_OX": 1e-4,
        "time_OX": 300.0,
        "delta_R": 0.001,
    }

    ox_time = stch.funFBOX.step(impossible_oxidation)
    ox_frac_data = impossible_oxidation.copy()
    ox_frac_data.pop("time_OX")
    ox_frac_data["frac_Ddelta_OX"] = 0.5
    ox_frac = stch.funFBOX.step(ox_frac_data)

    require(ox_time.get("status") == "infeasible", "impossible oxidation should be flagged")
    require(ox_time["delta_OX_min"] >= impossible_oxidation["delta_R"], "test case must be oxidation-infeasible")
    approx_equal(ox_time["delta_OX"], impossible_oxidation["delta_R"], CACHE_TOL, "infeasible oxidation delta_OX")
    approx_equal(ox_frac["delta_OX"], impossible_oxidation["delta_R"], CACHE_TOL, "infeasible frac_Ddelta_OX delta_OX")
    approx_equal(ox_time["frac_Ddelta_OX"], 0.0, CACHE_TOL, "infeasible oxidation frac_Ddelta_OX")
    approx_equal(ox_frac["frac_Ddelta_OX"], 0.0, CACHE_TOL, "infeasible frac_Ddelta_OX")

    return "  infeasible fixed-bed reduction and oxidation return no-op states"


def test_17():
    default_cache_dir = Path(DEFAULT_CACHE_DIR)
    require(Path(DEFAULT_CC_DB).parent == default_cache_dir, "DEFAULT_CC_DB should live in DEFAULT_CACHE_DIR")
    require(Path(DEFAULT_FB_DB).parent == default_cache_dir, "DEFAULT_FB_DB should live in DEFAULT_CACHE_DIR")
    require(Path(DEFAULT_CC_DB) != PACKAGE_CC_DB, "default CC cache should not be the packaged DB")
    require(Path(DEFAULT_FB_DB) != PACKAGE_FB_DB, "default FB cache should not be the packaged DB")

    cleanup_temp_dir(TMP_CACHE_DIR)
    nested_db = TMP_CACHE_DIR / "nested" / "counter_current_cache.db"
    try:
        data = {
            "delta_OX": 0.01,
            "delta_R": 0.04,
            "pO2_min": 10.0,
            "pO2_max": 1000.0,
            "x_min": 0.001,
            "x_max": 0.1,
        }
        stch.funCCSYS.system(data.copy(), use_db=True, db_name=str(nested_db))
        require(nested_db.exists(), "custom nested cache DB should be created")
        require(nested_db.parent.exists(), "custom nested cache directory should be created")
        rows = count_rows(nested_db)
        require(rows == 1, f"expected 1 cached nested row, found {rows}")
        require("inert" not in table_columns(nested_db), "inert must not be stored in custom system cache")
        return f"  cache defaults separated from packaged DBs; nested rows={rows}"
    finally:
        cleanup_temp_dir(TMP_CACHE_DIR)


run_test(1, "funR equilibrium reciprocity", test_1)
run_test(2, "funOX equilibrium reciprocity", test_2)
run_test(3, "funR thermo-electrochemical reciprocity", test_3)
run_test(4, "funOX thermo-electrochemical reciprocity", test_4)
run_test(5, "funCCR T_R_max / delta_R reciprocity", test_5)
run_test(6, "funCCOX T_OX_min / delta_OX reciprocity", test_6)
run_test(7, "funFBR time_R / frac_Ddelta_R consistency", test_7)
run_test(8, "funFBOX time_OX / frac_Ddelta_OX consistency", test_8)
run_test(9, "funCCSYS system-level coherence", test_9)
run_test(10, "funFBSYS system-level coherence", test_10)
run_test(11, "counter-current SQLite cache smoke test", test_11)
run_test(12, "fixed-bed SQLite cache smoke test", test_12)
run_test(13, "packaged counter-current DB example reproducibility", test_13)
run_test(14, "packaged fixed-bed DB example reproducibility", test_14)
run_test(15, "oxidant validation", test_15)
run_test(16, "fixed-bed infeasible direction handling", test_16)
run_test(17, "SQLite cache path handling", test_17)


print(f"Summary: {PASSED} passed, {len(FAILED)} failed")
if FAILED:
    print("Failed tests:")
    for index, title, message in FAILED:
        print(f"  {index}. {title}: {message}")
    sys.exit(1)

"""Full regression test with randomized inputs inside safe operating ranges.

Run from the repository root with:

    python tests/test_full.py
    python tests/test_full.py --cases 50 --seed 42
"""

import argparse
from pathlib import Path
import sqlite3
import sys

import numpy as np

from stchlab import STCHLab, CeO2


material = CeO2
stch = STCHLab(material)

EQ_TOL = 1e-6
REACTOR_TOL = 1e-3
SYSTEM_TOL = 1e-8

TMP_CC_DB = Path(__file__).resolve().parent / "__tmp_cc_test.db"
TMP_FB_DB = Path(__file__).resolve().parent / "__tmp_fb_test.db"

PASSED = 0
FAILED = []


def positive_int(value):
    """Parse a positive integer command-line argument."""
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return ivalue


def parse_args():
    """Parse command-line arguments for the randomized full test."""
    parser = argparse.ArgumentParser(description="Run the randomized STCHLab full test.")
    parser.add_argument(
        "--cases",
        type=positive_int,
        default=2,
        help="number of random cases per randomized test block (default: 2)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="random seed used to generate the cases (default: 42)",
    )
    return parser.parse_args()


ARGS = parse_args()
SEED = ARGS.seed
N_EQ_CASES = ARGS.cases
N_REACTOR_CASES = ARGS.cases
N_SYSTEM_CASES = ARGS.cases


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


def cleanup_temp_db(db_path):
    """Remove a temporary SQLite file if it exists."""
    try:
        db_path.unlink()
    except FileNotFoundError:
        pass


def run_test(index, title, callback):
    """Run one numbered test and track pass/fail status."""
    global PASSED

    print(f"Test {index}: {title}")
    try:
        detail = callback()
        if detail:
            print(detail, flush=True)
        print("PASSED", flush=True)
        print(flush=True)
        PASSED += 1
    except Exception as exc:
        print(f"FAILED: {exc}", flush=True)
        print(flush=True)
        FAILED.append((index, title, str(exc)))


def make_rng(offset):
    """Create a reproducible random generator for one test block."""
    return np.random.default_rng(SEED + offset)


def rand_uniform(rng, low, high):
    """Sample a uniform float."""
    return float(rng.uniform(low, high))


def rand_loguniform(rng, low, high):
    """Sample a log-uniform float."""
    return float(np.exp(rng.uniform(np.log(low), np.log(high))))


def random_delta_range(rng):
    """Generate a safe reduction extent range for CeO2."""
    delta_ox = rand_uniform(rng, 0.005, 0.02)
    delta_r = delta_ox + rand_uniform(rng, 0.01, 0.03)
    delta_r = min(delta_r, 0.05)
    require(delta_r > delta_ox, "random delta range must satisfy delta_R > delta_OX")
    return delta_ox, delta_r


def random_cc_data(rng):
    """Generate one random counter-current case."""
    delta_ox, delta_r = random_delta_range(rng)
    pO2_min = rand_loguniform(rng, 5.0, 50.0)
    pO2_max = pO2_min * rand_uniform(rng, 20.0, 150.0)
    pO2_max = min(pO2_max, 5.0e3)
    require(pO2_max > pO2_min, "pO2_max must be larger than pO2_min")

    x_min = rand_loguniform(rng, 1e-4, 5e-3)
    x_max = x_min + rand_uniform(rng, 0.02, 0.08)
    x_max = min(x_max, 0.12)
    require(x_max > x_min, "x_max must be larger than x_min")

    return {
        "delta_OX": delta_ox,
        "delta_R": delta_r,
        "pO2_min": pO2_min,
        "pO2_max": pO2_max,
        "x_min": x_min,
        "x_max": x_max,
    }


def random_fbr_data(rng):
    """Generate one random fixed-bed reduction-step case."""
    return {
        "pO2_min": rand_loguniform(rng, 0.08, 0.3),
        "T_R": rand_uniform(rng, 1100.0, 1250.0),
        "V_R": rand_uniform(rng, 0.2, 1.0),
        "ndot_inert": rand_loguniform(rng, 5e-3, 2e-2),
        "n_MOx": rand_loguniform(rng, 0.5, 1.5),
        "k_lim_R": rand_loguniform(rng, 1e-4, 5e-4),
        "time_R": rand_uniform(rng, 300.0, 1200.0),
        "delta_OX": rand_uniform(rng, 0.001, 0.01),
    }


def random_fbox_data(rng):
    """Generate one random fixed-bed oxidation-step case."""
    return {
        "x_min": rand_loguniform(rng, 1e-4, 2e-3),
        "T_OX": rand_uniform(rng, 1100.0, 1250.0),
        "ndot_inert": rand_loguniform(rng, 5e-3, 2e-2),
        "n_MOx": rand_loguniform(rng, 0.5, 1.5),
        "ndot_tot": rand_loguniform(rng, 5e-4, 2e-3),
        "time_OX": rand_uniform(rng, 120.0, 600.0),
        "delta_R": rand_uniform(rng, 0.04, 0.1),
    }


def random_fb_system_data(rng):
    """Generate one random fixed-bed system case."""
    tox = rand_uniform(rng, 1120.0, 1230.0)
    T_R_value = tox + rand_uniform(rng, 20.0, 180.0)
    return {
        "pO2_min": rand_loguniform(rng, 0.08, 0.3),
        "T_R": T_R_value,
        "V_R": rand_uniform(rng, 0.2, 0.9),
        "ndot_inert": rand_loguniform(rng, 2e-5, 1e-4),
        "n_MOx": rand_loguniform(rng, 3e-3, 8e-3),
        "x_min": rand_loguniform(rng, 1e-6, 1e-4),
        "T_OX": tox,
        "ndot_tot": rand_loguniform(rng, 2e-5, 1e-4),
        "k_lim_R": rand_loguniform(rng, 1e-5, 1e-4),
        "k_lim_OX": rand_loguniform(rng, 5e-5, 2e-4),
        "time_R": rand_uniform(rng, 300.0, 900.0),
        "time_OX": rand_uniform(rng, 150.0, 400.0),
    }


def test_1():
    rng = make_rng(1)
    for case in range(1, N_EQ_CASES + 1):
        delta = rand_uniform(rng, 0.005, 0.04)
        pO2 = rand_loguniform(rng, 10.0, 1000.0)
        temperature = stch.funR.T(delta, pO2)
        gibbs = stch.funR.G_R(delta, temperature, pO2)
        delta_check = stch.funR.delta(temperature, pO2)
        pO2_check = stch.funR.pO2(delta, temperature)
        approx_equal(gibbs, 0.0, EQ_TOL, f"G_R case {case}")
        approx_equal(delta_check, delta, EQ_TOL, f"delta case {case}")
        approx_equal(pO2_check, pO2, EQ_TOL, f"pO2 case {case}")
    return f"  {N_EQ_CASES} random reduction equilibrium cases (seed={SEED + 1})"


def test_2():
    rng = make_rng(2)
    for case in range(1, N_EQ_CASES + 1):
        delta = rand_uniform(rng, 0.005, 0.04)
        x = rand_loguniform(rng, 1e-3, 0.01)
        temperature = stch.funOX.T(delta, x)
        gibbs = stch.funOX.G_OX(delta, temperature, x)
        delta_check = stch.funOX.delta(temperature, x)
        x_check = stch.funOX.x(delta, temperature)
        approx_equal(gibbs, 0.0, EQ_TOL, f"G_OX case {case}")
        approx_equal(delta_check, delta, EQ_TOL, f"delta case {case}")
        approx_equal(x_check, x, EQ_TOL, f"x case {case}")
    return f"  {N_EQ_CASES} random oxidation equilibrium cases (seed={SEED + 2})"


def test_3():
    rng = make_rng(3)
    for case in range(1, N_EQ_CASES + 1):
        delta = rand_uniform(rng, 0.005, 0.04)
        pO2 = rand_loguniform(rng, 10.0, 1000.0)
        voltage = rand_uniform(rng, 0.2, 0.6)
        temperature = stch.funR.T(delta, pO2, voltage)
        gibbs = stch.funR.G_R(delta, temperature, pO2, voltage)
        delta_check = stch.funR.delta(temperature, pO2, voltage)
        pO2_check = stch.funR.pO2(delta, temperature, voltage)
        approx_equal(gibbs, 0.0, EQ_TOL, f"G_R case {case}")
        approx_equal(delta_check, delta, EQ_TOL, f"delta case {case}")
        approx_equal(pO2_check, pO2, EQ_TOL, f"pO2 case {case}")
    return f"  {N_EQ_CASES} random thermo-electrochemical reduction cases (seed={SEED + 3})"


def test_4():
    rng = make_rng(4)
    for case in range(1, N_EQ_CASES + 1):
        delta = rand_uniform(rng, 0.005, 0.04)
        x = rand_loguniform(rng, 1e-3, 0.01)
        voltage = rand_uniform(rng, -0.08, -0.01)
        temperature = stch.funOX.T(delta, x, voltage)
        gibbs = stch.funOX.G_OX(delta, temperature, x, voltage)
        delta_check = stch.funOX.delta(temperature, x, voltage)
        x_check = stch.funOX.x(delta, temperature, voltage)
        approx_equal(gibbs, 0.0, EQ_TOL, f"G_OX case {case}")
        approx_equal(delta_check, delta, EQ_TOL, f"delta case {case}")
        approx_equal(x_check, x, EQ_TOL, f"x case {case}")
    return f"  {N_EQ_CASES} random thermo-electrochemical oxidation cases (seed={SEED + 4})"


def test_5():
    rng = make_rng(5)
    for case in range(1, N_REACTOR_CASES + 1):
        delta_ox, delta_r = random_delta_range(rng)
        pO2_min = rand_loguniform(rng, 1.0, 100.0)
        pO2_max = min(pO2_min * rand_uniform(rng, 10.0, 500.0), 2.0e4)
        T_R_max = stch.funCCR.T_R_max(delta_ox, delta_r, pO2_min, pO2_max)
        lmbda = stch.funCCR.lmbda(pO2_min, pO2_max)
        delta_r_check = stch.funCCR.delta_R(delta_ox, pO2_min, lmbda, T_R_max, delta_r)
        approx_equal(delta_r_check, delta_r, REACTOR_TOL, f"delta_R case {case}")
    return f"  {N_REACTOR_CASES} random counter-current reduction cases (seed={SEED + 5})"


def test_6():
    rng = make_rng(6)
    for case in range(1, N_REACTOR_CASES + 1):
        delta_ox, delta_r = random_delta_range(rng)
        x_min = rand_loguniform(rng, 1e-5, 1e-2)
        x_max = min(x_min + rand_uniform(rng, 0.01, 0.15), 0.25)
        T_OX_min = stch.funCCOX.T_OX_min(delta_ox, delta_r, x_min, x_max)
        theta = stch.funCCOX.theta(x_min, x_max)
        delta_ox_check = stch.funCCOX.delta_OX(delta_r, x_min, theta, T_OX_min, delta_ox)
        approx_equal(delta_ox_check, delta_ox, REACTOR_TOL, f"delta_OX case {case}")
    return f"  {N_REACTOR_CASES} random counter-current oxidation cases (seed={SEED + 6})"


def test_7():
    rng = make_rng(7)
    for case in range(1, N_REACTOR_CASES + 1):
        data_time = random_fbr_data(rng)
        result_time = stch.funFBR.step(data_time)
        data_frac = {
            "pO2_min": data_time["pO2_min"],
            "T_R": data_time["T_R"],
            "V_R": data_time["V_R"],
            "ndot_inert": data_time["ndot_inert"],
            "n_MOx": data_time["n_MOx"],
            "k_lim_R": data_time["k_lim_R"],
            "frac_Ddelta_R": result_time["frac_Ddelta_R"],
            "delta_OX": data_time["delta_OX"],
        }
        result_frac = stch.funFBR.step(data_frac)
        approx_equal(result_frac["delta_R"], result_time["delta_R"], REACTOR_TOL, f"delta_R case {case}")
    return f"  {N_REACTOR_CASES} random fixed-bed reduction-step cases (seed={SEED + 7})"


def test_8():
    rng = make_rng(8)
    for case in range(1, N_REACTOR_CASES + 1):
        data_time = random_fbox_data(rng)
        result_time = stch.funFBOX.step(data_time)
        data_frac = {
            "x_min": data_time["x_min"],
            "T_OX": data_time["T_OX"],
            "ndot_inert": data_time["ndot_inert"],
            "n_MOx": data_time["n_MOx"],
            "ndot_tot": data_time["ndot_tot"],
            "frac_Ddelta_OX": result_time["frac_Ddelta_OX"],
            "delta_R": data_time["delta_R"],
        }
        result_frac = stch.funFBOX.step(data_frac)
        approx_equal(result_frac["delta_OX"], result_time["delta_OX"], REACTOR_TOL, f"delta_OX case {case}")
    return f"  {N_REACTOR_CASES} random fixed-bed oxidation-step cases (seed={SEED + 8})"


def test_9():
    rng = make_rng(9)
    ar_stch = STCHLab(material, inert="Ar")
    for case in range(1, N_SYSTEM_CASES + 1):
        data = random_cc_data(rng)
        state = stch.funCCSYS.system(data.copy(), use_db=False)
        result = stch.funCCSYS.STH(state, type="Classic")
        ar_state = ar_stch.funCCSYS.system(data.copy(), use_db=False)
        ar_result = ar_stch.funCCSYS.STH(ar_state, type="Classic")
        require(state["delta_R"] > state["delta_OX"], f"case {case}: delta_R must exceed delta_OX")
        require(state["T_R_max"] > 0.0, f"case {case}: T_R_max must be positive")
        require(state["T_OX_min"] > 0.0, f"case {case}: T_OX_min must be positive")
        require(state["lmbda"] > 0.0, f"case {case}: lmbda must be positive")
        require(state["theta"] > 0.0, f"case {case}: theta must be positive")
        require(result["Q_total"] > 0.0, f"case {case}: Q_total must be positive")
        require(ar_result["Q_total"] > 0.0, f"case {case}: Ar Q_total must be positive")
        require(0.0 < result["STH"] < 1.0, f"case {case}: STH must be between 0 and 1")
    return f"  {N_SYSTEM_CASES} random counter-current system cases (seed={SEED + 9})"


def test_10():
    rng = make_rng(10)
    ar_stch = STCHLab(material, inert="Ar")
    for case in range(1, N_SYSTEM_CASES + 1):
        data = random_fb_system_data(rng)
        state = stch.funFBSYS.system(data.copy(), use_db=False)
        result = stch.funFBSYS.STH(state)
        ar_state = ar_stch.funFBSYS.system(data.copy(), use_db=False)
        ar_result = ar_stch.funFBSYS.STH(ar_state)
        require(state["delta_R"] > state["delta_OX"], f"case {case}: delta_R must exceed delta_OX")
        require(state["Theta"] > 0.0, f"case {case}: Theta must be positive")
        require(state["Lmbda"] > 0.0, f"case {case}: Lmbda must be positive")
        require(result["Q_total"] > 0.0, f"case {case}: Q_total must be positive")
        require(ar_result["Q_total"] > 0.0, f"case {case}: Ar Q_total must be positive")
        require(result["STH_available"] > 0.0, f"case {case}: STH_available must be positive")
        require(0.0 < result["STH"] < 1.0, f"case {case}: STH must be between 0 and 1")
    return f"  {N_SYSTEM_CASES} random fixed-bed system cases (seed={SEED + 10})"


def test_11():
    cleanup_temp_db(TMP_CC_DB)
    try:
        rng = make_rng(11)
        data = random_cc_data(rng)
        first = stch.funCCSYS.system(data.copy(), use_db=True, db_name=str(TMP_CC_DB))
        second = stch.funCCSYS.system(data.copy(), use_db=True, db_name=str(TMP_CC_DB))
        compare_dict_values(first, second, ["delta_OX", "delta_R", "T_R_max", "T_OX_min", "lmbda", "theta"], SYSTEM_TOL)
        rows = count_rows(TMP_CC_DB)
        require(rows == 1, f"expected 1 cached counter-current row, found {rows}")
        return f"  random CC cache case reused successfully (seed={SEED + 11}, rows={rows})"
    finally:
        cleanup_temp_db(TMP_CC_DB)


def test_12():
    cleanup_temp_db(TMP_FB_DB)
    try:
        rng = make_rng(12)
        data = random_fb_system_data(rng)
        first = stch.funFBSYS.system(data.copy(), use_db=True, db_name=str(TMP_FB_DB))
        second = stch.funFBSYS.system(data.copy(), use_db=True, db_name=str(TMP_FB_DB))
        compare_dict_values(first, second, ["delta_R", "delta_OX", "Theta", "Lmbda"], SYSTEM_TOL)
        rows = count_rows(TMP_FB_DB)
        require(rows == 1, f"expected 1 cached fixed-bed row, found {rows}")
        return f"  random FB cache case reused successfully (seed={SEED + 12}, rows={rows})"
    finally:
        cleanup_temp_db(TMP_FB_DB)


run_test(1, "random funR equilibrium reciprocity", test_1)
run_test(2, "random funOX equilibrium reciprocity", test_2)
run_test(3, "random funR thermo-electrochemical reciprocity", test_3)
run_test(4, "random funOX thermo-electrochemical reciprocity", test_4)
run_test(5, "random funCCR T_R_max / delta_R reciprocity", test_5)
run_test(6, "random funCCOX T_OX_min / delta_OX reciprocity", test_6)
run_test(7, "random funFBR time_R / frac_Ddelta_R consistency", test_7)
run_test(8, "random funFBOX time_OX / frac_Ddelta_OX consistency", test_8)
run_test(9, "random funCCSYS system-level coherence", test_9)
run_test(10, "random funFBSYS system-level coherence", test_10)
run_test(11, "random counter-current SQLite cache smoke test", test_11)
run_test(12, "random fixed-bed SQLite cache smoke test", test_12)


print(f"Summary: {PASSED} passed, {len(FAILED)} failed")
if FAILED:
    print("Failed tests:")
    for index, title, message in FAILED:
        print(f"  {index}. {title}: {message}")
    sys.exit(1)

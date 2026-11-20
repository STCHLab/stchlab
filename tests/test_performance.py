"""Performance watchdog tests for representative STCHLab workflows.

Run from the repository root with:

    python tests/test_performance.py
    python tests/test_performance.py --cases 5 --seed 42 --timeout 20
"""

import argparse
import math
import multiprocessing as mp
import queue
import random
import sys
import time
import traceback


PASSED = 0
FAILED = []


def positive_int(value):
    """Parse a positive integer command-line argument."""
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return ivalue


def positive_float(value):
    """Parse a positive float command-line argument."""
    fvalue = float(value)
    if fvalue <= 0:
        raise argparse.ArgumentTypeError("value must be a positive number")
    return fvalue


def parse_args():
    """Parse command-line arguments for the performance watchdog."""
    parser = argparse.ArgumentParser(description="Run STCHLab performance watchdog tests.")
    parser.add_argument(
        "--cases",
        type=positive_int,
        default=3,
        help="number of representative cases per multi-case block (default: 3)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="random seed used to generate representative cases (default: 42)",
    )
    parser.add_argument(
        "--timeout",
        type=positive_float,
        default=15.0,
        help="maximum seconds allowed for each test block (default: 15.0)",
    )
    return parser.parse_args()


def require(condition, message):
    """Raise an assertion with a readable message when a check fails."""
    if not condition:
        raise AssertionError(message)


def rand_uniform(rng, low, high):
    """Sample a uniform float."""
    return float(rng.uniform(low, high))


def rand_loguniform(rng, low, high):
    """Sample a log-uniform float."""
    return float(math.exp(rng.uniform(math.log(low), math.log(high))))


def random_delta_range(rng):
    """Generate a safe reduction extent range for CeO2."""
    delta_ox = rand_uniform(rng, 0.005, 0.02)
    delta_r = delta_ox + rand_uniform(rng, 0.01, 0.03)
    delta_r = min(delta_r, 0.05)
    require(delta_r > delta_ox, "random delta range must satisfy delta_R > delta_OX")
    return delta_ox, delta_r


def random_cc_data(rng):
    """Generate one representative counter-current case."""
    delta_ox, delta_r = random_delta_range(rng)
    pO2_min = rand_loguniform(rng, 5.0, 50.0)
    pO2_max = min(pO2_min * rand_uniform(rng, 20.0, 150.0), 5.0e3)
    x_min = rand_loguniform(rng, 1e-4, 5e-3)
    x_max = min(x_min + rand_uniform(rng, 0.02, 0.08), 0.12)

    return {
        "delta_OX": delta_ox,
        "delta_R": delta_r,
        "pO2_min": pO2_min,
        "pO2_max": pO2_max,
        "x_min": x_min,
        "x_max": x_max,
    }


def random_fb_system_data(rng):
    """Generate one representative fixed-bed system case."""
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


def assert_positive_finite(value, label):
    """Check that a scalar result is positive and finite."""
    scalar = float(value)
    require(math.isfinite(scalar), f"{label} must be finite")
    require(scalar > 0.0, f"{label} must be positive")


def block_fixed_bed_infeasible_reduction(cases, seed):
    """Run the known infeasible fixed-bed reduction case."""
    del cases, seed
    from stchlab import CeO2, STCHLab

    stch = STCHLab(CeO2)
    data = {
        "pO2_min": 0.12923484877986363,
        "T_R": 1188.98258902546,
        "V_R": 0.5135601695713496,
        "ndot_inert": 0.011870647010743744,
        "n_MOx": 1.0277160178533251,
        "k_lim_R": 0.00010221064050146226,
        "time_R": 998.7020424640032,
        "delta_OX": 0.009961790369998526,
    }
    result = stch.funFBR.step(data)
    require(result.get("status") == "infeasible", "reduction case should be infeasible")
    require(float(result["frac_Ddelta_R"]) == 0.0, "infeasible reduction should have frac_Ddelta_R=0")


def block_fixed_bed_infeasible_oxidation(cases, seed):
    """Run the known infeasible fixed-bed oxidation case."""
    del cases, seed
    from stchlab import CeO2, STCHLab

    stch = STCHLab(CeO2)
    data = {
        "x_min": 0.1,
        "T_OX": 1400.0,
        "ndot_tot": 1e-3,
        "n_MOx": 1.0,
        "k_lim_OX": 1e-4,
        "time_OX": 300.0,
        "delta_R": 0.001,
    }
    result = stch.funFBOX.step(data)
    require(result.get("status") == "infeasible", "oxidation case should be infeasible")
    require(float(result["frac_Ddelta_OX"]) == 0.0, "infeasible oxidation should have frac_Ddelta_OX=0")


def block_counter_current_system(cases, seed):
    """Run representative uncached counter-current system cases."""
    from stchlab import CeO2, STCHLab

    rng = random.Random(seed + 101)
    stch = STCHLab(CeO2)
    for case in range(1, cases + 1):
        state = stch.funCCSYS.system(random_cc_data(rng), use_db=False)
        assert_positive_finite(state["T_R_max"], f"case {case} T_R_max")
        assert_positive_finite(state["T_OX_min"], f"case {case} T_OX_min")
        assert_positive_finite(state["lmbda"], f"case {case} lmbda")
        assert_positive_finite(state["theta"], f"case {case} theta")
        require(state["delta_R"] > state["delta_OX"], f"case {case}: delta_R must exceed delta_OX")


def block_fixed_bed_system(cases, seed):
    """Run representative uncached fixed-bed system cases."""
    from stchlab import CeO2, STCHLab

    rng = random.Random(seed + 202)
    stch = STCHLab(CeO2)
    for case in range(1, cases + 1):
        state = stch.funFBSYS.system(random_fb_system_data(rng), use_db=False)
        assert_positive_finite(state["Lmbda"], f"case {case} Lmbda")
        assert_positive_finite(state["Theta"], f"case {case} Theta")


BLOCKS = [
    ("fixed-bed infeasible reduction", block_fixed_bed_infeasible_reduction),
    ("fixed-bed infeasible oxidation", block_fixed_bed_infeasible_oxidation),
    ("counter-current system without cache", block_counter_current_system),
    ("fixed-bed system without cache", block_fixed_bed_system),
]


def run_block_worker(block_name, cases, seed, result_queue):
    """Execute one watchdog block inside a child process."""
    block_map = dict(BLOCKS)
    start = time.perf_counter()
    try:
        block_map[block_name](cases, seed)
        result_queue.put(("passed", time.perf_counter() - start, ""))
    except BaseException:
        result_queue.put(("failed", time.perf_counter() - start, traceback.format_exc()))


def run_test(index, title, block_name, cases, seed, timeout):
    """Run one watchdog block with a hard timeout."""
    global PASSED

    print(f"Test {index}: {title}")
    ctx = mp.get_context("spawn")
    result_queue = ctx.Queue()
    process = ctx.Process(target=run_block_worker, args=(block_name, cases, seed, result_queue))

    process.start()
    process.join(timeout)

    if process.is_alive():
        process.terminate()
        process.join(2.0)
        message = f"block exceeded {timeout:.1f}s timeout"
        print(f"FAILED: {message}")
        print()
        FAILED.append((index, title, message))
        return

    try:
        status, elapsed, detail = result_queue.get_nowait()
    except queue.Empty:
        message = f"block exited without reporting a result (exit code {process.exitcode})"
        print(f"FAILED: {message}")
        print()
        FAILED.append((index, title, message))
        return

    if status == "passed":
        print(f"  finished in {elapsed:.3f}s")
        print("PASSED")
        print()
        PASSED += 1
        return

    message = detail.strip().splitlines()[-1] if detail.strip() else "unknown failure"
    print(f"FAILED: {message}")
    print(detail)
    FAILED.append((index, title, message))


def main():
    """Run all performance watchdog blocks."""
    args = parse_args()
    for index, (title, _) in enumerate(BLOCKS, start=1):
        run_test(index, title, title, args.cases, args.seed, args.timeout)

    print(f"Summary: {PASSED} passed, {len(FAILED)} failed")
    if FAILED:
        print("Failed tests:")
        for index, title, message in FAILED:
            print(f"  {index}. {title}: {message}")
        sys.exit(1)


if __name__ == "__main__":
    main()

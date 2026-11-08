#!/usr/bin/env python3
"""Self-test for check_pbl_ysu.py (plain asserts, synthetic columns)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_pbl_ysu as cp  # noqa: E402


def diffuse(m, g, q, dt):
    """One backward-Euler closed-column diffusion step (reference for the tests)."""
    nz = len(m)
    a = np.diag(m.astype(float))
    for k in range(nz - 1):
        c = dt * g[k + 1]
        a[k, k] += c
        a[k + 1, k + 1] += c
        a[k, k + 1] -= c
        a[k + 1, k] -= c
    return np.linalg.solve(a, m * q)


def test_mass_drift_is_zero_for_identical_records():
    m = np.full((3, 4), 2.0)
    q = np.linspace(1.0, 2.0, 12).reshape(3, 4)
    assert cp.max_relative_drift([cp.gas_mass(m, q), cp.gas_mass(m, q)]) == 0.0


def test_mass_drift_detects_a_change():
    m = np.full((3, 4), 2.0)
    q = np.linspace(1.0, 2.0, 12).reshape(3, 4)
    q2 = q.copy()
    q2[1, 2] *= 1.0 + 1e-9
    assert cp.max_relative_drift([cp.gas_mass(m, q), cp.gas_mass(m, q2)]) > 1e-12


def test_diffusion_reduces_the_column_second_moment():
    m = np.array([120.0, 100.0, 80.0, 60.0])
    g = np.array([0.0, 0.5, 0.2, 0.1, 0.0])
    q0 = np.array([4.0, 1.0, 0.5, 0.2])
    q1 = diffuse(m, g, q0, 30.0)
    reduction = cp.second_moment_reduction(m[None, :], q1[None, :], q0[None, :])
    assert reduction.shape == (1,) and reduction[0] > 0.0


def test_uniform_column_has_no_reduction():
    m = np.array([120.0, 100.0, 80.0, 60.0])
    g = np.array([0.0, 0.5, 0.2, 0.1, 0.0])
    q0 = np.full(4, 3.0)
    q1 = diffuse(m, g, q0, 30.0)
    assert abs(cp.second_moment_reduction(m[None, :], q1[None, :], q0[None, :])[0]) <= 1e-12


def test_antidiffusion_is_not_a_reduction():
    m = np.array([120.0, 100.0, 80.0, 60.0])
    q_mixed = np.array([2.0, 1.8, 1.7, 1.6])
    q_sharp = np.array([3.0, 1.0, 1.5, 1.6])
    assert cp.second_moment_reduction(m[None, :], q_sharp[None, :], q_mixed[None, :])[0] < 0.0


def test_domain_fraction_counts_mixing_despite_an_inflow_column():
    # Column 0 is mixed; column 1 gains a little of column 0's difference by horizontal transport,
    # so its own variance rises. The domain still loses variance: the check is the domain total.
    m = np.array([[120.0, 100.0, 80.0, 60.0], [120.0, 100.0, 80.0, 60.0]])
    g = np.array([0.0, 0.5, 0.2, 0.1, 0.0])
    q0 = np.array([[4.0, 1.0, 0.5, 0.2], [4.0, 1.0, 0.5, 0.2]])
    q_off = q0.copy()
    q_on = np.array([diffuse(m[0], g, q0[0], 30.0), q0[1].copy()])
    q_on[1, 0] += 1e-3
    fraction, per_column = cp.variance_mixed_out(m, q_on, q_off, m, q0)
    assert per_column[0] > 0.0 and per_column[1] < 0.0
    assert fraction > 0.0


def test_domain_fraction_is_negative_for_antidiffusion():
    m = np.array([[120.0, 100.0, 80.0, 60.0]])
    q_mixed = np.array([[2.0, 1.8, 1.7, 1.6]])
    q_sharp = np.array([[3.0, 1.0, 1.5, 1.6]])
    fraction, _ = cp.variance_mixed_out(m, q_sharp, q_mixed, m, q_mixed)
    assert fraction < 0.0


def test_stronger_compares_only_the_named_gases():
    strong = {"run": "a", "gases": {"qNO": {"variance_mixed_out_fraction": 0.1},
                                    "qNO2": {"variance_mixed_out_fraction": 1e-5}}}
    weak = {"run": "b", "gases": {"qNO": {"variance_mixed_out_fraction": 1e-3},
                                  "qNO2": {"variance_mixed_out_fraction": 2e-5}}}
    assert cp.weaker_gases(strong, weak, ["qNO"]) == []
    assert cp.weaker_gases(strong, weak, ["qNO", "qNO2"]) == ["qNO2"]


def test_sign_violations_allow_roundoff_below_zero():
    q0 = np.array([[1e-9, 2e-9], [0.0, 0.0]])
    q = np.array([[-1e-22, 1e-9], [-5e-21, 0.0]])        # within -(1e-20 + 1e-12*max_column(q0))
    assert cp.sign_violations(q, q0) == 0


def test_sign_violations_count_negatives_and_nonfinite():
    q0 = np.array([[1e-9, 2e-9], [0.0, 0.0]])
    q = np.array([[-2e-20, 1e-9], [np.nan, -2e-20]])     # below -(1e-20 + 2e-21); NaN; below -1e-20
    assert cp.sign_violations(q, q0) == 3


def test_bitwise_difference_counts_one_ulp():
    a = np.linspace(0.0, 1.0, 10)
    b = a.copy()
    assert cp.count_bitwise_differences(a, b) == 0
    b[3] = np.nextafter(b[3], 2.0)
    assert cp.count_bitwise_differences(a, b) == 1


if __name__ == "__main__":
    tests = [value for name, value in sorted(globals().items()) if name.startswith("test_")]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"test_check_pbl_ysu: {len(tests)} of {len(tests)} tests passed")

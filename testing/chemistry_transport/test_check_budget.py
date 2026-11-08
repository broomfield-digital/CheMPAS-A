#!/usr/bin/env python3
"""Self-test for check_budget.py (plain asserts; run: python test_check_budget.py)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_budget as cb  # noqa: E402


def test_thresholds_are_frozen():
    assert cb.MASS_RESID_PER_CALL == 1e-12
    assert cb.MASS_RESID_1000_CALLS == 1e-10
    assert cb.FULL_RUN_RESID == 1e-8
    assert (cb.UNIFORM_ATOL, cb.UNIFORM_RTOL) == (1e-20, 1e-12)
    assert (cb.PARALLEL_ATOL, cb.PARALLEL_RTOL) == (1e-18, 1e-10)


def test_conservative_redistribution_has_no_residual():
    m_d = np.array([100.0, 160.0, 150.0])
    before = np.array([[4e-9, 2e-9, 1e-9]])
    moved = 0.5 * before[0, 0] * m_d[0]
    after = before.copy()
    after[0, 0] -= moved / m_d[0]
    after[0, 1] += moved / m_d[1]
    relative, absolute = cb.relative_residual(m_d, before, after)
    assert relative.shape == (1,) and relative[0] <= cb.MASS_RESID_PER_CALL
    assert absolute[0] <= 1e-12 * cb.column_mass(m_d, before)[0]


def test_source_is_measured_and_boundary_flux_is_credited():
    m_d = np.array([100.0, 160.0, 150.0])
    before = np.full((1, 3), 3e-9)
    after = before.copy()
    after[0, 2] *= 2.0
    relative, _ = cb.relative_residual(m_d, before, after)
    assert abs(relative[0] - 4.5e-7 / 1.23e-6) <= 1e-12
    relative, _ = cb.relative_residual(m_d, before, after, boundary_flux=1.5e-7, dt=3.0)
    assert relative[0] <= 1e-12


def test_zero_species_normalization_is_finite():
    m_d = np.array([100.0, 160.0, 150.0])
    zero = np.zeros((1, 3))
    relative, absolute = cb.relative_residual(m_d, zero, zero)
    assert relative[0] == 0.0 and absolute[0] == 0.0
    relative, _ = cb.relative_residual(m_d, zero, np.full((1, 3), 1e-40))
    assert np.isfinite(relative[0]) and abs(relative[0] - 1e-10) <= 1e-22


def test_nonnegativity_allows_only_roundoff():
    q0 = np.array([[1e-9, 2e-9, 3e-9]])
    assert cb.nonnegative(np.array([[-1e-22, 1e-9, 3e-9]]), q0)
    assert not cb.nonnegative(np.array([[-1e-19, 1e-9, 3e-9]]), q0)


def test_zero_species_stays_zero():
    q0 = np.array([[0.0, 0.0, 0.0], [1e-9, 1e-9, 1e-9]])
    assert cb.zero_species_stays_zero(np.array([[0.0, 0.0, 0.0], [2e-9, 0.0, 1e-9]]), q0)
    assert not cb.zero_species_stays_zero(np.array([[0.0, 1e-30, 0.0], [1e-9, 1e-9, 1e-9]]), q0)


def test_uniform_profile_tolerance():
    q = np.full((2, 5), 3e-9)
    assert cb.uniform_profile_unchanged(q * (1 + 1e-13), q)
    assert not cb.uniform_profile_unchanged(q * (1 + 1e-9), q)


def test_disabled_comparison_is_bitwise():
    base = {"qNO": np.array([1.0, 2.0, np.nan]), "qO3": np.array([0.0, 5.0])}
    same = {"qNO": np.array([1.0, 2.0, np.nan]), "qO3": np.array([0.0, 5.0])}
    assert cb.disabled_matches_baseline(same, base)
    ulp = {"qNO": np.array([1.0, np.nextafter(2.0, 3.0), np.nan]), "qO3": np.array([0.0, 5.0])}
    assert not cb.disabled_matches_baseline(ulp, base)
    signed = {"qNO": np.array([1.0, 2.0, np.nan]), "qO3": np.array([-0.0, 5.0])}
    assert not cb.disabled_matches_baseline(signed, base)
    assert not cb.disabled_matches_baseline({"qNO": base["qNO"]}, base)


if __name__ == "__main__":
    tests = [value for name, value in sorted(globals().items()) if name.startswith("test_")]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"test_check_budget: {len(tests)} of {len(tests)} tests passed")

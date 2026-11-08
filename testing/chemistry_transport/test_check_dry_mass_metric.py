#!/usr/bin/env python3
"""Self-test for check_dry_mass_metric.py (plain asserts)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_dry_mass_metric as cm  # noqa: E402


def mesh(terrain):
    zw = np.linspace(0.0, 20000.0, 11)                     # zeta at interfaces
    h = np.array([0.0, terrain, 0.5 * terrain])            # surface height per cell
    ztop = zw[-1]
    zgrid = h[:, None] + zw[None, :] * (1.0 - h[:, None] / ztop)
    zz = (zw[1:] - zw[:-1])[None, :] / np.diff(zgrid, axis=1)
    rdzw = 1.0 / (zw[1:] - zw[:-1])
    return zz, zgrid, rdzw


def test_consistent_metric_passes():
    zz, zgrid, rdzw = mesh(500.0)
    assert cm.metric_error(zz, zgrid, rdzw) <= 1e-12


def test_inconsistent_metric_fails():
    zz, zgrid, rdzw = mesh(500.0)
    zz = zz.copy()
    zz[1, 3] *= 1.0 + 1e-9
    assert cm.metric_error(zz, zgrid, rdzw) > 1e-12


if __name__ == "__main__":
    tests = [value for name, value in sorted(globals().items()) if name.startswith("test_")]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"test_check_dry_mass_metric: {len(tests)} of {len(tests)} tests passed")

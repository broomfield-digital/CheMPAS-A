#!/usr/bin/env python3
"""Check that the transport dry-mass weight equals the weight MPAS conserves.

The MICM gas transport kernels weight layers by m_d = zz*rho_zz*dz, with
dz = zgrid(k+1) - zgrid(k). The MPAS scalar advance applies the vertical flux
divergence with rdzw (mesh 'rdzw', used as rdnw in atm_advance_scalars), so
the dry tracer mass MPAS conserves per column is sum(rho_zz*q/rdzw). Both
weights agree when zz*dz*rdzw = 1. The initialization defines
zz = (zw(k+1)-zw(k))/(zgrid(k+1)-zgrid(k)) and rdzw = 1/(zw(k+1)-zw(k)),
so the identity holds to roundoff on any terrain-following mesh. This script
checks it on actual MPAS files (docs/plans/pbl-transport.md, shared item 4).

Usage: check_dry_mass_metric.py FILE [FILE ...]   (files with zz, zgrid, rdzw)
"""
import sys

import numpy as np

TOLERANCE = 1e-12


def metric_error(zz, zgrid, rdzw):
    """max |zz*dz*rdzw - 1| over all cells and layers."""
    dz = np.diff(np.asarray(zgrid), axis=-1)
    return float(np.max(np.abs(np.asarray(zz) * dz * np.asarray(rdzw)[None, :] - 1.0)))


def main(paths):
    from netCDF4 import Dataset

    status = 0
    for path in paths:
        with Dataset(path) as ds:
            zz, zgrid, rdzw = ds["zz"][:], ds["zgrid"][:], ds["rdzw"][:]
            terrain = float(np.max(zgrid[:, 0]) - np.min(zgrid[:, 0]))
        error = metric_error(zz, zgrid, rdzw)
        ok = error <= TOLERANCE
        status |= not ok
        print(f"{'PASS' if ok else 'FAIL'} {path}: {zz.shape[0]} cells, {zz.shape[1]} layers, "
              f"terrain relief {terrain:g} m, max|zz*dz*rdzw - 1| = {error:.3e}")
    return status


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1:]))

#!/usr/bin/env python3
"""Make a horizontally uniform, resting initial state from a supercell init file.

Every time-dependent cell field is replaced by its values in one reference
column far from the warm bubble. The state is then made dry and set equal to
the dynamics base state (qv = 0, rho = rho_base, theta = theta_base), with
zero wind. MPAS integrates perturbations from that base state, so the state
has no buoyancy or pressure-gradient forcing and stays at rest. A moist column
from the test-case initialization is not in the model's discrete hydrostatic
balance and drifts at a few mm/s. The MICM tracers get positive, vertically
varying (horizontally uniform) profiles. Mesh fields are not touched. Used by check_realized_increment.sh
to compare realized tracer increments with requested ones (P2 of
docs/plans/pbl-transport.md).

Usage: make_resting_init.py SOURCE_INIT OUTPUT_INIT [X_KM Y_KM]
"""
import shutil
import sys

import numpy as np
from netCDF4 import Dataset

TRACER_PROFILES = {
    # name: f(k), k = 0..nz-1, kg/kg
    "qO3": lambda k, nz: 8.275862068965518e-08 * (1.0 + 0.5 * np.sin(0.3 * k)),
    "qNO2": lambda k, nz: 1.0e-9 * (1.0 + k / nz),
    "qNO": lambda k, nz: np.where((k >= nz // 6) & (k < nz // 6 + 5), 2.0e-9, 1.0e-12),
}


def main(source, output, x_km=5.0, y_km=5.0):
    shutil.copyfile(source, output)
    with Dataset(output, "r+") as ds:
        x = ds["xCell"][:] / 1000.0
        y = ds["yCell"][:] / 1000.0
        ref = int(np.argmin((x - x_km) ** 2 + (y - y_km) ** 2))
        nz = len(ds.dimensions["nVertLevels"])
        changed = []
        for name, var in ds.variables.items():
            dims = var.dimensions
            if "Time" in dims and "nCells" in dims and var.dtype.kind == "f":
                data = var[:]
                cell_axis = dims.index("nCells")
                column = np.take(data, [ref], axis=cell_axis)
                var[:] = np.broadcast_to(column, data.shape)
                changed.append(name)
        ds["u"][:] = 0.0
        ds["w"][:] = 0.0
        ds["qv"][:] = 0.0
        ds["rho"][:] = ds["rho_base"][:]
        ds["theta"][:] = ds["theta_base"][:]
        k = np.arange(nz)
        for name, profile in TRACER_PROFILES.items():
            if name in ds.variables:
                ds[name][:] = np.broadcast_to(profile(k, nz), ds[name].shape)
    print(f"reference cell {ref} at ({x[ref]:.1f}, {y[ref]:.1f}) km; dry base state at rest; "
          f"uniform fields: {', '.join(changed)}")


if __name__ == "__main__":
    if len(sys.argv) not in (3, 5):
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], *map(float, sys.argv[3:]))

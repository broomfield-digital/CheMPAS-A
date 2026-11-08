#!/usr/bin/env python3
"""Give a supercell init file MICM gas profiles that PBL mixing can change.

The baseline supercell init has vertically uniform O3 and no NO or NO2, which
PBL mixing leaves unchanged. This copy keeps the supercell meteorology and
mesh, and sets horizontally uniform gas profiles with vertical structure in
the lowest kilometres (z is the layer-centre height above the surface):

  qNO   surface-confined:    2e-9 * exp(-z / 300 m) + 1e-12 kg/kg
  qNO2  a step: 1e-9 below 1000 m, 1e-11 above (on the supercell grid the
        step is the interface near 1.2 km, above both PBL tops, so it is
        mixed by YSU's free-atmosphere coefficient)
  qO3   depleted near the surface, increasing to the baseline value at 3 km:
        8.276e-8 * (0.5 + 0.5 * min(z / 3000 m, 1))

--skin-temperature T sets skintemp and sst to T everywhere. The baseline
surface (288 K) is colder than the lowest model layer (about 299 K), so the
surface layer is stable; a surface warmer than the air drives an upward heat
flux and a convective PBL. Used by check_pbl_ysu.sh (P3 of
docs/plans/pbl-transport.md).

Usage: make_pbl_init.py SOURCE_INIT OUTPUT_INIT [--skin-temperature T]
"""
import argparse
import shutil

import numpy as np
from netCDF4 import Dataset

O3_BASELINE = 8.275862068965518e-08


def gas_profiles(z):
    """Gas mixing ratios [kg/kg] at layer-centre heights z [m] above the surface."""
    return {
        "qNO": 2.0e-9 * np.exp(-z / 300.0) + 1.0e-12,
        "qNO2": np.where(z < 1000.0, 1.0e-9, 1.0e-11),
        "qO3": O3_BASELINE * (0.5 + 0.5 * np.minimum(z / 3000.0, 1.0)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("source")
    parser.add_argument("output")
    parser.add_argument("--skin-temperature", type=float, default=None)
    args = parser.parse_args()

    shutil.copyfile(args.source, args.output)
    with Dataset(args.output, "r+") as ds:
        zgrid = np.asarray(ds["zgrid"][:])
        z = 0.5 * (zgrid[:, :-1] + zgrid[:, 1:]) - zgrid[:, :1]
        for name, values in gas_profiles(z).items():
            var = ds[name]
            assert var.dimensions == ("Time", "nCells", "nVertLevels"), var.dimensions
            var[:] = np.broadcast_to(values, var.shape)
        if args.skin_temperature is not None:
            for name in ("skintemp", "sst"):
                ds[name][:] = args.skin_temperature
    surface = "baseline surface" if args.skin_temperature is None else f"skintemp = sst = {args.skin_temperature} K"
    print(f"{args.output}: gas profiles with vertical structure; {surface}")


if __name__ == "__main__":
    main()

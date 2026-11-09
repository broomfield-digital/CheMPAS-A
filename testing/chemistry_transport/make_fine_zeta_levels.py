#!/usr/bin/env python3
"""Write a specified-zeta-levels file for a supercell grid that resolves the
boundary layer (P6 of docs/plans/pbl-transport.md: "a finer-grid fixture").

The baseline supercell grid (supercell_zeta_levels.txt) has interfaces at 0,
299, 712 and 1182 m, so its PBL spans only one or two layers. This grid keeps
the baseline's 50-km top and puts uniform DZ_FINE layers up to FINE_TOP, then
lets the spacing grow geometrically by GROWTH per layer up to DZ_MAX, then
holds DZ_MAX, with the layers above FINE_TOP scaled so the top is exactly ZTOP.
init_atmosphere_model reads it through config_specified_zeta_levels, one
interface height [m] per line, from 0 to ZTOP (nVertLevels + 1 lines).

Usage: make_fine_zeta_levels.py OUT [--dz-fine 50] [--fine-top 1500]
           [--growth 1.1] [--dz-max 1000] [--ztop 50000]
       make_fine_zeta_levels.py --self-test
"""
import argparse
import sys

import numpy as np


def fine_levels(dz_fine=50.0, fine_top=1500.0, growth=1.1, dz_max=1000.0, ztop=50000.0):
    """Interface heights [m], 0 .. ztop."""
    n_fine = int(round(fine_top / dz_fine))
    z = list(np.arange(n_fine + 1) * dz_fine)
    upper = []
    dz = dz_fine
    while z[-1] + sum(upper) < ztop:
        dz = min(dz * growth, dz_max)
        upper.append(dz)
    upper = np.array(upper)
    # scale the stretched part so the last interface is exactly ztop
    upper *= (ztop - z[-1]) / upper.sum()
    return np.concatenate([z, z[-1] + np.cumsum(upper)])


def self_test():
    z = fine_levels()
    dz = np.diff(z)
    assert z[0] == 0.0 and abs(z[-1] - 50000.0) < 1e-6, "spans 0 to the 50-km top"
    assert np.all(dz > 0.0), "monotonic"
    assert np.all(dz[z[1:] <= 1500.0 + 1e-6] <= 50.0 + 1e-9), "at most 50 m below 1.5 km"
    assert np.count_nonzero(z <= 1500.0 + 1e-6) - 1 >= 30, "at least 30 layers below 1.5 km"
    assert np.all(dz[1:] / dz[:-1] <= 1.1 + 1e-9), "spacing grows by at most 10% per layer"
    assert dz.max() <= 1000.0 * (50000.0 - 1500.0) / (50000.0 - 1500.0) + 50.0, "no layer much above 1 km"
    print(f"make_fine_zeta_levels self-test: PASS ({len(z) - 1} layers, {dz.min():.0f}-{dz.max():.0f} m)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out", nargs="?")
    ap.add_argument("--dz-fine", type=float, default=50.0)
    ap.add_argument("--fine-top", type=float, default=1500.0)
    ap.add_argument("--growth", type=float, default=1.1)
    ap.add_argument("--dz-max", type=float, default=1000.0)
    ap.add_argument("--ztop", type=float, default=50000.0)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)
    if args.self_test:
        self_test()
        return 0
    if not args.out:
        ap.error("OUT is required")
    z = fine_levels(args.dz_fine, args.fine_top, args.growth, args.dz_max, args.ztop)
    with open(args.out, "w") as f:
        f.write("".join(f"{v:.6f}\n" for v in z))
    print(f"{args.out}: {len(z) - 1} layers, interfaces 0 .. {z[-1]:.0f} m")
    return 0


if __name__ == "__main__":
    sys.exit(main())

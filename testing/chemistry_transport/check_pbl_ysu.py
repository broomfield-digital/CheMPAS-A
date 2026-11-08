#!/usr/bin/env python3
"""Checks for the in-model YSU gas PBL qualification (check_pbl_ysu.sh; P3 of
docs/plans/pbl-transport.md). Reads MPAS history files (output.nc).

  mixing ON_RUN OFF_RUN --json PATH
      ON_RUN and OFF_RUN differ only in config_micm_pbl_transport. Checks that
      (1) every non-gas output field is bitwise identical (the gases are
      passive, so transport touches nothing else); (2) YSU diagnosed a PBL in
      every column (hpbl > 0) after the initial record; (3) each gas in
      ON_RUN is finite and no lower than -(1e-20 + 1e-12*max_column(q0))
      kg/kg, and its domain mass drifts by at most 1e-13 relative; (4) at the final record
      the domain total of each gas's dry-mass-weighted column variance (about
      the column mean) is lower with transport than without: every gradient
      gas is mixed down-gradient. Writes the metrics to PATH, including the
      number of columns whose own variance fell. The operator lowers every
      column's variance at every step, but resolved transport then moves the
      difference it made between columns, so near the storm a column can gain
      variance; the column reference tests check the operator per column.
  stronger STRONG_JSON WEAK_JSON [--gases GAS ...]
      Each named gas (default: all) has more variance mixed out (domain
      fraction) in STRONG than WEAK. Name the gases whose gradients lie inside
      both boundary layers: a gradient above the PBL top is mixed by the
      scheme's free-atmosphere coefficient in either regime.
  compare RUN_A RUN_B
      Every output field is bitwise identical.

Layer mass m = areaCell * dz * zz * rho_zz, the dry-air weight MPAS conserves.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

GASES = ("qNO", "qNO2", "qO3")
MASS_TOLERANCE = 1e-13


def gas_mass(m, q):
    return float(np.sum(m * q))


def max_relative_drift(masses):
    m0 = masses[0]
    return max(abs(m - m0) / abs(m0) for m in masses) if m0 != 0 else max(abs(m) for m in masses)


def column_variance(m, q):
    """Dry-mass-weighted variance about the column mean, summed over the column: (nCells,)."""
    mean = np.sum(m * q, axis=1) / np.sum(m, axis=1)
    return np.sum(m * (q - mean[:, None]) ** 2, axis=1)


def second_moment_reduction(m, q_mixed, q_reference):
    """Column variance removed by mixing: variance(q_reference) - variance(q_mixed), per column."""
    return column_variance(m, q_reference) - column_variance(m, q_mixed)


def variance_mixed_out(m, q_mixed, q_reference, m0, q0):
    """Domain fraction of the initial column variance (m0, q0) removed by mixing, and the per-column reduction."""
    per_column = second_moment_reduction(m, q_mixed, q_reference)
    return float(per_column.sum() / column_variance(m0, q0).sum()), per_column


def weaker_gases(strong, weak, gases):
    """Gases whose domain fraction of variance mixed out is not larger in STRONG than in WEAK."""
    return [gas for gas in gases
            if not strong["gases"][gas]["variance_mixed_out_fraction"] >
            weak["gases"][gas]["variance_mixed_out_fraction"]]


def sign_violations(q, q0):
    """Values (nCells, nVertLevels) that are non-finite or below the frozen roundoff floor
    -(1e-20 + 1e-12*max_column(q0)) kg/kg (docs/plans/pbl-transport.md thresholds)."""
    floor = -(1e-20 + 1e-12 * np.max(q0, axis=1))
    return int(np.count_nonzero(~np.isfinite(q)) + np.count_nonzero(q < floor[:, None]))


def count_bitwise_differences(a, b):
    a = np.ascontiguousarray(a)
    b = np.ascontiguousarray(b)
    if a.shape != b.shape or a.dtype != b.dtype:
        return max(a.size, b.size)
    if a.dtype.kind == "f":
        return int(np.count_nonzero(a.view(f"u{a.itemsize}") != b.view(f"u{b.itemsize}")))
    return int(np.count_nonzero(a != b))


class Run:
    def __init__(self, run_dir):
        self.path = Path(run_dir) / "output.nc"
        self.ds = Dataset(self.path)
        self.ds.set_auto_mask(False)
        self.times = ["".join(c.decode() for c in row).strip() for row in self.ds["xtime"][:]]
        area = self.ds["areaCell"][:]
        dz = np.diff(self.ds["zgrid"][:], axis=1)
        self.geometry = area[:, None] * dz * self.ds["zz"][:]

    def layer_mass(self, record):
        return self.geometry * self.ds["rho_zz"][record]

    def gas(self, name, record):
        return self.ds[name][record]

    def numeric_variables(self):
        return {name for name, var in self.ds.variables.items() if var.dtype.kind in "iuf"}


def check(condition, message, failures):
    print(("PASS " if condition else "FAIL ") + message)
    if not condition:
        failures.append(message)


def mass_drifts(run):
    out = {}
    for gas in GASES:
        masses = [gas_mass(run.layer_mass(t), run.gas(gas, t)) for t in range(len(run.times))]
        out[gas] = max_relative_drift(masses)
    return out


def compare_runs(a, b, exclude=()):
    names_a, names_b = a.numeric_variables() - set(exclude), b.numeric_variables() - set(exclude)
    differing = {}
    for name in sorted(names_a & names_b):
        n = count_bitwise_differences(a.ds[name][:], b.ds[name][:])
        if n:
            differing[name] = n
    return names_a ^ names_b, differing


def cmd_mixing(args):
    failures = []
    on, off = Run(args.on), Run(args.off)
    check(on.times == off.times and len(on.times) > 1, f"same output times ({len(on.times)} records)", failures)
    if failures:
        return 1
    unmatched, differing = compare_runs(on, off, exclude=GASES)
    check(not unmatched and not differing,
          f"non-gas fields bitwise identical with and without gas transport (differing: {differing or 'none'})",
          failures)
    hpbl = on.ds["hpbl"][1:]
    check(bool(np.all(hpbl > 0)), f"hpbl > 0 in every column after the initial record (min {hpbl.min():.1f} m)",
          failures)
    for gas in GASES:
        q0 = on.gas(gas, 0)
        bad = sum(sign_violations(on.gas(gas, t), q0) for t in range(len(on.times)))
        check(bad == 0, f"{gas} finite and nonnegative to roundoff at every record with transport "
                        f"({bad} violations)", failures)
    drifts = mass_drifts(on)
    for gas, drift in drifts.items():
        check(drift <= MASS_TOLERANCE, f"{gas} domain mass conserved with transport (max drift {drift:.2e})",
              failures)

    final = len(on.times) - 1
    m0, m = on.layer_mass(0), on.layer_mass(final)
    metrics = {"run": str(on.path.parent), "times": on.times,
               "hpbl_final_m": {"min": float(on.ds["hpbl"][final].min()),
                                "median": float(np.median(on.ds["hpbl"][final])),
                                "max": float(on.ds["hpbl"][final].max())},
               "hfx_final_W_m2": {"min": float(on.ds["hfx"][final].min()),
                                  "median": float(np.median(on.ds["hfx"][final])),
                                  "max": float(on.ds["hfx"][final].max())},
               "mass_drift": drifts, "gases": {}}
    for gas in GASES:
        q_on, q_off = on.gas(gas, final), off.gas(gas, final)
        fraction, reduction = variance_mixed_out(m, q_on, q_off, m0, on.gas(gas, 0))
        initial = column_variance(m0, on.gas(gas, 0))
        lowest_change = (q_on[:, 0] - q_off[:, 0]) / np.maximum(np.abs(q_off[:, 0]), 1e-30)
        difference = m * (q_on - q_off)
        burden = difference.sum(axis=1) / np.maximum(np.abs(difference).sum(axis=1), 1e-300)
        not_reduced = reduction <= 0
        metrics["gases"][gas] = {
            "variance_mixed_out_fraction": fraction,
            "column_reduction_min_relative": float(np.min(reduction / initial)),
            "columns_reduced": int(np.count_nonzero(~not_reduced)),
            "columns": int(reduction.size),
            # sum(m*d)/sum(m*|d|) over the columns whose variance did not fall: +1 or -1 when the
            # transport difference there arrived or left by resolved transport, 0 for local mixing
            "not_reduced_burden_ratio": [float(v) for v in burden[not_reduced]],
            "lowest_level_relative_change_median": float(np.median(lowest_change)),
        }
        check(fraction > 0,
              f"{gas}: mixed down-gradient (domain fraction of column variance mixed out {fraction:.3e}; "
              f"{metrics['gases'][gas]['columns_reduced']}/{reduction.size} columns reduced)", failures)
    metrics["failures"] = failures
    Path(args.json).write_text(json.dumps(metrics, indent=2) + "\n")
    return 1 if failures else 0


def cmd_stronger(args):
    failures = []
    strong = json.loads(Path(args.strong).read_text())
    weak = json.loads(Path(args.weak).read_text())
    weaker = weaker_gases(strong, weak, args.gases)
    for gas in args.gases:
        s = strong["gases"][gas]["variance_mixed_out_fraction"]
        w = weak["gases"][gas]["variance_mixed_out_fraction"]
        check(gas not in weaker, f"{gas}: {Path(strong['run']).name} mixes more than {Path(weak['run']).name} "
                                 f"({s:.3e} > {w:.3e})", failures)
    return 1 if failures else 0


def cmd_compare(args):
    failures = []
    a, b = Run(args.a), Run(args.b)
    unmatched, differing = compare_runs(a, b)
    check(a.times == b.times, f"same output times ({len(a.times)} records)", failures)
    check(not unmatched, f"same output variables (unmatched: {sorted(unmatched) or 'none'})", failures)
    check(not differing, f"every common output field bitwise identical (differing: {differing or 'none'})",
          failures)
    return 1 if failures else 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("mixing")
    p.add_argument("on")
    p.add_argument("off")
    p.add_argument("--json", required=True)
    p.set_defaults(func=cmd_mixing)
    p = sub.add_parser("stronger")
    p.add_argument("strong")
    p.add_argument("weak")
    p.add_argument("--gases", nargs="+", choices=GASES, default=list(GASES))
    p.set_defaults(func=cmd_stronger)
    p = sub.add_parser("compare")
    p.add_argument("a")
    p.add_argument("b")
    p.set_defaults(func=cmd_compare)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

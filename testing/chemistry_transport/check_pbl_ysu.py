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
  compare-at RUN_A RUN_B TIME [--atol A --rtol R] [--gases-only | --fields F ...] [--exclude F ...]
             [--json PATH]
      The records at xtime TIME, with cells aligned by indexToCellID: bitwise
      identical by default, else |a - b| <= A + R*|b| in every element (restart
      against continuous, and different MPI decompositions; P5). Reports the
      largest differences of every differing field.

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


def _xtimes(ds):
    return ["".join(c.decode() for c in row).strip() for row in ds["xtime"][:]]


def compare_at(path_a, path_b, time, atol=0.0, rtol=0.0, fields=None, exclude=()):
    """Compare the record at TIME (xtime) of two MPAS history files, with cells aligned by indexToCellID.

    With atol = rtol = 0 an element matches only if it is bitwise identical; otherwise it matches when
    |a - b| <= atol + rtol*|b|. fields restricts the comparison to the named variables (default: every
    numeric variable in both files), less the names in exclude. Returns {"passed", "time", "missing", "fields": {name: {"n_outside",
    "size", "max_abs", "max_rel"}}}."""
    result = {"time": time, "fields": {}, "missing": [], "passed": False}
    with Dataset(path_a) as da, Dataset(path_b) as db:
        da.set_auto_mask(False)
        db.set_auto_mask(False)
        ta, tb = _xtimes(da), _xtimes(db)
        if time not in ta or time not in tb:
            result["missing"].append(f"record {time} (a: {time in ta}, b: {time in tb})")
            return result
        ra, rb = ta.index(time), tb.index(time)
        order_a = np.argsort(da["indexToCellID"][:])
        order_b = np.argsort(db["indexToCellID"][:])
        if not np.array_equal(da["indexToCellID"][:][order_a], db["indexToCellID"][:][order_b]):
            result["missing"].append("different cell sets")
            return result
        numeric = lambda d: {n for n, v in d.variables.items() if v.dtype.kind in "iuf"}
        names = sorted(numeric(da) & numeric(db)) if fields is None else list(fields)
        names = [n for n in names if n not in set(exclude)]
        for name in names:
            if name not in da.variables or name not in db.variables:
                result["missing"].append(name)
                continue
            va, vb = da[name], db[name]
            xa, xb = va[:], vb[:]
            if va.dimensions and va.dimensions[0] == "Time":
                xa, xb = xa[ra], xb[rb]
            dims = va.dimensions[1:] if va.dimensions and va.dimensions[0] == "Time" else va.dimensions
            if "nCells" in dims:
                axis = dims.index("nCells")
                xa, xb = np.take(xa, order_a, axis=axis), np.take(xb, order_b, axis=axis)
            xa, xb = np.asarray(xa), np.asarray(xb)
            if atol == 0.0 and rtol == 0.0:
                n_out = count_bitwise_differences(xa, xb)
            else:
                n_out = int(np.count_nonzero(~(np.abs(xa.astype(float) - xb.astype(float))
                                               <= atol + rtol * np.abs(xb.astype(float)))))
            diff = np.abs(xa.astype(float) - xb.astype(float)) if xa.shape == xb.shape else np.array([np.inf])
            scale = np.maximum(np.abs(xb.astype(float)), np.finfo(float).tiny) if xa.shape == xb.shape else 1.0
            result["fields"][name] = {"n_outside": n_out, "size": int(xa.size),
                                      "max_abs": float(diff.max()) if diff.size else 0.0,
                                      "max_rel": float((diff / scale).max()) if diff.size else 0.0}
    result["passed"] = not result["missing"] and all(r["n_outside"] == 0 for r in result["fields"].values())
    return result


def cmd_compare_at(args):
    fields = list(GASES) if args.gases_only else args.fields
    result = compare_at(Path(args.a) / "output.nc", Path(args.b) / "output.nc", args.time,
                        atol=args.atol, rtol=args.rtol, fields=fields, exclude=args.exclude or ())
    mode = "bitwise" if args.atol == 0.0 and args.rtol == 0.0 else f"atol={args.atol:g}, rtol={args.rtol:g}"
    differing = {n: r for n, r in result["fields"].items() if r["n_outside"]}
    for name, r in sorted(differing.items()):
        print(f"  {name}: {r['n_outside']} of {r['size']} outside; max |diff| {r['max_abs']:.3e}, "
              f"max rel {r['max_rel']:.3e}")
    if args.json:
        Path(args.json).write_text(json.dumps(result, indent=1) + "\n")
    failures = []
    check(not result["missing"], f"record {args.time} and fields present (missing: {result['missing'] or 'none'})",
          failures)
    check(not differing, f"{len(result['fields'])} fields at {args.time} agree ({mode}); differing: "
          f"{sorted(differing) or 'none'}", failures)
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
    p = sub.add_parser("compare-at")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("time", help="xtime of the record, e.g. 0000-01-01_18:02:24")
    p.add_argument("--atol", type=float, default=0.0)
    p.add_argument("--rtol", type=float, default=0.0)
    p.add_argument("--gases-only", action="store_true", help=f"compare only {', '.join(GASES)}")
    p.add_argument("--fields", nargs="+")
    p.add_argument("--exclude", nargs="+", help="leave these fields out (e.g. the gases)")
    p.add_argument("--json")
    p.set_defaults(func=cmd_compare_at)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

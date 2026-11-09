#!/usr/bin/env python3
"""Reactive regression for MICM gas PBL transport (docs/plans/pbl-transport.md,
"Coupled/regression"; P6).

  invariants RUN --invariant NAME=SPECIES:MOLAR_MASS,... [...] [--max-drift 1e-12]
             [--nonnegative SPECIES ...] [--json PATH]
      For each invariant, the domain total of sum(q_s/M_s) weighted by the dry
      air mass areaCell*dz*zz*rho_zz, per output record, and its largest drift
      relative to the first record. Linear invariants of the mechanism (LNOx:
      nitrogen NO + NO2 without lightning, odd oxygen NO2 + O3) are conserved
      by the reactions, by closed-column transport and by resolved advection
      on the periodic supercell domain.
      Each named species must stay above -(1e-20 + 1e-12*max(q at record 0)).
  --self-test
      Synthetic checks of the invariant and sign arithmetic.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np


def invariant_totals(mass, fields, terms):
    """Domain total of sum(q/M) over terms [(species, molar mass)], per record: mass (T, cells, levels) [kg],
    fields {species: (T, cells, levels)} [kg/kg]."""
    total = np.zeros(mass.shape[0])
    for name, molar in terms:
        total += (mass * fields[name]).sum(axis=(1, 2)) / molar
    return total


def max_drift(totals):
    return float(np.max(np.abs(totals - totals[0])) / max(abs(totals[0]), np.finfo(float).tiny))


def floor_violations(q):
    floor = -(1.0e-20 + 1.0e-12 * float(np.max(np.abs(q[0]))))
    return int(np.count_nonzero(~np.isfinite(q))) + int(np.count_nonzero(q < floor))


def parse_invariant(text):
    name, terms = text.split("=", 1)
    return name, [(s, float(m)) for s, m in (t.split(":") for t in terms.split(","))]


def cmd_invariants(args):
    from netCDF4 import Dataset
    with Dataset(Path(args.run) / "output.nc") as d:
        d.set_auto_mask(False)
        dz = np.diff(d["zgrid"][:], axis=1)
        mass = d["areaCell"][:][None, :, None] * dz[None] * d["zz"][:][None] * d["rho_zz"][:]
        invariants = [parse_invariant(t) for t in args.invariant]
        names = {s for _, terms in invariants for s, _ in terms} | set(args.nonnegative)
        fields = {s: np.asarray(d[f"q{s}"][:], dtype=float) for s in names}
    failures, report = [], {"invariants": {}, "sign_violations": {}}
    for name, terms in invariants:
        totals = invariant_totals(mass, fields, terms)
        drift = max_drift(totals)
        report["invariants"][name] = {"drift": drift, "first": float(totals[0]), "last": float(totals[-1])}
        ok = drift <= args.max_drift
        print(f"{'PASS' if ok else 'FAIL'} invariant {name}: largest relative drift {drift:.3e} "
              f"(limit {args.max_drift:g}) over {len(totals)} records")
        if not ok:
            failures.append(name)
    for s in args.nonnegative:
        n = floor_violations(fields[s])
        report["sign_violations"][s] = n
        print(f"{'PASS' if n == 0 else 'FAIL'} q{s}: {n} values non-finite or below the floor")
        if n:
            failures.append(f"q{s}")
    if args.json:
        Path(args.json).write_text(json.dumps(report, indent=1) + "\n")
    return 1 if failures else 0


def self_test():
    rng = np.random.default_rng(1)
    mass = np.broadcast_to(rng.uniform(1.0, 2.0, (5, 4)), (3, 5, 4))
    no = np.broadcast_to(rng.uniform(0.5, 1.0, (5, 4)), (3, 5, 4)).copy()
    o3 = np.broadcast_to(rng.uniform(1.0, 2.0, (5, 4)), (3, 5, 4)).copy()
    # NO + O3 -> NO2 + O2 in place, mole for mole: NO2 gains what NO and O3 each lose
    moles = 0.1 * no / 0.030
    moles[0] = 0.0
    fields = {"NO": no - moles * 0.030, "O3": o3 - moles * 0.048, "NO2": moles * 0.046}
    nitrogen = [("NO", 0.030), ("NO2", 0.046)]
    odd_oxygen = [("NO2", 0.046), ("O3", 0.048)]
    assert max_drift(invariant_totals(mass, fields, nitrogen)) < 1e-14, "NO + O3 -> NO2 conserves NO + NO2"
    assert max_drift(invariant_totals(mass, fields, odd_oxygen)) < 1e-14, "NO + O3 -> NO2 conserves NO2 + O3"
    fields["NO"] = fields["NO"] * 1.001
    assert max_drift(invariant_totals(mass, fields, nitrogen)) > 1e-6, "a source is caught"
    q = np.full((2, 3, 4), 1.0)
    assert floor_violations(q) == 0
    q[1, 0, 0] = -2.0e-12
    assert floor_violations(q) == 1, "below -(1e-20 + 1e-12*max) is a violation"
    q[1, 0, 0] = -5.0e-13
    assert floor_violations(q) == 0, "roundoff below zero is allowed"
    q[1, 1, 1] = np.nan
    assert floor_violations(q) == 1, "non-finite is a violation"
    assert parse_invariant("Ox=NO2:0.046,O3:0.048") == ("Ox", [("NO2", 0.046), ("O3", 0.048)])
    print("check_pbl_reactive self-test: PASS")


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if argv == ["--self-test"]:
        self_test()
        return 0
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="command", required=True)
    p = sub.add_parser("invariants")
    p.add_argument("run")
    p.add_argument("--invariant", action="append", default=[])
    p.add_argument("--max-drift", type=float, default=1.0e-12)
    p.add_argument("--nonnegative", nargs="+", default=[])
    p.add_argument("--json")
    p.set_defaults(func=cmd_invariants)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

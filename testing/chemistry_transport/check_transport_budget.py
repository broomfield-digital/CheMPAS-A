#!/usr/bin/env python3
"""Check the per-step MICM gas transport budget that the model logs with
config_micm_transport_budget = .true. (docs/plans/pbl-transport.md, P5).

Rank 0 writes one line per process, gas and model step, with the domain totals
reduced over ranks:

  [transport] budget PBL t=<s> gas=<n> scalar=<idx> before=<kg> after=<kg>
      increment=<kg> rel=<increment/before> residual=<largest column relative
      residual> qmin=<kg/kg> qmax=<kg/kg> columns=<active columns> refresh=<s>
      dt=<s>

t is the model time since the start of the run, refresh the model time of the
PBL call that built the operator, dt the time step the tendency used.

Checks on the PBL lines:
  - |rel| <= --max-rel and residual <= --max-residual (closed columns conserve);
  - with --expect-gases, every step has exactly that many gases, each once;
  - with --expect-steps, that many distinct steps appear;
  - with --dt, every line used that time step;
  - with --interval, the operator was refreshed less than one interval ago
    (0 <= t - refresh < interval), at a model step time (with --dt), and no
    more often than the interval allows (the pbl alarm rings at the first step
    at or after each interval, so the interval need not be a multiple of dt);
  - with --min-q, qmin >= that floor.

Usage: check_transport_budget.py LOG [--expect-gases N] [--expect-steps N]
           [--dt S] [--interval S] [--max-rel 1e-12] [--max-residual 1e-12]
           [--min-q Q] [--json OUT]
"""
import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

LINE = re.compile(
    r"\[transport\] budget (?P<process>PBL|CU) t=(?P<t>\S+) gas=(?P<gas>\d+) scalar=(?P<scalar>\d+) "
    r"before=(?P<before>\S+) after=(?P<after>\S+) increment=(?P<increment>\S+) rel=(?P<rel>\S+) "
    r"residual=(?P<residual>\S+) qmin=(?P<qmin>\S+) qmax=(?P<qmax>\S+) columns=(?P<columns>\d+) "
    r"refresh=(?P<refresh>\S+) dt=(?P<dt>\S+)")
TIME_TOL = 1.0e-6   # [s]; model times are exact multiples of the step


def parse(text):
    """Budget records (dicts) of every budget line in the log text, in order."""
    rows = []
    for m in LINE.finditer(text):
        row = {k: float(v) for k, v in m.groupdict().items() if k not in ("process", "gas", "scalar", "columns")}
        row.update(process=m["process"], gas=int(m["gas"]), scalar=int(m["scalar"]), columns=int(m["columns"]))
        rows.append(row)
    return rows


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("log", type=Path)
    ap.add_argument("--process", default="PBL", choices=("PBL", "CU"))
    ap.add_argument("--expect-gases", type=int)
    ap.add_argument("--expect-steps", type=int)
    ap.add_argument("--dt", type=float)
    ap.add_argument("--interval", type=float)
    ap.add_argument("--max-rel", type=float, default=1.0e-12)
    ap.add_argument("--max-residual", type=float, default=1.0e-12)
    ap.add_argument("--min-q", type=float)
    ap.add_argument("--json", type=Path)
    return ap


def check(rows, args):
    """(failures, summary) for the budget rows of args.process."""
    rows = [r for r in rows if r["process"] == args.process]
    failures = []
    if not rows:
        return [f"no {args.process} budget lines in the log"], {"lines": 0, "steps": 0}
    by_step = defaultdict(list)
    for r in rows:
        by_step[r["t"]].append(r["gas"])
        where = f"t={r['t']:g} gas={r['gas']}"
        if abs(r["rel"]) > args.max_rel:
            failures.append(f"{where}: rel {r['rel']:.3e} exceeds {args.max_rel:g}")
        if r["residual"] > args.max_residual:
            failures.append(f"{where}: column residual {r['residual']:.3e} exceeds {args.max_residual:g}")
        if args.dt is not None and abs(r["dt"] - args.dt) > TIME_TOL:
            failures.append(f"{where}: dt {r['dt']:g} is not {args.dt:g}")
        if args.min_q is not None and r["qmin"] < args.min_q:
            failures.append(f"{where}: qmin {r['qmin']:.3e} is below {args.min_q:g}")
        if args.interval is not None:
            if r["refresh"] < 0:
                failures.append(f"{where}: the operator was never refreshed")
            elif args.dt is not None:
                k = r["refresh"] / args.dt
                if abs(k - round(k)) * args.dt > TIME_TOL:
                    failures.append(f"{where}: refresh {r['refresh']:g} is not a model step time")
            age = r["t"] - r["refresh"]
            if age < -TIME_TOL or age > args.interval - TIME_TOL:
                failures.append(f"{where}: the operator from {r['refresh']:g} s is older than the "
                                f"{args.interval:g} s interval")
    if args.expect_gases is not None:
        for t, gases in sorted(by_step.items()):
            if sorted(gases) != list(range(1, args.expect_gases + 1)):
                failures.append(f"t={t:g}: gas lines {sorted(gases)}, expected 1..{args.expect_gases} once each")
    if args.interval is not None:
        # the pbl alarm rings at the first step at or after each interval, so a run spanning [t0, t1] has at
        # most floor((t1-t0)/interval)+1 refreshes; more means the PBL ran more often than its interval
        refreshes = {r["refresh"] for r in rows}
        times = [r["t"] for r in rows]
        allowed = int((max(times) - min(times) + TIME_TOL) // args.interval) + 1
        if len(refreshes) > allowed:
            failures.append(f"{len(refreshes)} operator refreshes in {min(times):g}-{max(times):g} s, at most "
                            f"{allowed} expected with a {args.interval:g} s interval")
    if args.expect_steps is not None and len(by_step) != args.expect_steps:
        failures.append(f"{len(by_step)} steps logged, expected {args.expect_steps}")
    summary = {
        "lines": len(rows), "steps": len(by_step),
        "max_abs_rel": max(abs(r["rel"]) for r in rows),
        "max_residual": max(r["residual"] for r in rows),
        "min_qmin": min(r["qmin"] for r in rows),
        "columns": sorted({r["columns"] for r in rows}),
        "dt": sorted({r["dt"] for r in rows}),
        "refresh_times": sorted({r["refresh"] for r in rows}),
    }
    return failures, summary


def main(argv=None):
    args = build_parser().parse_args(argv)
    failures, summary = check(parse(args.log.read_text()), args)
    summary["failures"] = failures
    if args.json:
        args.json.write_text(json.dumps(summary, indent=1) + "\n")
    for f in failures[:20]:
        print(f"FAIL {f}")
    if len(failures) > 20:
        print(f"FAIL ... {len(failures) - 20} more")
    print(f"{'PASS' if not failures else 'FAIL'}: transport budget, {summary['lines']} lines over "
          f"{summary['steps']} steps, max |rel| {summary.get('max_abs_rel', 0):.2e}, "
          f"max column residual {summary.get('max_residual', 0):.2e}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())

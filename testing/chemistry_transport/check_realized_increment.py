#!/usr/bin/env python3
"""Compare realized MICM tracer increments with the ones the transport requested.

P2 of docs/plans/pbl-transport.md (shared items 4 and 7; code document §7.3).
The run uses the PBL gas test double (model_doubles/mpas_atmphys_pbl_gas.F) on a
horizontally uniform resting state with an inert mechanism, so after each
model step the tracer change in every cell must equal what the double asked
for:

  source     q(t+dt) - q(t) = dt*S(n,k), S(n,k) = S0*n*(1 + 0.5*sin(k))
  diffusion  q(t+dt) = backward-Euler closed-column diffusion of q(t) with
             G = G0 at interior interfaces and m_d = zz*rho_zz(t)*dz, solved
             here as a dense numpy system, independently of the model kernel
  off        q(t+dt) = q(t)

Per-layer comparisons use the frozen parallel tolerances
(|realized - requested| <= 1e-18 + 1e-10*|requested|) and need a state at
rest. In every run the column tracer mass that MPAS conserves,
M = sum(rho_zz*q/rdzw), must change by exactly the requested amount:
dt*sum(S*rho_zz(t)/rdzw) for source, 0 for diffusion and off, within 1e-12 of
M. With --column-only (a moving state, where density changes during the step)
only this budget is checked. It fixes the dry-mass weight of the tendency to
rho_zz at the start of the step, whatever the advection mode.

With --fresh-against OFF_RUN (diffusion mode; photolysis active, state at
rest) only the first step is checked: the transport must act on the
post-chemistry state, so q(t1) = D(q_off(t1)), where q_off is the matching run
with transport off. A stale implementation that diffused the pre-chemistry
state would give q_off(t1) + D(q(t0)) - q(t0); its error is reported and must
exceed the tolerance, so the check discriminates.

Usage: check_realized_increment.py RUN_DIR --mode {source,diffusion,off} --init INIT_FILE
                                   [--column-only | --fresh-against OFF_RUN]
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_budget import PARALLEL_ATOL, PARALLEL_RTOL  # noqa: E402

S0 = 1.0e-13      # must match model_doubles/mpas_atmphys_pbl_gas.F
G0 = 0.5
MASS_TOL = 1e-12


def transport_order(log_path):
    """MICM species names in transport order, from the chemistry log."""
    names = re.findall(r"\[MUSICA\] transport: (\S+) -> scalars\(", Path(log_path).read_text())
    if not names:
        raise SystemExit(f"no transport map in {log_path}")
    return ["q" + name for name in names]


def seconds(xtime):
    text = b"".join(xtime).decode().strip()
    hms = text.split("_")[1].split(":")
    return int(text.split("_")[0].split("-")[2]) * 86400 + int(hms[0]) * 3600 + int(hms[1]) * 60 + float(hms[2])


def diffuse(q, m_d, dt):
    """Dense backward-Euler flux-form solve; q, m_d: (nCells, nz)."""
    ncells, nz = q.shape
    out = np.empty_like(q)
    k = np.arange(nz)
    for start in range(0, ncells, 2000):
        m = m_d[start:start + 2000]
        n = m.shape[0]
        A = np.zeros((n, nz, nz))
        g_below = np.where(k > 0, G0, 0.0)
        g_above = np.where(k < nz - 1, G0, 0.0)
        A[:, k, k] = m / dt + g_below + g_above
        A[:, k[1:], k[1:] - 1] = -G0
        A[:, k[:-1], k[:-1] + 1] = -G0
        out[start:start + n] = np.linalg.solve(A, (m / dt * q[start:start + n])[..., None])[..., 0]
    return out


def check_fresh(args):
    if args.mode != "diffusion":
        raise SystemExit("--fresh-against needs --mode diffusion")
    order = transport_order(args.run_dir / "log.atmosphere.0000.out")
    with Dataset(args.run_dir / "output.nc") as ds, Dataset(args.fresh_against / "output.nc") as off:
        times = np.array([seconds(x) for x in ds["xtime"][:2]])
        dt = times[1] - times[0]
        m_d = ds["zz"][:] * ds["rho_zz"][0] * np.diff(ds["zgrid"][:], axis=1)
        report = {"mode": "fresh", "against": str(args.fresh_against), "dt": dt, "species": {}}
        failures = []
        for name in order:
            q0, q1 = ds[name][0], ds[name][1]
            post = off[name][1]                        # chemistry applied, no transport
            chem = float(np.max(np.abs(post - off[name][0])))
            fresh = diffuse(post, m_d, dt) - post
            stale = diffuse(q0, m_d, dt) - q0
            tol = PARALLEL_ATOL + PARALLEL_RTOL * np.abs(fresh)
            error = float(np.max(np.abs((q1 - post) - fresh) - tol))
            stale_error = float(np.max(np.abs((q1 - post) - stale) - tol))
            report["species"][name] = {"chemistry_change": chem, "excess_fresh": error, "excess_if_stale": stale_error}
            if error > 0.0:
                failures.append(f"{name}: post-chemistry diffusion misses by {error:.3e} over tolerance")
            if chem > 0.0 and stale_error <= 0.0:
                failures.append(f"{name}: a stale (pre-chemistry) tendency would also pass; check not discriminating")
        if max(e["chemistry_change"] for e in report["species"].values()) <= 1e-14:
            failures.append("chemistry changed nothing; the freshness check needs active photolysis")
    report["result"] = "PASS" if not failures else "FAIL"
    report["failures"] = failures
    (args.run_dir / "realized_increment.json").write_text(json.dumps(report, indent=1))
    worst = max(e["excess_fresh"] for e in report["species"].values())
    stale = min(e["excess_if_stale"] for e in report["species"].values() if e["chemistry_change"] > 0)
    print(f"{report['result']} {args.run_dir.name}: first step uses the post-chemistry state "
          f"(excess over tolerance {worst:.2e}; a stale tendency would exceed it by {stale:.2e})")
    for line in failures:
        print("  " + line)
    return 0 if not failures else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--mode", choices=("source", "diffusion", "off"), required=True)
    ap.add_argument("--init", type=Path, required=True)
    ap.add_argument("--column-only", action="store_true",
                    help="check only the column tracer-mass budget (state not at rest)")
    ap.add_argument("--fresh-against", type=Path,
                    help="matching run with transport off: check the first step uses the post-chemistry state")
    args = ap.parse_args()
    if args.fresh_against:
        return check_fresh(args)

    order = transport_order(args.run_dir / "log.atmosphere.0000.out")
    with Dataset(args.init) as init:
        rdzw = init["rdzw"][:]
    with Dataset(args.run_dir / "output.nc") as ds:
        times = np.array([seconds(x) for x in ds["xtime"][:]])
        zz = ds["zz"][:]
        dz = np.diff(ds["zgrid"][:], axis=1)
        rho_zz = ds["rho_zz"][:]
        q = {name: ds[name][:] for name in order}
        w_max = float(np.max(np.abs(ds["w"][:])))
        u_max = float(np.max(np.abs(ds["u"][:])))

    nz = zz.shape[1]
    k = np.arange(1, nz + 1)
    failures = []
    report = {"mode": args.mode, "column_only": args.column_only, "species": order, "times_s": times.tolist(),
              "max_abs_w": w_max, "max_abs_u": u_max, "steps": []}
    for t in range(len(times) - 1):
        dt = times[t + 1] - times[t]
        m_d = zz * rho_zz[t] * dz
        step = {"t": times[t], "dt": dt, "species": {}}
        for n, name in enumerate(order, start=1):
            before, after = q[name][t], q[name][t + 1]
            source = np.broadcast_to(S0 * n * (1.0 + 0.5 * np.sin(k)), before.shape)
            mass_before = np.sum(rho_zz[t] * before / rdzw, axis=1)
            mass_after = np.sum(rho_zz[t + 1] * after / rdzw, axis=1)
            mass_requested = (np.sum(dt * source * rho_zz[t] / rdzw, axis=1) if args.mode == "source"
                              else np.zeros_like(mass_before))
            mass_error = float(np.max(np.abs(mass_after - mass_before - mass_requested) / mass_before))
            entry = {"column_mass_requested_relative": float(np.max(np.abs(mass_requested) / mass_before)),
                     "column_mass_error_relative": mass_error}
            if mass_error > MASS_TOL:
                failures.append(f"t={times[t]:g}s {name}: column mass budget misses by {mass_error:.3e}")
            if not args.column_only:
                if args.mode == "source":
                    requested = dt * source
                elif args.mode == "diffusion":
                    requested = diffuse(before, m_d, dt) - before
                else:
                    requested = np.zeros_like(before)
                realized = after - before
                excess = np.abs(realized - requested) - (PARALLEL_ATOL + PARALLEL_RTOL * np.abs(requested))
                worst = float(np.max(excess))
                entry.update({"max_abs_requested": float(np.max(np.abs(requested))),
                              "max_abs_error": float(np.max(np.abs(realized - requested))),
                              "max_excess_over_tolerance": worst})
                if worst > 0.0:
                    failures.append(f"t={times[t]:g}s {name}: realized - requested exceeds tolerance by {worst:.3e}")
            step["species"][name] = entry
        report["steps"].append(step)

    report["result"] = "PASS" if not failures else "FAIL"
    report["failures"] = failures
    (args.run_dir / "realized_increment.json").write_text(json.dumps(report, indent=1))
    entries = [e for s in report["steps"] for e in s["species"].values()]
    worst_mass = max(e["column_mass_error_relative"] for e in entries)
    layer = ("column budget only" if args.column_only else
             f"max|realized-requested| {max(e['max_abs_error'] for e in entries):.3e}")
    print(f"{report['result']} {args.run_dir.name}: {len(report['steps'])} steps, {layer}, "
          f"column-mass budget error {worst_mass:.3e}, max|w| {w_max:.2e}, max|u| {u_max:.2e}")
    for line in failures[:6]:
        print("  " + line)
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())

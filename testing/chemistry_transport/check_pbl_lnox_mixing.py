#!/usr/bin/env python3
"""Relative size of the MICM gas PBL transport differences in the supercell LNOx case.

Compares the transport-on and transport-off runs of run_supercell_lnox_pbl.sh as
ON - OFF, per output snapshot, in percent of the OFF values:

  nox_moved      half the dry-mass-weighted sum of |dNOx| over the domain NOx:
                 the share of the NOx that transport put in another cell
  plume          |dNOx|/NOx and |dO3|/O3 in cells with NOx >= 10 ppbv
  pbl            |dNOx|/NOx (NOx >= 0.1 ppbv) and |dO3|/O3 below the PBL height
  o3_hole        |dO3| over the O3 deficit, 50 ppbv (the initial O3) - O3,
                 where the deficit is >= 1 ppbv
  height bands   each band's share of the summed |dNOx|, and the plume's
                 NOx-weighted mean height and spread
Distributions are given as median, 90th and 99th percentiles and maximum.

With --kzh-run (a run_supercell_kzh_restart.sh directory) it also checks the size
against YSU's own diffusivity. It reports kzh on the bottom faces of the plume
cells and the layer mixing time dz^2/kzh, and diffuses the OFF NOx of the
restart output offline (backward Euler, dry-mass flux form) for 3 and 10
minutes. The predicted increments are compared with the accumulated ON - OFF
difference at the nearest output time. Above the PBL top kzh is the coefficient
the gases see (bl_ysu.F90: xkzq = xkzh); inside the PBL they use YSU's xkzq, so
the prediction is for the free atmosphere, where the differences are.

These are consistency checks against the scheme's own coefficients, not a
physical validation of YSU. Mixing ratios are converted to ppbv as in
plot_pbl_transport_diff.py. Nothing is smoothed.

Usage: check_pbl_lnox_mixing.py ON_RUN OFF_RUN [--kzh-run DIR] [--minutes 30 60 90 120] [--json OUT]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, chartostring

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_pbl_transport_diff import M_AIR, M_NO, M_NO2, NLEV_KEEP, load, minutes  # noqa: E402

O3_INITIAL = 50.0      # ppbv, the standard case's initial O3
PLUME_NOX = 10.0       # ppbv
PBL_NOX = 0.1          # ppbv
BANDS = ((0, 1), (1, 3), (3, 5), (5, 8), (8, 12), (12, 20))   # km


def pct_stats(ratio):
    """Median, 90th and 99th percentiles and maximum of a ratio, in percent."""
    if ratio.size == 0:
        return {"n": 0}
    return {"n": int(ratio.size), "median": 100 * float(np.median(ratio)),
            "p90": 100 * float(np.percentile(ratio, 90)), "p99": 100 * float(np.percentile(ratio, 99)),
            "max": 100 * float(ratio.max())}


def snapshot(on, off, it):
    zc, m = off["zc"], off["f"]["mass"][it]
    nox_on, nox_off = on["f"]["nox"][it], off["f"]["nox"][it]
    o3_on, o3_off = on["f"]["o3"][it], off["f"]["o3"][it]
    dn, do = np.abs(nox_on - nox_off), np.abs(o3_on - o3_off)
    plume = nox_off >= PLUME_NOX
    below = zc < on["f"]["hpbl"][it][:, None] / 1000.0
    deficit = O3_INITIAL - o3_off
    hole = deficit >= 1.0

    def moments(q):
        w = m * q
        zbar = (w * zc).sum() / w.sum()
        return float(zbar), float(np.sqrt((w * (zc - zbar) ** 2).sum() / w.sum()))

    (z_on, s_on), (z_off, s_off) = moments(nox_on), moments(nox_off)
    return {
        "nox_moved_pct": 100 * 0.5 * float((m * dn).sum() / (m * nox_off).sum()),
        "plume_nox_pct": pct_stats(dn[plume] / nox_off[plume]),
        "plume_o3_pct": pct_stats(do[plume] / o3_off[plume]),
        "pbl_nox_pct": pct_stats(dn[below & (nox_off >= PBL_NOX)] / nox_off[below & (nox_off >= PBL_NOX)]),
        "pbl_o3_pct": pct_stats(do[below] / o3_off[below]),
        "o3_hole_pct": pct_stats(do[hole] / deficit[hole]),
        "max_abs_diff_ppbv": {"NOx": float(dn.max()), "O3": float(do.max())},
        "abs_dnox_share_by_height_km": {f"{lo}-{hi}": float(dn[(zc >= lo) & (zc < hi)].sum() / dn.sum())
                                        for lo, hi in BANDS},
        "plume_mean_height_km": {"off": z_off, "on": z_on},
        "plume_spread_m": {"off": 1000 * s_off, "on": 1000 * s_on},
    }


def diffuse(q, rho, dz, zc, kface, dt):
    """Backward-Euler vertical diffusion of q in each column. Layer dry mass rho*dz; on each
    interior face the conductance rho_face*K/(zc_k - zc_{k-1}), with kface[:, k] on the bottom
    face of layer k; no flux through the bottom and top."""
    ncell, nlev = q.shape
    g = np.zeros((ncell, nlev + 1))
    g[:, 1:nlev] = 0.5 * (rho[:, 1:] + rho[:, :-1]) * kface[:, 1:] / (zc[:, 1:] - zc[:, :-1])
    mass = rho * dz
    a, c = -dt * g[:, :nlev], -dt * g[:, 1:]
    b, r = mass - a - c, mass * q
    cp, dp = np.zeros_like(q), np.zeros_like(q)
    cp[:, 0], dp[:, 0] = c[:, 0] / b[:, 0], r[:, 0] / b[:, 0]
    for k in range(1, nlev):
        den = b[:, k] - a[:, k] * cp[:, k - 1]
        cp[:, k] = c[:, k] / den
        dp[:, k] = (r[:, k] - a[:, k] * dp[:, k - 1]) / den
    x = np.zeros_like(q)
    x[:, -1] = dp[:, -1]
    for k in range(nlev - 2, -1, -1):
        x[:, k] = dp[:, k] - cp[:, k] * x[:, k + 1]
    return x


def kzh_check(kzh_run, on, off, tmin):
    with Dataset(kzh_run / "output_kzh.nc") as d:
        time = str(chartostring(d["xtime"][-1])).strip()
        if not np.array_equal(np.asarray(d["indexToCellID"][:]), off["ids"]):
            raise SystemExit(f"{kzh_run}: cell order differs from the runs")

        def get(name):
            return np.asarray(d[name][-1][:, :NLEV_KEEP], dtype=np.float64)

        kface = get("kzh")
        theta = get("theta")
        nox = (get("qNO") * M_AIR / M_NO + get("qNO2") * M_AIR / M_NO2) * 1.0e9
        rho = get("rho_zz") * np.asarray(d["zz"][:, :NLEV_KEEP], dtype=np.float64)
        zgrid = np.asarray(d["zgrid"][:, :NLEV_KEEP + 1], dtype=np.float64)
        area = np.asarray(d["areaCell"][:], dtype=np.float64)[:, None]
    dz = np.diff(zgrid, axis=1)
    zc = 0.5 * (zgrid[:, 1:] + zgrid[:, :-1])
    # the run nearest the restart output, which the restart must continue
    t_kzh = minutes([off["times"][0], time])[1]
    it = int(np.argmin(np.abs(tmin - t_kzh)))
    corr = float(np.corrcoef(theta.ravel(), off["f"]["theta"][it].ravel())[0, 1])
    plume = nox >= PLUME_NOX
    k_plume = kface[plume]
    tau = (zc[:, 1:] - zc[:, :-1]) ** 2 / np.maximum(kface[:, 1:], 1.0e-12)
    tau_plume = tau[plume[:, 1:]]
    out = {"kzh_time": time, "kzh_minutes": float(t_kzh), "compared_record_minutes": float(tmin[it]),
           "theta_correlation_with_off_run": corr,
           "plume_kzh_m2_s": {"median": float(np.median(k_plume)), "p90": float(np.percentile(k_plume, 90)),
                              "p99": float(np.percentile(k_plume, 99)), "max": float(k_plume.max())},
           "plume_layer_mixing_time_min": {"median": float(np.median(tau_plume) / 60),
                                           "p1": float(np.percentile(tau_plume, 1) / 60)},
           "plume_faces_mixing_within_120_min_pct": 100 * float((tau_plume < 7200.0).mean()),
           "predicted": {}}
    m = rho * dz * area
    for dt in (180.0, 600.0):
        dq = np.abs(diffuse(nox, rho, dz, zc, kface, dt) - nox)
        out["predicted"][f"{dt / 60:.0f} min"] = {
            "nox_moved_pct": 100 * 0.5 * float((m * dq).sum() / (m * nox).sum()),
            "plume_nox_pct": pct_stats(dq[plume] / nox[plume]),
            "max_abs_diff_ppbv": float(dq.max())}
    out["actual"] = {key: value for key, value in snapshot(on, off, it).items()
                     if key in ("nox_moved_pct", "plume_nox_pct", "max_abs_diff_ppbv")}
    return out


def fmt(s):
    if s.get("n", 0) == 0:
        return "no cells"
    return f"median {s['median']:.3g}%  p90 {s['p90']:.3g}%  p99 {s['p99']:.3g}%  max {s['max']:.3g}%  (n={s['n']})"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("on_run", type=Path)
    ap.add_argument("off_run", type=Path)
    ap.add_argument("--kzh-run", type=Path)
    ap.add_argument("--minutes", type=float, nargs="+", default=[30.0, 60.0, 90.0, 120.0])
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()

    on, off = load(args.on_run), load(args.off_run)
    if on["times"] != off["times"] or not np.array_equal(on["ids"], off["ids"]):
        raise SystemExit("the two runs differ in output times or cells")
    if on["f"]["hpbl"] is None:
        raise SystemExit("hpbl is not in the output; run_supercell_lnox_pbl.sh adds it")
    tmin = minutes(on["times"])
    report = {"on_run": str(args.on_run.resolve()), "off_run": str(args.off_run.resolve()),
              "theta_bitwise_identical": bool(np.array_equal(on["f"]["theta"], off["f"]["theta"])),
              "snapshots": {}}
    for want in args.minutes:
        it = int(np.argmin(np.abs(tmin - want)))
        s = snapshot(on, off, it)
        report["snapshots"][f"{tmin[it]:.0f} min"] = s
        print(f"--- {tmin[it]:.0f} min")
        print(f"  NOx moved to another cell: {s['nox_moved_pct']:.3g}% of the NOx")
        print(f"  plume |dNOx|/NOx: {fmt(s['plume_nox_pct'])}")
        print(f"  plume |dO3|/O3:   {fmt(s['plume_o3_pct'])}")
        print(f"  PBL   |dNOx|/NOx: {fmt(s['pbl_nox_pct'])}")
        print(f"  PBL   |dO3|/O3:   {fmt(s['pbl_o3_pct'])}")
        print(f"  O3 hole |dO3|/deficit: {fmt(s['o3_hole_pct'])}")
        bands = ", ".join(f"{k} {100 * v:.1f}%" for k, v in s["abs_dnox_share_by_height_km"].items())
        print(f"  share of |dNOx| by height (km): {bands}")
        print(f"  plume spread {s['plume_spread_m']['off']:.0f} -> {s['plume_spread_m']['on']:.0f} m")
    if args.kzh_run:
        k = kzh_check(args.kzh_run, on, off, tmin)
        report["kzh_check"] = k
        print(f"--- YSU diffusivity at {k['kzh_minutes']:.1f} min "
              f"(theta correlation with the off run at {k['compared_record_minutes']:.0f} min: "
              f"{k['theta_correlation_with_off_run']:.6f})")
        pk = k["plume_kzh_m2_s"]
        print(f"  plume kzh: median {pk['median']:.3g}  p90 {pk['p90']:.3g}  p99 {pk['p99']:.3g}  "
              f"max {pk['max']:.3g} m2/s; faces mixing a layer within 120 min: "
              f"{k['plume_faces_mixing_within_120_min_pct']:.2g}%")
        for label, p in k["predicted"].items():
            print(f"  predicted over {label}: moved {p['nox_moved_pct']:.3g}%; plume {fmt(p['plume_nox_pct'])}")
        a = k["actual"]
        print(f"  actual at {k['compared_record_minutes']:.0f} min: moved {a['nox_moved_pct']:.3g}%; "
              f"plume {fmt(a['plume_nox_pct'])}")
    if args.json:
        args.json.write_text(json.dumps(report, indent=1) + "\n")


if __name__ == "__main__":
    main()

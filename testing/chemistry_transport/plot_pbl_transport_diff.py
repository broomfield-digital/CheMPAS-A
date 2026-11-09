#!/usr/bin/env python3
"""Difference diagnostics for MICM gas PBL transport in the supercell LNOx case.

Compares two runs of the standard LNOx case that differ only in
config_micm_pbl_transport (run_supercell_lnox_pbl.sh), as ON - OFF:

  pbl_diff_lowest_layer  per snapshot, lowest model layer: NOx with transport,
                         and the NOx and O3 differences
  pbl_diff_xsection      per snapshot, x-z section through the storm (lowest
                         ZTOP km): the NOx and O3 differences, with the PBL
                         height of the transport run
  pbl_diff_xsection_full the same sections up to FULL_ZTOP km, with the 10 and
                         100 ppbv NOx contours of the transport run outlining
                         the lightning plume. Above the PBL top YSU still
                         mixes, with its local Richardson-number diffusivity,
                         and the gases now share it, so the largest
                         differences lie in the storm aloft
  pbl_diff_profiles      domain-mean NOx and O3 differences against height:
                         the whole column, and the lowest 3 km
  pbl_diff_timeseries    mass-weighted mean NOx and O3 below 1 km, with and
                         without transport, and their difference

Mixing ratios are converted to ppbv with the constants of the established
LNOx plots (scripts/plot_lnox_o3.py: M_air 0.029, M_NO 0.030, M_NO2 0.046,
M_O3 0.048 kg/mol). Nothing is smoothed or clipped. The difference maps and
sections use symmetric logarithmic contour levels spanning three decades
below the largest |difference| among the plotted snapshots (one scale per
species and figure); smaller differences are left white. A JSON summary
records the largest differences, where they lie by height, and checks that
the two runs carry the same total nitrogen and the same meteorology.

Usage: plot_pbl_transport_diff.py ON_RUN OFF_RUN --outdir OUT [--minutes 30 60 90 120]
           [--ztop-km 4] [--full-ztop-km 16] [--label YSU]
Each run directory holds output.nc and supercell_init.nc (for xCell, yCell).
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import BoundaryNorm  # noqa: E402
from matplotlib.tri import Triangulation  # noqa: E402
from netCDF4 import Dataset, chartostring  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import style  # noqa: E402

M_AIR, M_NO, M_NO2, M_O3 = 0.029, 0.030, 0.046, 0.048
NLEV_KEEP = 30          # cell centres below about 20 km on this mesh


def ppbv(q, m_species):
    return q * (M_AIR / m_species) * 1.0e9


def load(run):
    """Fields of one run, levels below about 20 km."""
    run = Path(run)
    with Dataset(run / "output.nc") as d:
        times = [str(t).strip() for t in chartostring(d["xtime"][:])]
        f = {name: np.asarray(d[name][:, :, :NLEV_KEEP], dtype=np.float64)
             for name in ("qNO", "qNO2", "qO3", "rho_zz", "theta")}
        f["w"] = np.asarray(d["w"][:, :, :NLEV_KEEP + 1], dtype=np.float64)
        f["hpbl"] = np.asarray(d["hpbl"][:], dtype=np.float64) if "hpbl" in d.variables else None
        zgrid = np.asarray(d["zgrid"][:, :NLEV_KEEP + 1], dtype=np.float64)
        zz = np.asarray(d["zz"][:, :NLEV_KEEP], dtype=np.float64)
        area = np.asarray(d["areaCell"][:], dtype=np.float64)
        ids = np.asarray(d["indexToCellID"][:])
    with Dataset(run / "supercell_init.nc") as d:
        if not np.array_equal(np.asarray(d["indexToCellID"][:]), ids):
            raise SystemExit(f"{run}: supercell_init.nc cell order differs from output.nc")
        x = np.asarray(d["xCell"][:], dtype=np.float64) / 1000.0
        y = np.asarray(d["yCell"][:], dtype=np.float64) / 1000.0
    dz = np.diff(zgrid, axis=1)
    f["mass"] = f["rho_zz"] * zz[None] * dz[None] * area[None, :, None]   # dry air [kg]
    f["nox"] = ppbv(f["qNO"], M_NO) + ppbv(f["qNO2"], M_NO2)
    f["o3"] = ppbv(f["qO3"], M_O3)
    f["nitrogen"] = (f["qNO"] / M_NO + f["qNO2"] / M_NO2) * f["mass"]  # mol N per cell layer
    zc = 0.5 * (zgrid[:, :-1] + zgrid[:, 1:]) / 1000.0
    return dict(times=times, f=f, x=x, y=y, zc=zc, zface=zgrid / 1000.0, ids=ids)


def minutes(times):
    def sec(t):
        hh, mm, ss = t.split("_")[1].split(":")
        return int(hh) * 3600 + int(mm) * 60 + float(ss)
    t0 = sec(times[0])
    return np.array([(sec(t) - t0) / 60.0 for t in times])


def symlog_levels(vmax, decades=3):
    """Symmetric levels: +-vmax and three decades of round values below it."""
    if not np.isfinite(vmax) or vmax <= 0.0:
        return np.array([-1.0, 1.0])
    top = 10.0 ** np.ceil(np.log10(vmax))
    pos = []
    for e in range(int(np.log10(top)) - decades, int(np.log10(top)) + 1):
        for m in (1.0, 2.0, 5.0):
            v = m * 10.0 ** e
            if v <= top * (1 + 1e-9):
                pos.append(v)
    pos = np.array(sorted(set(pos)))
    pos = pos[pos >= top / 10.0 ** decades]
    return np.concatenate([-pos[::-1], pos])


def nice_levels(vmax, target=10):
    """Levels 0..vmax at round steps (1, 2 or 5 times a power of ten)."""
    if not np.isfinite(vmax) or vmax <= 0.0:
        return np.linspace(0.0, 1.0, 11)
    raw = vmax / target
    step = min((m * 10.0 ** np.floor(np.log10(raw)) for m in (1.0, 2.0, 5.0, 10.0)
                if m * 10.0 ** np.floor(np.log10(raw)) >= raw))
    return np.arange(0.0, vmax + step, step)


def decade_ticks(levels):
    """Colour-bar ticks at +-10^k among the levels."""
    mant = np.abs(levels) / 10.0 ** np.floor(np.log10(np.abs(levels)))
    return levels[np.isclose(mant, 1.0)]


def diff_contour(ax, xx, zz_or_tri, field, levels, tri=False):
    cmap = plt.get_cmap(style.get_bias_cmap())
    norm = BoundaryNorm(levels, cmap.N, extend="neither")
    # values inside the smallest level band (|d| < levels[n/2]) show as white
    if tri:
        cs = ax.tricontourf(zz_or_tri, field, levels=levels, cmap=cmap, norm=norm)
    else:
        cs = ax.contourf(xx, zz_or_tri, field, levels=levels, cmap=cmap, norm=norm)
    return cs


def storm_row(run, it):
    """Row of cells through the strongest updraft near 5 km (face-averaged w)."""
    f, x, y = run["f"], run["x"], run["y"]
    zc_mean = run["zc"].mean(axis=0)
    k5 = int(np.argmin(np.abs(zc_mean - 5.0)))
    wmid = 0.5 * (f["w"][it, :, k5] + f["w"][it, :, k5 + 1])
    y0 = y[int(np.argmax(wmid))]
    dx = (x.max() - x.min()) / 200
    bins = np.digitize(x, np.arange(x.min(), x.max() + dx, dx))
    sel = []
    for b in np.unique(bins):
        cells = np.where(bins == b)[0]
        sel.append(cells[np.argmin(np.abs(y[cells] - y0))])
    sel = np.array(sel)
    return sel[np.argsort(x[sel])], y0


def xsection_figure(on, d_nox, d_o3, snaps, rows, tmin, ztop, outstem, title, subtitle, nox_lines=()):
    """x-z sections of the NOx and O3 differences through each snapshot's storm row, up to ZTOP km.
    NOX_LINES: NOx levels (ppbv) of the transport run drawn as contours."""
    kz = int(np.searchsorted(on["zc"].mean(axis=0), ztop)) + 1
    lv_nox = symlog_levels(max(np.abs(d_nox[i][rows[r][0], :kz]).max() for r, i in enumerate(snaps)))
    lv_o3 = symlog_levels(max(np.abs(d_o3[i][rows[r][0], :kz]).max() for r, i in enumerate(snaps)))
    fig, axes = plt.subplots(len(snaps), 2, figsize=(14, 3.6 * len(snaps)), squeeze=False)
    for r, i in enumerate(snaps):
        sel, y0 = rows[r]
        xs = on["x"][sel]
        zs = on["zc"][sel, :kz]
        X = np.repeat(xs[:, None], kz, axis=1)
        for c, (dd, lv, name) in enumerate(((d_nox, lv_nox, "NOx"), (d_o3, lv_o3, "O3"))):
            ax = axes[r, c]
            cs = diff_contour(ax, X, zs, dd[i][sel, :kz], lv)
            cb = fig.colorbar(cs, ax=ax, label="ppbv", ticks=decade_ticks(lv), format="%g")
            cb.ax.tick_params(labelsize=7)
            handles = []
            if nox_lines:
                styles = ["--", "-", ":"][:len(nox_lines)]
                cl = ax.contour(X, zs, on["f"]["nox"][i][sel, :kz], levels=list(nox_lines),
                                colors="0.35", linewidths=0.7, linestyles=styles)
                ax.clabel(cl, fmt="%g", fontsize=6)
                handles += [plt.Line2D([], [], color="0.35", lw=0.7, ls=ls, label=f"NOx = {v:g} ppbv")
                            for v, ls in zip(nox_lines, styles)]
            if on["f"]["hpbl"] is not None:
                handles += ax.plot(xs, on["f"]["hpbl"][i, sel] / 1000.0, color="0.15", lw=1.2,
                                   label="PBL height (on)")
            if handles:
                ax.legend(handles=handles, loc="upper right", fontsize=8, frameon=False)
            ax.set_ylim(0, ztop)
            ax.set_xlabel("x (km)")
            ax.set_ylabel("z (km)")
            ax.set_title(style.panel_heading("abcdefghijklmnopqrstuvwxyz"[2 * r + c],
                                             f"Δ{name} at y = {y0:.0f} km, {tmin[i]:.0f} min"))
    style.apply_figure_heading(fig, title, subtitle)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    for ext in ("png", "pdf"):
        fig.savefig(f"{outstem}.{ext}")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("on_run")
    ap.add_argument("off_run")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--minutes", type=float, nargs="+", default=[30.0, 60.0, 90.0, 120.0])
    ap.add_argument("--ztop-km", type=float, default=4.0, help="top of the low-level cross-sections")
    ap.add_argument("--full-ztop-km", type=float, default=16.0, help="top of the full-height cross-sections")
    ap.add_argument("--label", default="YSU", help="PBL scheme named in the titles")
    args = ap.parse_args()

    on, off = load(args.on_run), load(args.off_run)
    if on["times"] != off["times"] or not np.array_equal(on["ids"], off["ids"]):
        raise SystemExit("the two runs differ in output times or cells")
    tmin = minutes(on["times"])
    snaps = [int(np.argmin(np.abs(tmin - m))) for m in args.minutes]
    args.outdir.mkdir(parents=True, exist_ok=True)
    style.setup()
    d_nox = on["f"]["nox"] - off["f"]["nox"]
    d_o3 = on["f"]["o3"] - off["f"]["o3"]
    summary = {"on_run": str(Path(args.on_run).resolve()), "off_run": str(Path(args.off_run).resolve()),
               "times": on["times"], "snapshots_min": [float(tmin[i]) for i in snaps]}

    # same meteorology and the same total nitrogen in both runs
    summary["theta_bitwise_identical"] = bool(np.array_equal(on["f"]["theta"], off["f"]["theta"]))
    n_on = on["f"]["nitrogen"].sum(axis=(1, 2))
    n_off = off["f"]["nitrogen"].sum(axis=(1, 2))
    summary["total_nitrogen_rel_diff_max"] = float(np.max(np.abs(n_on - n_off) / np.maximum(n_off, 1e-300)))
    summary["max_abs_diff_ppbv"] = {
        f"{tmin[i]:.0f} min": {"NOx": float(np.abs(d_nox[i]).max()), "O3": float(np.abs(d_o3[i]).max())}
        for i in snaps}
    # where the differences lie: height of the largest |dNOx|, and each band's share of the summed |dNOx|
    bands = ((0, 1), (1, 3), (3, 5), (5, 8), (8, 12), (12, 20))
    summary["max_abs_diff_nox_height_km"] = {}
    summary["abs_diff_nox_share_by_height_km"] = {}
    for i in snaps:
        a = np.abs(d_nox[i])
        cell, k = np.unravel_index(int(np.argmax(a)), a.shape)
        summary["max_abs_diff_nox_height_km"][f"{tmin[i]:.0f} min"] = float(on["zc"][cell, k])
        summary["abs_diff_nox_share_by_height_km"][f"{tmin[i]:.0f} min"] = {
            f"{lo}-{hi}": float(a[(on["zc"] >= lo) & (on["zc"] < hi)].sum() / a.sum()) for lo, hi in bands}

    # (1) lowest layer maps
    tri = Triangulation(on["x"], on["y"])
    lv_nox = symlog_levels(max(np.abs(d_nox[i, :, 0]).max() for i in snaps))
    lv_o3 = symlog_levels(max(np.abs(d_o3[i, :, 0]).max() for i in snaps))
    nox_max = max(on["f"]["nox"][i, :, 0].max() for i in snaps)
    fig, axes = plt.subplots(len(snaps), 3, figsize=(15, 5.0 * len(snaps)), squeeze=False)
    for r, i in enumerate(snaps):
        ax = axes[r, 0]
        cs = ax.tricontourf(tri, on["f"]["nox"][i, :, 0], levels=nice_levels(nox_max),
                            cmap=style.get_sequential_cmap())
        fig.colorbar(cs, ax=ax, label="ppbv")
        ax.set_title(style.panel_heading("abcdefghijklmnopqrstuvwxyz"[3 * r], f"NOx with transport, {tmin[i]:.0f} min"))
        for c, (dd, lv, name) in enumerate(((d_nox, lv_nox, "NOx"), (d_o3, lv_o3, "O3")), start=1):
            ax = axes[r, c]
            cs = diff_contour(ax, None, tri, dd[i, :, 0], lv, tri=True)
            cb = fig.colorbar(cs, ax=ax, label="ppbv", ticks=decade_ticks(lv), format="%g")
            cb.ax.tick_params(labelsize=7)
            ax.set_title(style.panel_heading("abcdefghijklmnopqrstuvwxyz"[3 * r + c],
                                             f"Δ{name} (on − off), {tmin[i]:.0f} min"))
        for ax in axes[r]:
            ax.set_xlabel("x (km)")
            ax.set_ylabel("y (km)")
            ax.set_aspect("equal")
    z0 = float(on["zc"][:, 0].mean() * 1000.0)
    style.apply_figure_heading(fig, f"Gas PBL Transport, Lowest Layer ({z0:.0f} m)",
                               f"Supercell LNOx, {args.label}: transport on − off; symmetric log levels, 3 decades")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.subplots_adjust(hspace=0.35)
    for ext in ("png", "pdf"):
        fig.savefig(args.outdir / f"pbl_diff_lowest_layer.{ext}")
    plt.close(fig)

    # (2) cross-sections through the storm
    rows = [storm_row(on, i) for i in snaps]
    xsection_figure(on, d_nox, d_o3, snaps, rows, tmin, args.ztop_km, args.outdir / "pbl_diff_xsection",
                    "Gas PBL Transport Through the Storm",
                    f"Supercell LNOx, {args.label}: transport on − off; section through the 5 km updraft maximum")
    xsection_figure(on, d_nox, d_o3, snaps, rows, tmin, args.full_ztop_km, args.outdir / "pbl_diff_xsection_full",
                    "Gas PBL Transport Through the Storm, Full Height",
                    f"Supercell LNOx, {args.label}: colours, transport on − off; lines, NOx = 10 and 100 ppbv "
                    "outlining the lightning plume",
                    nox_lines=(10.0, 100.0))

    # (3) domain-mean difference profiles (mass-weighted per level): whole column, and the lowest 3 km
    zmean = on["zc"].mean(axis=0)
    m = on["f"]["mass"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 11))
    colors = plt.get_cmap("viridis")(np.linspace(0.1, 0.85, len(snaps)))
    for row, ztop in enumerate((args.full_ztop_km, 3.0)):
        kp = int(np.searchsorted(zmean, ztop)) + 1
        for col_ax, (dd, name) in enumerate(((d_nox, "NOx"), (d_o3, "O3"))):
            ax = axes[row, col_ax]
            for col, i in zip(colors, snaps):
                prof = (dd[i] * m[i]).sum(axis=0) / m[i].sum(axis=0)
                ax.plot(prof[:kp], zmean[:kp], color=col, lw=1.6, marker="o", ms=3, label=f"{tmin[i]:.0f} min")
            ax.axvline(0.0, color="0.6", lw=0.8)
            ax.set_ylim(0, ztop)
            ax.locator_params(axis="x", nbins=5)
            ax.ticklabel_format(axis="x", style="sci", scilimits=(-2, 3))
            ax.set_xlabel(f"domain-mean Δ{name} (ppbv)")
            ax.set_ylabel("z (km)")
            ax.legend(frameon=False, fontsize=9)
            where = "whole column" if row == 0 else "lowest 3 km"
            ax.set_title(style.panel_heading("abcd"[2 * row + col_ax], f"Δ{name} (on − off), {where}"))
    style.apply_figure_heading(fig, "Gas PBL Transport, Domain-Mean Profiles",
                               f"Supercell LNOx, {args.label}: transport on − off, mass-weighted per level (markers: layer centres)")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    for ext in ("png", "pdf"):
        fig.savefig(args.outdir / f"pbl_diff_profiles.{ext}")
    plt.close(fig)

    # (4) mass-weighted means below 1 km against time
    below = on["zc"] < 1.0
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.8))
    series = {}
    for run, tag in ((on, "on"), (off, "off")):
        mm = run["f"]["mass"] * below[None]
        series[tag] = {name: (run["f"][name] * mm).sum(axis=(1, 2)) / mm.sum(axis=(1, 2)) for name in ("nox", "o3")}
    for ax, name, label in ((axes[0], "nox", "NOx"), (axes[1], "o3", "O3")):
        ax.plot(tmin, series["on"][name], color=style.get_palette(2)[0], lw=1.8, label="transport on")
        ax.plot(tmin, series["off"][name], color=style.get_palette(2)[1], lw=1.8, ls="--", label="transport off")
        ax.set_xlabel("time (min)")
        ax.set_ylabel(f"{label} below 1 km (ppbv)")
        ax.ticklabel_format(axis="y", useOffset=False)
        ax.legend(frameon=False)
    for name, label, col in (("nox", "NOx", style.get_palette(2)[0]), ("o3", "O3", style.get_palette(2)[1])):
        axes[2].plot(tmin, series["on"][name] - series["off"][name], color=col, lw=1.8, label=f"Δ{label}")
    axes[2].axhline(0.0, color="0.6", lw=0.8)
    axes[2].set_xlabel("time (min)")
    axes[2].set_ylabel("on − off below 1 km (ppbv)")
    axes[2].legend(frameon=False)
    for k, (ax, t) in enumerate(zip(axes, ("NOx below 1 km", "O3 below 1 km", "Differences"))):
        ax.set_title(style.panel_heading("abc"[k], t))
    style.apply_figure_heading(fig, "Gas PBL Transport, Lowest Kilometre",
                               f"Supercell LNOx, {args.label}: mass-weighted mean below 1 km, transport on and off")
    fig.tight_layout(rect=(0, 0, 1, 0.84))
    for ext in ("png", "pdf"):
        fig.savefig(args.outdir / f"pbl_diff_timeseries.{ext}")
    plt.close(fig)
    summary["below_1km_mean_diff_ppbv"] = {
        f"{tmin[i]:.0f} min": {"NOx": float(series["on"]["nox"][i] - series["off"]["nox"][i]),
                               "O3": float(series["on"]["o3"][i] - series["off"]["o3"][i])} for i in snaps}
    (args.outdir / "pbl_diff_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()

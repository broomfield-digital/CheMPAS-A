#!/usr/bin/env python
"""LNOx supercell pages of the CheMPAS-A paper figure gallery.

The function bodies below are copied unchanged from
NCAR/CheMPAS-A-Dev scripts/plot_paper_gallery.py at
5766b208cf160ae685a656eacbad015a7c2ae9b4: ppbv, zc_km, species_stream,
finish, and the LNOx supercell section (lnox_common, lnox_nox_xsection,
lnox_species_timeseries, lnox_profiles). The other gallery pages, their data
paths and the cartopy map helpers are not included.

Set DATA["lnox"] to the run's output.nc before calling a page function.
The pages expect the 2-h run with 10-min output: record 4 is 40 min.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import netCDF4
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import style  # noqa: E402
import style_gallery as sg  # noqa: E402

DATA: dict[str, Path | None] = {"lnox": None}

SINGLE = (5.6, 4.2)
WIDE = (7.2, 4.4)


def ppbv(q, species):
    return sg.mmr_to_vmr_ppbv(np.asarray(q), species)


def zc_km(zgrid):
    z = np.asarray(zgrid)
    return 0.5 * (z[:, :-1] + z[:, 1:]) / 1000.0


def species_stream(key):
    """High-frequency species stream written by the figure-source reruns
    (species.nc beside the case output); None when not yet produced."""
    path = DATA[key].parent / "species.nc"
    if not path.exists():
        return None
    ds = netCDF4.Dataset(path)
    xt = ds.variables["xtime"]
    times = ["".join(c.decode() for c in row).strip() for row in xt[:]]
    return ds, times


def finish(fig, ax, stem, title, caption, outdir, manifest, *, map_plot=False,
           subtitle=None):
    sub_kw = dict(ha="center", va="bottom",
                  fontsize=plt.rcParams["legend.fontsize"],
                  color=style.NCAR_COLORS["gray"])
    if map_plot:
        # The cartopy gridliner reports infinite label extents, which drives
        # matplotlib's auto title position to y=inf (and breaks tight bboxes);
        # pin the title y explicitly on GeoAxes.
        ax.set_title(sg.title_case(title), y=1.10 if subtitle else 1.04)
        if subtitle:
            ax.text(0.5, 1.045, sg.title_case(subtitle),
                    transform=ax.transAxes, **sub_kw)
    else:
        if subtitle:
            ax.set_title(sg.title_case(title), pad=24)
            ax.text(0.5, 1.015, sg.title_case(subtitle),
                    transform=ax.transAxes, **sub_kw)
        else:
            ax.set_title(sg.title_case(title))
    if not map_plot:
        # tight_layout fights the cartopy gridliner label layout on GeoAxes
        # ("Singular matrix"), so it is skipped for maps.
        fig.tight_layout()
    # Tight bounding box on every page: no wide canvas border around the
    # figure. (Titles are pinned to finite y on GeoAxes, so the historical
    # gridliner-driven degenerate crop no longer applies.)
    pad = 0.05
    for ext in ("pdf", "png"):
        fig.savefig(Path(outdir) / f"{stem}.{ext}", bbox_inches="tight",
                    pad_inches=pad)
    # Record the cropped width: the gallery scales pages by canvas width so
    # typography matches across figures.
    width = float(fig.get_size_inches()[0])
    try:
        tw = float(fig.get_tightbbox().width) + 2 * pad
        if np.isfinite(tw) and 0 < tw <= width:
            width = tw
    except Exception:
        pass
    plt.close(fig)
    manifest.append({"stem": stem, "title": title, "caption": caption,
                     "width_in": round(width, 2)})
    print(f"  rendered {stem}")


# ---------------------------------------------------------------------------
# LNOx supercell
# ---------------------------------------------------------------------------

def lnox_common(ds):
    x = ds["xCell"][:] / 1000.0
    y = ds["yCell"][:] / 1000.0
    zkm = zc_km(ds["zgrid"][:])
    return x, y, zkm


def lnox_nox_xsection(ds, mode, outdir, manifest):
    t = 4  # 18:00 + 40 min at 10-min output
    x, y, zkm = lnox_common(ds)
    y0 = y[np.argmin(np.abs(y - 46.0))]  # a single mesh row, not two staggered
    row = np.abs(y - y0) < 0.1
    order = np.argsort(x[row])
    xs = x[row][order]
    z = zkm[row][order]
    nox = np.ma.filled(
        (ppbv(ds["qNO"][t], "qNO") + ppbv(ds["qNO2"][t], "qNO2"))[row][order], np.nan
    )
    lwc = np.ma.filled(((ds["qc"][t] + ds["qr"][t]) * 1.0e3)[row][order], np.nan)
    do3 = np.ma.filled(
        (ppbv(ds["qO3"][t], "qO3") - ppbv(ds["qO3"][0], "qO3"))[row][order], np.nan
    )
    X = np.repeat(xs[:, None], z.shape[1], axis=1)

    fig, ax = plt.subplots(figsize=WIDE)
    ntarget = 8 if mode == sg.COLOR else 5
    levels = sg.nice_levels(0.0, float(np.nanpercentile(nox, 99.5)), ntarget)
    sg.field_on_axes(ax, X, z, nox, mode=mode, levels=levels[1:], fig=fig,
                     cbar_label="NOx [ppbv]")
    lwc_levels, lws = [0.01, 0.5, 2.0], [0.7, 1.4, 2.2]
    lwc_color = style.NCAR_COLORS["aqua"] if mode == sg.COLOR else "0.6"
    cs = ax.contour(X, z, lwc, levels=lwc_levels, colors=[lwc_color],
                    linewidths=lws)
    handles = [plt.Line2D([], [], color=lwc_color, lw=w) for w in lws]
    do3_color = style.NCAR_COLORS["gray"] if mode == sg.COLOR else "black"
    span = float(np.nanmax(np.abs(do3)))
    step = sg.nice_levels(0.0, span, 3)[1]
    dlev = [-2 * step, -step, step, 2 * step]
    dlw = [1.6 if abs(lv) > 1.5 * step else 0.8 for lv in dlev]
    cs2 = ax.contour(X, z, do3, levels=dlev, colors=[do3_color], linewidths=dlw,
                     linestyles=["dotted" if lv > 0 else "dashed" for lv in dlev])
    leg = ax.legend(handles, [f"{v:g}" for v in lwc_levels],
                    title=sg.title_case("LWC [g kg$^{-1}$]"), loc="upper left")
    leg.get_title().set_fontsize(plt.rcParams["legend.fontsize"])
    ax.add_artist(leg)
    d_handles = [plt.Line2D([], [], color=do3_color, lw=w,
                            linestyle="dashed" if lv < 0 else "dotted")
                 for lv, w in zip(dlev, dlw)]
    leg2 = ax.legend(d_handles, [f"{lv:+g}" for lv in dlev],
                     title=style.format_title("\u0394O3 [ppbv]"),
                     loc="upper right")
    leg2.get_title().set_fontsize(plt.rcParams["legend.fontsize"])
    ax.set_xlabel(sg.title_case("Distance [km]"))
    ax.set_ylabel(sg.title_case("Height [km]"))
    ax.set_ylim(0, 20)
    finish(fig, ax, "lnox_nox_xsection",
           "LNOx Supercell NOx Cross Section at 40 min",
           "Lightning-NOx supercell, x-z section at y = 46 km, t = 40 min. "
           "NOx (contours/fill, ppbv), liquid+rain water content (line-weight "
           "legend, g kg-1), and ozone change from the initial profile "
           "(legend, ppbv).",
           outdir, manifest, subtitle="Idealized Supercell, Minute 40 of 120")


def lnox_species_timeseries(ds, mode, outdir, manifest):
    x, y, zkm = lnox_common(ds)
    lwc_col = np.max((ds["qc"][4] + ds["qr"][4]) * 1e3, axis=1)
    region = lwc_col > 0.1
    hf = species_stream("lnox")
    src, cadence = (hf[0], 2.0) if hf else (ds, 10.0)
    nt = src.dimensions["Time"].size
    minutes = np.arange(nt) * cadence
    series = {}
    below = zkm[region] < 15.0
    for sp in ("qNO", "qNO2"):
        vals = [np.mean(ppbv(src[sp][t], sp)[region][below]) for t in range(nt)]
        series[sp] = np.array(vals)
    series["qNOx"] = series["qNO"] + series["qNO2"]
    o3 = np.array([np.mean(ppbv(src["qO3"][t], "qO3")[region][below]) for t in range(nt)])
    series["dO3"] = o3 - o3[0]

    fig, ax = plt.subplots(figsize=SINGLE)
    labels = {"qNO": "NO", "qNO2": "NO2", "qNOx": "NOx", "dO3": "\u0394O3"}
    for i, (sp, v) in enumerate(series.items()):
        kw = sg.species_style("qO3" if sp == "dO3" else sp, mode, i)
        ax.plot(minutes, v, label=style.format_title(labels[sp]), **kw)
    ax.axhline(0.0, color="0.7", lw=0.6)
    ax.set_xlabel(sg.title_case("Time [min]"))
    ax.set_ylabel(sg.title_case("Mixing Ratio [ppbv]"))
    ax.legend()
    finish(fig, ax, "lnox_species_timeseries",
           "Storm-Region Species Evolution",
           "Storm-region mean NO, NO2, NOx, and ozone change (ppbv) below "
           "15 km over the 2-h LNOx supercell run. The storm region is "
           "fixed at t = 40 min: the columns whose peak cloud-plus-rain "
           "water exceeds 0.1 g kg-1.",
           outdir, manifest, subtitle="Idealized 2-h Supercell Simulation")


def lnox_profiles(ds, mode, outdir, manifest):
    t = 4
    x, y, zkm = lnox_common(ds)
    lwc_col = np.max((ds["qc"][t] + ds["qr"][t]) * 1e3, axis=1)
    region = lwc_col > 0.1
    zmean = zkm[region].mean(axis=0)
    fig, ax = plt.subplots(figsize=(4.6, 4.6))
    for i, sp in enumerate(("qNO", "qNO2")):
        prof = ppbv(ds[sp][t], sp)[region].mean(axis=0)
        ax.plot(prof, zmean, label=style.format_title(sp[1:]),
                **sg.species_style(sp, mode, i))
    nox = (ppbv(ds["qNO"][t], "qNO") + ppbv(ds["qNO2"][t], "qNO2"))[region].mean(axis=0)
    ax.plot(nox, zmean, label=style.format_title("NOx"),
            **sg.species_style("qNOx", mode, 2))
    ax.set_xlabel(sg.title_case("Mixing Ratio [ppbv]"))
    ax.set_ylabel(sg.title_case("Height [km]"))
    ax.set_ylim(0, 20)
    ax.legend()
    finish(fig, ax, "lnox_profiles", "Storm-Region NOx Profiles at 40 min",
           "Storm-region mean NO, NO2, and NOx profiles (ppbv) at t = 40 min "
           "of the lightning-NOx supercell simulation. Storm region as "
           "defined for the species time series.",
           outdir, manifest, subtitle="Idealized Supercell, Minute 40 of 120")

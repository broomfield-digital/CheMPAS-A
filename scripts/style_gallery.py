"""Interim style extensions for the paper figure gallery.

Black-and-white labeled-contour vocabulary layered on top of ``style.py``
WITHOUT modifying it: ``style.py`` is byte-pinned in three committed figure
manifests, so the mono-mode vocabulary lives here during the co-author review
phase and merges into ``style.py`` (protocol v2) once the gallery selection is
made. See FABLE_PLAN_PLOTS.md.

Directives encoded here (maintainer, 2026-08-16):
- single-panel bias; Title Case labels; annotations minimal (detail in captions)
- literature units (ppbv, DU, 10^15 molec cm^-2, ng N m^-2 s^-1, ...)
- "nice" contour levels: 1-2 significant decimals, from a shared picker
- mono mode: grayscale ramp or labeled contours; negative levels dashed;
  species distinguished by line style (stable per species), not color
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import style  # noqa: E402

# --------------------------------------------------------------------------
# Modes
# --------------------------------------------------------------------------

COLOR = "color"
MONO = "bw"

#: Stable line-style vocabulary (mono mode). Order matters; assignment is by
#: species via SPECIES_LINESTYLES, never by plot-local cycling.
LINE_STYLES = ["-", "--", "-.", ":", (0, (5, 1, 1, 1, 1, 1))]

#: One stable line style per species within the whole gallery (mono analogue
#: of style.species_color). Extend here, never ad hoc in a plot.
SPECIES_LINESTYLES = {
    "qO3": "-",
    "qNO": "--",
    "qNO2": "-.",
    "qNOx": ":",
    "qCO": (0, (5, 1, 1, 1, 1, 1)),
    "qCH4": (0, (3, 1)),
    "qO": "--",
    "qO1D": ":",
    "qA": "-",
    "qB": "--",
    "qAB": "-.",
}

MONO_LINE_COLOR = "black"
MONO_SEQUENTIAL_CMAP = "Greys"

#: Contour-label font size (pt); set per mode by setup() from the active
#: preset (legend_small). Merges into style.FontSizes as ``contour_label``
#: at protocol-v2 time.
CONTOUR_LABEL_SIZE = 9


def setup(mode: str) -> None:
    """Initialize matplotlib for a gallery mode.

    color -> the standard NCAR language (style.setup default context).
    bw    -> publication-context sizes, grayscale-safe rcParams.
    """
    global CONTOUR_LABEL_SIZE
    if mode == COLOR:
        style.setup(context="default")
        CONTOUR_LABEL_SIZE = style.FONT_SIZES_DEFAULT.legend_small
    elif mode == MONO:
        style.setup(context="publication")
        import matplotlib as mpl

        mpl.rcParams["axes.prop_cycle"] = mpl.cycler(
            color=[MONO_LINE_COLOR] * len(LINE_STYLES), linestyle=LINE_STYLES
        )
        mpl.rcParams["image.cmap"] = MONO_SEQUENTIAL_CMAP
        CONTOUR_LABEL_SIZE = style.FONT_SIZES_PUBLICATION.legend_small
    else:
        raise ValueError(f"unknown gallery mode: {mode!r}")


def species_style(name: str, mode: str, index: int = 0) -> dict:
    """Return plot kwargs giving a species its stable identity in a mode."""
    if mode == MONO:
        ls = SPECIES_LINESTYLES.get(name, LINE_STYLES[index % len(LINE_STYLES)])
        return {"color": MONO_LINE_COLOR, "linestyle": ls}
    return {"color": style.species_color(name, index), "linestyle": "-"}


# --------------------------------------------------------------------------
# Nice levels
# --------------------------------------------------------------------------

_NICE_STEPS = np.array([1.0, 2.0, 2.5, 5.0, 10.0])


def nice_levels(vmin: float, vmax: float, target: int = 7) -> np.ndarray:
    """Contour levels on a 1/2/2.5/5 mantissa grid, 1-2 decimals preferred.

    Returns levels covering [vmin, vmax] with roughly ``target`` intervals.
    """
    if not np.isfinite([vmin, vmax]).all() or vmax <= vmin:
        raise ValueError(f"bad range for levels: [{vmin}, {vmax}]")
    raw = (vmax - vmin) / max(target, 1)
    mag = 10.0 ** np.floor(np.log10(raw))
    step = _NICE_STEPS[np.searchsorted(_NICE_STEPS, raw / mag)] * mag
    lo = np.floor(vmin / step) * step
    hi = np.ceil(vmax / step) * step
    levels = np.arange(lo, hi + 0.5 * step, step)
    # Snap tiny float dust so labels print clean (<= 2 decimals where possible)
    return np.round(levels, 10)


def level_format(levels: np.ndarray) -> str:
    """A printf format giving the fewest decimals that distinguish levels."""
    step = np.min(np.diff(levels)) if len(levels) > 1 else abs(levels[0]) or 1.0
    for dec in range(0, 4):
        if np.isclose(round(step, dec), step, rtol=0, atol=1e-12 * max(1, step)):
            return f"%.{dec}f"
    return "%.3g"


# --------------------------------------------------------------------------
# Field rendering: filled color vs labeled mono contours
# --------------------------------------------------------------------------

def field_on_axes(
    ax,
    x,
    y,
    z,
    *,
    mode: str,
    levels: np.ndarray,
    tri: bool = False,
    cbar_label: str | None = None,
    fig=None,
    cmap: str | None = None,
):
    """Draw a scalar field: filled contours + colorbar (color mode) or
    labeled line contours (mono mode). Returns the contour set."""
    fmt = level_format(np.asarray(levels))
    args = (x, z) if (tri and y is None) else (x, y, z)
    if mode == COLOR:
        f = ax.tricontourf if tri else ax.contourf
        cs = f(*args, levels=levels, cmap=cmap or "ncar_sunset", extend="max")
        if fig is not None:
            cb = fig.colorbar(cs, ax=ax, shrink=0.9, pad=0.02)
            if cbar_label:
                cb.set_label(style.format_title(cbar_label))
        return cs
    c = ax.tricontour if tri else ax.contour
    neg = [lv for lv in levels if lv < 0]
    pos = [lv for lv in levels if lv >= 0]
    css = []
    if pos:
        cs = c(*args, levels=pos, colors=MONO_LINE_COLOR, linewidths=0.9)
        ax.clabel(cs, fmt=fmt, fontsize=CONTOUR_LABEL_SIZE, inline=True)
        css.append(cs)
    if neg:
        cs = c(*args, levels=neg, colors=MONO_LINE_COLOR, linewidths=0.9,
               linestyles="dashed")
        ax.clabel(cs, fmt=fmt, fontsize=CONTOUR_LABEL_SIZE, inline=True)
        css.append(cs)
    return css[0] if css else None


# --------------------------------------------------------------------------
# Literature units
# --------------------------------------------------------------------------

M_AIR = 28.9647e-3  # kg mol-1
AVOGADRO = 6.02214076e23
G0 = 9.80665

MOLAR_MASS = {  # kg mol-1
    "qO3": 47.997e-3,
    "qNO": 30.006e-3,
    "qNO2": 46.006e-3,
    "qCO": 28.010e-3,
    "qCH4": 16.043e-3,
    "qO": 15.999e-3,
    "qO1D": 15.999e-3,
    "qA": 28.0e-3,   # idealized tracers: report as mmr-scaled, not converted
    "qB": 28.0e-3,
    "qAB": 28.0e-3,
}


def mmr_to_vmr_ppbv(q: np.ndarray, species: str) -> np.ndarray:
    """Mass mixing ratio (kg/kg) -> volume mixing ratio in ppbv."""
    return q * (M_AIR / MOLAR_MASS[species]) * 1.0e9


def column_molec_cm2(q, rho, zgrid, species: str) -> np.ndarray:
    """Vertical column (molec cm-2) from mmr, dry density (kg m-3), and the
    staggered geometric heights zgrid (nCells, nVertLevels+1)."""
    dz = np.diff(zgrid, axis=-1)  # (nCells, nVertLevels)
    mass_col = np.sum(q * rho * dz, axis=-1)  # kg m-2
    return mass_col / MOLAR_MASS[species] * AVOGADRO * 1.0e-4


DU_MOLEC_CM2 = 2.6867e16


def column_du(q, rho, zgrid, species: str = "qO3") -> np.ndarray:
    """Vertical column in Dobson Units."""
    return column_molec_cm2(q, rho, zgrid, species) / DU_MOLEC_CM2


# --------------------------------------------------------------------------
# Labels
# --------------------------------------------------------------------------

UNIT_TOKENS = {"min", "km", "kg", "g", "m", "s", "h", "hPa", "ppbv", "ppmv",
               "pptv", "molec", "cm", "mol", "DU", "UTC"}


def title_case(text: str) -> str:
    """Title Case for label prose; unit tokens and bracketed unit strings keep
    their case exactly (never 'Kg' or 'Min'). Chemical formatting applied last."""
    minor = {"a", "an", "and", "at", "by", "for", "in", "of", "on", "or",
             "the", "to", "vs", "with", "above", "below", "over", "under",
             "from", "into", "onto", "upon", "since", "between", "through",
             "during", "without", "within", "toward", "towards", "against"}
    out = []
    depth = 0
    for i, w in enumerate(text.split()):
        bare = w.strip("[](),.;:")
        opens = w.count("[") + w.count("(")
        closes = w.count("]") + w.count(")")
        in_unit = depth > 0 or w.startswith(("[", "("))
        depth = max(0, depth + opens - closes)
        if in_unit or bare in UNIT_TOKENS or bare.lower() in UNIT_TOKENS \
                or any(ch.isdigit() for ch in bare) or "$" in w:
            out.append(w)
        elif i > 0 and w.lower() in minor:
            out.append(w.lower())
        else:
            out.append(w[:1].upper() + w[1:])
    return style.format_title(" ".join(out))

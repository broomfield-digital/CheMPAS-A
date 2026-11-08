#!/usr/bin/env python3
"""Figures for P2 of docs/plans/pbl-transport.md (column diffusion kernel and tendency proof).

Normally run through make_p2_plots.sh, which first writes the kernel data with
dump_pbl_diffusion.F90 (the production kernel on the reference fixtures and a
stiffness sweep).

Usage: plot_p2_evidence.py KERNEL_DATA_DIR REALIZED_RUN_ROOT OUTDIR
  KERNEL_DATA_DIR    output of dump_pbl_diffusion
  REALIZED_RUN_ROOT  a RUN_ROOT of check_realized_increment.sh (rest, moving and fresh cases)
  OUTDIR             receives p2_1..p2_4 as PNG and PDF
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from netCDF4 import Dataset

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "testing/chemistry_transport"))
import style  # noqa: E402
import check_realized_increment as cri  # noqa: E402

COLUMNS = REPO / "testing/chemistry_transport/columns"
BLUE, AQUA, ORANGE, GREEN, GRAY = (style.NCAR_COLORS[k] for k in ("ncar_blue", "aqua", "orange", "green", "gray"))
SPECIES_STYLE = {"qO3": (BLUE, "o"), "qNO2": (ORANGE, "s"), "qNO": (GREEN, "^")}
SPECIES_TEXT = {"qO3": "O3", "qNO2": "NO2", "qNO": "NO"}
FS = style.FONT_SIZES_DEFAULT


def save(fig, outdir, stem):
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {stem}.png/.pdf")


def read_fixture(name):
    lines = (COLUMNS / f"diffusion_{name}.txt").read_text().split("\n")
    nz, nt = map(int, lines[1].split())
    dt = float(lines[2])
    m_d = np.array(lines[3].split(), float)
    G = np.array(lines[4].split(), float)
    q0 = np.array([l.split() for l in lines[5:5 + nt]], float)
    ref = np.array([l.split() for l in lines[5 + nt:5 + 2 * nt]], float)
    return dt, m_d, G, q0, ref


# ---------------------------------------------------------------------------------------------------------------
# Figure 1: kernel profiles against the exact solution
# ---------------------------------------------------------------------------------------------------------------
def fig_kernel_profiles(kdir, outdir):
    panels = [("pulses", "Isolated Pulses", ["bottom pulse", "middle pulse", "top pulse"]),
              ("sharp_gradient", "Sharp Gradient", ["step profile"]),
              ("large_timestep", "Large Time Step (dt·G/m = 960)", ["smooth", "step", "middle pulse"])]
    colors = [BLUE, AQUA, ORANGE]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), sharey=True)
    for ax, (name, title, labels) in zip(axes, panels):
        dt, m_d, G, q0, ref = read_fixture(name)
        kern = np.loadtxt(kdir / f"kernel_{name}.txt", ndmin=2)
        k = np.arange(1, m_d.size + 1)
        for n, label in enumerate(labels):
            c = colors[n]
            ax.plot(q0[n] * 1e9, k, color=c, lw=1.0, ls="--", alpha=0.8)
            ax.plot(kern[n] * 1e9, k, color=c, lw=2.0, label=label)
            ax.plot(ref[n] * 1e9, k, ls="none", marker="o", ms=8, mfc="none", mec=c, mew=1.2)
        ax.set_title(title)
        ax.set_xlabel("Mixing Ratio (10$^{-9}$ kg kg$^{-1}$)")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=len(labels), fontsize=FS.legend_small,
                  frameon=False)
    axes[0].set_ylabel("Layer (1 = Ground)")
    style_handles = [Line2D([], [], color=GRAY, lw=1.0, ls="--", label="initial"),
                     Line2D([], [], color=GRAY, lw=2.0, label="kernel (Fortran)"),
                     Line2D([], [], color=GRAY, ls="none", marker="o", mfc="none", ms=8, label="exact solution")]
    fig.legend(handles=style_handles, loc="upper center", ncol=3, fontsize=FS.legend_small,
               bbox_to_anchor=(0.5, 0.93), frameon=False)
    fig.suptitle("Closed-Column Implicit Diffusion: Kernel vs. Exact Solution", fontsize=style.SUPTITLE_SIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, outdir, "p2_1_kernel_profiles")


# ---------------------------------------------------------------------------------------------------------------
# Figure 2: conservation and uniform-profile drift against stiffness for the three forms
# ---------------------------------------------------------------------------------------------------------------
def numpy_forms(kappa):
    """Earlier kernel forms, reproduced in NumPy on the test_pbl_diffusion column."""
    nz = 60
    k = np.arange(1, nz + 1)
    m = 50.0 + 150.0 * np.sin(0.37 * k) ** 2
    q0 = np.vstack([1e-9 * (1 + np.sin(k.astype(float))), np.where(k <= 5, 5e-9, 0.0)])
    dt = 600.0
    G = np.full(nz + 1, kappa * m.mean() / dt)
    a = np.array([-dt * (G[i] if i > 0 else 0) / m[i] for i in range(nz)])
    c = np.array([-dt * (G[i + 1] if i < nz - 1 else 0) / m[i] for i in range(nz)])
    b = 1 - a - c

    def thomas(r):
        cp = np.zeros(nz); dp = np.zeros((r.shape[0], nz))
        cp[0] = c[0] / b[0]; dp[:, 0] = r[:, 0] / b[0]
        for i in range(1, nz):
            w = b[i] - a[i] * cp[i - 1]; cp[i] = c[i] / w; dp[:, i] = (r[:, i] - a[i] * dp[:, i - 1]) / w
        x = np.zeros_like(r); x[:, -1] = dp[:, -1]
        for i in range(nz - 2, -1, -1):
            x[:, i] = dp[:, i] - cp[i] * x[:, i + 1]
        return x

    def flux(q):
        F = np.zeros((q.shape[0], nz + 1))
        F[:, 1:nz] = G[1:nz] * (q[:, :-1] - q[:, 1:])
        return F

    state = thomas
    def state_flux(q):
        F = flux(thomas(q)); return q + dt * (F[:, :-1] - F[:, 1:]) / m

    out = {}
    for name, f in (("state", state), ("state_flux", state_flux)):
        q = q0.copy()
        for _ in range(1000):
            q = f(q)
        mass = max(abs((m * q[n]).sum() - (m * q0[n]).sum()) / (m * q0[n]).sum() for n in range(2))
        u = f(np.full((2, nz), 3e-9))
        out[name] = (mass, np.abs(u - 3e-9).max() / 3e-9)
    return out


def fig_stiffness(kdir, outdir):
    data = np.loadtxt(kdir / "kernel_stiffness.txt")
    kap = data[:, 0]
    forms = {"state": [], "state_flux": []}
    for kv in kap:
        res = numpy_forms(kv)
        for name in forms:
            forms[name].append(res[name])
    floor = 1e-18
    series = [("Thomas solve for the state (NumPy)", BLUE, "o", "--", np.array(forms["state"])),
              ("state + flux-form finish (NumPy)", AQUA, "s", "-.", np.array(forms["state_flux"])),
              ("increment form, final (Fortran)", ORANGE, "^", "-", data[:, 1:3])]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    panels = [(0, "Mass Residual After 1000 Calls", 1e-10, "frozen limit 1e-10"),
              (1, "Uniform Profile: Max Relative Change per Call", 1e-12, "frozen limit 1e-12")]
    for ax, (col, title, limit, limit_label) in zip(axes, panels):
        for label, color, marker, ls, vals in series:
            y = np.maximum(vals[:, col], floor)
            ax.loglog(kap, y, color=color, ls=ls, lw=2.0, marker=marker, ms=6, label=label)
        ax.axhline(limit, color=style.NCAR_COLORS["red"], lw=1.2, ls=":")
        ax.text(kap[0], limit * 1.6, limit_label, color=style.NCAR_COLORS["red"], fontsize=FS.legend_small)
        ax.set_ylim(floor / 3, 1e-6)
        ax.set_xlabel("Stiffness κ = dt·G / m$_d$")
        ax.set_title(title)
        ax.grid(True, which="major", alpha=0.3)
    axes[1].text(kap[len(kap) // 2], floor * 6, "increment form: exactly 0 (plotted at 1e-18)",
                 color=GRAY, fontsize=FS.legend_small, ha="center")
    axes[0].set_ylabel("Relative Error")
    axes[0].legend(loc="upper left", fontsize=FS.legend_small)
    fig.suptitle("Kernel Forms vs. Stiffness (60-Layer Column)", fontsize=style.SUPTITLE_SIZE)
    fig.tight_layout()
    save(fig, outdir, "p2_2_kernel_stiffness")


# ---------------------------------------------------------------------------------------------------------------
# Figure 3: realized vs requested increments in the model
# ---------------------------------------------------------------------------------------------------------------
def load_run(run):
    order = cri.transport_order(run / "log.atmosphere.0000.out")
    ds = Dataset(run / "output.nc")
    zk = 0.5 * (ds["zgrid"][0, :-1] + ds["zgrid"][0, 1:]) / 1000.0
    return order, ds, zk


def fig_realized(root, outdir):
    fig, axes = plt.subplots(2, 2, figsize=(15, 13))

    # (a) source at rest, step 1, cell 0 (every cell identical)
    ax = axes[0, 0]
    order, ds, zk = load_run(root / "rest-mono-source")
    k = np.arange(1, zk.size + 1)
    dt = 3.0
    for n, name in enumerate(order, start=1):
        c, mk = SPECIES_STYLE[name]
        req = dt * cri.S0 * n * (1 + 0.5 * np.sin(k))
        real = ds[name][1, 0] - ds[name][0, 0]
        ax.plot(req * 1e13, zk, color=c, lw=1.6, label=f"{style.species_label(SPECIES_TEXT[name])} requested")
        ax.plot(real * 1e13, zk, ls="none", marker=mk, ms=5, mfc="none", mec=c,
                label=f"{style.species_label(SPECIES_TEXT[name])} realized")
    ax.set_xlabel("Increment in One Step (10$^{-13}$ kg kg$^{-1}$)")
    ax.set_ylabel("Height (km)")
    ax.set_title("(a) Source Test at Rest, Monotonic")
    ax.grid(True, alpha=0.3)
    ds.close()

    # (b) diffusion at rest, step 1
    ax = axes[0, 1]
    order, ds, zk = load_run(root / "rest-mono-diffusion")
    m_d = ds["zz"][:] * ds["rho_zz"][0] * np.diff(ds["zgrid"][:], axis=1)
    for name in order:
        c, mk = SPECIES_STYLE[name]
        q0 = ds[name][0][:1]
        req = (cri.diffuse(q0, m_d[:1], dt) - q0)[0]
        real = ds[name][1, 0] - ds[name][0, 0]
        scale = np.max(np.abs(req))
        ax.plot(req / scale, zk, color=c, lw=1.6,
                label=f"{style.species_label(SPECIES_TEXT[name])} requested (dense solve)")
        ax.plot(real / scale, zk, ls="none", marker=mk, ms=5, mfc="none", mec=c,
                label=f"{style.species_label(SPECIES_TEXT[name])} realized")
    ax.set_xlabel("Increment / Species' Max |Requested|")
    ax.set_title("(b) Diffusion Test at Rest, Monotonic")
    ax.grid(True, alpha=0.3)
    ds.close()

    # (c) column-budget error for every case
    ax = axes[1, 0]
    advs = ["mono", "ordinary", "posdef", "mono_split", "ordinary_split", "mono_chem2"]
    adv_text = {"mono": "monotonic", "ordinary": "ordinary", "posdef": "positive-definite*",
                "mono_split": "monotonic, split", "ordinary_split": "ordinary, split",
                "mono_chem2": "monotonic, chem/2 steps"}
    rows, ylabels = [], []
    for regime, color in (("rest", BLUE), ("moving", ORANGE)):
        for adv in advs:
            for mode, mk in (("source", "o"), ("diffusion", "s")):
                rows.append((regime, adv, mode, color, mk))
        rows.append((regime, "mono", "off", color, "^"))
    floor = 1e-18                                   # exact zeros plotted here
    for i, (regime, adv, mode, color, mk) in enumerate(rows):
        rep = json.loads((root / f"{regime}-{adv}-{mode}" / "realized_increment.json").read_text())
        err = max(e["column_mass_error_relative"] for s in rep["steps"] for e in s["species"].values())
        ax.semilogx(max(err, floor), i, ls="none", marker=mk, ms=7, color=color)
        ylabels.append(f"{regime}: {adv_text[adv]}, {mode}")
    ax.axvline(1e-12, color=style.NCAR_COLORS["red"], lw=1.2, ls=":")
    ax.axvline(3.9e-10, color=GRAY, lw=1.2, ls="--")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(ylabels, fontsize=FS.legend_small)
    ax.invert_yaxis()
    ax.set_xlim(1e-17, 1e-8)
    ax.set_xlabel("Column Tracer-Mass Budget Error (relative)")
    ax.set_title("(c) Column Budget, Every Case")
    ax.legend(handles=[Line2D([], [], color=style.NCAR_COLORS["red"], ls=":", label="threshold 1e-12"),
                       Line2D([], [], color=GRAY, ls="--", label="error if weighted by ρ(t+dt)"),
                       Line2D([], [], color=BLUE, marker="o", ls="none", label="rest regime"),
                       Line2D([], [], color=ORANGE, marker="o", ls="none", label="moving regime"),
                       Line2D([], [], color=GRAY, marker="o", ls="none", label="source"),
                       Line2D([], [], color=GRAY, marker="s", ls="none", label="diffusion"),
                       Line2D([], [], color=GRAY, marker="^", ls="none", label="off")],
              fontsize=FS.legend_small, loc="lower right")
    ax.grid(True, axis="x", alpha=0.3)

    # (d) freshness: post-chemistry vs stale prediction
    ax = axes[1, 1]
    width = 0.38
    labels, fresh_err, stale_err, tol = [], [], [], []
    for adv in ("mono", "ordinary"):
        run, off = root / f"fresh-{adv}-diffusion", root / f"fresh-{adv}-off"
        order, ds, _ = load_run(run)
        offd = Dataset(off / "output.nc")
        m_d = ds["zz"][:] * ds["rho_zz"][0] * np.diff(ds["zgrid"][:], axis=1)
        for name in order:
            q0, q1, post = ds[name][0], ds[name][1], offd[name][1]
            fresh = cri.diffuse(post, m_d, 3.0) - post
            stale = cri.diffuse(q0, m_d, 3.0) - q0
            realized = q1 - post
            fresh_err.append(float(np.max(np.abs(realized - fresh))))
            stale_err.append(float(np.max(np.abs(realized - stale))))
            tol.append(float(np.max(cri.PARALLEL_ATOL + cri.PARALLEL_RTOL * np.abs(fresh))))
            labels.append(f"{style.species_label(SPECIES_TEXT[name])}\n{'mono.' if adv == 'mono' else 'ord.'}")
        ds.close(); offd.close()
    x = np.arange(len(labels))
    zfloor = 1e-26                                  # below every value: nothing is clipped
    ax.bar(x - width / 2, np.maximum(fresh_err, zfloor), width, color=BLUE, label="vs. post-chemistry diffusion (implemented)")
    ax.bar(x + width / 2, np.maximum(stale_err, zfloor), width, color=ORANGE, hatch="//", edgecolor="white",
           label="vs. stale (pre-chemistry) diffusion")
    ax.plot(x, tol, ls="none", marker="_", ms=22, mew=2, color=style.NCAR_COLORS["red"], label="tolerance")
    ax.set_yscale("log")
    ax.set_ylim(zfloor / 3, 1e-8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=FS.legend_small)
    ax.set_ylabel("Max |realized − predicted| (kg kg$^{-1}$)")
    ax.set_title("(d) Photolysis On: Post-Chemistry State Used")
    ax.legend(fontsize=FS.legend_small, loc="upper left")
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle("Realized vs. Requested Gas Increments in the Model (PBL Test Double, 8 Ranks)",
                 fontsize=style.SUPTITLE_SIZE)
    fig.text(0.5, -0.005, "Horizontally uniform dry base state; panels (a, b) show cell 0, and every cell is identical. "
             "*MPAS runs positive-definite as monotonic. In (d) the post-chemistry errors are roundoff, 4e-25 to 3e-23.",
             ha="center", fontsize=FS.legend_small, color=GRAY)
    fig.tight_layout(rect=(0, 0.01, 1, 1))
    fig.subplots_adjust(hspace=0.34)
    shared = [Line2D([], [], color=SPECIES_STYLE[n][0], marker=SPECIES_STYLE[n][1], mfc="none", lw=1.6,
                     label=style.species_label(SPECIES_TEXT[n])) for n in ("qO3", "qNO2", "qNO")]
    shared += [Line2D([], [], color=GRAY, lw=1.6, label="requested: (a) S·dt, (b) dense solve"),
               Line2D([], [], color=GRAY, ls="none", marker="o", mfc="none", label="realized in the model")]
    top = axes[0, 0].get_position().y0
    fig.legend(handles=shared, loc="upper center", ncol=5, fontsize=FS.legend_small, frameon=False,
               bbox_to_anchor=(0.5, top - 0.045))
    save(fig, outdir, "p2_3_realized_increment")


# ---------------------------------------------------------------------------------------------------------------
# Figure 4: the moving regime
# ---------------------------------------------------------------------------------------------------------------
def fig_moving(root, outdir):
    order, off, zk = load_run(root / "moving-mono-off")
    shades = [plt.cm.Blues(v) for v in (0.45, 0.7, 0.95)]
    fig, axes = plt.subplots(2, 2, figsize=(13, 10.5))

    ax = axes[0, 0]
    zw = off["zgrid"][0] / 1000.0
    for t in range(1, 4):
        ax.plot(off["w"][t, 0] * 1000.0, zw, color=shades[t - 1], lw=2.0, label=f"t = {3 * t} s")
    ax.axvline(0, color=GRAY, lw=0.8)
    ax.set_xlabel("w (mm s$^{-1}$)")
    ax.set_ylabel("Height (km)")
    ax.set_title("(a) Vertical Velocity (Uniform in Every Cell)")
    ax.legend(fontsize=FS.legend_small)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    rho = off["rho_zz"][:, 0, :]
    for t in range(1, 4):
        ax.plot(np.abs(rho[t] - rho[t - 1]) / rho[t - 1] * 1e6, zk, color=shades[t - 1], lw=2.0,
                label=f"step ending {3 * t} s")
    ax.set_xlabel("|Δρ| / ρ per Step (10$^{-6}$)")
    ax.set_title("(b) Dry Density Change per Step")
    ax.legend(fontsize=FS.legend_small)
    ax.grid(True, alpha=0.3)

    _, mono, _ = load_run(root / "moving-mono-source")
    _, ordn, _ = load_run(root / "moving-ordinary-source")
    k = np.arange(1, zk.size + 1)
    ax = axes[1, 0]
    n = order.index("qO3") + 1
    req = 3.0 * cri.S0 * n * (1 + 0.5 * np.sin(k))
    dm = mono["qO3"][1, 0] - off["qO3"][1, 0]
    do = ordn["qO3"][1, 0] - off["qO3"][1, 0]
    ax.plot(req * 1e13, zk, color=GRAY, lw=1.2, ls="--", label="requested S·dt")
    ax.plot(dm * 1e13, zk, color=BLUE, lw=2.0, label="monotonic path")
    ax.plot(do * 1e13, zk, color=ORANGE, lw=1.4, ls="-.", label="ordinary path")
    ax.set_xlabel("O$_3$ Source Increment After One Step (10$^{-13}$ kg kg$^{-1}$)")
    ax.set_ylabel("Height (km)")
    ax.set_title("(c) Source Increment vs. Transport Off (Curves Coincide)")
    ax.legend(fontsize=FS.legend_small)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    for nn, name in enumerate(order, start=1):
        c, mk = SPECIES_STYLE[name]
        req = 3.0 * cri.S0 * nn * (1 + 0.5 * np.sin(k))
        dm = mono[name][1, 0] - off[name][1, 0]
        do = ordn[name][1, 0] - off[name][1, 0]
        rel = np.maximum(np.abs(dm - do) / np.abs(req), 1e-18)
        ax.semilogx(rel, zk, color=c, marker=mk, ms=4, lw=1.2, label=style.species_label(SPECIES_TEXT[name]))
    ax.set_xlabel("|Monotonic − Ordinary| / Requested, per Layer")
    ax.set_title("(d) Per-Layer Difference Between Advection Paths")
    ax.legend(fontsize=FS.legend_small, title="Species", title_fontsize=FS.legend_small)
    ax.grid(True, alpha=0.3)

    for a in axes.flat:
        a.set_ylim(0, 30)
    fig.suptitle("Moving Regime: Uniform Drift Driven by the Case's Eddy Viscosity", fontsize=style.SUPTITLE_SIZE)
    fig.text(0.5, -0.005, "Monotonic stage 3 adds the source before transport, so it moves with the air; the ordinary "
             "path adds it in place. Both close the column budget to 4e-16. Exact zeros in (d) are plotted at 1e-18. "
             "Heights shown to 30 km of 50 km.",
             ha="center", fontsize=FS.legend_small, color=GRAY)
    fig.tight_layout(rect=(0, 0.01, 1, 1))
    save(fig, outdir, "p2_4_moving_regime")
    for d in (off, mono, ordn):
        d.close()


def main():
    kdir, root, outdir = (Path(a) for a in sys.argv[1:4])
    outdir.mkdir(parents=True, exist_ok=True)
    style.setup()
    fig_kernel_profiles(kdir, outdir)
    fig_stiffness(kdir, outdir)
    fig_realized(root, outdir)
    fig_moving(root, outdir)


if __name__ == "__main__":
    main()

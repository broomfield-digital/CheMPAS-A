#!/usr/bin/env python3
"""Render the established CheMPAS-A-Dev LNOx diagnostics for a supercell run.

Runs these plotting functions unchanged from this repository's scripts/:

  scripts/lnox_gallery.py    lnox_nox_xsection, lnox_species_timeseries,
                             lnox_profiles
  scripts/plot_lnox_o3.py    vertical, evolution, horizontal, comparison,
                             profiles, source, photolysis

scripts/plot_lnox_o3.py, style.py and style_gallery.py are byte-identical
copies from NCAR/CheMPAS-A-Dev at 5766b208cf160ae685a656eacbad015a7c2ae9b4.
scripts/lnox_gallery.py holds the LNOx pages of that revision's
plot_paper_gallery.py, also unchanged.

Usage:
  python testing/chemistry_transport/render_lnox_diagnostics.py RUN_DIR \\
      --outdir OUT_DIR

RUN_DIR must contain output.nc and supercell_init.nc. Figures go to
OUT_DIR/diagnostics and OUT_DIR/gallery as PNG and PDF.

The adaptations this run needs are applied at explicit seams instead of
by editing the copied scripts:

1. Mesh coordinates. The output stream omits xCell/yCell, so they are
   served from the run's supercell_init.nc after checking that its
   indexToCellID, zgrid and areaCell match the output exactly.
2. Colormap import. plot_lnox_o3.py builds a white-plateau colormap with
   two stops at position 0.10. Matplotlib >= 3.11 rejects repeated
   positions, so a repeated position is moved down by 1e-9. The rendered
   result is the same hard white-to-ramp step the script intends.
3. Vertical-velocity staggering. plot_lnox_o3.load_data keeps the lower
   face of w (w[:, :, :-1]). The model gates the source on the layer
   midpoint average 0.5 * (w[k] + w[k+1]) (mpas_lightning_nox.F), so the
   face average replaces it.
4. Vertical extent. This mesh has 60 levels up to about 50 km. The
   plot_lnox_o3 functions draw the whole column, so the data is limited
   to cell centres below --ztop-km (20 km by default, matching the gallery's
   y-limit and the 20-km design grid). The run aborts if NOx above the cut
   exceeds --max-nox-above-ztop-ppbv.
5. Updraft search level. find_updraft_y searches level nVertLevels // 4,
   which is about 5 km on the 40-level, 20-km grid it was written for. On this mesh
   that index is a different height, so the search uses the level nearest
   --updraft-search-km (5 km by default).
6. Source parameters come from the output file's global attributes. The
   renderer aborts unless config_lnox_gating_mode is 'altitude', because
   the plot_lnox_o3 source reconstructions implement altitude gating.
   The --nox-single panel, which reconstructs isotherm gating, is not
   rendered.
7. Updraft-profile height window. The source-gating z_min/z_max also
   clip the updraft-profile plot. The profiles use --profile-ztop-km
   instead, so NOx carried above z_max stays visible.
8. Gallery NOx line style. style.species_color has no qNOx entry, so the
   gallery's colour mode draws NOx in the same orange as NO2. NOx is drawn
   grey and dashed instead, following plot_lnox_o3.py.
"""
from __future__ import annotations

import argparse
import json
import sys
from contextlib import contextmanager
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from netCDF4 import Dataset  # noqa: E402

REPO_SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
M_AIR = 0.0289644  # kg/mol, host CHEM_M_AIR
M_NO, M_NO2 = 0.030, 0.046


class MeshAugmentedDataset:
    """Read-only output view that adds mesh coordinates from the init file."""

    MESH_VARS = ("xCell", "yCell")

    def __init__(self, output: Path, mesh: Path):
        self._out = Dataset(output)
        self._mesh = Dataset(mesh)
        for name in ("indexToCellID", "zgrid", "areaCell"):
            if not np.array_equal(self._out[name][:], self._mesh[name][:]):
                raise SystemExit(f"{mesh} does not match {output} in {name}")
        self.variables = {name: self._mesh.variables[name]
                          for name in self.MESH_VARS}
        self.variables.update(self._out.variables)
        self.dimensions = self._out.dimensions

    def __getitem__(self, name):
        return self.variables[name]

    def ncattrs(self):
        return self._out.ncattrs()

    def getncattr(self, name):
        return self._out.getncattr(name)

    def close(self):
        self._out.close()
        self._mesh.close()


@contextmanager
def strict_from_list():
    """Move repeated colormap stop positions to strictly increasing values."""
    original = LinearSegmentedColormap.from_list

    def from_list(name, colors, *args, **kwargs):
        stops = list(colors)
        pairs = stops and all(
            isinstance(c, (tuple, list)) and len(c) == 2
            and np.isscalar(c[0]) for c in stops)
        if pairs:
            pos = [float(c[0]) for c in stops]
            for i in range(len(pos) - 2, -1, -1):
                if pos[i] >= pos[i + 1]:
                    pos[i] = pos[i + 1] - 1.0e-9
            stops = [(p, c[1]) for p, c in zip(pos, stops)]
        return original(name, stops, *args, **kwargs)

    LinearSegmentedColormap.from_list = staticmethod(from_list)
    try:
        yield
    finally:
        LinearSegmentedColormap.from_list = original


def time_labels(ds) -> list[str]:
    return [b"".join(row).decode().strip() for row in ds["xtime"][:]]


def lnox_config(ds) -> dict:
    cfg = {k: ds.getncattr(f"config_lnox_{k}")
           for k in ("gating_mode", "source_rate", "w_threshold", "w_ref",
                     "z_min", "z_max", "t_min", "t_max", "nox_tau")}
    cfg["gating_mode"] = str(cfg["gating_mode"]).strip()
    return cfg


def render_gallery(output: Path, mesh: Path, outdir: Path) -> list[dict]:
    import lnox_gallery as pg
    import style_gallery as sg

    pg.DATA["lnox"] = output  # species_stream() then looks beside the run
    sg.setup(sg.COLOR)

    # Seam 8: species_color has no qNOx entry, so colour mode falls back to
    # palette index 2, which is the same orange as qNO2, and NO2 and NOx
    # become indistinguishable. Use plot_lnox_o3.py's NOx convention instead:
    # grey, dashed.
    import style
    species_style = sg.species_style

    def nox_distinct_style(name, mode, index=0):
        if name == "qNOx" and mode == sg.COLOR:
            return {"color": style.NCAR_COLORS["gray"], "linestyle": "--"}
        return species_style(name, mode, index)

    sg.species_style = nox_distinct_style
    ds = MeshAugmentedDataset(output, mesh)
    labels = time_labels(ds)
    # The gallery functions hard-code record 4 as 40 min at 10-min output.
    if labels[4] != "0000-01-01_18:40:00" or len(labels) != 13:
        raise SystemExit(f"gallery record layout mismatch: {labels}")
    manifest: list[dict] = []
    outdir.mkdir(parents=True, exist_ok=True)
    for fn in (pg.lnox_nox_xsection, pg.lnox_species_timeseries,
               pg.lnox_profiles):
        fn(ds, sg.COLOR, outdir, manifest)
    ds.close()
    (outdir / "gallery_manifest.json").write_text(json.dumps(manifest, indent=1))
    return manifest


def render_diagnostics(output: Path, mesh: Path, outdir: Path, cfg: dict,
                       every_min: float, level_km: float, ztop_km: float,
                       updraft_km: float, profile_ztop_km: float,
                       max_nox_above_ztop_ppbv: float) -> dict:
    with strict_from_list():
        import plot_lnox_o3 as lo3
    import style

    def dataset_factory(filename, mode="r"):
        if Path(filename).resolve() == output.resolve():
            return MeshAugmentedDataset(output, mesh)
        return Dataset(filename, mode)

    lo3.Dataset = dataset_factory
    style.setup()

    data = lo3.load_data(str(output))

    # Seam 3: layer-midpoint w, as in lightning_nox_inject.
    with Dataset(output) as ds:
        w_faces = ds["w"][:]
    data["w"] = 0.5 * (w_faces[:, :, :-1] + w_faces[:, :, 1:])

    # Seam 4: keep cell centres below ztop_km after checking NOx aloft.
    zc_km = 0.5 * (data["zgrid"][:, :-1] + data["zgrid"][:, 1:])
    keep = int((zc_km.mean(axis=0) < ztop_km).sum())
    nox_above = float((lo3.to_ppbv(data["qNO"][:, :, keep:], M_NO)
                       + lo3.to_ppbv(data["qNO2"][:, :, keep:], M_NO2)).max())
    if nox_above > max_nox_above_ztop_ppbv:
        raise SystemExit(f"NOx above {ztop_km} km reaches {nox_above:.3g} ppbv")
    for name in ("qNO", "qNO2", "qO3", "w", "qc", "theta", "pressure",
                 "j_no2"):
        if name in data:
            data[name] = data[name][:, :, :keep]
    data["zgrid"] = data["zgrid"][:, :keep + 1]
    data["nVertLevels"] = keep

    heights = np.array([lo3.get_level_height(data, k) for k in range(keep)])

    # Seam 5: updraft search at the design height.
    search_level = int(np.argmin(np.abs(heights - updraft_km)))

    def find_updraft_y(d, time_idx):
        return d["yCell"][np.argmax(d["w"][time_idx, :, search_level])]

    lo3.find_updraft_y = find_updraft_y

    tm = data["time_minutes"]
    snapshots = [int(np.argmin(np.abs(tm - t)))
                 for t in np.arange(every_min, tm[-1] + 1e-6, every_min)]
    level = int(np.argmin(np.abs(heights - level_km)))
    src = dict(w_threshold=cfg["w_threshold"], w_ref=cfg["w_ref"],
               z_min_m=cfg["z_min"], z_max_m=cfg["z_max"],
               source_rate_ppbv=cfg["source_rate"])
    outdir.mkdir(parents=True, exist_ok=True)
    base = str(outdir / "lightning_nox")

    def suffix(i):
        return f"_t{int(round(tm[i])):03d}min"

    no_vmax = max(lo3.to_ppbv(data["qNO"][i].max(), M_NO) for i in snapshots)
    no2_vmax = max(lo3.to_ppbv(data["qNO2"][i].max(), M_NO2) for i in snapshots)
    o3_anom = []
    for i in snapshots:
        o3 = lo3.to_ppbv(data["qO3"][i], 0.048)
        o3_anom.append(float(np.abs(o3 - np.median(o3)).max()))
    for i in snapshots:
        lo3.plot_vertical_cross_section(
            data, i, f"{base}_vertical{suffix(i)}.png",
            no_vmax=no_vmax, no2_vmax=no2_vmax,
            o3_anom_max=max(o3_anom + [1e-6]), **src)
        lo3.plot_species_comparison(
            data, level, i, f"{base}_comparison{suffix(i)}.png", **src)
        lo3.plot_updraft_profiles(
            data, i, f"{base}_profiles{suffix(i)}.png",
            z_min_m=0.0, z_max_m=profile_ztop_km * 1000.0)
    lo3.plot_time_evolution(data, f"{base}_evolution.png")
    lo3.plot_horizontal_evolution(data, level, f"{base}_horizontal.png")
    lo3.plot_source_profiles(data, f"{base}_source.png", **src)
    lo3.plot_photolysis(data, f"{base}_photolysis.png")
    return dict(snapshot_minutes=[float(tm[i]) for i in snapshots],
                horizontal_level=level,
                horizontal_level_km=float(heights[level]),
                updraft_search_level=search_level,
                updraft_search_level_km=float(heights[search_level]),
                levels_kept=keep, ztop_km=ztop_km,
                max_nox_above_ztop_ppbv=nox_above)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--every-min", type=float, default=30.0)
    ap.add_argument("--level-km", type=float, default=10.0)
    ap.add_argument("--ztop-km", type=float, default=20.0)
    ap.add_argument("--updraft-search-km", type=float, default=5.0)
    ap.add_argument("--profile-ztop-km", type=float, default=16.0)
    ap.add_argument("--max-nox-above-ztop-ppbv", type=float, default=0.01)
    ap.add_argument("--part", choices=("all", "diagnostics", "gallery"),
                    default="all")
    args = ap.parse_args()

    if args.part == "all":
        # style.setup() registers colormaps once per process, so each script
        # family renders in its own interpreter.
        import subprocess
        argv = [a for a in sys.argv[1:] if not a.startswith("--part")]
        for part in ("diagnostics", "gallery"):
            subprocess.run([sys.executable, __file__, *argv, "--part", part],
                           check=True)
        return 0

    output = args.run_dir / "output.nc"
    mesh = args.run_dir / "supercell_init.nc"
    sys.path.insert(0, str(REPO_SCRIPTS))

    with Dataset(output) as ds:
        cfg = {k: (v if isinstance(v, str) else float(v))
               for k, v in lnox_config(ds).items()}
        labels = time_labels(ds)
    if cfg["gating_mode"] != "altitude":
        raise SystemExit(f"gating mode {cfg['gating_mode']!r} is not altitude")
    print(f"LNOx config: {cfg}")
    print(f"records: {len(labels)} ({labels[0]} .. {labels[-1]})")

    if args.part == "diagnostics":
        info = render_diagnostics(
            output, mesh, args.outdir / "diagnostics", cfg, args.every_min,
            args.level_km, args.ztop_km, args.updraft_search_km,
            args.profile_ztop_km, args.max_nox_above_ztop_ppbv)
    else:
        info = dict(stems=[e["stem"] for e in
                           render_gallery(output, mesh, args.outdir / "gallery")])

    record = dict(run_dir=str(args.run_dir), lnox=cfg, records=labels,
                  part=args.part, info=info)
    (args.outdir / f"render_info_{args.part}.json").write_text(
        json.dumps(record, indent=1))
    print(json.dumps(info, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

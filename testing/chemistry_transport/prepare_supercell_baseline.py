#!/usr/bin/env python3
"""Stage the resolved supercell LNOx regression case without changing old runs."""
import argparse
from pathlib import Path
import shutil
import subprocess
import tarfile

MUSICA_REV = "1403e3d22717bc87f3bf9d0aa591caf039c92bbc"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--musica-source", type=Path, required=True)
    parser.add_argument("--duration-seconds", type=int, choices=(60, 600, 7200), default=7200,
                        help="60/600 for development smoke tests; 7200 for the baseline")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    run = args.run_dir.resolve()
    if any((run / name).exists() for name in ("namelist.atmosphere", "supercell_init.nc", "output.nc")):
        raise SystemExit("Use a fresh run directory; existing case inputs/results are preserved.")
    run.mkdir(parents=True, exist_ok=True)
    (run / "provenance").mkdir(exist_ok=True)
    for name in ("supercell_grid.nc", "supercell.graph.info.part.8"):
        shutil.copy2(args.input_dir / name, run / name)
    for path in (repo / "test_cases/supercell").iterdir():
        if path.is_file():
            shutil.copy2(path, run / path.name)
    for name in ("lnox_o3.yaml", "tuvx_no2.json", "tuvx_upper_atm.csv"):
        shutil.copy2(repo / "micm_configs" / name, run / name)
    for name in ("namelist.init_atmosphere", "namelist.atmosphere"):
        path = run / name
        text = path.read_text().replace("0000-01-01_00:00:00", "0000-01-01_18:00:00")
        if name == "namelist.atmosphere":
            text = text.split("! Enable chemistry", 1)[0]
            duration = args.duration_seconds
            duration_text = f"{duration // 3600:02d}:{duration // 60 % 60:02d}:{duration % 60:02d}"
            text = text.replace("config_run_duration = '02:00:00'", f"config_run_duration = '{duration_text}'")
            text = text.replace("    config_microp_scheme", "    config_pbl_scheme = 'off'\n    config_convection_scheme = 'off'\n    config_microp_scheme")
            text += """&chemistry
    config_micm_file = 'lnox_o3.yaml'
    config_chemistry_interval = 0.0
    config_chem_substeps = 1
    config_micm_relative_tolerance = 1.0e-6
/
&photolysis
    config_tuvx_config_file = 'tuvx_no2.json'
    config_tuvx_top_extension = .true.
    config_tuvx_upper_column_mode = 'legacy_static'
    config_tuvx_extension_file = 'tuvx_upper_atm.csv'
    config_tuvx_update_interval = 600.0
    config_j_no2_max = 0.01
    config_chemistry_latitude = 35.86
    config_chemistry_longitude = -97.93
/
&lnox
    config_lnox_gating_mode = 'altitude'
    config_lnox_source_rate = 0.5
    config_lnox_w_threshold = 5.0
    config_lnox_w_ref = 10.0
    config_lnox_z_min = 5000.0
    config_lnox_z_max = 12000.0
    config_lnox_nox_tau = 0.0
/
"""
        path.write_text(text)
    if args.duration_seconds < 7200:
        path = run / "streams.atmosphere"
        interval = args.duration_seconds // 2
        path.write_text(path.read_text().replace('output_interval="00:10:00"',
                        f'output_interval="00:{interval // 60:02d}:{interval % 60:02d}"'))
    fields = (run / "stream_list.atmosphere.output").read_text().splitlines()
    fields += ["j_jNO2", "areaCell", "indexToCellID", "zz", "rho_zz", "rainnc", "rainc"]
    (run / "stream_list.atmosphere.output").write_text("\n".join(dict.fromkeys(fields)) + "\n")
    archive = run / "provenance/tuvx-data.tar"
    subprocess.run(["git", "-C", str(args.musica_source), "archive", "--format=tar",
                    f"--output={archive}", MUSICA_REV, "configs/tuvx/data"], check=True)
    with tarfile.open(archive) as tar:
        tar.extractall(run / "provenance", filter="data")
    shutil.move(run / "provenance/configs/tuvx/data", run / "data")
    (run / "provenance/configs/tuvx").rmdir()
    (run / "provenance/configs").rmdir()
    for path in (repo / "src/core_atmosphere/physics/physics_wrf/files").iterdir():
        if path.is_file():
            shutil.copy2(path, run / path.name)
    for name in ("init_atmosphere_model", "atmosphere_model"):
        (run / name).symlink_to(Path("bin") / name)
    print(run)


if __name__ == "__main__":
    main()

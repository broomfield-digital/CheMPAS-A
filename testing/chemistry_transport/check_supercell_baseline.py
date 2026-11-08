#!/usr/bin/env python3
"""Validate and summarize the two-hour, eight-rank resolved supercell baseline.

This checks an existing chemistry/advection regression case. It does not qualify
PBL or cumulus gas transport, and reports reactive inventories separately from
transport conservation tests. All numeric output records are inspected.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from netCDF4 import Dataset, chartostring

MOLAR_MASS = {"qNO": 0.030, "qNO2": 0.046, "qO3": 0.048}
HOST_M_AIR = 0.0289644  # mpas_chem_constants.F; initializer historically uses 0.029.


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(var, record=None):
    if "Time" in var.dimensions:
        require(var.dimensions[0] == "Time", f"Unexpected axes: {var.name}")
        value = var[record]
    else:
        value = var[:]
    require(not np.any(np.ma.getmaskarray(value)), f"Masked output: {var.name}, record {record}")
    value = np.asarray(value)
    require(np.all(np.isfinite(value)), f"Nonfinite output: {var.name}, record {record}")
    return value


def summarize(run, duration=7200):
    require((run / "provenance/run-exit-status").read_text().strip() == "0", "MPI run failed")
    log = (run / "log.atmosphere.0000.out").read_text()
    require("Logging complete." in log, "Missing final log footer")
    require(re.search(r"MPI task count:\s+8\b", log), "Expected eight MPI ranks")
    for kind in ("Error", "Critical error"):
        require(re.search(rf"^\s*{kind} messages\s*=\s*0\s*$", log, re.M), f"Nonzero/missing {kind} count")
    require(not list(run.glob("*.abort")), "MPAS abort files found")
    for path in run.glob("log.atmosphere.*.err"):
        require(path.stat().st_size == 0, f"Rank error log: {path.name}")
    for scheme in ("pbl", "convection"):
        require(re.search(rf"config_{scheme}_scheme\s*=\s*off\b", log), f"{scheme} must be off")
    require("[LNOx] Lightning NOx source enabled:" in log, "LNOx is disabled")
    steps = len(re.findall(r"^\s*Begin timestep ", log, re.M))
    require(steps == duration // 3, f"Expected {duration // 3} model steps, found {steps}")
    recovery = [line.strip() for line in log.splitlines()
                if re.search(r"\[.*(?:MICM|MUSICA|CHEM).*\].*(?:failed|failure|recover|rollback)", line, re.I)]
    require(not recovery, f"Chemistry failure/recovery messages: {recovery[:3]}")
    summary = {
        "case": "supercell_lnox", "purpose": "resolved-advection and coupled-chemistry regression",
        "pbl_scheme": "off", "convection_scheme": "off", "mpi_ranks": 8,
        "dt_seconds": 3.0, "simulated_seconds": duration, "expected_steps": duration // 3,
        "all_numeric_records_finite_and_unmasked": True,
        "completed_steps": steps,
        "ppbv_dry_air_molar_mass_kg_mol": HOST_M_AIR,
        "initial_O3_kg_kg": 50.0 * 1e-9 * (0.048 / 0.029),
        "initial_O3_note": "Nominal 50 ppbv uses init_lnox_o3.py's 0.029 kg/mol air convention; reported ppbv uses the host's 0.0289644 kg/mol.",
        "startup_warnings": [line.strip() for line in log.splitlines() if line.startswith("WARNING:")],
        "rank_error_files": [], "chemistry_recovery_messages": recovery,
        "variables": {}, "records": [],
    }
    with Dataset(run / "output.nc") as ds:
        times = [str(t).strip() for t in chartostring(ds.variables["xtime"][:])]
        interval = 600 if duration == 7200 else duration // 2
        expected = [f"0000-01-01_{18 + seconds // 3600:02d}:{seconds // 60 % 60:02d}:{seconds % 60:02d}"
                    for seconds in range(0, duration + 1, interval)]
        require(times == expected, f"Incomplete/unexpected output times: {times}")
        summary["times"] = times
        summary["n_cells"] = len(ds.dimensions["nCells"])
        summary["n_levels"] = len(ds.dimensions["nVertLevels"])
        require((summary["n_cells"], summary["n_levels"]) == (28080, 60), "Unexpected mesh")
        required = set(MOLAR_MASS) | {"j_jNO2", "rho", "rho_zz", "zz", "zgrid", "areaCell", "indexToCellID", "w"}
        require(required <= ds.variables.keys(), f"Missing fields: {required - ds.variables.keys()}")
        for name, var in ds.variables.items():
            if var.dtype.kind not in "iuf":
                continue
            digest = hashlib.sha256()
            vmin, vmax = float("inf"), float("-inf")
            for record in range(len(times)) if "Time" in var.dimensions else [None]:
                value = read(var, record)
                digest.update(np.ascontiguousarray(value).tobytes())
                vmin = min(vmin, float(value.min()))
                vmax = max(vmax, float(value.max()))
            summary["variables"][name] = {"min": vmin, "max": vmax, "dtype": str(var.dtype),
                                           "shape": list(var.shape), "values_sha256": digest.hexdigest()}
        area = read(ds.variables["areaCell"], 0)
        zgrid = read(ds.variables["zgrid"], 0)
        dz = np.diff(zgrid, axis=-1)
        require(np.all(area > 0) and np.all(dz > 0), "Nonpositive area/layer thickness")
        ids = read(ds.variables["indexToCellID"], 0)
        require(len(np.unique(ids)) == len(ids), "Repeated cell IDs/halo cells in output")
        zmid = 0.5 * (zgrid[:, :-1] + zgrid[:, 1:])
        zz = read(ds.variables["zz"], 0)
        density_difference = 0.0
        for i, time in enumerate(times):
            rho = read(ds.variables["rho"], i)
            require(np.all(rho > 0), "Nonpositive dry density")
            rho_state = zz * read(ds.variables["rho_zz"], i)
            density_difference = max(density_difference, float(np.max(np.abs(rho - rho_state) / rho)))
            dry_mass = area[:, None] * dz * rho
            w = read(ds.variables["w"], i)
            wmid = 0.5 * (w[:, :-1] + w[:, 1:])
            record = {"xtime": time, "dry_air_mass_kg": float(dry_mass.sum()),
                      "w_max_m_s": float(w.max()), "w_min_m_s": float(w.min()),
                      "lightning_gate_cell_levels_at_output": int(np.count_nonzero(
                          (wmid > 5) & (zmid >= 5000) & (zmid <= 12000))),
                      "species": {}}
            moles = {}
            for name, molar_mass in MOLAR_MASS.items():
                q = read(ds.variables[name], i)
                floor = -(1e-20 + 1e-12 * summary["variables"][name]["max"])
                require(float(q.min()) >= floor, f"Negative {name} at {time}: {q.min()}")
                if i == 0:
                    value = 50.0 * 1e-9 * (0.048 / 0.029) if name == "qO3" else 0.0
                    require(np.all(q == value), f"Unexpected initial {name}")
                mass = float(np.sum(dry_mass * q))
                moles[name] = mass / molar_mass
                record["species"][name] = {"mass_kg": mass, "min_kg_kg": float(q.min()),
                                            "max_kg_kg": float(q.max()),
                                            "max_ppbv": float(q.max() * HOST_M_AIR / molar_mass * 1e9)}
            record["NOx_mol"] = moles["qNO"] + moles["qNO2"]
            record["Ox_mol"] = moles["qO3"] + moles["qNO2"]
            summary["records"].append(record)
        require(density_difference < 1e-12, f"Dry density mismatch: {density_difference}")
        summary["max_rho_vs_zz_rho_zz_relative_difference"] = density_difference
        require(summary["variables"]["j_jNO2"]["max"] > 0, "No photolysis activity")
        require(summary["variables"]["j_jNO2"]["min"] >= 0, "Negative photolysis rate")
        # MPAS omits inactive-package diagnostics when cumulus is off.
        if "rainc" in summary["variables"]:
            require(summary["variables"]["rainc"]["min"] == summary["variables"]["rainc"]["max"] == 0,
                    "Unexpected parameterized convective precipitation")
        if duration == 7200:
            for name in ("qNO", "qNO2"):
                require(summary["records"][-1]["species"][name]["mass_kg"] > 0, f"No final {name}")
            require(any(r["lightning_gate_cell_levels_at_output"] > 0 for r in summary["records"]),
                    "No active lightning source cells at output times")
    for quantity in ("dry_air_mass_kg", "Ox_mol"):
        initial = summary["records"][0][quantity]
        summary[f"{quantity}_max_relative_drift"] = max(abs(r[quantity] - initial) / initial for r in summary["records"])
    summary["inventory_notes"] = (
        "Inventories use areaCell * layer thickness * dry rho. NOx = NO + NO2 in moles includes lightning production; "
        "Ox = O3 + NO2 in moles is invariant for these two chemical reactions and the NO-only source. "
        "Reported drifts characterize this coupled baseline; no isolated PBL/convection budget is claimed. "
        "Lightning activation is sampled at output times, not integrated over all steps. "
        "Optimized MPAS writes a regular log only for rank 0; all rank error/abort files and the MPI exit status are checked."
    )
    summary["status"] = "pass"
    summary["qualification"] = "coupled regression baseline" if duration == 7200 else "development smoke only"
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--duration-seconds", type=int, choices=(60, 600, 7200), default=7200)
    parser.add_argument("--reference", type=Path,
                        help="Require exact saved-field agreement with a same-toolchain summary")
    args = parser.parse_args()
    summary = summarize(args.run_dir, args.duration_seconds)
    if args.reference:
        reference = json.loads(args.reference.read_text())
        require(summary["times"] == reference["times"], "Reference output times differ")
        differences = [name for name, value in reference["variables"].items()
                       if summary["variables"].get(name) != value]
        require(not differences, f"Reference fields differ: {', '.join(differences)}")
        summary["reference_comparison"] = {"path": str(args.reference), "result": "exact saved-field match",
                                            "fields_compared": len(reference["variables"])}
    path = args.output or args.run_dir / "summary.json"
    path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(f"PASS: {len(summary['times'])} records, {args.duration_seconds} simulated seconds, eight ranks; {path}")
    print(f"Maximum relative drift: dry air={summary['dry_air_mass_kg_max_relative_drift']:.3e}, Ox={summary['Ox_mol_max_relative_drift']:.3e}")


if __name__ == "__main__":
    main()

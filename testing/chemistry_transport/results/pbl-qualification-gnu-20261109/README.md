# MICM gas PBL transport qualification, GNU Fortran, 2026-11-09

Status: **pass**. All 15 checks passed and none skipped a part. Neither the
source nor any executable changed during the run.

This is the final qualification of MICM gas PBL transport (issue #40), P5 and
P6 of [the plan](../../../../docs/plans/pbl-transport.md).
`testing/chemistry_transport/run_pbl_qualification.sh` ran it on Orion from
06:37 to 10:29 UTC (3 h 52 min).

## Build and inputs

- Source: `chempas-a-pbl-complete` at `bf3f1abf`, with no tracked changes. No
  file under `src/` changed after `72fb7399`, the last model change.
- Build: `MUSICA=true`, double precision (`-fdefault-real-8`), `OPENMP=false`.
  Model SHA-256 `2eaacb4753ea…` (full hash in `manifest.json`).
- Toolchain: GNU Fortran 15.2.0 (conda-forge), Open MPI 5.0.10, NetCDF-C
  4.10.1, NetCDF-Fortran 4.6.3, PIO 2.6.9, MUSICA-Fortran 0.16.5 (`1403e3d2`),
  MICM `bb57684a`.
- Ranks: eight, except where a check compares other counts.
- Baseline: `supercell_lnox/20261107T201901Z-2c6ff111` (supercell init
  SHA-256 `b79a002abe68…`).
- Fine grid: 100 layers, 50 m up to 1.5 km, prepared by `prepare_fine_grid.sh`
  at `f62f5423`. The init hashes are in `manifest.json`.

The "transport off is unchanged" comparisons use these reference executables:

| Name | Built from | SHA-256 |
| --- | --- | --- |
| pre_p3 | `50f2666c`, before YSU transport | `2d34e9c905a3…` |
| pre_p4 | `b6373312`, before MYNN transport, with the MYNN snow-content fix | `9e69a9ad341b…` |
| pre_p5 | `872e9255`, the last P4 model change | `4d959f018971…` |

## Results

| Check | Minutes | Result |
| --- | --- | --- |
| python-unit | 0 | 42 checker unit tests pass. The reactive, fine-grid, level and record tools pass their self-tests |
| fortran-unit | 1 | 12 Fortran unit programs pass: kernels, operators, adapters and wrapper exports |
| mmm-patches | 0 | The MMM-physics patches apply, repeat as a no-op, and revert cleanly |
| wiring | 0.3 | The LNOx gases reach physics; host-bound H2O is excluded |
| realized-increment | 10 | In every advection mode the realized increment matches the requested one (largest difference 5.3e-23; column-mass error ≤ 6.2e-16) |
| ysu | 21 | The gates stop the model. With transport, non-gas fields are bitwise identical to transport off, and gas mass drifts ≤ 4.1e-16. Convective mixes more than stable. Transport off is bitwise identical to pre_p3 |
| mynn | 37 | The same for local diffusion, updrafts, and updrafts with downdrafts. Transport off is bitwise identical to pre_p4 |
| cadence | 3 | The PBL and surface schemes run 20 times in 20 steps with `'none'` or an interval of `config_dt`, and 7 times with 9 s (DT_PBL = 9 s). An interval of `config_dt` is bitwise identical to `'none'`. The default cadence is bitwise identical to pre_p5 |
| intervals | 4 | YSU and MYNN at 9 s, YSU at 10 s: non-gas fields bitwise, gas mass conserved, budget max \|rel\| 3.3e-16. PBL off with transport requested runs and logs transport as off |
| restart | 15 | Five restarts on a PBL-interval boundary are bitwise identical to the continuous runs (YSU and MYNN, `'none'` and 9 s, and native MPAS). The restart off the boundary differs, as documented |
| chemistry-intervals | 8 | YSU and MYNN with chemistry every 3, 9 and 21 s: non-gas fields bitwise, nitrogen and odd-oxygen drift ≤ 4.9e-16, budget max \|rel\| 1.1e-15 |
| mpi | 14 | YSU and MYNN on 1, 2 and 4 ranks against 8: every output field bitwise identical |
| reactive | 50 | LNOx with NO2 photolysis, YSU and MYNN: nitrogen and odd-oxygen drift ≤ 4.9e-16. With lightning (YSU, 30 min): odd-oxygen drift 1.2e-16, budget over 600 steps |
| fine-grid | 68 | The PBL spans 15 (YSU) and 18 (MYNN) layers. NO mixes at 0.81 (YSU) and 1.04 (MYNN) of an offline replay with the scheme's own `kzh`. MYNN updrafts are active in 98.9% of columns |
| smoke-60s | 0.5 | 60 s of the standard LNOx supercell on eight ranks; dry-air drift 0, odd-oxygen drift 1.2e-16 |

Notes:

- `logs/restart.log` has one `FAIL` line. It belongs to the restart off the
  PBL-interval boundary (`ysu-9-offb`), which is reported, not required to
  match. The next line, a `NOTE`, marks it as the expected difference.
- The MPI check allows atol 1e-18 and rtol 1e-10. Every field was bitwise
  identical anyway.
- These are numerical checks of conservation, mixing and reproducibility. They
  do not validate the chemistry or the lightning source.
- Ten minutes on the fine grid is shorter than the PBL's mixing time, so NO is
  not yet well mixed there. The plan explains the replay criterion.

## Files

- `manifest.json`: source, executables and their hashes, toolchain, baseline,
  fine-grid inputs, start and end times, and each check's result.
- `summary.json`: each check's result, elapsed seconds and log path in the run
  directory.
- `logs/`: the 15 check logs, copied from the run directory.

The model output stays in `~/Data/CheMPAS-A/p6-dev/qual2` on Orion.

## Rerun

```bash
CHEMPAS_FINE_GRID=FINE_GRID_DIR \
CHEMPAS_PRE_P3_MODEL=PRE_P3_MODEL \
CHEMPAS_PRE_P4_MODEL=PRE_P4_MODEL \
CHEMPAS_PRE_P5_MODEL=PRE_P5_MODEL \
testing/chemistry_transport/run_pbl_qualification.sh BASELINE_DIR NEW_RUN_ROOT [MODEL]
```

Without `CHEMPAS_FINE_GRID` the runner prepares the fine grid itself. Without a
reference executable, that comparison is skipped and the status is
"incomplete".

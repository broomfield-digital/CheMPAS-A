# PBL transport of MICM gases

CheMPAS-A can mix the gases of a MICM mechanism vertically with the planetary
boundary layer (PBL) scheme. The gases then get the same turbulent mixing that
the scheme gives heat and moisture. The option is off by default. The
development plan, with every test and result, is in
[docs/plans/pbl-transport.md](../../plans/pbl-transport.md).

## What it does

When the PBL scheme runs, it computes how strongly each column mixes. CheMPAS-A
keeps that operator and applies it to every transported MICM gas once per model
step, after chemistry, as a closed-column, mass-conserving update:

- **YSU**: diffusion with YSU's scalar exchange coefficient, including its
  entrainment at the PBL top and its local mixing above the PBL.
- **MYNN**: diffusion with MYNN's scalar exchange coefficient, plus its
  mass-flux updrafts (`config_mynn_edmf = 1`) and downdrafts
  (`config_mynn_edmf_dd = 1`) when those are on.

Water vapour and the other moisture fields keep their own path; only MICM gases
are affected. Host-bound species (H2O) and prescribed species are left alone.

## Turning it on

In `&physics` of `namelist.atmosphere`:

```text
config_pbl_scheme         = 'bl_ysu'      ! or 'bl_mynn'
config_sfclayer_scheme    = 'sf_monin_obukhov'   ! sf_mynn with MYNN
config_convection_scheme  = 'off'
config_micm_pbl_transport = .true.
```

The model log then reports, during initialization:

```text
[transport] PBL: gas operator refreshed at each PBL call (config_pbl_interval = none) and applied every model step
[transport] nMicmTransport = 3; PBL gas transport enabled = T; ...
```

## Supported configurations

| Setting | Supported |
| --- | --- |
| Build | `MUSICA=true`, double precision, GNU Fortran; one OpenMP thread |
| PBL scheme | `bl_ysu` (with or without `config_ysu_pblmix`), `bl_mynn` (local diffusion; updrafts; updrafts and downdrafts), or `off` (no transport) |
| Surface layer | Required with an active PBL scheme (MPAS's own rule) |
| `config_pbl_interval` | `'none'`, or at least `config_dt` (3, 6, 9 and 10 s tested with 3-s steps). Between PBL calls the last operator is applied every step. A shorter interval is not supported: MPAS then runs the PBL and surface schemes every step with the shorter interval as their time step |
| `config_convection_scheme` | `'off'` only (convection-resolving configurations) |
| `config_chemistry_interval` | Any |
| Restart | Supported; the operator is rebuilt at the restart's first step. A restart on a PBL-interval boundary continues bitwise. Off a boundary, the PBL calls shift to the restart time, as the meteorology's do |
| MPI | Any rank count |

The qualification (`testing/chemistry_transport/run_pbl_qualification.sh`) runs
the idealized supercell on eight ranks with GNU Fortran. A Clang/Flang build
has not been qualified yet.

## What stops the model

With `config_micm_pbl_transport = .true.`, initialization stops with a
`[transport]` message for:

- a PBL scheme without a gas adapter;
- an active convection scheme, including Grell-Freitas, which the
  `convection_permitting` suite selects. Parameterized convection does not yet
  move MICM gases (issue #41), so the gases would stay trapped below cloud
  base. Set `config_convection_scheme = 'off'` for a convection-resolving run;
- more than one OpenMP thread.

During the run, the model stops if the MYNN mass-flux operator needs more than
10,000 substeps (implausible plume data). It also stops if the operator was not
refreshed within `config_pbl_interval`, which would mean a missed PBL call.

Without MICM gases (a `MUSICA=false` build, or a mechanism with no transported
gas) the option does nothing and the log says so.

## Budget diagnostics

`config_micm_transport_budget = .true.` (in `&physics`) writes one line per gas
and model step to `log.atmosphere.0000.out`, reduced over all ranks:

```text
[transport] budget PBL t=<s> gas=<n> scalar=<index> before=<kg> after=<kg> increment=<kg> rel=<increment/before>
    residual=<largest column residual> qmin=<kg/kg> qmax=<kg/kg> columns=<active columns> refresh=<s> dt=<s>
```

Each line is written on one line in the log. The masses are dry-air weighted
totals. `t` and `refresh` are seconds since this run's start, so they count
from the restart time in a restarted run. `refresh` is the time of the PBL call
whose operator was applied. `rel` and `residual` should be at round-off (about
1e-15). The checker below also fails any value that is not finite, and any
increment that does not match `after - before`. Check a log with:

```bash
python testing/chemistry_transport/check_transport_budget.py log.atmosphere.0000.out --dt 3 --interval 3
```

The option adds MPI reductions every step; leave it off for production runs.

## What to expect

On the supercell grid (interfaces at 0, 299, 712 and 1182 m), the PBL spans
only one or two layers. Most of the effect in a deep storm then comes from
YSU's mixing above the PBL top. YSU applies its Richardson-number diffusivity
there, as it does for heat and moisture. In the standard LNOx case most of the
NOx difference made by transport lies at 8–12 km, in the lightning plume, and
under 1% lies below 1 km. A grid with finer layers near the surface resolves
the mixed layer itself; see the fine-grid results in the plan.

## Finer vertical grids

The supercell namelist sets an explicit vertical eddy viscosity,
`config_v_mom_eddy_visc2 = config_v_theta_eddy_visc2 = 500.0` (m²/s), in
`&nhyd_model`. MPAS applies it as a tendency computed once per step, which is
stable only while K·dt/dz² stays below 0.5. With 50-m layers and
`config_dt = 3` that number is 0.6, and the run blows up within two minutes,
with or without gas transport. With a PBL scheme on such a grid, set both
viscosities to 0.0 and let the PBL scheme mix the boundary layer, or shorten
`config_dt`. `testing/chemistry_transport/prepare_fine_grid.sh` builds the
qualification's fine grid, with 50-m layers up to 1.5 km.

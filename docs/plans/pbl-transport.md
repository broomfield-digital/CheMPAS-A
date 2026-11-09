# Implementation plan: PBL transport of MICM tracers

Status: the P1 shared foundation (`0e2f15f6`), P2, the column diffusion kernel
and tendency proof (`7ae5f39d` shared, `0eb2c8ad`), P3, YSU gas transport at the
default PBL interval (see [P3](#p3-ysu-gas-pbl-transport-implemented)), and P4,
MYNN gas transport with its local diffusion, updrafts and downdrafts (see
[P4](#p4-mynn-gas-pbl-transport-implemented)), are implemented. So are P5, PBL
intervals, restart, chemistry intervals, MPI and the gate against parameterized
convection (see the [P5 plan](#p5-plan-pbl-intervals-restart-and-mpi)), and P6,
the GNU Fortran eight-rank qualification (see
[P6](#p6-gnu-fortran-eight-rank-qualification-implemented)). The work completes
for convection-resolving configurations, with parameterized convection off, the
first science target being the resolved supercell, when the final qualification
run is recorded in P6. Qualification with
parameterized convection belongs to the convective transport plan (#41); the
Clang/Flang build is deferred. Users start from the
[user guide](../chempas/guides/pbl-gas-transport.md).

- Tracking issue: [CheMPAS-A #40](https://github.com/NCAR/CheMPAS-A/issues/40).
- Branch: `chempas-a-pbl`, from `chempas-a` at
  `89008ed50ce8f9ffd10876deed7c9a48ec573a38`. P5 Tasks 3–7 and P6 are on
  `chempas-a-pbl-complete`, from `chempas-a-pbl` at `a267582b` (tag
  `pbl-checkpoint-2026-11-08`).
- Companion: [convective transport plan](https://github.com/NCAR/CheMPAS-A/blob/chempas-a-convect/docs/plans/convective-transport.md).
- Coordination: [MUSICA #1073](https://github.com/NCAR/musica/issues/1073).
- Source/issue audit: 2026-10-18. These are development plans, not supported
  capabilities. Standalone MICM contribution #39 is a separate effort.

Implement vertical gas transport for every prognostic MICM tracer discovered at
runtime, retaining the existing native MPAS moisture transport. Emissions,
lightning sources, deposition, scavenging, and aerosol transport are excluded
from this qualification. Keep those sources disabled in the transport fixtures.

## Findings at the branch base

| Location | Finding and consequence |
| --- | --- |
| `chemistry/musica/mpas_musica.F`: `species_entry_t`, `resolve_mpas_indices`, `host_binding` | `chem_species(:)%mpas_index` is the authoritative mapping. H2O is bound to native `qv`; prescribed O2/N2 can also use `qv` as an index anchor. An index alone does not identify an independently transported gas. |
| `physics/physics_wrf/module_bl_ysu.F`: `ysu` | The wrapper declares `nmix = 0`; its existing pack/unpack loops cannot transport MICM species. |
| MMM-physics `bl_ysu.F90`: `bl_ysu_run` | A generic `qmix`/`qmixtnp` loop already exists. It uses the moisture diffusion matrix, with its own `xkzq` and entrainment treatment. It does not copy the heat/water nonlocal source terms to passive gases. The solve uses `dt2 = 2*dt`, with tendencies divided by `dt2`. |
| `physics/mpas_atmphys_driver_pbl.F`: `pbl_to_MPAS` | The exported `kzq` is assigned `kzh_p`. It is not sufficient evidence of the scalar operator actually used inside YSU. |
| `physics/physics_wrf/module_bl_mynn.F` and MMM-physics `bl_mynn.F90` | Chemistry arguments are behind `WRF_CHEM`; the MPAS call supplies none. The dormant wrapper has inconsistent chemistry identifiers. The external chemistry writeback is also guarded by `WRF`, with a positive concentration floor. Enabling a macro is not a complete adapter. |
| MMM-physics `bl_mynn_subroutines.F90`: `mynn_mix_chem`, `DMP_mf` | Existing code combines diffusion and plume tracer fluxes. Its chemistry solve includes deposition arguments and a prescribed top value. Closed-column gas transport requires explicit zero-flux boundaries and no artificial tracer floor. |
| `physics/mpas_atmphys_interface.F`: `MPAS_to_physics` | `rho_p = zz*rho_zz*(1+qv)` includes water vapor; chemistry obtains dry density as `zz*rho_zz`. Moist pressure/density weights cannot be assumed to conserve dry-air tracer mass. |
| `physics/mpas_atmphys_manager.F` and `mpas_atmphys_driver.F` | `config_pbl_interval` creates an alarm and sets `dt_pbl`, but the driver is called every model step when PBL and surface layer are enabled. No PBL alarm check gates that call at this base. |
| `mpas_atm_core.F`: `atm_do_timestep` | Order is physics, chemistry on scalar time level 1, then dynamics. Tracer tendencies computed inside the earlier physics call can already be stale when dynamics uses them. |
| `physics/mpas_atmphys_todynamics.F`: `physics_get_tend_work` | Only fixed native scalar indices receive PBL tendencies, multiplied by `rho_zz` from time level 2. |
| `dynamics/mpas_atm_time_integration.F` | Ordinary and monotonic scalar advances handle source tendencies differently; the monotonic path applies them to old state and clears them. Verify both the mass time level and single application before enabling a transport mode. |

Paths in this table are relative to `src/core_atmosphere/`. The external source
audited was [NCAR/MMM-physics at a4baf7f3](https://github.com/NCAR/MMM-physics/tree/a4baf7f3243d1db0dbc5f63473f895bdbdc05c30),
selected by `Externals.cfg` tag `20250616-MPASv8.3`.

## Shared interface with the convection branch

Implement this foundation first in a separate commit on `chempas-a-pbl`; reuse
that commit on `chempas-a-convect`. Do not create two independent species maps,
budget implementations, or tendency application loops. The following names are
proposed interfaces, not existing routines or namelist options.

1. Build a compact `transport_to_scalar(:)` map after MICM resolves its indices.
   Include every non-host-bound, writable prognostic MICM species, regardless of
   mechanism name, ordering, or scalar contiguity. Validate range and uniqueness;
   fail on an unexpected collision with a native field. Log each species,
   scalar index, and exclusion reason. Do not filter species by initial value,
   chemical lifetime, or a fixed list. An empty map is a valid no-op.
2. Store the map in block-owned host data. Physics accepts indices and ordinary
   Fortran arrays, without importing MICM solver objects or `mpas_musica`.
   Populate the host data from chemistry initialization and supply an empty map
   without MUSICA. Preserve the current one-block-per-rank chemistry restriction.
   This avoids a physics/chemistry build cycle: these libraries currently build
   as independent prerequisites.
3. Add a physics-side helper, provisionally `mpas_atmphys_tracer_transport.F`,
   with compact buffers `q(nTransport,nVertLevels,nOwnedCells)` and separate PBL
   and convection tendency buffers. Register only needed persistent cache and
   diagnostic fields through existing pools/runtime registration. Allocate outside
   threaded loops; keep tile workspaces private and clear every inactive column.
   Transport all host vertical levels, including any outside a restricted MICM
   solve domain. Bound workspace memory with tiles/species chunks; do not rerun
   the meteorological scheme for each species.
4. Units are `q = kg species / kg dry air`, `R = dq/dt` in `s-1`, dry layer mass
   per area `m_d = rho_dry*dz`, and interface tracer flux in
   `kg species m-2 s-1`. For an upward-positive flux,
   `m_d(k)*R(k) = F(k-1/2) - F(k+1/2)`. The domain species budget is
   `sum_owned(areaCell*m_d*q)`. Establish the exact equivalence to the MPAS
   `rho_zz`, `zz`, and vertical metrics at the application time level.
5. Cache meteorological transport coefficients/plume properties during physics.
   Construct gas tendencies later, at the once-per-model-step tendency assembly
   after chemistry, using current authoritative `scalars` and actual model `dt`.
   Carry `dt` explicitly to the helper. Cache coefficients, not tracer-dependent
   increments that would be replayed after chemistry or another transport update.
6. Use a scratch column for the deterministic sequence PBL then convection:
   `q1 = PBL(q0,dt)`, `q2 = CU(q1,dt)`,
   `R_pbl = (q1-q0)/dt`, `R_cu = (q2-q1)/dt`.
   Add the sum once at each mapped scalar index in `physics_get_tend`, after
   tendency clearing, using the proven MPAS mass conversion. Do not also write
   these trial updates to the prognostic state. The normal scalar advance owns
   state mutation; native moisture tendencies retain their existing path.
7. Prove the resulting full-step increment through ordinary/monotonic transport,
   RK stages, and dynamics substeps. Include fresh post-chemistry state, chemistry
   skipped intervals, and the existing scalar halo exchanges. If the advection
   source treatment cannot preserve the required positivity/conservation, fix
   that treatment for mapped gas sources before declaring that mode supported;
   do not hide failures with post-update clipping.
8. Persist any held operator coefficients, active masks, and validity times needed
   across restart. Per-step gas tendencies are reconstructed and are not replayed.
   Fail explicitly on missing/incompatible enabled cache state unless a tested
   deterministic reconstruction is available. Free map/cache/workspace data at
   finalize and test repeated initialization.

Introduce `config_micm_pbl_transport = .false.` as an independent opt-in control.
Its disabled path must preserve the baseline. `MUSICA=false` remains a no-op.
An empty `config_micm_file` is not one: in MUSICA builds it already aborts in
`musica_init` at the branch base, independently of transport. Resolve `suite`
before checking support; log the
resolved scheme and transport capability. An enabled, unsupported combination
must stop during initialization with an actionable message. `pbl_scheme='off'`
is a supported zero-transport case; an active PBL with surface layer off is not.

## Scheme implementation stages

| Stage/configuration | Planned treatment | Gate before claiming support |
| --- | --- | --- |
| YSU, default PBL interval | Replace the zero tracer interface with runtime-sized inputs/outputs. Expose the actual scalar exchange/entrainment operator so it can be evaluated on fresh gas state, independent of the earlier meteorology solve. Use conservative dry-air weights and closed boundaries. | Resolve the `2*dt` convention, top/interface indexing, and moist-pressure versus dry-mass weighting with column reference tests. |
| YSU stable/unstable and top-down mixing | Retain the selected scheme's scalar diffusivity and entrainment response. Document precisely which passive-gas nonlocal terms exist; zero surface tracer flux does not justify copying thermal/moisture source terms. | Compare normal and `config_ysu_pblmix` configurations against the exposed scheme operator; any omitted gas entrainment/nonlocal term must be addressed or the configuration excluded. |
| MYNN with `config_mynn_edmf=0`, `config_mynn_edmf_dd=0` | Add a clean MPAS passive-gas interface using the scheme's scalar diffusivity and an implicit dry-mass flux-divergence solve. | Uniform, zero, pulse, sharp-gradient, and nonuniform-density columns; zero flux at both external boundaries. |
| MYNN with updraft EDMF | Export/reuse plume geometry and entrainment coefficients; recompute plume concentrations from current gas profiles, include plume tracer flux and compensating environment transport together with diffusion. | Updraft-only and combined diffusion/EDMF budgets and nonnegativity. `edmf_a`, `edmf_w`, or `exch_h` alone are not the complete operator. |
| MYNN downdraft EDMF | Extend the adapter through the downdraft pathway and its environmental compensation. | Remains explicitly unsupported until downdraft tracer budgets pass; do not silently use the updraft-only operator. |
| Nondefault PBL intervals | Diagnose and repair the alarm/`dt_pbl` contract as a separately reviewable change before qualification. Specify coefficient refresh times and per-step gas application. | Initially reject nondefault intervals when gas PBL transport is enabled. Lift the restriction only after interval/restart tests and native-meteorology regressions establish the intended cadence. Disabled gas transport must not acquire an incidental cadence change. Lifted in P5 (`d2ac0898`): the cadence is repaired (`d90b9fd7`), and the operator from the last PBL call is applied every step and stamped, with interval, restart, chemistry-interval and MPI checks. |
| Parameterized convection active | Not qualified here. The native convection schemes do not carry MICM gases, so PBL mixing alone would leave them trapped below cloud base. The PBL-then-convection tendency sequence of the shared interface stays ready for #41. | Reject gas PBL transport with any active convection scheme, including scale-aware Grell-Freitas in the `convection_permitting` suite. The convective transport plan (#41) lifts the restriction after its coupled PBL + convection + MICM qualification. |

Do not enable WRF chemistry/emissions macros globally. Adapt passive-gas
algorithms with no deposition velocity, fire enhancement, or imposed positive
floor. Use a monotone implicit diffusion solve and conservative flux limiting
or substeps for nonlocal terms. Any roundoff repair must be conservative and
reported; a zero species must remain zero.

Required MMM-physics exports are an external dependency. Prefer a versioned
MMM-physics change with a recorded immutable revision in the external manifest.
If a temporary patch is needed, check it into CheMPAS-A with source hash checks,
deterministic application, and retained license notices. An edited ignored
`physics_mmm` checkout is not a reproducible implementation deliverable.

## Commit-sized work and tracking

- [x] Audit #40, branch base, native ownership, scheme hooks, and time ordering.
- [x] Record this branch plan and link the companion effort.
- [x] P1: shared map, empty/no-MUSICA path, block lifecycle, capability checks,
  dry-mass budget helper, and diagnostic schema. Commit `0e2f15f6`; its API,
  Makefile dependencies and evidence are under
  [P1 shared foundation](#p1-shared-foundation-implemented).
- [x] P2: standalone conservative column kernels and independent reference
  fixtures; prove `rho_zz` conversion and tendency consumption for each supported
  advection mode before connecting active mixing. Commits `7ae5f39d` (shared) and
  `0eb2c8ad`; see [P2](#p2-column-diffusion-kernel-and-tendency-proof-implemented).
- [x] P3: YSU runtime interface/operator exports, generic tendency assembly,
  default-interval integration, and stable/unstable/top-down qualification.
  Commits `85f68008` (MMM-physics export patch) and `c1bd210a`; see
  [P3](#p3-ysu-gas-pbl-transport-implemented).
- [x] P4: MYNN local diffusion, then updraft EDMF, then downdraft EDMF in separate
  commits. Keep unqualified combinations rejected at initialization. Commits
  `a68d7598`, `c5e144da`, `9eca3f99`, `0c21a646`, `2181961a`, `872e9255` and
  `fd9ccbe8`; see
  [P4](#p4-mynn-gas-pbl-transport-implemented).
- [x] P5: nondefault interval contract, cache restart behavior, chemistry
  intervals, MPI/halo checks, and the gate against parameterized convection.
  Commits `d90b9fd7` (native PBL cadence repair), `25d09c38`, `f62f5423`,
  `d2ac0898`, `72fb7399`, `a9a7d2c7`, `d4e7c954` and `bdf8aea2`; see the
  [P5 plan](#p5-plan-pbl-intervals-restart-and-mpi).
- [ ] P6: GCC/gfortran eight-rank results, supported-configuration table
  (convection-resolving configurations), reproducible fixtures, and user
  configuration documentation. The Clang/Flang build is deferred (decided
  2026-11-08): development and the first science runs are on Linux with
  gfortran. Flang qualification, including the absent-optional pattern, follows
  when a Flang platform is in use.

Coupled PBL + parameterized convection + MICM qualification is not a completion
criterion of this plan. It belongs to the convective transport plan (#41), after
its stages, because that is where the convective gas operator is added.

Likely files: `Registry.xml`; chemistry initialization/runtime registration in
`chemistry/musica/mpas_atm_chemistry.F`; the new physics helper; physics driver,
interface, manager, control, initialization/finalization and tendency assembly;
YSU/MYNN wrappers; necessary MMM exports; and affected Makefiles. Modify core
time integration or halo groups only where the verified application contract
requires it. Keep shared foundation commits separate from scheme commits.

## P1 shared foundation (implemented)

P1 is commit `0e2f15f6` on `chempas-a-pbl`. Reuse it unchanged as C1 on
`chempas-a-convect`; each branch then replaces only the body of its own process
module. With no process enabled, P1 changes no model output.

| File (under `src/core_atmosphere/`) | Role |
| --- | --- |
| `chemistry/musica/mpas_chem_transport_map.F` | Dependency-free map kernel: every non-host-bound, write-back MICM species in mechanism order. Rejects an index out of range, outside the runtime `chemistry` constituent group (`chemistry_start:chemistry_end`), or repeated. |
| `chemistry/musica/mpas_musica.F` | `musica_build_transport_map` applies the kernel to `chem_species` and logs each species' scalar index or exclusion reason. |
| `chemistry/musica/mpas_atm_chemistry.F` | `chemistry_init` publishes `nMicmTransport` and `micmTransportIndex` to `block % dimensions` after `resolve_mpas_indices`. |
| `physics/mpas_atmphys_tracer_column.F` | Pure closed-column kernels (dry mass `zz*rho_zz*dz`, flux divergence, relative residual), the `transport_budget_t` schema, and pure capability and cache checks. Two wrappers log and stop the model: `tracer_transport_check_capability(process, enabled, scheme, supported, reason)` and `transport_cache_validate(meta, signature, restart_has_cache, enabled)`. |
| `physics/mpas_atmphys_pbl_gas.F`, `physics/mpas_atmphys_cu_gas.F` | Identity stubs: `*_init(block, nTransport)`, `*_enabled()`, `*_active(iCell)`, `*_apply(iCell, dt, m_d, q)`, `*_finalize()`. |
| `physics/mpas_atmphys_tracer_transport.F` | Map and per-rank budgets. Runs `q1 = PBL(q0)`, `q2 = CU(q1)` on scratch columns, and `tracer_transport_add_tend` adds `(R_pbl + R_cu)*rho_zz` once per mapped scalar. Initializes the processes before sizing buffers; the `(nTransport, nVertLevels, nCellsSolve)` tendency buffers exist only when a process is enabled. |
| `physics/mpas_atmphys_todynamics.F` | `physics_get_tend(..., dt)` calls assembly and `add_tend` only when `tracer_transport_enabled()`. |
| `dynamics/mpas_atm_time_integration.F` | `atm_srk3` passes the model `dt`. |
| `mpas_atm_core.F` | `tracer_transport_init` after `chemistry_init`; `tracer_transport_finalize` before the physics finalize. |

Makefile dependencies added in P1:

```makefile
# src/core_atmosphere/physics/Makefile (OBJS += the four new objects)
mpas_atmphys_todynamics.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_tracer_transport.o \
	mpas_atmphys_vars.o

mpas_atmphys_tracer_column.o: \
	mpas_atmphys_utilities.o

mpas_atmphys_pbl_gas.o: \
	mpas_atmphys_tracer_column.o

mpas_atmphys_cu_gas.o: \
	mpas_atmphys_tracer_column.o

mpas_atmphys_tracer_transport.o: \
	mpas_atmphys_cu_gas.o \
	mpas_atmphys_pbl_gas.o \
	mpas_atmphys_tracer_column.o

# src/core_atmosphere/chemistry/musica/Makefile (OBJS += mpas_chem_transport_map.o, always built)
mpas_musica.o: mpas_chem_constants.o mpas_chem_transport_map.o
```

libphys still does not depend on libchem: the map crosses only through
`block % dimensions`.

How P1 differs from the [code sketches](pbl-transport-code.md) (§5):

- The native-collision check uses the runtime `chemistry` constituent group
  range instead of a list of native indices. Runtime registration places every
  MICM tracer in that group, so any index outside it aliases another field.
- Process init takes `(block, nTransport)` instead of `(configs)`. A process
  can then reach configs, structs and dimensions, and can log a no-op for an
  empty map.
- Each process exposes `*_enabled()` as a function. When neither process is
  enabled, physics skips all transport work, leaves `tend_scalars` untouched,
  and allocates no tendency buffers.
- The map and column kernels are pure, so they are unit tested without MICM or
  an initialized MPAS log. Only the two wrappers named above log.
- Budgets are weighted by `areaCell`, so masses are in kg. P1 defines the schema
  and the per-rank accumulation over active columns. Rank reduction, output and
  OpenMP directives come with the first active process, when they can be
  exercised.
- `config_micm_pbl_transport` is not part of P1; it belongs to the PBL-only
  stages (code document §6.1).

Verification (GNU gfortran 15.2, double precision, eight ranks, 2026-11-07):

- `testing/chemistry_transport/unit/run_unit_tests.sh`: 4 suites, 66 checks,
  built with `-O0 -g -fcheck=all`.
  - Map: 0, 1, 3 and 64 gases; reordered and noncontiguous indices; host-bound
    H2O, O2 and N2; rejected indices.
  - Column kernels, the budget schema, and the capability and cache checks.
  - Transport core with the real stubs: a bitwise no-op with no buffers.
  - Transport core with link-time test-double processes: PBL-then-CU
    ordering, single application at mapped indices only, inactive columns,
    budgets, and repeated setup.
- `testing/chemistry_transport/test_check_budget.py` (frozen thresholds): 8/8.
- `testing/chemistry_transport/check_transport_wiring.sh`:
  - LNOx: NO, NO2 and O3 map to scalars 6, 5 and 4, and physics reports
    `nMicmTransport = 3` with tendency buffer size 0.
  - LNOx plus H2O: H2O is excluded as host-bound.
  - A `MUSICA=false` build: `nMicmTransport = 0` and the run completes.
- Disabled path: the 60-s and 600-s supercell smoke runs match the pre-P1
  output exactly on all 20 saved fields at every output time. The 60-s
  reference is the tracked `smoke-60s-summary.json`; the 600-s reference came
  from the pre-P1 baseline executable.
- A fresh whole-branch review found no critical issues. Its two important
  findings (buffers allocated while disabled, and a shared interface too narrow
  for both branches) were fixed before the squash. The deferred minor items are
  listed below.

Deferred to the stage that first runs the pool path (P2/P3):

- Compute `m_d` from the same `rho_zz` time level as the tendency mass. Under
  `MPAS_OPENACC` the host copy of time level 1 can be stale.
- Add a test where `nCellsSolve` is smaller than the array extent.
- Rename the budget masses, or document that they cover active columns only
  and are not the domain inventory.
- Add range checks to the budget accessor.
- Tighten the wiring script's `nMicmTransport` match.
- Publish the map through the `inout` dimensions argument of
  `chemistry_init`.

For the convection branch: the shared comments cite this plan, and
`check_transport_wiring.sh` needs the supercell baseline tooling from
`chempas-a-pbl`. The unit tests are self-contained.

## P2 column diffusion kernel and tendency proof (implemented)

P2 has two commits on `chempas-a-pbl`:

- `7ae5f39d`, a shared commit. Cherry-pick it onto `chempas-a-convect` as well;
  together with P1 it applies cleanly there.
  - Assembly takes the dry-mass weight from `rho_zz` time level 2, the array
    that also scales the tendency. At that call it equals level 1, and under
    `MPAS_OPENACC` it is the level refreshed on the host.
  - It adds an owned-cell (`nCellsSolve`) test guard and an exact wiring match.
  - It adds `check_dry_mass_metric.py`.
- `0eb2c8ad`, the PBL-only kernel and the proofs.

No PBL scheme is connected; the production `pbl_gas` stub stays disabled.

### Kernel

`physics/mpas_atmphys_pbl_diffusion.F` provides
`pbl_diffuse_column(dt, m_d, G, q)`, a closed-column backward-Euler diffusion.
`G(k) = rho*Kz/dz` at the interface below layer k, and the ground and top
interfaces are always closed.

The kernel solves for the increment, whose right-hand side is the flux
divergence of `q0`, using one Thomas factorization for every species. It then
applies the interface fluxes of `q0 + increment` to `q0` in flux form. The
results, to roundoff at any stiffness:

- the column mass telescopes exactly;
- a uniform profile has zero right-hand side and stays exactly unchanged;
- values stay within their initial range;
- a zero species stays exactly zero (there is no floor).

The forms considered along the way:

- Solving for the state alone loses mass in proportion to `dt*G/m_d` times
  machine precision: 6e-10 over 1,000 calls at a ratio of 1e5, which fails the
  frozen threshold.
- A flux-form finish on the state conserves mass, but a uniform profile drifts
  past the frozen 1e-12 limit above a ratio of about 2e3.
- The increment form keeps both. Its error relative to the exact backward-Euler
  solution grows like `kappa*eps`, where `kappa = max dt*(G(k)+G(k+1))/m_d(k)`
  (1.7e-12 at a kappa of 960). That error is negligible next to the
  truncation error.

Preconditions, which the P3 caller must check: `m_d > 0`, `G >= 0` and finite,
and consistent sizes.

Verification (`unit/test_pbl_diffusion.F90`, 53 checks):

- agreement with six reference fixtures (`columns/make_diffusion_fixtures.py`),
  each an exact rational dense solve of the flux-form system, within
  `max(1e-12, 32*eps*kappa)`. The fixtures cover nonuniform `m_d` and `G`;
  zero, near-zero and step species; pulses at the bottom, middle and top; a
  sharp gradient; zero coefficients; a large timestep (kappa = 960); and a
  uniform profile;
- mass residual ≤ 1e-12 per call and ≤ 1e-10 over 1,000 calls, up to
  `dt*G/m_d = 1e6`;
- a uniform profile unchanged, and the maximum principle, at 1e4 and 1e6;
- nonnegativity and zero species;
- bitwise invariance to species order and to the values supplied at the closed
  interfaces;
- the analytic backward-Euler decay of cosine modes.

### Dry-mass weight

MPAS applies the vertical scalar flux divergence with `rdzw`, so the tracer
mass it conserves per column is `sum(rho_zz*q/rdzw)`. The kernels use
`m_d = zz*rho_zz*dz`. Initialization defines `zz` and `rdzw` from the same
interface heights, so `zz*dz*rdzw = 1`. `check_dry_mass_metric.py` measures it
on real meshes: 1.1e-16 on the flat supercell and 2.2e-16 on the mountain
wave (250 m relief). `rho_zz` is the dry density (`mpas_init_atm_cases.F`
divides by `1 + qv`).

### Realized-increment proof

`check_realized_increment.sh` builds the model from tracked files with
`model_doubles/mpas_atmphys_pbl_gas.F`, a test double whose mode comes from
`CHEMPAS_TEST_PBL_GAS`:

- source: `q += dt*S`, with S in closed form;
- diffusion: the kernel with constant G;
- off.

Each case runs three 3-s steps on a horizontally uniform dry base state
(`make_resting_init.py`) with an inert mechanism. `check_realized_increment.py`
compares the results with independent predictions.

It covers six advection settings:

- monotonic;
- ordinary;
- positive-definite (MPAS warns that this is unimplemented and runs the
  monotonic path; the outputs are bitwise identical);
- monotonic with split transport (3 dynamics steps per transport step);
- ordinary with split transport;
- monotonic with chemistry every 2 steps.

There are 28 checks on 30 runs, all passing (GNU, eight ranks):

| Regime | Checks | Result |
| --- | --- | --- |
| At rest (eddy viscosity off; \|w\| ≤ 4e-14 m/s), all six settings × source, diffusion, plus off | every layer and cell: realized vs. requested increment within 1e-18 + 1e-10·\|requested\| (source: `S*dt`; diffusion: a numpy dense solve with `m_d` from the start of the step); column budget | per-layer error ≤ 6e-23 kg/kg; budget ≤ 3.4e-16 |
| Moving (the case's eddy viscosity drives a uniform drift, \|Δρ\|/ρ ≈ 9e-6 per step), same matrix | column budget: realized change = `dt*Σ R*m_d(t)` (0 for diffusion and off), within 1e-12 of the column mass | ≤ 4.2e-16; a weight taken at `t + dt` would miss by 3.9e-10 |
| Fresh (NO₂ photolysis on, at rest), monotonic and ordinary | first step: `q(t1) = D(q_off(t1))`, i.e. diffusion of the post-chemistry state from the matching off run | exact within tolerance; a stale, pre-chemistry tendency would miss by 1.3e-13 |

So in the four distinct advection paths (ordinary or monotonic, each split or
not), the gas tendency is applied exactly once and is not replayed. It acts on
the post-chemistry state, and its dry-mass weight is `rho_zz` at the start of
the step. Diffusion conserves the column in every mode.

### Findings

- Under motion, the per-layer placement of a source differs by advection path.
  The monotonic stage 3 adds `dt*R` to the old state before transport, so the
  increment moves with the air. The ordinary path adds it in flux form, in
  place. In the moving case the difference reached 8e-5 of the increment. Both
  conserve the column exactly and agree exactly at rest. This is ordinary
  operator splitting, not a defect, so per-layer agreement is asserted only at
  rest.
- Even from its dry base state, the supercell case's eddy viscosity (500 m² s⁻¹
  on momentum and theta) drives a horizontally uniform vertical drift of a few
  mm/s.
- MICM rejects a mechanism with no rate parameters, so the inert test mechanism
  keeps its reactions, with Arrhenius `A = 0` and photolysis `j = 0`. Turning
  TUV-x off also requires upper-column mode `none`, no top extension, and an
  empty extension file.
- A fresh whole-branch review found one critical issue and two important ones,
  all fixed before the squash:
  - the reference fixtures were ignored by the repository's `*.txt` rule;
  - the earlier state-form finish drifted on uniform profiles when stiff;
  - the coverage claims for the positive-definite and chemistry-interval cases
    were overstated.

Deferred:

- No moist or terrain end-to-end column has been run: the test state has
  `zz = 1` and `qv = 0`. P3's YSU columns will be the first.
- Kernel input checks belong to the P3 caller.
- Script robustness (Open MPI's `mpiexec -x`, provenance when reusing a build)
  and an entry point for the Python self-tests.
- Halo and MPI coverage with horizontally varying sources (P5).
- P1's remaining minor items: budget naming, accessor range checks, and
  `chemistry_init` aliasing.

## P3 YSU gas PBL transport (implemented)

P3 connects YSU to the P2 kernel. `config_micm_pbl_transport = .true.` mixes every
independently transported MICM gas with YSU's scalar exchange at the default
PBL interval. It has two commits on `chempas-a-pbl`:

- `85f68008`, the MMM-physics export patch and its application;
- `c1bd210a`, the YSU gas transport, with its unit tests and in-model
  qualification.

### Operator

YSU mixes its moisture tracers (`qmix`) with the water-vapor matrix. That matrix
is built from `xkzq`, YSU's scalar exchange coefficient after its
entrainment-zone adjustment, which sits on the interface between layers k and
k+1. MPAS's `exch_h` is a different quantity, the heat coefficient `xkzh`, and
YSU stores it on the bottom face of layer k+1. MMM-physics does not export
`xkzq`, so P3 carries a temporary patch that does, as this plan allows:

- `physics/patches/mmm-physics/bl_ysu-export-xkzq.patch` adds an optional
  `xkzqx` output to `bl_ysu_run`. It is set just before the tracer solve, with
  `xkzqx(kte) = 0`. The file keeps its BSD-3 notice.
- `physics/apply_mmm_patches.sh` applies the patch from the physics Makefile.
  It applies only to the pristine file of the pinned tag (`20250616-MPASv8.3`,
  `a4baf7f3`) or leaves an already patched file alone, both identified by
  SHA-256 in `MANIFEST`. It refuses any other contents.
- `checkout_externals` skips every external, with only a warning, while any
  external has local changes. Leaving the patch in place would therefore hold
  MMM-physics and UGWP at their old revisions after an `Externals.cfg` update.
  The Makefile runs `apply_mmm_patches.sh --revert` before `checkout_externals`
  and applies the patch again afterwards. An unchanged external keeps the
  patched file's modification time, so it is not recompiled. An updated
  external is checked out and then stops the build until the patch is
  refreshed. `test_apply_mmm_patches.sh` covers all of these cases (13 checks).
- The MMM-physics BSD-3 text is beside the patch, since the patch contains
  portions of that source.

The MPAS YSU wrapper (`physics_wrf/module_bl_ysu.F`) passes `xkzqx` out as the
optional `exch_q`, on the same k/k+1 interface. The PBL driver allocates the
array it passes only when transport is enabled; an unallocated allocatable
array is an absent argument (Fortran 2008).

After YSU, the driver stores the interface exchange for each owned column
(`pbl_gas_cache_store_moist_column`, then `pbl_gas_exchange`):
`G(k+1) = 0.5*(m_d(k)+m_d(k+1))*xkzq(k)/dzc**2`, with
`dzc = 0.5*(dz(k)+dz(k+1))` and `m_d = rho_p/(1+qv)*dz = zz*rho_zz*dz`. This is
the dry-air form of YSU's `dsig/(g*dza)*xkzq/dza`. A negative or non-finite K,
or a nonpositive mass, stops the model. Later in the same step, transport
assembly applies P2's closed-column backward-Euler solve, with the model `dt`,
to the post-chemistry gases.

Only the operator is held, never an increment. With `config_pbl_interval =
'none'` every step refreshes it before it is used, so there are no restart
fields (P5).

YSU's `qmix` tracer solve has no nonlocal or entrainment source term; those
appear only in the heat and moisture equations. Entrainment reaches tracers
through `xkzq`, and so does top-down mixing (`config_ysu_pblmix`).

### Column reference tests

`unit/test_ysu_gas_operator.F90` (28 checks) calls `bl_ysu_run` itself on
stable, convective and cloud-topped columns, with top-down mixing off and on:

- Exactness of the export. Rebuilt with YSU's moist pressure weights and
  `dt2 = 2*dt`, the export reproduces YSU's own tracer tendency `qmixtnp` to
  ≤ 1e-10. `xkzq` at the k/k+1 interface is the operator YSU applies.
- The CheMPAS-A operator conserves dry mass to roundoff. It converges to YSU's
  as dt → 0: at `dt = 1e-3` s they differ by ≤ 9.3e-4, the moist versus dry
  weighting, of order qv.
- At the model step the two differ by design. Per mode, YSU's tracer update is
  `(1+λ)/(1+2λ)`, so it never removes more than half of a mode in one step.
  Backward Euler gives `1/(1+λ)`. Against a near-exact solution from 1,000
  substeps, the CheMPAS-A step is closer in every column:

  | Column | CheMPAS-A step error | YSU step error |
  | --- | --- | --- |
  | Convective | 0.33 | 0.76 |
  | Cloud-topped, top-down mixing on | 0.28 | 0.65 |
  | Stable | 4.7e-3 | 1.4e-2 |

  So the gas step does not reproduce YSU's moisture step: the operator matches,
  the time discretization deliberately does not.
- In the cloud-topped column, top-down mixing raises the export's maximum K
  from 50 to 110 m² s⁻¹. That branch needs the cloud to reach the PBL top and
  radiative cooling below `kpbl`.

The rest of the chain from YSU to the kernel is unit-tested:

- `unit/test_ysu_wrapper_export.F90` (9 checks) runs the MPAS wrapper with the
  driver's `(ims:ime, kms:kme, jms:jme)` arrays. `exch_q(k)` must equal
  `bl_ysu_run`'s `xkzqx(k)` bitwise, and level `kme` must not be written. With
  `exch_q` absent, or given an unallocated array, every other output must be
  bitwise unchanged. Shifting the copy-out up one level fails three of these
  checks.
- `unit/test_pbl_gas_exchange.F90` (9 checks) covers the construction and input
  validation of G.
- `unit/test_pbl_gas.F90` (19 checks) covers:
  - the gate table;
  - the operator cache and its mask;
  - the dry layer mass from moist density;
  - `apply`, which is bitwise the kernel with the cached G.

### Gate

`pbl_gas_init` runs once, after `chemistry_init`.

| Configuration with `config_micm_pbl_transport = .true.` | Result |
| --- | --- |
| No transported MICM gas | logged no-op |
| `config_pbl_scheme = 'off'` | supported zero transport, logged |
| `bl_ysu` with `sf_monin_obukhov`/`sf_monin_obukhov_rev`, `config_pbl_interval = 'none'`, `config_ysu_pblmix` either | enabled |
| `bl_mynn` | stopped until P4; enabled since P4 (see the P4 gate) |
| a nondefault `config_pbl_interval` | stops: until P5 |
| another scheme, or a surface layer of `off` | stops (MPAS's `physics_namelist_check` already rejects YSU without a surface layer) |

### In-model qualification

These runs are evidence that the operator is integrated correctly and conserves
mass. They are not a physical qualification of YSU gas mixing; that rests on the
column tests above. The supercell grid has 60 stretched levels up to 50 km, with
interfaces at 0, 299, 712, 1182 and 1694 m, so its PBL spans only one or two
layers. A finer-grid fixture belongs before P6.

`check_pbl_ysu.sh BASELINE_DIR RUN_ROOT [MODEL [PRE_P3_MODEL]]` runs the
supercell with YSU and the Monin-Obukhov surface layer. Each run is 300 s on
eight ranks, with output every 60 s. The mechanism is inert (`A = 0`, `j = 0`,
TUV-x off) and lightning is off.

The baseline initial gases are vertically uniform, and mixing cannot change a
uniform profile. `make_pbl_init.py` therefore keeps the supercell meteorology
and sets horizontally uniform profiles:

- NO confined near the surface;
- NO2 1e-9 kg/kg below 1 km. The step falls at the 1182 m interface, above
  both PBL tops, so YSU's free-atmosphere coefficient mixes it in both
  regimes;
- O3 depleted toward the surface.

The stable case uses the baseline surface (288 K, about 11 K colder than the
lowest layer). The convective case sets `skintemp = sst = 305` K. Each case is
run with and without gas transport. `check_pbl_ysu.py` checks:

- the non-gas fields are bitwise identical with and without transport;
- `hpbl > 0` in every column;
- each gas stays finite and no lower than
  `-(1e-20 + 1e-12*max_column(q0))` at every record;
- each gas's domain mass is conserved (≤ 1e-13 relative);
- every gradient gas is mixed down-gradient: the domain total of its dry-mass
  column variance is lower with transport than without.

Results (GNU, eight ranks, 300 s), on final-record data:

| Case | hpbl min / median (m) | hfx median (W m⁻²) | Column variance mixed out: NO / NO2 / O3 | Columns reduced (NO2) | Max gas mass drift |
| --- | --- | --- | --- | --- | --- |
| Stable | 150 / 716 | −56 | 4.9e-4 / 4.5e-5 / 2.9e-5 | 28,047 of 28,080 | 4.1e-16 |
| Convective | 487 / 561 | +76 | 1.2e-1 / 3.8e-5 / 1.6e-3 | 28,046 | 1.9e-16 |
| Convective, top-down on | 487 / 561 | +76 | 1.2e-1 / 3.8e-5 / 1.6e-3 | 28,046 | 3.8e-16 |

Every check passes, and in every pair the non-gas fields are bitwise identical.
In the convective case the domain-median lowest-layer NO falls 11% and O3 rises
1.8%, against 0.04% and 0.006% in the stable case. Convection mixes NO 250 times
and O3 55 times more strongly than the stable case does. NO2's step lies above
both PBL tops, so it mixes about equally in both regimes (3.8e-5 against
4.5e-5) and is left out of that comparison. Every column reduced its NO and O3
variance.

The operator lowers every column's variance at every step, but resolved
transport moves the difference it made between columns. Near the warm bubble a
few columns therefore gain variance: in the stable case, 33 of 28,080 for NO2,
in a ring 5–7 km from the bubble centre, plus one far-field column at
−3.5e-7. In each of them the transport difference is almost entirely one-signed
(`sum(m*d)/sum(m*|d|)` ≥ 0.9955), so the column gained gas, which is horizontal
inflow, rather than being redistributed in place. The required check is
therefore the domain total; the per-column counts are reported, and the column
tests check the operator column by column.

YSU's PBL is not deeper in the convective case. Its stable height uses a
critical bulk Richardson number of 0.25 with the supercell's low-level shear;
the unstable height uses 0 with a thermal excess. The convective case does mix
more strongly, and the check uses that.

With `config_ysu_pblmix` set, YSU also raises the PBL height using the
liquid-ice potential temperature. That changes the meteorology slightly, but
the gas metrics agree with the plain convective case to the digits shown. The
radiatively driven cloud-top branch needs cooling, and the supercell has no
radiation, so only the column tests exercise that branch.

Gate cases: transport requested with MYNN or with a nondefault interval stops
the model with the gate's message. YSU without a surface layer stops in
`physics_namelist_check`. PBL off runs and logs the transport as not enabled.

Transport off, with YSU on, is bitwise identical to the pre-P3 executable in
every output field. The PBL-off 60-s and 600-s smoke runs match the P1
references exactly (`run_supercell_smoke.sh` with `CHEMPAS_SMOKE_REFERENCE`).
The P1 wiring check passes. The final runs were at `c1bd210a`, with GNU on eight
ranks. They are kept outside the repository, under
`~/Data/CheMPAS-A/p3-dev/`: `pbl-ysu-final`, `smoke-60s-final`,
`smoke-600s-final` and `wiring-final`. The pre-P3 executable for comparison is
`pre-p3-atmosphere_model` (revision `50f2666c`).

### Standard LNOx case with YSU gas transport (2026-11-08)

The standard 2-h supercell LNOx case was run with YSU and the Monin-Obukhov
surface layer, once with gas transport off and once with it on. The case has
lightning on, the full mechanism and 10-min output. The runs were at
`872e9255`, on eight ranks. The baseline has the PBL off, so the meteorology of
these runs differs from it.

Tools, in `testing/chemistry_transport/`:

- `run_supercell_lnox_pbl.sh BASELINE_DIR RUN_DIR off|on [MODEL [ysu|mynn]]`
  makes the runs.
- `plot_pbl_transport_diff.py` draws on − off:
  - the lowest layer;
  - sections through the storm to 4 and to 16 km, the 16-km set with the 10
    and 100 ppbv NOx contours;
  - domain-mean profiles;
  - the mean below 1 km against time.

  The colour levels are symmetric logarithmic: three bands per decade, white
  more than three decades below the largest difference.
- `check_pbl_lnox_mixing.py` gives the differences in percent of the off run.
- `run_supercell_kzh_restart.sh` restarts a run for 30 s to write YSU's `kzh`,
  which the output stream does not carry.

Both runs finished cleanly. θ is bitwise identical between them, and the total
nitrogen agrees to 8e-13.

70–86% of the summed |ΔNOx| lies at 8–12 km, in the lightning plume; under
0.5% lies below 1 km. Above the PBL top, YSU's scalar coefficient is its
free-atmosphere heat diffusivity (`xkzq = xkzh`, `bl_ysu.F90:1126`). That is a
local Richardson-number closure driven by shear, with a moist correction in
cloud (`bl_ysu.F90:989–1030`). Heat, moisture and condensate are mixed with it
in both runs; transport extends it to the gases. In a deep-convection case,
"PBL transport" is therefore mostly turbulent mixing in the free troposphere
of the storm.

The sections show dipoles inside the plume. The NOx maximum near 10 km is
lowered and its flanks gain NOx. O3 mirrors NOx: mixing partly fills the
titration minimum, and the redistributed NO titrates O3 elsewhere.

On − off, in percent of the off run:

| | 30 min | 60 min | 90 min | 120 min |
| --- | --- | --- | --- | --- |
| NOx moved to another cell | 0.89 | 0.20 | 0.13 | 0.13 |
| Plume NOx (≥ 10 ppbv): median / p99 / max | 1.3 / 21 / 136 | 0.25 / 5.0 / 52 | 0.16 / 3.7 / 22 | 0.15 / 3.4 / 23 |
| Plume O3: median / p99 / max | 0.28 / 2.2 / 3.4 | 0.05 / 0.94 / 1.9 | 0.03 / 0.66 / 1.9 | 0.03 / 0.56 / 2.0 |
| PBL NOx (≥ 0.1 ppbv): median / p99 / max | no cells | 1.1 / 13 / 18 | 0.59 / 4.7 / 7.8 | 0.89 / 5.7 / 14 |
| PBL O3: max | < 1e-5 | 0.19 | 0.04 | 0.07 |
| Largest \|ΔNOx\| / \|ΔO3\| (ppbv) | 17.6 / 1.66 | 9.2 / 0.88 | 8.2 / 0.52 | 7.1 / 0.76 |

The PBL percentages are of small amounts. The mean below 1 km differs by at
most 6e-4 ppbv of NOx, against about 0.2 ppbv. The lowest-layer differences are
1e-3 to 1e-2 ppbv. O3 changes less than NOx in percent because it stays near
37–50 ppbv in the plume, while NOx rises from 0 to about 300 ppbv within a few
kilometres.

The sizes match YSU's own coefficients. The 30-s restart from the off run's
2-h state continues that run: the θ correlation is 0.999998, and the scripted
restart is bitwise identical to one made by hand. In the plume:

- the median `kzh` is YSU's background, 0.01 m²/s;
- 1.4% of the faces have `kzh` large enough to mix adjacent layers within 2 h
  (p99 112, max 372 m²/s, under the 1000 m²/s cap).

Diffusing the off run's NOx offline with this `kzh` for 10 min moves 0.039% of
the NOx, with a plume p99 of 2.9%. The accumulated difference at 120 min is
0.13% and 3.4%, about three intervals of mixing. That is expected while
advection carries air through the few turbulent layers and lightning adds
unmixed NOx.

This is consistency with the scheme's coefficients, not a physical validation
of YSU. The decline from 30 to 120 min fits a young storm with stronger
turbulent layers, but `kzh` was written only at 120 min. The plots are in the
shared plot folder under `supercell-lnox/2026-11-08-pbl-ysu/`. The runs are
under `~/Data/CheMPAS-A/p4-dev/`: `lnox-ysu-off`, `lnox-ysu-on` and
`lnox-ysu-off-kzh`.

### Deferred

- The MMM-physics export upstream, in a versioned release, to replace the
  patch. Opening that request is an outward-facing step for the maintainers.
- `exch_q(k)` is on the k/k+1 interface, while YSU's `exch_h(k+1)` is on the
  same interface by the bottom-face convention. MYNN (P4) must not copy one
  convention into the other.
- The absent-optional pattern is verified only with gfortran. Flang is deferred
  beyond P6 (see the P6 checklist item).
- Restart: by construction, every step at the default interval refreshes the
  operator before it is used. A restart-versus-continuous run, and a stamp on
  the cache's step, belong with the P5 interval contract.
- Physical qualification on a finer grid, and a configuration with radiation
  for the cloud-top branch in the model.
- `run_unit_tests.sh` applies the MMM patch to the build tree as a side effect.
  The name of `test_tracer_transport_stubs` no longer fits, and repeated
  `pbl_gas_init` has no unit test.
- With `flag_qi = .false.`, the YSU wrapper copies out an undefined ice
  tendency. This is pre-existing and unused by the driver.

## P4 MYNN gas PBL transport (implemented)

P4 connects MYNN to the gas transport in every mass-flux mode. It has seven
commits on `chempas-a-pbl`, besides this record:

- `a68d7598` (P4a), local diffusion;
- `c5e144da`, the MMM-physics patch exporting MYNN's updraft plumes;
- `9eca3f99` (P4b), updraft EDMF;
- `0c21a646`, a fix to MPAS's MYNN wrapper found while testing (below);
- `2181961a` (P4c), downdraft EDMF, with the patch extended to downdrafts;
- `872e9255`, the wrapper allocating the plume export only when asked for it;
- `fd9ccbe8`, the in-model qualification for every MYNN mode.

### Operator

**Local diffusion.** MYNN mixes every scalar with `khdz(k) = rhoz(k)*dfh(k)`
on the bottom face of layer k, and exports `exch_h(k) = dfh(k)*dzk` on that face
(`exch_h(1) = 0`). The gas operator is MYNN's own, in dry form, with no patch:
`G(k) = rhoz_d(k)*exch_h(k)/dzk`, `dzk = 0.5*(dz(k)+dz(k-1))`, with MYNN's face
interpolation `rhoz_d(k) = (rho_d(k)*dz(k-1)+rho_d(k-1)*dz(k))/(dz(k-1)+dz(k))`
and `rho_d = rho_p/(1+qv)` (`pbl_gas_exchange_faces`). The P2 kernel applies it
with the model `dt`. MYNN raises `khdz` to at least `0.5*s_aw` to stabilize its
centered compensating subsidence; that floor is not copied, because the gas
operator's upwind compensation supplies the same order of diffusion and is
monotone.

**Mass flux.** MYNN does not export its plumes, so two hash-checked patches
(`bl_mynn_subroutines-export-plumes.patch`, `bl_mynn-export-plumes.patch`, in
the P3 patch mechanism) add per-plume outputs:

- updrafts (`DMP_mf`, 8 plumes): `aw_up(k,i) = UPA*UPW*Psig_w` on the top face of
  layer k after the flux limiter (summed with `rhoz`, MYNN's `s_aw(k+1)`), and
  the mixing fraction `ent_up = ENT*dz` where the plume leaves layer k;
- downdrafts (`DDMF_JPL`, 5 downdrafts): `aw_dn(k,j) = DOWNA*(-DOWNW)` on the
  bottom face of layer k, where MYNN's `sd_aw(k)` stores it, and `ent_dn`.

MYNN mixes a scalar into a plume with `eps = ENT*dz` whatever the plume's mass
flux does. Where a plume accelerates, that drains a layer of gas it never had,
so the gas operator makes each plume mass-consistent when it stores it:

    E(k) = max(eps*M_in/(1-eps), M_out - M_in),   D(k) = M_in + E(k) - M_out >= 0,
    c(k) = (M_in*c_in + E(k)*q(k))/(M_in + E(k)).

Where MYNN's fraction already keeps `D >= 0`, the rebuilt plume is exactly
MYNN's `(1-eps)*c + eps*q`; elsewhere just enough environmental air is
entrained. An updraft is rooted in layer 1, a downdraft in its start layer, and
a plume detrains wholly in the layer it ends in. Dry mass fluxes use the dry
density at the face the flux crosses.

After the diffusion, `mpas_atmphys_pbl_edmf` rebuilds the plume gas from the
post-chemistry profiles at each substep and applies the plume fluxes with
compensating environmental motion from the upwind layer: downward from above
for updrafts, upward from below for downdrafts. In flux form,

    m_d(k)*dq(k)/dt = sum over plumes [D*c(k) - E*q(k)]
                      - (Mtot(k) + Ntot(k+1))*q(k) + Mtot(k+1)*q(k+1) + Ntot(k)*q(k-1),

so an explicit step of length `dts` is a convex combination of `q(k)`, its
neighbours and the plume concentrations when `dts*(Eu(k) + Mtot(k) + Ed(k) +
Ntot(k+1)) <= m_d(k)`. The step uses the fewest equal substeps that meet this in
every layer: monotone, nonnegative, conservative to roundoff, a uniform gas
unchanged and a zero gas zero. A bound above 1e4 substeps, a non-finite
operator, or invalid plume data stops the model.

Only operators are held, refreshed every step at the default interval, as in P3.

### Column tests

- `unit/test_mynn_gas_operator.F90` runs MYNN through the MPAS wrapper on
  synthetic columns (`unit/mynn_column.F90`):
  - local diffusion reproduces MYNN's own mixing of a passive scalar
    (water-friendly aerosol) to 2.5e-16 in stable and convective dry columns.
    The pinned MYNN mixes ozone but reports a zero ozone tendency, so ozone
    cannot serve;
  - with updrafts, a moist convective column conserves, bounds and leaves
    uniform every gas; the plumes lift a surface-rich aerosol 0.91 times as much
    as MYNN's own aerosol mixing does beyond local diffusion (the
    discretizations differ, centered implicit against upwind explicit);
  - with updrafts and downdrafts, a stratocumulus-topped column conserves and
    bounds every gas, and the downdrafts carry cloud-top gas downward.
- `unit/test_mynn_wrapper_export.F90` checks that the exported updraft flux,
  summed over plumes, is MYNN's `edmf_a*edmf_w` to 2.4e-16; that the export is
  confined to the faces it should occupy; that level `kme` is not written; and
  that every other output is bitwise unchanged when the export is absent or an
  unallocated array. It also checks that repeated MYNN calls on a cloudy column
  without snow microphysics are bitwise identical (see below).
- `unit/test_pbl_edmf.F90` checks the kernels on synthetic plumes: mass flux on
  faces, mass-consistent entrainment, MYNN's mixing wherever its fraction is
  consistent (97 updraft and 31 downdraft layers), an independent layer-budget
  form (agreement to 4e-16), conservation per call and over 1,000 calls, bounds,
  zero and uniform gases, bitwise substep equivalence, and rejection of invalid
  input.
- `unit/test_pbl_gas.F90` and `unit/test_pbl_gas_exchange.F90` cover the gate,
  the face exchange, the plume caches and their mask, and `apply` (bitwise the
  diffusion followed by the plume step).

### Gate

| Configuration with `config_micm_pbl_transport = .true.` | Result |
| --- | --- |
| `bl_mynn` (MPAS selects `sf_mynn`), `config_pbl_interval = 'none'` | enabled; logs the operator: local diffusion, plus updraft EDMF with `config_mynn_edmf = 1`, plus downdraft EDMF with `config_mynn_edmf_dd = 1` as well |
| `bl_mynn` with a nondefault `config_pbl_interval` | stops: until P5 |

MYNN runs its downdrafts only inside its updraft scheme, so
`config_mynn_edmf_dd` has no effect with `config_mynn_edmf = 0`.

### MYNN defects found

- MPAS's `bl_mynn_pre_run` zeroed the specific cloud water and ice but not the
  snow content, which it sets only when the microphysics carries snow. With
  Kessler, MYNN's cloud PDF read an uninitialized stack array, so its results in
  cloudy columns depended on earlier calls (8.6 m² s⁻¹ in `exch_h` at a cloud
  top between two identical calls). `0c21a646` zeroes it; runs with snow are
  unchanged. Transport-off MYNN runs are therefore compared with the pre-P4
  executable plus that fix.
- Not fixed, pinned MMM-physics:
  - `DDMF_JPL` limits the downdraft velocity with `MIN(1.25*dz/200, -2)` where
    a -2 m/s bound is meant, so every downdraft slows by `1.25*dz/200` m/s per
    level and survives MYNN's two-level minimum only when it starts faster
    than that. Realistic stratocumulus downdrafts barely transport; the column
    test uses strong cloud-top cooling.
  - `DDMF_JPL`'s diagnostic averaging reads `DOWNA(k-1)` and `ENT(k-1)` below
    their bounds. The unit tests compile that file without bounds checking.
  - `bl_mynn_pre_run` converts snow with `qs/(1+qs)` instead of `qs/(1+qv)`.
    Fixing it would change runs with snow microphysics.

### In-model qualification

`check_pbl_mynn.sh BASELINE_DIR RUN_ROOT [MODEL [PRE_P4_MODEL]]` runs the
supercell with MYNN (MPAS selects `sf_mynn`) and the P3 setup: inert mechanism,
lightning off, `make_pbl_init.py` gas profiles, 300 s on eight ranks with output
every 60 s, each case with and without gas transport, checked by
`check_pbl_ysu.py` (non-gas fields bitwise, `hpbl > 0`, nonnegativity, mass
drift ≤ 1e-13, down-gradient mixing). Results (GNU, eight ranks, final record):

| Case | hpbl min / median (m) | hfx median (W m⁻²) | Column variance mixed out: NO / NO2 / O3 | Max gas mass drift |
| --- | --- | --- | --- | --- |
| Local, stable (`config_mynn_edmf = 0`) | 299 / 299 | −32 | 7.1e-4 / 1.5e-6 / 1.3e-5 | 4.1e-16 |
| Local, convective (305 K) | 299 / 401 | +73 | 9.4e-3 / 5.1e-6 / 1.3e-4 | 4.1e-16 |
| Updraft EDMF (MPAS default), convective | 299 / 401 | +73 | 9.4e-3 / 5.1e-6 / 1.3e-4 | 4.1e-16 |
| Updraft and downdraft EDMF, convective | 299 / 401 | +73 | 9.4e-3 / 5.1e-6 / 1.3e-4 | 4.1e-16 |

Every check passes; every on/off pair has bitwise identical non-gas fields, and
convective mixing exceeds stable mixing (NO by 13 times, O3 by 10). Transport
off with downdrafts on, and local diffusion off, are bitwise identical in every
output field to the executable before P4 with the snow fix
(`pre-p4-fixed-atmosphere_model`, `b6373312`: `8a727d94` plus `0c21a646`).

MYNN's mass flux does not activate on this grid, so the EDMF rows repeat the
local ones. MYNN switches its plumes off unless the first plume rises a full
model level, and its PBL top sits at the first interface (299 m) under the
413 m second layer. Scale awareness (the mesh spacing is 540 m) would suppress
the plumes as well, but without it a run at 315 K still stays plume-free for
15 simulated minutes, its PBL below 432 m, and a model column replayed through the unit
harness confirms the plume deactivation. No stratocumulus top triggered the
downdraft scheme either. The EDMF pairs are therefore integration evidence: the
mass-flux code paths run with gas transport on, conserve, and leave the
meteorology bitwise unchanged. The mass-flux operator itself is qualified by the
column tests on MYNN's own plumes. A finer-grid fixture belongs before P6, as
for P3.

The P3 qualifications rerun on the final P4 executable pass unchanged:
`check_pbl_ysu.sh` (45 checks; the YSU mixing fractions equal P3's to the digits
shown in its table, and transport off is bitwise identical to the pre-P3
executable), the PBL-off 60-s and 600-s smoke runs against the P1 references,
and the P1 wiring check. The final runs were at `872e9255` (scripts at
`fd9ccbe8`; the MYNN run's provenance records `9e2e5f88`, the same script
before a comment correction), with GNU on eight ranks, kept outside the repository under
`~/Data/CheMPAS-A/p4-dev/`: `pbl-mynn-final`, `pbl-ysu-final`,
`smoke-60s-final`, `smoke-600s-final` and `wiring-final`, with the executables
`p4-872e9255-atmosphere_model` and `pre-p4-fixed-atmosphere_model`.

### Deferred (P4)

- Report upstream (outward-facing; for the maintainers): the export patches;
  `DDMF_JPL`'s velocity limiter and out-of-bounds diagnostic reads; the snow
  conversion in `bl_mynn_pre_run`.
- MYNN's downdraft mixing fraction is not capped (the `0.9/dz` cap is commented
  out upstream), so `eps >= 1` would stop the model as invalid plume data. The
  limiter defect makes such downdrafts unreachable today.
- The plume caches hold `(2*nlev+1)*13` values per column, 26 times the diffusion
  cache, and each substep rebuilds all plumes over the full column. Storing only
  up to the highest plume top would cut both (with P6 performance work).
- The unit wrapper test does not drive MYNN's flux limiter; a review probe
  confirmed the export equals `edmf_a*edmf_w` when it fires.
- `pbl_edmf_column` relies on `M(1,:) = 0` and `N(nz+1,:) = 0` from its
  producers, and its error message names only the substep bound. The kernel
  test's layer-budget reference reuses the plume routines it checks.
- `test_apply_mmm_patches.sh` stops without a FAIL line if its mtime loop fails.

## P5 plan: PBL intervals, restart and MPI

P5 covers what the PBL transport can qualify alone: nondefault PBL intervals,
restart, chemistry intervals, MPI decomposition, and a gate against
parameterized convection.

The first science target is the resolved supercell: about 540 m spacing, with
`config_convection_scheme = 'off'`. Resolved advection carries the gases through
the storm and no convection scheme runs, so the PBL work completes without the
convective transport of #41. Coupling with parameterized convection moves to
#41's plan (decided 2026-11-08). An earlier draft split P5 into P5a (PBL-only)
and P5b (coupled, after #41); the commits made under that draft keep the P5a
label.

### Audit (2026-11-08)

- `config_pbl_interval` does not set a call interval in MPAS. The physics
  manager creates the `pbl` alarm and sets `dt_pbl` to the interval, but nothing
  tests the alarm (`mpas_is_alarm_ringing` is called for radiation, convection
  and the other physics alarms, never for `pbl`). `physics_driver` calls the surface
  layer, the land surface and the PBL on every physics step, and the PBL schemes
  (YSU `dt`, MYNN `delt`), Noah (`dt`), the gravity-wave drag, the ocean mixed
  layer and sea ice all integrate with `dt_pbl`. With an interval longer than
  `config_dt`, each step therefore advances those schemes by the interval: the
  land surface runs fast by interval/dt, and the PBL tendencies come from an
  implicit step of the wrong length. Upstream MPAS (`NSF-NCAR/MPAS`) has the same
  code. Only `'none'` (`dt_pbl = dt`) is consistent.
- The PBL tendencies (`rublten`, `rvblten`, `rthblten`, `rqvblten`, …) live in
  `tend_physics`, are applied at every dynamics step and are in the restart
  stream, so MPAS already supports holding them between PBL calls (as WRF does
  with `bldt`).
- The gas operator cache is module memory, refreshed by every PBL call and not
  restarted. At the default interval every step refreshes it before use.
- `config_chemistry_interval` (seconds, 0 = every step) sets the MICM cadence.
  Gas PBL transport acts on the post-chemistry gases every step.
- Partition files for 2, 4, 8, 12, 16, 24 and 32 ranks exist in
  `~/Data/MPAS/supercell/`.

### Decision before Task 2 (made 2026-11-08: option 1, repair)

Nondefault PBL intervals are broken in native MPAS, independent of chemistry.
Gas transport at such intervals cannot be qualified until the native cadence is
right. The options:

1. **Repair the cadence in CheMPAS-A (recommended).** Run the surface layer,
   land surface and PBL group only when the `pbl` alarm rings, with `dt_pbl` the
   interval, and hold their tendencies in between: WRF's `bldt` semantics, which
   the `dt_pbl` arguments already assume. `'none'` stays bitwise unchanged. This
   changes results only for runs that set a nondefault interval, which are wrong
   today. The defect would be reported upstream (outward-facing: your call).
2. **Keep rejecting nondefault intervals.** Document the native defect, report
   it upstream, and close P5a's interval work as "not supported". Restart, chemistry
   intervals and MPI still proceed.

### Tasks

1. **Reproduce the native defect** (no code change): a supercell YSU run with
   `config_pbl_interval` three times `config_dt`, transport off, logging PBL
   calls and `DT_PBL`, against `'none'`. Expected: the PBL runs every step with
   `dt_pbl` = 3·dt. Gives the evidence for the decision and the upstream report.
2. **Native cadence repair** (option 1 only; its own commit). Gate the
   surface-layer/LSM/PBL group (and the schemes that take `dt_pbl`: GWDO, OML,
   sea ice) on the `pbl` alarm, reset it like the convection alarm, and keep the
   held tendencies. Tests: `'none'` and an interval equal to `config_dt` are
   bitwise identical to the pre-P5 executable; at 3·dt the group runs every third
   step (log), tendencies are held between calls, and restart files carry them.
   A native-meteorology comparison against `'none'` is reported, not required to
   match.
3. **Gas transport at a nondefault interval.** The held operator is applied every
   step with the model `dt` between refreshes, consistent with the held
   meteorological tendencies. Add a refresh stamp to the cache (the step or time of
   the last PBL call) and stop the model if the operator is applied before any
   refresh or after its interval has expired. Lift the interval gate for YSU and
   MYNN. Tests: unit tests for the stamp and gate; in-model on/off pairs at 3·dt
   for YSU and MYNN (non-gas fields bitwise, conservation, down-gradient mixing).
4. **Restart.** Continuous against restarted runs, gas transport on, at `'none'`
   and at 3·dt, with the restart on an interval boundary: every field bitwise
   identical (the alarm rings at the restart time, as radiation's does, so the
   operator is rebuilt before use). A restart off the interval boundary shifts
   the PBL cadence as it does for radiation; it is reported and documented, not
   required to match. Add the repeated-`pbl_gas_init` unit test (P3 deferred).
5. **Chemistry intervals.** Transport on with `config_chemistry_interval` of 1, 3
   and 7 model steps (YSU and MYNN): conservation, nonnegativity, and the P2
   realized-increment check that the PBL increment is applied exactly once per
   step whether chemistry ran that step or not.
6. **MPI decomposition.** The same 60-s case, transport on, on 1, 2, 4 and 8
   ranks: gas fields aligned by `indexToCellID` agree within `atol = 1e-18 kg/kg`,
   `rtol = 1e-10`, with observed differences reported (the spec's thresholds; the
   meteorology is not expected to be bitwise across decompositions). Include a
   decomposition with ranks whose columns are all inactive (PBL off in a region is
   not available, so a rank count high enough to leave a rank with only
   zero-exchange columns, or a unit test of an all-inactive rank).
7. **Parameterized convection gate.** Stop at initialization when
   `config_micm_pbl_transport` is set and the resolved `config_convection_scheme`
   is not `'off'`, with a message naming #41. This includes the scale-aware
   Grell-Freitas of the `convection_permitting` suite. Tests: a unit test of
   `pbl_gas_capability` with each convection scheme and `'off'`; an in-model gate
   case (YSU with a convection scheme stops with the message); the resolved
   supercell (`'off'`) still runs.
8. **Documentation**: the P5 section, the gate table, the supported
   configurations (convection-resolving, with parameterized convection off), and
   the upstream report text for the native defect (filing it is your call).

### Progress

- Tasks 1 and 2 (2026-11-08). `check_pbl_cadence.sh` runs the supercell for 20
  steps of 3 s and counts scheme calls in the timer summary. Before the repair,
  with `config_pbl_interval = 00:00:09` YSU, MYNN and their surface layers ran on
  all 20 steps while integrating with `DT_PBL = 9`. The repair adds `l_pbl`, set
  from the `pbl` alarm in `physics_timetracker` like `l_conv`, and runs the
  surface layer, land surface (Noah, and Noah-MP, whose step is now `dt_pbl`),
  ocean mixed layer, sea ice, surface diagnostics, PBL and gravity-wave drag only
  when it is set; their tendencies are held in between. Gravity-wave drag has to be
  in the group: it adds its drag to the held PBL tendency, so it would otherwise add
  it again every step. After the repair the 9 s interval gives 7 calls (steps 1, 4,
  …, 19); an interval equal to `config_dt` is bitwise identical to `'none'`, and
  `'none'` is bitwise identical to the executable before the repair. After 60 s the
  9 s run differs from the every-step run by at most 1.1e-3 K in θ and 2 m in the
  PBL height. Gas transport still rejects nondefault intervals until Task 3.
- Task 3 (`25d09c38`, `d2ac0898`). After each PBL call the driver stamps the gas
  operator cache with the model time and `dt_pbl`
  (`pbl_gas_cache_refreshed`). Transport applies the held operator every step
  with the model `dt`. Before each use, `pbl_gas_begin_step` stops the model if
  the operator was never built, or if `age*dt >= interval`, which means a refresh
  was missed. The pbl alarm rings at the first step at or after each interval,
  so the oldest operator used is `ceiling(interval/dt) - 1` steps old. The
  interval gate is lifted for YSU and MYNN.

  At the default interval, 60-s YSU and MYNN runs with transport on are bitwise
  identical to the pre-P5 executable. `check_pbl_intervals.sh` covers 9-s
  intervals for YSU and MYNN:
  - the non-gas fields are bitwise identical with transport off and on;
  - gas mass is conserved (drift ≤ 4e-16), and the gases mix down-gradient;
  - the budget log shows refreshes at 0, 9, …, 54 s.

  With a 10-s interval and 3-s steps the refreshes are at 0, 12, 21, 30, 42 and
  51 s. PBL off, with transport requested and an interval, runs as a
  zero-transport case. `check_pbl_ysu.sh` and `check_pbl_mynn.sh` pass, with
  their interval gate cases turned into runs that must succeed.
- Budget diagnostics (`f62f5423`), for the completion criteria.
  `config_micm_transport_budget` logs, each step, per process and gas, the
  domain budget reduced over ranks:
  - the dry-air-weighted mass before and after, the increment and its relative
    size;
  - the largest column residual and the extrema;
  - the active columns, the refresh time and `dt`.

  `check_transport_budget.py` checks it. A 60-s YSU run has |rel| ≤ 1.6e-16 and
  column residuals ≤ 1e-15, and the output is bitwise identical with the option
  off.
- Task 4 (`a9a7d2c7`). The operator is not restarted. The pbl alarm starts at
  the restart time, so the first step rebuilds the operator, a tested
  deterministic reconstruction (shared interface item 8).
  `check_pbl_restart.sh` runs 144 s continuously with a restart at 72 s, and
  72 s again from that restart, with transport on. YSU and MYNN, at `'none'`
  and at 9 s: every output field at 144 s is bitwise identical, as in a
  transport-off pair. A restart at 60 s, off the 9-s boundary, runs to
  completion. The PBL cadence then shifts, as radiation's does: after 60 s, θ
  differs by up to 1.2e-5 K and NO by up to 3e-5 relative. That is reported,
  not required to vanish. A unit test covers repeated initialization.
- Task 7 (`72fb7399`). Gas PBL transport stops at initialization with any
  active convection scheme, with a message naming #41. That includes
  Grell-Freitas of the `convection_permitting` suite; `physics_namelist_check`
  resolves the suite before transport starts. It also stops with more than one
  OpenMP thread. PBL off stays a zero-transport case. In-model, YSU with new
  Tiedtke stops with the message, while transport with convection off runs, and
  so does convection with transport off. The OpenMP gate is unit-tested only:
  the builds use `OPENMP=false`.
- Task 5 (`d4e7c954`). `check_pbl_chemistry_intervals.sh` runs YSU and MYNN
  with chemistry every 3, 9 and 21 s, with the reactive LNOx mechanism and NO2
  photolysis, transport off and on, for 63 s. In each pair:
  - the non-gas fields at 63 s are bitwise identical;
  - nitrogen (NO + NO2) and odd oxygen (NO2 + O3) drift by at most 4.9e-16;
  - no gas falls below the nonnegativity floor;
  - the budget log is conservative on all 21 steps (|rel| ≤ 3.3e-16, column
    residuals ≤ 1.3e-15).

  `check_realized_increment.sh` adds chemistry every 9 and 21 s to its
  advection modes; the realized increments match the requested ones to 5e-23.
- Task 6 (`bdf8aea2`). `check_pbl_mpi.sh` runs YSU and MYNN at a 9-s PBL
  interval on 1, 2, 4 and 8 ranks. It requires the gases at 60 s to agree with
  the 8-rank run within atol 1e-18 and rtol 1e-10; in fact every output field
  is bitwise identical for each rank count. A unit case with every transport
  process inactive leaves the gases untouched.

Out of scope here: coupling with parameterized convection (#41's plan), OpenMP
thread-private buffers (the builds use `OPENMP=false`), the deferred Clang/Flang
build, and the P3/P4 items deferred to P6 (see [P6](#p6-gnu-fortran-eight-rank-qualification-implemented)).

### Upstream report (draft, not filed)

Filing it is an outward-facing step for the maintainers.

> **`config_pbl_interval` sets `dt_pbl`, but the surface and PBL schemes still
> run every step**
>
> In `mpas_atmphys_manager.F`, `physics_run_init` creates the `pbl` alarm from
> `config_pbl_interval` and sets `dt_pbl` to that interval. Nothing tests the
> alarm, though: `physics_timetracker` checks the radiation, convection and
> other physics alarms, never `pbl`. `physics_driver` (`mpas_atmphys_driver.F`)
> calls `driver_sfclayer`, `driver_lsm`, `driver_pbl` and `driver_gwdo` on every
> physics step. Those schemes, the ocean mixed layer and sea ice integrate with
> `dt_pbl`. With an interval longer than `config_dt`, each step therefore
> advances them by the interval. The land surface runs interval/dt times too
> fast, and the PBL tendencies come from an implicit step of the wrong length.
> Only `'none'` (`dt_pbl = dt`) is consistent.
>
> To reproduce, run the idealized supercell with `config_dt = 3` and
> `config_pbl_interval = '00:00:09'` for 20 steps. The log reports
> `DT_PBL = 9`, and the timer summary shows the surface layer and PBL called on
> all 20 steps.
>
> A fix: set a flag from the `pbl` alarm in `physics_timetracker`, as `l_conv`
> is set from the convection alarm. Run the surface layer, land surface, ocean
> mixed layer, sea ice, surface diagnostics, PBL and gravity-wave drag only when
> it is set, and hold their tendencies in between. The tendencies already live
> in `tend_physics` and the restart stream. Gravity-wave drag must be in the
> group, because it adds its drag to the held PBL tendency. With the fix, a 9-s
> interval gives 7 calls in 20 steps. An interval equal to `config_dt` is
> bitwise identical to `'none'`, and `'none'` is unchanged. The fix is in
> CheMPAS-A commit `d90b9fd7`.

## P6 GNU Fortran eight-rank qualification (implemented)

P6 qualifies MICM gas PBL transport for convection-resolving configurations,
built with GNU Fortran 15.2.0 (conda-forge), Open MPI 5.0.10, NetCDF-C 4.10.1,
NetCDF-Fortran 4.6.3, PIO 2.6.9 and the pinned MUSICA 0.16.5, in double
precision without OpenMP, on eight ranks.

### Supported configurations

| Setting | Supported | Stops at initialization |
| --- | --- | --- |
| PBL scheme | `bl_ysu`, with or without `config_ysu_pblmix`; `bl_mynn` with local diffusion, with updrafts (`config_mynn_edmf = 1`), and with updrafts and downdrafts (`config_mynn_edmf_dd = 1`); `off` (zero transport) | any other scheme |
| `config_pbl_interval` | `'none'` or at least `config_dt` (3, 6, 9 and 10 s tested with 3-s steps); the operator from the last PBL call is applied every step. Shorter intervals are not supported: MPAS then integrates the PBL and surface schemes every step with the shorter interval | — |
| `config_convection_scheme` | `'off'` | any other, including the Grell-Freitas of `convection_permitting` (#41) |
| `config_chemistry_interval` | any (3, 9 and 21 s tested) | — |
| Restart | yes; the operator is rebuilt at the restart's first step; bitwise on a PBL-interval boundary | — |
| MPI | any rank count (1, 2, 4 and 8 tested) | — |
| OpenMP | one thread | more than one thread |
| Toolchain | GNU Fortran; Clang/Flang not qualified | — |

### Fixtures and results

`testing/chemistry_transport/run_pbl_qualification.sh` runs every check in
order and writes a summary and a manifest. The final run is recorded in
⟨TBD: results link⟩.

⟨TBD: results table⟩

### Finer vertical grid

`prepare_fine_grid.sh` builds a supercell grid with 100 layers: 50 m up to
1.5 km, then stretched by at most 10% per layer to the baseline's 50-km top
(`make_fine_zeta_levels.py`). On the baseline grid, with interfaces at 0, 299,
712 and 1182 m, the PBL spans only one or two layers and MYNN's plumes never
start. `check_pbl_fine.sh` runs YSU, MYNN with updrafts, and MYNN at a 9-s
interval for 10 min on the convective init.

Each pair runs with transport off and on. In every pair, every non-gas field
is bitwise identical between the two runs. Gas mass is conserved within
4.4e-16, every gas stays nonnegative, and the budget log passes for all
200 steps.

| Case | PBL layers (median hpbl) | NO range: on / off / replay | Mixing, fraction of the replay | Updraft columns | Budget max \|rel\| |
| --- | --- | --- | --- | --- | --- |
| YSU | 15 (727 m) | 0.87 / 1.79 / 0.74 | 0.81 | — | 4.8e-16 |
| MYNN with updrafts | 18 (901 m) | 1.64 / 2.10 / 1.66 | 1.04 | 98.9% | 6.2e-16 |
| MYNN, 9-s interval | 18 (902 m) | 1.64 / 2.10 / 1.66 | 1.04 | 98.9% | 6.2e-16 |

The NO range is the domain median of `(max - min)/mean` over the layers below
0.8·hpbl at 10 min; the initial NO falls off as exp(-z/300 m). The replay
diffuses the initial NO offline with the scheme's own `kzh` from the output.
The mixing column is ln(off/on) divided by ln(off/replay); the check requires
0.5 to 2. A run without transport scores 0 and fails.

Ten minutes is shorter than the PBL's mixing time, so NO is not yet well mixed.
For YSU, the mean `kzh` in the PBL is 37–52 m²/s and the slowest mode's time
scale, (hpbl/π)²/K, is about 16 min. MYNN starts with little TKE: its mean
`kzh` grows from 2 to 7 m²/s over the run. The plan's first criterion, a range
5 times smaller with transport, was out of reach for a correct operator; the
replay itself reaches only 2.4 for YSU and 1.27 for MYNN. The criterion is now
agreement with the replay. YSU mixes slightly less than the replay, as
expected. Its scalar coefficient `xkzq` does not exceed the heat coefficient
`kzh` inside the PBL, and the replay holds each minute's latest, growing K.

MYNN's updrafts are active in 98.9% of the columns on this grid. On the
baseline grid they never started. The mass-flux path therefore now runs in
the model with gas transport on, conserves, and leaves the meteorology
bitwise unchanged. Over these 10 min the updrafts add little to the mixing of
NO measured this way. The mass-flux rate itself rests on the column tests.

The supercell namelist's explicit vertical eddy viscosity
(`config_v_mom_eddy_visc2 = config_v_theta_eddy_visc2 = 500 m2/s`) is off in
these cases. It is a tendency computed once per step, so it is stable only for
K·dt/dz² ≤ 0.5. On 50-m layers at dt = 3 s that number is 0.6, against 0.017
on the baseline grid, and every fine-grid run blew up after about 100 s,
transport off included (w rose from 1 to 13 m/s within a few steps, then went
to NaN). With the term off at dt = 3 s, or kept at dt = 2 s, 4-min runs stay
stable. Off is the choice here: it keeps the 3-s step of every other check and
leaves the resolved boundary layer to the PBL scheme.

### Reactive regression

`check_pbl_reactive.sh` uses the baseline LNOx mechanism with NO2 photolysis
(`config_j_no2_max = 0.01`). Each case runs with transport off and on. Without
lightning, nitrogen (NO + NO2) and odd oxygen (NO2 + O3) are conserved by the
reactions and by transport. Lightning adds NO only, so odd oxygen stays
conserved with it. Invariants are in moles, weighted by dry-air mass.

| Case | Non-gas fields | Nitrogen drift | Odd-oxygen drift | Budget |
| --- | --- | --- | --- | --- |
| Stable init, YSU, 5 min | 22 bitwise | 1.6e-16 | 1.3e-16 | 100 steps, max \|rel\| 7.3e-16, residual 2.9e-15 |
| Stable init, MYNN, 5 min | 24 bitwise | 4.9e-16 | 1.3e-16 | 100 steps, max \|rel\| 1.7e-16, residual 1.5e-15 |
| Supercell init with lightning, YSU, 30 min | 22 bitwise | (source) | 1.2e-16 | 600 steps, max \|rel\| 1.9e-16, residual 1.1e-15 |

The drift limit is 1e-12. Every value of NO, NO2 and O3 stayed above the
nonnegativity floor. These are numerical checks of conservation; they do not
validate the chemistry or the lightning source.

### Rulings

- The Chapman+NOx case runs on the global JW mesh without a land surface. At
  that spacing it needs parameterized convection, which is outside the
  convection-resolving configurations qualified here.
- The reactive regression uses the LNOx mechanism's linear invariants. The ABBA
  mechanism was a development case and is being retired.
- The LNOx case with lightning runs 30 min, because lightning starts after about
  10 min of storm spin-up.
- MPI decomposition: the check requires the gases to agree within atol 1e-18 and
  rtol 1e-10; every output field turned out bitwise identical.
- A restart off the PBL interval boundary shifts the PBL cadence, as it does for
  radiation. It is reported, not required to match.
- The fine-grid cases switch off the explicit vertical eddy viscosity (above).
- The fine-grid mixing check compares the model with an offline replay using
  the scheme's own diffusivity, not with a fixed factor. In 10 min, a correct
  operator cannot mix the PBL well (above). The band of a factor of 2 accepts
  a rate error of up to 2 in the model. The column tests pin the rate exactly.

### Deferred

- The Clang/Flang build, including the absent-optional pattern of the YSU
  wrapper, until a Flang platform is in use.
- OpenMP: thread-private buffers for the operator caches. Until then,
  initialization stops with more than one thread.
- MYNN plume cache sizing: the caches hold `(2*nlev+1)*13` values per column,
  and each substep rebuilds the plumes over the full column. Storing only up to
  the highest plume top would cut both.
- A configuration with radiation, for YSU's cloud-top branch in the model.
- The budget option's per-step MPI reductions; it is a diagnostic, off by
  default.
- Upstream reports, outward-facing and for the maintainers: the native PBL
  cadence defect (draft above), the MMM-physics export patches, `DDMF_JPL`'s
  velocity limiter and out-of-bounds diagnostic reads, and the snow conversion
  in `bl_mynn_pre_run`.

## Verification and completion criteria

Create version-controlled generators, namelists, mechanism variants, checks,
and small column reference data under `testing/chemistry_transport/`. Record
external input URLs, checksums, exact revisions, compiler/dependency versions,
build commands, MPI command, and elapsed simulated time. Large model inputs and
NetCDF outputs belong in the run-data tree, with a manifest and summarized
results tracked alongside this plan.

| Test group | Required evidence |
| --- | --- |
| Mapping/ownership | Empty/inactive chemistry; 1, 3, and 64 independently transported gases; reordered mechanisms; synthetic noncontiguous indices interleaved with native/passive fields; bound H2O and prescribed O2/N2. Every eligible species participates exactly once; excluded fields and unrelated scalars do not receive gas tendencies. |
| Column physics | Nonuniform dry density and layer spacing; uniform profiles; isolated lower/middle/top pulses; sharp gradients; zero and near-zero species; diffusion/nonlocal components separately; zero coefficients and large stable timesteps. Check every species, not just a sum over gases. Compare diffusion against an independently assembled matrix or analytic constant-coefficient mode. |
| Active PBL | Stable and convective boundary-layer fixtures with YSU and each qualified MYNN mode; show nonzero exchange and redistribution of every gradient tracer. Keep convection and chemical reactions off initially, with an inert MICM mechanism still loaded. |
| Time/parallel/restart | Repeat many updates with chemistry intervals of 1, 3, and 7 model steps; test qualified physics intervals separately, including no-mixing transitions. Compare one and eight ranks by `indexToCellID`, include ranks with no active columns, and restart between refresh times. Check thread-private buffers if OpenMP is enabled. |
| Coupled/regression | Reactive cases with linear invariants (P6: LNOx nitrogen and odd oxygen; the ABBA development mechanism is being retired, and Chapman/NOx is outside the convection-resolving subset), then the resolved supercell LNOx case with PBL transport. Check transport-only budgets before chemistry and reaction invariants afterwards. Compare disabled transport to the unmodified base, including native meteorology. Exercise each supported scalar advection setting. PBL + parameterized convection + MICM is qualified under #41, not here. |
| Toolchains/output | A double-precision GCC/gfortran build of the full pinned CheMPAS-A/MUSICA stack, with eight-rank NetCDF runs. A build of plain MPAS or a chemistry run from before the transport work does not satisfy this gate. A Clang/Flang build is deferred until a Flang platform is in use. |

Initial double-precision acceptance thresholds (freeze them with fixtures before
implementation; justify any later changes):

- Isolated closed-column relative species-mass residual <= `1e-12` per call and
  <= `1e-10` over 1,000 calls; full integration accumulated transport-budget
  residual <= `1e-8`. Use both absolute and relative errors. Normalize using
  `max(initial species mass, sum(m_d)*1e-30)` for per-area columns, and include
  area for domain budgets. Exactly zero species must remain zero.
- Uniform-profile and permutation comparisons: absolute error <= `1e-20 kg/kg`
  plus `1e-12` times the reference magnitude. Nonnegativity allows roundoff no
  smaller than `-(1e-20 + 1e-12*max_column(q_initial)) kg/kg`; larger violations
  fail even if a later chemistry clamp would remove them.
- Serial/eight-rank, restart/continuous, and short compiler-comparison fixtures:
  `atol=1e-18 kg/kg`, `rtol=1e-10` after cell-ID alignment. Report observed
  differences separately; long nonlinear meteorology is not a bitwise test.
- Inactive/disabled tendency is exactly zero. Disabled transport reproduces the
  base fields on the same build/settings, preferably bitwise, excluding logs and
  timing metadata. All output is finite; expected duration and restart events
  are reached; no rank reports an error.

Write per-species PBL tendency, integrated increment, min/max, before/after dry
mass, boundary flux, conservation residual, active-column count, coefficient
refresh time, and applied `dt` to selectable diagnostics/NetCDF or a compact
budget stream. Check the realized state increment as well as the computed
tendency so omitted or duplicate application is visible.

The resolved supercell LNOx case is the initial regression baseline, with PBL
and parameterized convection explicitly off. Use
`testing/chemistry_transport/prepare_supercell_baseline.py` to stage fresh inputs
and pinned TUV-x data. The full case runs for two simulated hours with a 3-second
timestep on eight MPI ranks, starting at `0000-01-01_18:00:00`; output is saved
every ten minutes. Initialize NO and NO2 to zero and O3 to 50 ppbv with
`scripts/init_lnox_o3.py` (nominal 50 ppbv using that initializer's 0.029 kg/mol
air convention, or `8.275862068965518e-8 kg/kg`). Reported ppbv diagnostics use
the host's 0.0289644 kg/mol air constant. This establishes existing resolved
transport and reactive chemistry behavior; active-PBL qualification still needs
separate fixtures.

The GNU baseline on 2026-11-07 at `2c6ff111` passed all 2,400 steps and 13 output
records with zero model errors. Maximum relative dry-air mass and Ox inventory
drifts were `2.181e-16` and `2.413e-16`, respectively. The 60-second developer
smoke also passed; the 600-second option is available but was not run in this
session. See the recorded [manifest](../../testing/chemistry_transport/results/supercell-lnox-gnu-20261107/manifest.json)
and [summary](../../testing/chemistry_transport/results/supercell-lnox-gnu-20261107/summary.json).
The manifest records existing startup warnings and the local run-data location;
large NetCDF outputs remain outside the repository. A final restart was written,
but restart/continuous equivalence has not been tested.

For development, reuse the frozen baseline inputs with a fresh executable:

```bash
bash testing/chemistry_transport/run_supercell_smoke.sh BASELINE_DIR NEW_RUN_DIR 60
# Optional longer smoke test, using a different new run directory:
bash testing/chemistry_transport/run_supercell_smoke.sh BASELINE_DIR NEW_RUN_DIR 600
```

The runner defaults to 60 simulated seconds, copies the current checkout's
`atmosphere_model` (or an explicit fourth argument), and refuses to overwrite
an existing run directory. It records the executable/input hashes, source
revision and model-source diff, MPI exit status, timing, and a checked JSON
summary. It uses the dependency environment recorded with the baseline.
Optionally set `CHEMPAS_SMOKE_REFERENCE` to a saved summary for the same duration
and toolchain to require exact agreement of every reference output field.

`check_supercell_baseline.py` checks the exact duration and step/output counts,
all numeric fields at every saved time, positivity of the chemical species,
photolysis activity, dry-density consistency, and rank error/abort files. It
records per-variable value hashes and dry-air/chemical inventories. The full
baseline also requires lightning activation and nonzero final NO and NO2.
Smoke tests do not require lightning production before storm spin-up and do
not qualify PBL transport. Chemical inventory drift is reported separately
from the isolated transport-budget thresholds above.

Complete #40 only after the supported subset is explicitly listed, all required
tests for that subset pass with GCC/gfortran, remaining configurations fail
clearly, and the result manifest and implementation commit IDs are recorded
here. No transport test results are claimed by this planning commit.

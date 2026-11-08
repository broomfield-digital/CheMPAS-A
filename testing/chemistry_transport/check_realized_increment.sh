#!/usr/bin/env bash
# Prove that a MICM gas tendency built by the transport core is realized once,
# with the dry-mass weight MPAS conserves, in every scalar-advection mode (P2
# of docs/plans/pbl-transport.md; code document §7.3).
#
# Usage: check_realized_increment.sh BASELINE_DIR RUN_ROOT [DOUBLE_MODEL]
#
# 1. Builds the atmosphere model (GNU, MUSICA) from this checkout's tracked
#    files in RUN_ROOT/build, with model_doubles/mpas_atmphys_pbl_gas.F in
#    place of the production stub. The checkout's own build is untouched.
# 2. Makes a horizontally uniform resting state from the baseline
#    supercell_init.nc (make_resting_init.py) and an inert mechanism.
# 3. Runs three 3-s steps for each advection mode and each double mode, and
#    compares realized with requested increments (check_realized_increment.py),
#    in two regimes: "rest" (eddy viscosity off, so the uniform state stays at
#    rest and every layer is compared) and "moving" (the case's eddy viscosity
#    drives a uniform vertical drift with density changes of about 1e-5 per
#    step, and the column tracer-mass budget is compared). A third regime,
#    "fresh", turns NO2 photolysis on at rest and checks that transport acts on
#    the post-chemistry state. Note: MPAS runs config_positive_definite with the
#    monotonic path (it warns that positive-definite advection is
#    unimplemented), so the posdef cases repeat the monotonic path.
# RUN_ROOT must not exist. DOUBLE_MODEL reuses an executable built by an
# earlier run of this script (RUN_ROOT/build/atmosphere_model) instead of step 1.
set -euo pipefail
if (( $# < 2 || $# > 3 )); then
    echo "Usage: $0 BASELINE_DIR RUN_ROOT [DOUBLE_MODEL]" >&2
    exit 2
fi
baseline=$(cd "$1" && pwd)
run_root=$2
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
if [[ -e "$run_root" ]]; then
    echo "Use a new run root; existing runs are preserved: $run_root" >&2
    exit 2
fi
mkdir -p "$run_root/provenance"
run_root=$(cd "$run_root" && pwd)
source "$baseline/environment.sh"

if (( $# == 3 )); then
    model=$(realpath "$3")
    echo "== reusing the test-double model $model"
else
    echo "== build with the PBL gas test double"
    git -C "$repo" rev-parse HEAD > "$run_root/provenance/chempas_revision"
    git -C "$repo" status --short > "$run_root/provenance/source-status"
    git -C "$repo" diff --binary -- src Makefile > "$run_root/provenance/source.patch"
    build=$run_root/build
    mkdir "$build"
    (cd "$repo" && git ls-files -z | xargs -0 cp --parents -t "$build")
    for dir in physics_mmm physics_noaa physics_wrf/files; do
        mkdir -p "$build/src/core_atmosphere/physics/$(dirname "$dir")"
        cp -a "$repo/src/core_atmosphere/physics/$dir" "$build/src/core_atmosphere/physics/$dir"
    done
    cp "$here/model_doubles/mpas_atmphys_pbl_gas.F" "$build/src/core_atmosphere/physics/mpas_atmphys_pbl_gas.F"
    (cd "$build" && make clean CORE=atmosphere > "$run_root/build-clean.log" 2>&1 && \
        make -j8 gfortran CORE=atmosphere PRECISION=double OPENMP=false PIO="$PIO" MUSICA=true \
        AUTOCLEAN=true > "$run_root/build.log" 2>&1)
    model=$build/atmosphere_model
fi
test -x "$model"
sha256sum "$model" > "$run_root/provenance/model-sha256.txt"

echo "== resting initial state and inert mechanism"
python "$here/make_resting_init.py" "$baseline/supercell_init.nc" "$run_root/resting_init.nc"
# MICM needs at least one rate parameter, so the reactions stay but do nothing:
# the Arrhenius titration gets A = 0, and photolysis has j = 0 (TUV-x off and
# config_j_no2_max = 0 in the namelist).
python - "$baseline/lnox_o3.yaml" "$run_root/inert.yaml" <<'PY'
import re
import sys
from pathlib import Path
text = Path(sys.argv[1]).read_text()
text, count = re.subn(r'^(\s*A:\s*)\S+', r'\g<1>0.0', text, flags=re.M)
assert count == 1, 'expected one Arrhenius A'
assert 'type: PHOTOLYSIS' in text
Path(sys.argv[2]).write_text(text)
PY

stage() {
    # stage CASE_DIR ADVECTION REGIME
    local case_dir=$1 adv=$2 regime=$3
    mkdir -p "$case_dir"
    for name in namelist.atmosphere streams.atmosphere stream_list.atmosphere.output \
                supercell.graph.info.part.8 tuvx_no2.json tuvx_upper_atm.csv; do
        cp "$baseline/$name" "$case_dir/$name"
    done
    cp "$run_root/inert.yaml" "$case_dir/lnox_o3.yaml"
    ln -s "$run_root/resting_init.nc" "$case_dir/supercell_init.nc"
    ln -s "$model" "$case_dir/atmosphere_model"
    ln -s "$baseline/data" "$case_dir/data"
    shopt -s nullglob
    for path in "$baseline"/*.TBL "$baseline"/*DATA*; do
        ln -s "$path" "$case_dir/$(basename "$path")"
    done
    shopt -u nullglob
    python - "$case_dir" "$adv" "$regime" <<'PY'
import re
import sys
from pathlib import Path
run, adv, regime = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
nml = run / 'namelist.atmosphere'
text = nml.read_text()

def setopt(name, value, group='nhyd_model'):
    global text
    pattern = rf"^(\s*){name}\s*=.*$"
    if re.search(pattern, text, flags=re.M):
        text = re.sub(pattern, rf"\g<1>{name} = {value}", text, count=1, flags=re.M)
    else:
        text, count = re.subn(rf"^&{group}\s*$", f"&{group}\n    {name} = {value}", text, count=1, flags=re.M)
        assert count == 1, f'no &{group} group'

modes = {
    'mono':           dict(monotonic='.true.',  posdef='.false.', split='.false.', steps='1'),
    'ordinary':       dict(monotonic='.false.', posdef='.false.', split='.false.', steps='1'),
    'posdef':         dict(monotonic='.false.', posdef='.true.',  split='.false.', steps='1'),
    'mono_split':     dict(monotonic='.true.',  posdef='.false.', split='.true.',  steps='3'),
    'ordinary_split': dict(monotonic='.false.', posdef='.false.', split='.true.',  steps='3'),
    'mono_chem2':     dict(monotonic='.true.',  posdef='.false.', split='.false.', steps='1'),
}
m = modes[adv]
setopt('config_run_duration', "'00:00:09'")
setopt('config_monotonic', m['monotonic'])
setopt('config_positive_definite', m['posdef'])
setopt('config_split_dynamics_transport', m['split'])
setopt('config_dynamics_split_steps', m['steps'])
setopt('config_lnox_source_rate', '0.0', group='lnox')
setopt('config_tuvx_config_file', "''", group='photolysis')
setopt('config_tuvx_top_extension', '.false.', group='photolysis')
setopt('config_tuvx_upper_column_mode', "'none'", group='photolysis')
setopt('config_tuvx_extension_file', "''", group='photolysis')
setopt('config_chemistry_interval', '6.0' if adv == 'mono_chem2' else '0.0', group='chemistry')
setopt('config_j_no2_max', '0.01' if regime == 'fresh' else '0.0', group='photolysis')
if regime in ('rest', 'fresh'):
    for name in ('config_h_mom_eddy_visc2', 'config_h_mom_eddy_visc4', 'config_v_mom_eddy_visc2',
                 'config_h_theta_eddy_visc2', 'config_h_theta_eddy_visc4', 'config_v_theta_eddy_visc2'):
        setopt(name, '0.0')
nml.write_text(text)

streams = run / 'streams.atmosphere'
s = streams.read_text()
assert 'output_interval="00:10:00"' in s
streams.write_text(s.replace('output_interval="00:10:00"', 'output_interval="00:00:03"'))
listing = run / 'stream_list.atmosphere.output'
names = [n for n in listing.read_text().split() if n != 'j_jNO2']
listing.write_text('\n'.join(names) + '\n')
PY
}

status=0
run_case() {
    # run_case REGIME ADVECTION MODE [EXTRA CHECK ARGS...]
    local regime=$1 adv=$2 mode=$3 case_dir=$run_root/$1-$2-$3
    local column_only=("${@:4}")
    [[ $regime == moving ]] && column_only=(--column-only)
    stage "$case_dir" "$adv" "$regime"
    set +e
    (cd "$case_dir" && CHEMPAS_TEST_PBL_GAS=$mode mpiexec -x CHEMPAS_TEST_PBL_GAS -n 8 ./atmosphere_model > run.out 2>&1)
    local run_status=$?
    set -e
    if (( run_status != 0 )); then
        echo "FAIL $regime-$adv-$mode: model exited with status $run_status; see $case_dir"
        status=1
        return
    fi
    grep -q "PBL gas test double: mode = $mode" "$case_dir/log.atmosphere.0000.out" || \
        { [[ $mode == off ]] || { echo "FAIL $regime-$adv-$mode: test double did not report mode $mode"; status=1; return; }; }
    [[ $regime == fresh && $mode == off ]] && return
    python "$here/check_realized_increment.py" "$case_dir" --mode "$mode" --init "$run_root/resting_init.nc" \
        "${column_only[@]}" || status=1
}

echo "== realized increments"
for regime in rest moving; do
    for adv in mono ordinary posdef mono_split ordinary_split mono_chem2; do
        for mode in source diffusion; do
            run_case "$regime" "$adv" "$mode"
        done
    done
    run_case "$regime" mono off
done
for adv in mono ordinary; do
    run_case fresh "$adv" off
    run_case fresh "$adv" diffusion --fresh-against "$run_root/fresh-$adv-off"
done
if (( status == 0 )); then
    echo "PASS: realized increments match requested increments in every advection mode"
else
    echo "FAIL: realized increments"
fi
exit "$status"

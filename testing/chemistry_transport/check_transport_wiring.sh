#!/usr/bin/env bash
# Check that the MICM gas transport map reaches physics (P1 of
# docs/plans/pbl-transport.md). Runs short (two-step) supercell cases from a
# frozen baseline and checks the map logged by chemistry and physics:
#   lnox          NO, NO2 and O3 are transported; nMicmTransport = 3
#   lnox_h2o      the mechanism plus an unreacting H2O species: H2O is bound to
#                 qv and excluded; nMicmTransport = 3
#   no_musica     MODEL built with MUSICA=false: no map is published;
#                 nMicmTransport = 0 and the run completes
# Usage: check_transport_wiring.sh BASELINE_DIR RUN_ROOT [MODEL [CASE ...]]
# Default cases: lnox lnox_h2o (MODEL built with MUSICA=true). RUN_ROOT must
# not exist; each case runs in RUN_ROOT/<case>.
# An empty config_micm_file is not a case: MUSICA builds abort on it in
# musica_init, independently of transport.
set -euo pipefail
if (( $# < 2 )); then
    echo "Usage: $0 BASELINE_DIR RUN_ROOT [MODEL [CASE ...]]" >&2
    exit 2
fi
baseline=$(cd "$1" && pwd)
run_root=$2
script_dir=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$script_dir/../.." && pwd)
model=$(realpath "${3:-$repo/atmosphere_model}")
if (( $# > 3 )); then
    cases=("${@:4}")
else
    cases=(lnox lnox_h2o)
fi
test -x "$model"
if [[ -e "$run_root" ]]; then
    echo "Use a new run root; existing runs are preserved: $run_root" >&2
    exit 2
fi
mkdir -p "$run_root"
run_root=$(cd "$run_root" && pwd)
source "$baseline/environment.sh"

stage() {
    local case_dir=$1
    mkdir -p "$case_dir"
    for name in environment.sh namelist.atmosphere streams.atmosphere stream_list.atmosphere.output \
                supercell.graph.info.part.8 lnox_o3.yaml tuvx_no2.json tuvx_upper_atm.csv; do
        cp "$baseline/$name" "$case_dir/$name"
    done
    cp --reflink=auto "$baseline/supercell_init.nc" "$case_dir/"
    ln -s "$model" "$case_dir/atmosphere_model"
    ln -s "$baseline/data" "$case_dir/data"
    shopt -s nullglob
    for path in "$baseline"/*.TBL "$baseline"/*DATA*; do
        ln -s "$path" "$case_dir/$(basename "$path")"
    done
    shopt -u nullglob
}

edit() {
    # edit CASE_DIR CASE: two 3-s steps, no history output, and case-specific inputs.
    python - "$1" "$2" <<'PY'
from pathlib import Path
import re
import sys
run, case = Path(sys.argv[1]), sys.argv[2]
nml = run / 'namelist.atmosphere'
text, count = re.subn(r"config_run_duration\s*=\s*'[^']*'", "config_run_duration = '00:00:06'", nml.read_text())
assert count == 1, 'expected one run duration'
nml.write_text(text)
if case == 'no_musica':
    # j_jNO2 is a runtime chemistry diagnostic; it does not exist without MUSICA.
    listing = run / 'stream_list.atmosphere.output'
    names = [n for n in listing.read_text().split() if n != 'j_jNO2']
    listing.write_text('\n'.join(names) + '\n')
streams = run / 'streams.atmosphere'
text = streams.read_text()
assert 'output_interval="00:10:00"' in text
streams.write_text(text.replace('output_interval="00:10:00"', 'output_interval="none"'))
if case == 'lnox_h2o':
    mech = run / 'lnox_o3.yaml'
    text = mech.read_text()
    species = ('  - name: H2O\n'
               '    __absolute tolerance: 1e-12\n'
               '    __long name: water_vapor\n'
               '    __initial concentration: 0.0\n'
               '    __molar mass: 0.018\n\n'
               'phases:')
    assert text.count('\nphases:') == 1
    text = text.replace('\nphases:', '\n' + species, 1)
    text, count = re.subn(r'(      - O3\n)', r'\1      - H2O\n', text, count=1)
    assert count == 1, 'expected the gas phase species list'
    mech.write_text(text)
PY
}

expect() {
    # expect LOG PATTERN DESCRIPTION
    if grep -qF -- "$2" "$1"; then
        echo "  PASS $3"
    else
        echo "  FAIL $3 (missing: $2)"
        status=1
    fi
}

status=0
for case in "${cases[@]}"; do
    case_dir=$run_root/$case
    stage "$case_dir"
    edit "$case_dir" "$case"
    echo "== $case"
    set +e
    (cd "$case_dir" && mpiexec -n 8 ./atmosphere_model > run.out 2>&1)
    run_status=$?
    set -e
    log=$case_dir/log.atmosphere.0000.out
    if (( run_status != 0 )) || [[ ! -f "$log" ]]; then
        echo "  FAIL model exited with status $run_status; see $case_dir"
        status=1
        continue
    fi
    expect "$log" "Critical error messages =            0" "no critical errors"
    case $case in
        lnox)
            for sp in NO NO2 O3; do
                expect "$log" "[MUSICA] transport: $sp -> scalars(" "$sp is transported"
            done
            expect "$log" "[transport] nMicmTransport = 3" "physics receives three gases"
            ;;
        lnox_h2o)
            expect "$log" "[MUSICA] transport: exclude H2O (host-bound to qv)" "H2O is excluded"
            expect "$log" "[transport] nMicmTransport = 3" "physics receives three gases"
            ;;
        no_musica)
            expect "$log" "[transport] nMicmTransport = 0" "physics receives no gases"
            ;;
        *)
            echo "  FAIL unknown case $case"
            status=1
            ;;
    esac
done
if (( status == 0 )); then
    echo "PASS: transport wiring"
else
    echo "FAIL: transport wiring"
fi
exit "$status"

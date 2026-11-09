#!/usr/bin/env bash
# MICM gas PBL transport at nondefault PBL intervals (P5 Task 3 of
# docs/plans/pbl-transport.md). The PBL scheme refreshes the gas operator when
# the pbl alarm rings; transport applies the held operator every model step
# with the model dt, as MPAS applies the held meteorological PBL tendencies.
# Supercell, 60 s (20 steps of 3 s), stable init, inert mechanism (pbl_cases.sh),
# budget logging on.
#
# Usage: check_pbl_intervals.sh BASELINE_DIR RUN_ROOT [MODEL]
#
# Cases and checks:
#   ysu-9-off/on, mynn-9-off/on   config_pbl_interval = 9 s (3 steps), transport
#       off and on: check_pbl_ysu.py mixing (non-gas fields bitwise, gases finite
#       and nonnegative, domain mass <= 1e-13, down-gradient mixing); the budget
#       log of the transport run (check_transport_budget.py) has one PBL line per
#       gas and step, closed-column conservation, dt = 3 s, and an operator
#       refreshed every 9 s and never older than that.
#   ysu-10-on   10-s interval, not a multiple of dt: the alarm rings at 0, 12, 21,
#       30, ... s, and the held operator stays within its interval.
#   pbloff-9-on   PBL off with transport requested and a 9-s interval: a
#       zero-transport case that runs and logs transport as not enabled.
# RUN_ROOT must not exist.
set -euo pipefail
if (( $# < 2 || $# > 3 )); then
    echo "Usage: $0 BASELINE_DIR RUN_ROOT [MODEL]" >&2
    exit 2
fi
baseline=$(cd "$1" && pwd)
run_root=$2
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
model=$(realpath "${3:-$repo/atmosphere_model}")
test -x "$model"
if [[ -e "$run_root" ]]; then
    echo "Use a new run root; existing runs are preserved: $run_root" >&2
    exit 2
fi
mkdir -p "$run_root/provenance"
run_root=$(cd "$run_root" && pwd)
source "$baseline/environment.sh"
git -C "$repo" rev-parse HEAD > "$run_root/provenance/chempas_revision"
git -C "$repo" status --short > "$run_root/provenance/source-status"
sha256sum "$model" > "$run_root/provenance/model-sha256.txt"

PBL_DEFAULT=bl_ysu
SFCLAYER_DEFAULT=sf_monin_obukhov
# shellcheck source=pbl_cases.sh
source "$here/pbl_cases.sh"
prepare_inputs

init=$run_root/stable_init.nc
mynn=(pbl=bl_mynn sfclayer=sf_mynn)
stage ysu-9-off "$init" "$model" 00:01:00 transport=.false. pbl_interval=00:00:09
stage ysu-9-on "$init" "$model" 00:01:00 transport=.true. pbl_interval=00:00:09 budget=.true.
stage mynn-9-off "$init" "$model" 00:01:00 "${mynn[@]}" transport=.false. pbl_interval=00:00:09
stage mynn-9-on "$init" "$model" 00:01:00 "${mynn[@]}" transport=.true. pbl_interval=00:00:09 budget=.true.
stage ysu-10-on "$init" "$model" 00:01:00 transport=.true. pbl_interval=00:00:10 budget=.true.
stage pbloff-9-on "$init" "$model" 00:01:00 pbl=off sfclayer=off transport=.true. pbl_interval=00:00:09
cases=(ysu-9-off ysu-9-on mynn-9-off mynn-9-on ysu-10-on pbloff-9-on)
ran=1
for (( n = 0; n < ${#cases[@]}; n += 2 )); do
    pids=()
    for c in "${cases[@]:n:2}"; do
        run_model "$c" & pids+=("$!")
    done
    m=0
    for c in "${cases[@]:n:2}"; do
        wait "${pids[m]}" || { fail "$c: the model exited with an error; see $run_root/$c"; ran=0; }
        m=$((m + 1))
    done
done

budget() {
    # budget CASE INTERVAL
    python "$here/check_transport_budget.py" "$run_root/$1/log.atmosphere.0000.out" --expect-gases 3 \
        --expect-steps 20 --dt 3 --interval "$2" --json "$run_root/$1-budget.json" || status=1
}

if (( ran )); then
    for c in ysu-9-on mynn-9-on ysu-10-on; do
        has_log "$c" "PBL gas transport enabled = T" || fail "$c: PBL gas transport not logged as enabled"
    done
    echo "== YSU, 9-s PBL interval"
    python "$here/check_pbl_ysu.py" mixing "$run_root/ysu-9-on" "$run_root/ysu-9-off" \
        --json "$run_root/ysu-9.json" || status=1
    budget ysu-9-on 9
    echo "== MYNN, 9-s PBL interval"
    python "$here/check_pbl_ysu.py" mixing "$run_root/mynn-9-on" "$run_root/mynn-9-off" \
        --json "$run_root/mynn-9.json" || status=1
    budget mynn-9-on 9
    echo "== YSU, 10-s PBL interval (not a multiple of the 3-s step)"
    budget ysu-10-on 10
    echo "== PBL off, transport requested, 9-s interval"
    if has_log pbloff-9-on "PBL gas transport enabled = F"; then
        echo "PASS pbloff-9-on: runs; PBL gas transport logged as not enabled"
    else
        fail "pbloff-9-on: transport not logged as disabled"
    fi
fi

if (( status == 0 )); then
    echo "PASS: PBL gas transport at nondefault intervals"
else
    echo "FAIL: PBL gas transport at nondefault intervals"
fi
exit "$status"

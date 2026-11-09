#!/usr/bin/env bash
# Restart equivalence for MICM gas PBL transport (P5 Task 4 of
# docs/plans/pbl-transport.md). The gas operator is not in the restart file:
# the pbl alarm rings at a restart's first step, so the PBL scheme rebuilds the
# operator before transport uses it (pbl_gas_begin_step stops the model
# otherwise). Supercell, stable init, inert mechanism (pbl_cases.sh), transport
# and budget logging on, history and restart every 72 s.
#
# Usage: check_pbl_restart.sh BASELINE_DIR RUN_ROOT [MODEL]
#
# Pairs (continuous 0-144 s, writing a restart at 72 s; restarted 72-144 s from
# it), every output field at 0000-01-01_18:02:24 bitwise identical
# (check_pbl_ysu.py compare-at):
#   ysu-none, ysu-9, mynn-none, mynn-9     PBL interval 'none' or 9 s (72 s is on
#                                          the 9-s boundary)
#   native-ysu-none                        transport off: MPAS's own restart, the
#                                          reference for the pairs above
# Reported, not required to match: ysu-9-offb, restarted at 60 s, off the 9-s
# boundary (the pbl alarm starts again at the restart time, so the PBL cadence
# shifts, as radiation's does); it must run to completion.
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

run_pairs() {
    # run_pairs CASE...: run the cases two at a time; returns 1 if any model failed
    local rc=0 n m
    local cases=("$@")
    for (( n = 0; n < ${#cases[@]}; n += 2 )); do
        local pids=()
        for c in "${cases[@]:n:2}"; do
            run_model "$c" & pids+=("$!")
        done
        m=0
        for c in "${cases[@]:n:2}"; do
            wait "${pids[m]}" || { fail "$c: the model exited with an error; see $run_root/$c"; rc=1; }
            m=$((m + 1))
        done
    done
    return $rc
}

declare -A opts=(
    [ysu-none]="transport=.true. budget=.true."
    [ysu-9]="transport=.true. budget=.true. pbl_interval=00:00:09"
    [mynn-none]="pbl=bl_mynn sfclayer=sf_mynn transport=.true. budget=.true."
    [mynn-9]="pbl=bl_mynn sfclayer=sf_mynn transport=.true. budget=.true. pbl_interval=00:00:09"
    [native-ysu-none]="transport=.false."
)
pairs=(ysu-none ysu-9 mynn-none mynn-9 native-ysu-none)

echo "== continuous runs"
continuous=()
for p in "${pairs[@]}"; do
    # shellcheck disable=SC2086
    stage "$p-cont" "$init" "$model" 00:02:24 ${opts[$p]} output_interval=00:01:12 restart_interval=00:01:12
    continuous+=("$p-cont")
done
# shellcheck disable=SC2086
stage ysu-9-offb-cont "$init" "$model" 00:02:00 ${opts[ysu-9]} output_interval=00:01:00 restart_interval=00:01:00
continuous+=(ysu-9-offb-cont)
ran=1
run_pairs "${continuous[@]}" || ran=0

if (( ran )); then
    echo "== restarted runs"
    restarted=()
    for p in "${pairs[@]}"; do
        # shellcheck disable=SC2086
        stage "$p-rest" "$init" "$model" 00:01:12 ${opts[$p]} output_interval=00:01:12 \
            restart_interval=00:01:12 restart_from="$run_root/$p-cont/restart.0000-01-01_18.01.12.nc"
        restarted+=("$p-rest")
    done
    # shellcheck disable=SC2086
    stage ysu-9-offb-rest "$init" "$model" 00:01:00 ${opts[ysu-9]} output_interval=00:01:00 \
        restart_interval=00:01:00 restart_from="$run_root/ysu-9-offb-cont/restart.0000-01-01_18.01.00.nc"
    restarted+=(ysu-9-offb-rest)
    run_pairs "${restarted[@]}" || ran=0
fi

if (( ran )); then
    for p in "${pairs[@]}"; do
        echo "== $p: restarted against continuous at 0000-01-01_18:02:24"
        python "$here/check_pbl_ysu.py" compare-at "$run_root/$p-cont" "$run_root/$p-rest" 0000-01-01_18:02:24 \
            --json "$run_root/$p.json" || status=1
    done
    for p in ysu-none ysu-9 mynn-none mynn-9; do
        has_log "$p-rest" "PBL gas transport enabled = T" || fail "$p-rest: transport not logged as enabled"
        python "$here/check_transport_budget.py" "$run_root/$p-rest/log.atmosphere.0000.out" --expect-gases 3 \
            --expect-steps 24 --dt 3 --interval "$([[ $p == *-9 ]] && echo 9 || echo 3)" \
            --json "$run_root/$p-rest-budget.json" > /dev/null || fail "$p-rest: budget log"
    done
    echo "== ysu-9-offb: restarted off the 9-s boundary, at 60 s (reported only)"
    python "$here/check_pbl_ysu.py" compare-at "$run_root/ysu-9-offb-cont" "$run_root/ysu-9-offb-rest" \
        0000-01-01_18:02:00 --json "$run_root/ysu-9-offb.json" || \
        echo "NOTE ysu-9-offb differs, as expected when the PBL cadence shifts; see $run_root/ysu-9-offb.json"
fi

if (( status == 0 )); then
    echo "PASS: restart"
else
    echo "FAIL: restart"
fi
exit "$status"

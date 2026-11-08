#!/usr/bin/env bash
# PBL call cadence (P5a of docs/plans/pbl-transport.md). config_pbl_interval
# sets dt_pbl, the time step the surface layer, land surface and PBL schemes
# integrate with; those schemes must then run only when the interval has
# elapsed (WRF's bldt), holding their tendencies in between. Supercell, 60 s
# (20 steps of 3 s), gas transport off, inert mechanism (pbl_cases.sh).
#
# Usage: check_pbl_cadence.sh BASELINE_DIR RUN_ROOT [MODEL [PRE_P5_MODEL]]
#
# Cases and checks (call counts from the timer summary of rank 0's log):
#   ysu-none   config_pbl_interval = 'none': bl_ysu and sf_monin_obukhov run on
#              all 20 steps
#   ysu-dt     interval = config_dt (3 s): 20 calls, all output bitwise identical
#              to ysu-none
#   ysu-3dt    interval = 9 s: 7 calls (steps 1, 4, ..., 19), DT_PBL = 9 logged
#   mynn-3dt   MYNN with interval = 9 s: bl_mynn and sf_mynn run 7 times
#   ysu-none-pre-p5  ysu-none with PRE_P5_MODEL: all output bitwise identical
# RUN_ROOT must not exist.
set -euo pipefail
if (( $# < 2 || $# > 4 )); then
    echo "Usage: $0 BASELINE_DIR RUN_ROOT [MODEL [PRE_P5_MODEL]]" >&2
    exit 2
fi
baseline=$(cd "$1" && pwd)
run_root=$2
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
model=$(realpath "${3:-$repo/atmosphere_model}")
pre_p5_model=${4:+$(realpath "$4")}
test -x "$model"
[[ -z "$pre_p5_model" ]] || test -x "$pre_p5_model"
if [[ -e "$run_root" ]]; then
    echo "Use a new run root; existing runs are preserved: $run_root" >&2
    exit 2
fi
mkdir -p "$run_root/provenance"
run_root=$(cd "$run_root" && pwd)
source "$baseline/environment.sh"
git -C "$repo" rev-parse HEAD > "$run_root/provenance/chempas_revision"
git -C "$repo" status --short > "$run_root/provenance/source-status"
sha256sum "$model" ${pre_p5_model:+"$pre_p5_model"} > "$run_root/provenance/model-sha256.txt"

PBL_DEFAULT=bl_ysu
SFCLAYER_DEFAULT=sf_monin_obukhov
# shellcheck source=pbl_cases.sh
source "$here/pbl_cases.sh"
prepare_inputs

stage ysu-none "$run_root/stable_init.nc" "$model" 00:01:00
stage ysu-dt "$run_root/stable_init.nc" "$model" 00:01:00 pbl_interval=00:00:03
stage ysu-3dt "$run_root/stable_init.nc" "$model" 00:01:00 pbl_interval=00:00:09
stage mynn-3dt "$run_root/stable_init.nc" "$model" 00:01:00 pbl=bl_mynn sfclayer=sf_mynn pbl_interval=00:00:09
cases=(ysu-none ysu-dt ysu-3dt mynn-3dt)
if [[ -n "$pre_p5_model" ]]; then
    stage ysu-none-pre-p5 "$run_root/stable_init.nc" "$pre_p5_model" 00:01:00
    cases+=(ysu-none-pre-p5)
fi
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

calls() {
    # calls CASE TIMER: the call count of TIMER in the case's timer summary
    awk -v t="$2" '$2 == t {print $4; exit}' "$run_root/$1/log.atmosphere.0000.out"
}
expect_calls() {
    # expect_calls CASE TIMER N
    local got
    got=$(calls "$1" "$2")
    if [[ "$got" == "$3" ]]; then
        echo "PASS $1: $2 ran $3 times"
    else
        fail "$1: $2 ran ${got:-0} times, expected $3"
    fi
}

if (( ran )); then
    echo "== call cadence"
    expect_calls ysu-none bl_ysu 20
    expect_calls ysu-none sf_monin_obukhov 20
    expect_calls ysu-dt bl_ysu 20
    expect_calls ysu-dt sf_monin_obukhov 20
    expect_calls ysu-3dt bl_ysu 7
    expect_calls ysu-3dt sf_monin_obukhov 7
    expect_calls mynn-3dt bl_mynn 7
    expect_calls mynn-3dt sf_mynn 7
    if has_log ysu-3dt "DT_PBL    = 9.00000"; then
        echo "PASS ysu-3dt: DT_PBL = 9 s"
    else
        fail "ysu-3dt: DT_PBL = 9 s not logged"
    fi
    echo "== an interval equal to config_dt is the default cadence"
    python "$here/check_pbl_ysu.py" compare "$run_root/ysu-dt" "$run_root/ysu-none" || status=1
    echo "== the default cadence is unchanged from the pre-P5 executable"
    if [[ -n "$pre_p5_model" ]]; then
        python "$here/check_pbl_ysu.py" compare "$run_root/ysu-none" "$run_root/ysu-none-pre-p5" || status=1
    else
        echo "SKIP no PRE_P5_MODEL given"
    fi
fi

if (( status == 0 )); then
    echo "PASS: PBL cadence"
else
    echo "FAIL: PBL cadence"
fi
exit "$status"

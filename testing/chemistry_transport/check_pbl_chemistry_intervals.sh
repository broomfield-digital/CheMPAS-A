#!/usr/bin/env bash
# MICM gas PBL transport with chemistry intervals (P5 Task 5 of
# docs/plans/pbl-transport.md). Chemistry runs every config_chemistry_interval
# seconds; PBL transport is assembled every model step from the current gases,
# whether chemistry ran in that step or not. (The realized increment, applied
# exactly once per step, is proved with the test double in
# check_realized_increment.sh, modes mono_chem2, mono_chem3 and mono_chem7.)
#
# Usage: check_pbl_chemistry_intervals.sh BASELINE_DIR RUN_ROOT [MODEL]
#
# Supercell, stable init, the reactive LNOx mechanism (NO + O3 -> NO2 with the
# baseline rate, NO2 photolysis with config_j_no2_max = 0.01, TUV-x and
# lightning off), 63 s (21 steps of 3 s). YSU and MYNN, each with chemistry
# every 3, 9 and 21 s (1, 3 and 7 steps), transport off and on, budget logging
# on. Checks for each pair:
#   - every non-gas output field bitwise identical at the final record;
#   - with transport, nitrogen (NO + NO2) and odd oxygen (NO2 + O3), in moles
#     weighted by the dry air mass, conserved within 1e-12 relative (both are
#     invariants of the mechanism, and transport conserves each gas), and every
#     gas finite and above the nonnegativity floor;
#   - the budget log has one PBL assembly per gas and step (21 steps), each
#     conservative, with dt = 3 s.
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
sha256sum "$model" "$baseline/lnox_o3.yaml" > "$run_root/provenance/sha256.txt"

PBL_DEFAULT=bl_ysu
SFCLAYER_DEFAULT=sf_monin_obukhov
# shellcheck source=pbl_cases.sh
source "$here/pbl_cases.sh"
prepare_inputs
init=$run_root/stable_init.nc
reactive=(mechanism="$baseline/lnox_o3.yaml" j_no2_max=0.01 output_interval=00:00:21)

cases=()
for scheme in ysu mynn; do
    pbl=()
    [[ $scheme == mynn ]] && pbl=(pbl=bl_mynn sfclayer=sf_mynn)
    for ci in 3 9 21; do
        stage "$scheme-chem$ci-off" "$init" "$model" 00:01:03 "${pbl[@]}" "${reactive[@]}" \
            chemistry_interval="$ci.0" transport=.false.
        stage "$scheme-chem$ci-on" "$init" "$model" 00:01:03 "${pbl[@]}" "${reactive[@]}" \
            chemistry_interval="$ci.0" transport=.true. budget=.true.
        cases+=("$scheme-chem$ci-off" "$scheme-chem$ci-on")
    done
done
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

if (( ran )); then
    for scheme in ysu mynn; do
        for ci in 3 9 21; do
            c=$scheme-chem$ci
            echo "== $c: chemistry every $ci s"
            has_log "$c-on" "PBL gas transport enabled = T" || fail "$c-on: transport not logged as enabled"
            python "$here/check_pbl_ysu.py" compare-at "$run_root/$c-off" "$run_root/$c-on" 0000-01-01_18:01:03 \
                --exclude qNO qNO2 qO3 --json "$run_root/$c-nongas.json" || status=1
            python "$here/check_pbl_reactive.py" invariants "$run_root/$c-on" \
                --invariant N=NO:0.030,NO2:0.046 --invariant Ox=NO2:0.046,O3:0.048 \
                --nonnegative NO NO2 O3 --json "$run_root/$c-invariants.json" || status=1
            python "$here/check_transport_budget.py" "$run_root/$c-on/log.atmosphere.0000.out" --expect-gases 3 \
                --expect-steps 21 --dt 3 --interval 3 --json "$run_root/$c-budget.json" || status=1
        done
    done
fi

if (( status == 0 )); then
    echo "PASS: chemistry intervals"
else
    echo "FAIL: chemistry intervals"
fi
exit "$status"

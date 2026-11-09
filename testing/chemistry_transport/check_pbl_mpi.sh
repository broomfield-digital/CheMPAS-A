#!/usr/bin/env bash
# MPI decomposition check for MICM gas PBL transport (P5 Task 6 of
# docs/plans/pbl-transport.md). The transport tendencies are built per owned
# column and added once per mapped scalar; the scalar advance owns the halos.
#
# Usage: check_pbl_mpi.sh BASELINE_DIR RUN_ROOT [MODEL]
#
# YSU and MYNN at a 9-s PBL interval, transport on, stable init, inert
# mechanism, 60 s, on 1, 2, 4 and 8 ranks (the four rank counts of a scheme run
# together, 15 ranks). Against 8 ranks at 0000-01-01_18:01:00, with cells
# aligned by indexToCellID:
#   - the gases agree within atol = 1e-18 kg/kg, rtol = 1e-10 (the spec's
#     serial/eight-rank thresholds): required;
#   - every output field: reported (largest absolute and relative difference),
#     in RUN_ROOT/<scheme>-<N>-all.json.
# A rank whose columns are all inactive cannot be made on this mesh (YSU and
# MYNN mix every column); unit/test_tracer_transport_ops.F90 covers it.
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

ran=1
for scheme in ysu mynn; do
    pbl=()
    [[ $scheme == mynn ]] && pbl=(pbl=bl_mynn sfclayer=sf_mynn)
    pids=()
    for n in 1 2 4 8; do
        stage "$scheme-$n" "$init" "$model" 00:01:00 "${pbl[@]}" pbl_interval=00:00:09 transport=.true. ranks="$n"
        run_model "$scheme-$n" & pids+=("$!")
    done
    m=0
    for n in 1 2 4 8; do
        wait "${pids[m]}" || { fail "$scheme-$n: the model exited with an error; see $run_root/$scheme-$n"; ran=0; }
        m=$((m + 1))
    done
done

if (( ran )); then
    for scheme in ysu mynn; do
        for n in 1 2 4; do
            echo "== $scheme: $n ranks against 8"
            python "$here/check_pbl_ysu.py" compare-at "$run_root/$scheme-8" "$run_root/$scheme-$n" \
                0000-01-01_18:01:00 --gases-only --atol 1e-18 --rtol 1e-10 \
                --json "$run_root/$scheme-$n-gases.json" || status=1
            python "$here/check_pbl_ysu.py" compare-at "$run_root/$scheme-8" "$run_root/$scheme-$n" \
                0000-01-01_18:01:00 --json "$run_root/$scheme-$n-all.json" > "$run_root/$scheme-$n-all.txt" || true
            if grep -q 'differing: none' "$run_root/$scheme-$n-all.txt"; then
                echo "NOTE every output field is bitwise identical"
            else
                echo "NOTE fields that differ (reported, not required):"
                grep -E '^  ' "$run_root/$scheme-$n-all.txt" || true
            fi
        done
    done
fi

if (( status == 0 )); then
    echo "PASS: decomposition"
else
    echo "FAIL: decomposition"
fi
exit "$status"

#!/usr/bin/env bash
# Reactive regression with MICM gas PBL transport (docs/plans/pbl-transport.md,
# "Coupled/regression"; P6): the transport-only budget each step, then the
# reaction invariants. Supercell, 8 ranks, the baseline LNOx mechanism
# (lnox_o3.yaml) with NO2 photolysis (config_j_no2_max = 0.01, TUV-x off).
#
# Usage: check_pbl_reactive.sh BASELINE_DIR RUN_ROOT [MODEL]
#
# Without lightning, the stable init with vertical gas gradients (make_pbl_init.py:
# surface NO, NO2 below 1 km, O3 depleted near the surface), YSU and MYNN,
# transport off and on, 5 min:
#   - every non-gas output field bitwise identical at the final record;
#   - with transport, nitrogen NO + NO2 and odd oxygen NO2 + O3 (moles, dry-air
#     weighted) conserved within 1e-12 relative, each gas above the
#     nonnegativity floor;
#   - the budget log: conservative transport every step (100 steps).
# With lightning, the baseline supercell init (NO = NO2 = 0, O3 = 50 ppbv), YSU,
# transport off and on, 30 min (lightning starts after about 10 min):
#   - non-gas fields bitwise identical at the final record;
#   - odd oxygen NO2 + O3 conserved within 1e-12 (lightning adds NO only; the
#     reactions and transport conserve it), every gas above the floor;
#   - the budget log: conservative transport every step (600 steps).
# (The Chapman+NOx case runs on the global JW mesh without a land surface; at
# that spacing it needs parameterized convection, outside the convection-
# resolving configurations qualified here.)
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
lnox=(mechanism="$baseline/lnox_o3.yaml" j_no2_max=0.01)
cases=()
for scheme in ysu mynn; do
    pbl=()
    [[ $scheme == mynn ]] && pbl=(pbl=bl_mynn sfclayer=sf_mynn)
    stage "nolightning-$scheme-off" "$run_root/stable_init.nc" "$model" 00:05:00 "${pbl[@]}" "${lnox[@]}" \
        transport=.false.
    stage "nolightning-$scheme-on" "$run_root/stable_init.nc" "$model" 00:05:00 "${pbl[@]}" "${lnox[@]}" \
        transport=.true. budget=.true.
    cases+=("nolightning-$scheme-off" "nolightning-$scheme-on")
done
lightning=("${lnox[@]}" lnox_source_rate=0.5 output_interval=00:10:00)
stage lightning-ysu-off "$baseline/supercell_init.nc" "$model" 00:30:00 "${lightning[@]}" transport=.false.
stage lightning-ysu-on "$baseline/supercell_init.nc" "$model" 00:30:00 "${lightning[@]}" transport=.true. budget=.true.
cases+=(lightning-ysu-off lightning-ysu-on)
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
        c=nolightning-$scheme
        echo "== $c"
        python "$here/check_pbl_ysu.py" compare-at "$run_root/$c-off" "$run_root/$c-on" 0000-01-01_18:05:00 \
            --exclude qNO qNO2 qO3 j_jNO2 --json "$run_root/$c-nongas.json" || status=1
        python "$here/check_pbl_reactive.py" invariants "$run_root/$c-on" \
            --invariant N=NO:0.030,NO2:0.046 --invariant Ox=NO2:0.046,O3:0.048 \
            --nonnegative NO NO2 O3 --json "$run_root/$c-invariants.json" || status=1
        python "$here/check_transport_budget.py" "$run_root/$c-on/log.atmosphere.0000.out" --expect-gases 3 \
            --expect-steps 100 --dt 3 --interval 3 --json "$run_root/$c-budget.json" || status=1
    done
    echo "== lightning-ysu"
    python "$here/check_pbl_ysu.py" compare-at "$run_root/lightning-ysu-off" "$run_root/lightning-ysu-on" \
        0000-01-01_18:30:00 --exclude qNO qNO2 qO3 j_jNO2 --json "$run_root/lightning-ysu-nongas.json" || status=1
    python "$here/check_pbl_reactive.py" invariants "$run_root/lightning-ysu-on" --invariant Ox=NO2:0.046,O3:0.048 \
        --nonnegative NO NO2 O3 --json "$run_root/lightning-ysu-invariants.json" || status=1
    python "$here/check_transport_budget.py" "$run_root/lightning-ysu-on/log.atmosphere.0000.out" --expect-gases 3 \
        --expect-steps 600 --dt 3 --interval 3 --json "$run_root/lightning-ysu-budget.json" || status=1
fi

if (( status == 0 )); then
    echo "PASS: reactive regression with PBL transport"
else
    echo "FAIL: reactive regression with PBL transport"
fi
exit "$status"

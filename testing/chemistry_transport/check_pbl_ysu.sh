#!/usr/bin/env bash
# In-model qualification of YSU PBL mixing of MICM gases (P3 of
# docs/plans/pbl-transport.md), on the supercell with YSU and the
# Monin-Obukhov surface layer, an inert mechanism (titration A = 0, j = 0,
# TUV-x off) and lightning off. Gas profiles with vertical structure come from
# make_pbl_init.py; the gases are passive.
#
# Usage: check_pbl_ysu.sh BASELINE_DIR RUN_ROOT [MODEL [PRE_P3_MODEL]]
#
# Cases (MYNN's gate cases are in check_pbl_mynn.sh):
#   gate-interval              transport requested with a nondefault PBL
#                              interval stops the model with the gate's message
#   gate-no-sfclayer           YSU without a surface layer stops in MPAS's own
#                              physics_namelist_check, before transport init
#                              (the gas gate's message is unit-tested)
#   gate-pbl-off               transport requested with PBL off runs as a
#                              logged zero-transport case
#   stable-on / stable-off     baseline surface (288 K, colder than the air):
#                              with and without gas transport, 300 s
#   convective-on / -off       skintemp = sst = 305 K: upward heat flux
#   convective-topdown-on/-off the convective case with config_ysu_pblmix =
#                              true. That mode raises YSU's PBL height with the
#                              liquid-ice potential temperature, changing the
#                              meteorology; its radiatively driven cloud-top
#                              entrainment needs cooling, and the supercell has
#                              no radiation (the column reference tests cover
#                              that branch)
#   stable-off-pre-p3          stable-off with PRE_P3_MODEL (the executable
#                              before P3): all output bitwise identical
# Checks are in check_pbl_ysu.py: each on/off pair has bitwise identical
# non-gas fields, hpbl > 0, conserved gas mass and down-gradient mixing of
# every gas; the convective case mixes the gases whose gradients lie inside
# the boundary layer (NO, O3) more strongly than the stable case. RUN_ROOT must
# not exist.
set -euo pipefail
if (( $# < 2 || $# > 4 )); then
    echo "Usage: $0 BASELINE_DIR RUN_ROOT [MODEL [PRE_P3_MODEL]]" >&2
    exit 2
fi
baseline=$(cd "$1" && pwd)
run_root=$2
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
model=$(realpath "${3:-$repo/atmosphere_model}")
pre_p3_model=${4:+$(realpath "$4")}
test -x "$model"
[[ -z "$pre_p3_model" ]] || test -x "$pre_p3_model"
if [[ -e "$run_root" ]]; then
    echo "Use a new run root; existing runs are preserved: $run_root" >&2
    exit 2
fi
mkdir -p "$run_root/provenance"
run_root=$(cd "$run_root" && pwd)
source "$baseline/environment.sh"
git -C "$repo" rev-parse HEAD > "$run_root/provenance/chempas_revision"
git -C "$repo" status --short > "$run_root/provenance/source-status"
sha256sum "$model" ${pre_p3_model:+"$pre_p3_model"} > "$run_root/provenance/model-sha256.txt"

PBL_DEFAULT=bl_ysu
SFCLAYER_DEFAULT=sf_monin_obukhov
# shellcheck source=pbl_cases.sh
source "$here/pbl_cases.sh"
prepare_inputs

echo "== gate cases"
stage gate-interval "$run_root/stable_init.nc" "$model" 00:00:06 transport=.true. pbl_interval=00:00:06
expect_stop gate-interval "config_pbl_interval must be 'none'"
stage gate-no-sfclayer "$run_root/stable_init.nc" "$model" 00:00:06 sfclayer=off transport=.true.
expect_stop gate-no-sfclayer "wrong choice for surface layer scheme with YSU PBL"
stage gate-pbl-off "$run_root/stable_init.nc" "$model" 00:00:06 pbl=off sfclayer=off transport=.true.
if ! run_model gate-pbl-off; then
    fail "gate-pbl-off: the model stopped; see $run_root/gate-pbl-off"
elif has_log gate-pbl-off "PBL: resolved scheme = off" && \
     has_log gate-pbl-off "PBL gas transport enabled = F"; then
    echo "PASS gate-pbl-off: runs; PBL gas transport logged as not enabled"
else
    fail "gate-pbl-off: missing transport log lines"
fi

echo "== YSU runs (300 s, output every 60 s)"
stage stable-off "$run_root/stable_init.nc" "$model" 00:05:00
stage stable-on "$run_root/stable_init.nc" "$model" 00:05:00 transport=.true.
stage convective-off "$run_root/convective_init.nc" "$model" 00:05:00 transport=.false.
stage convective-on "$run_root/convective_init.nc" "$model" 00:05:00 transport=.true.
stage convective-topdown-off "$run_root/convective_init.nc" "$model" 00:05:00 transport=.false. pblmix=.true.
stage convective-topdown-on "$run_root/convective_init.nc" "$model" 00:05:00 transport=.true. pblmix=.true.
cases=(stable-off stable-on convective-off convective-on convective-topdown-off convective-topdown-on)
if [[ -n "$pre_p3_model" ]]; then
    stage stable-off-pre-p3 "$run_root/stable_init.nc" "$pre_p3_model" 00:05:00
    cases+=(stable-off-pre-p3)
fi
ran=1
for c in "${cases[@]}"; do
    run_model "$c" || { fail "$c: the model exited with an error; see $run_root/$c"; ran=0; }
done
for c in stable-on convective-on convective-topdown-on; do
    has_log "$c" "PBL gas transport enabled = T" || fail "$c: PBL gas transport not logged as enabled"
done
has_log stable-off "PBL gas transport enabled = F" || fail "stable-off: transport logged as enabled"

if (( ran )); then
    echo "== stable surface layer"
    python "$here/check_pbl_ysu.py" mixing "$run_root/stable-on" "$run_root/stable-off" \
        --json "$run_root/stable.json" || status=1
    echo "== convective surface layer"
    python "$here/check_pbl_ysu.py" mixing "$run_root/convective-on" "$run_root/convective-off" \
        --json "$run_root/convective.json" || status=1
    # NO2's step (about 1.2 km) lies above both PBL tops: both regimes mix it with YSU's
    # free-atmosphere coefficient, so only the in-PBL gradients are compared.
    python "$here/check_pbl_ysu.py" stronger "$run_root/convective.json" "$run_root/stable.json" \
        --gases qNO qO3 || status=1
    echo "== convective, top-down mixing on (config_ysu_pblmix)"
    python "$here/check_pbl_ysu.py" mixing "$run_root/convective-topdown-on" "$run_root/convective-topdown-off" \
        --json "$run_root/convective-topdown.json" || status=1
    echo "== transport off is unchanged from the pre-P3 executable"
    if [[ -n "$pre_p3_model" ]]; then
        python "$here/check_pbl_ysu.py" compare "$run_root/stable-off" "$run_root/stable-off-pre-p3" || status=1
    else
        echo "SKIP no PRE_P3_MODEL given"
    fi
fi

if (( status == 0 )); then
    echo "PASS: YSU gas PBL qualification"
else
    echo "FAIL: YSU gas PBL qualification"
fi
exit "$status"

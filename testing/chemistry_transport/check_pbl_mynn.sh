#!/usr/bin/env bash
# In-model qualification of MYNN PBL mixing of MICM gases (P4 of
# docs/plans/pbl-transport.md), on the supercell with MYNN and its surface
# layer (MPAS selects sf_mynn for bl_mynn), an inert mechanism and lightning
# off (pbl_cases.sh). The gases are passive.
#
# Usage: check_pbl_mynn.sh BASELINE_DIR RUN_ROOT [MODEL [PRE_P4_MODEL]]
#
# Cases:
#   gate-dd                    transport requested with MYNN's downdrafts
#                              (config_mynn_edmf = config_mynn_edmf_dd = 1)
#                              stops the model with the gate's message
#   gate-interval              MYNN with a nondefault PBL interval stops
#   local-stable-on / -off     config_mynn_edmf = 0, baseline surface (288 K):
#                              with and without gas transport, 300 s
#   local-convective-on / -off config_mynn_edmf = 0, skintemp = sst = 305 K
#   edmf-convective-on / -off  config_mynn_edmf = 1 (updraft plumes, the MPAS
#                              default), 305 K; edmf_a is written and must show
#                              plumes in some column
#   edmf-convective-off-pre-p4 edmf-convective-off with PRE_P4_MODEL: bitwise
#   local-stable-off-pre-p4    local-stable-off with PRE_P4_MODEL (the
#                              executable before P4): all output bitwise
#                              identical
# Checks are in check_pbl_ysu.py (scheme-independent): each on/off pair has
# bitwise identical non-gas fields, hpbl > 0, conserved gas mass and
# down-gradient mixing of every gas; the convective case mixes NO and O3 more
# strongly than the stable case. RUN_ROOT must not exist.
set -euo pipefail
if (( $# < 2 || $# > 4 )); then
    echo "Usage: $0 BASELINE_DIR RUN_ROOT [MODEL [PRE_P4_MODEL]]" >&2
    exit 2
fi
baseline=$(cd "$1" && pwd)
run_root=$2
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
model=$(realpath "${3:-$repo/atmosphere_model}")
pre_p4_model=${4:+$(realpath "$4")}
test -x "$model"
[[ -z "$pre_p4_model" ]] || test -x "$pre_p4_model"
if [[ -e "$run_root" ]]; then
    echo "Use a new run root; existing runs are preserved: $run_root" >&2
    exit 2
fi
mkdir -p "$run_root/provenance"
run_root=$(cd "$run_root" && pwd)
source "$baseline/environment.sh"
git -C "$repo" rev-parse HEAD > "$run_root/provenance/chempas_revision"
git -C "$repo" status --short > "$run_root/provenance/source-status"
sha256sum "$model" ${pre_p4_model:+"$pre_p4_model"} > "$run_root/provenance/model-sha256.txt"

PBL_DEFAULT=bl_mynn
SFCLAYER_DEFAULT=sf_mynn
# shellcheck source=pbl_cases.sh
source "$here/pbl_cases.sh"
prepare_inputs

echo "== gate cases"
stage gate-dd "$run_root/stable_init.nc" "$model" 00:00:06 transport=.true. mynn_edmf=1 mynn_edmf_dd=1
expect_stop gate-dd "MYNN downdraft EDMF gas transport is not yet qualified"
stage gate-interval "$run_root/stable_init.nc" "$model" 00:00:06 transport=.true. mynn_edmf=0 pbl_interval=00:00:06
expect_stop gate-interval "config_pbl_interval must be 'none'"

echo "== MYNN runs (300 s, output every 60 s)"
stage local-stable-off "$run_root/stable_init.nc" "$model" 00:05:00 mynn_edmf=0
stage local-stable-on "$run_root/stable_init.nc" "$model" 00:05:00 mynn_edmf=0 transport=.true.
stage local-convective-off "$run_root/convective_init.nc" "$model" 00:05:00 mynn_edmf=0 transport=.false.
stage local-convective-on "$run_root/convective_init.nc" "$model" 00:05:00 mynn_edmf=0 transport=.true.
edmf=(mynn_edmf=1 mynn_edmf_dd=0 mynn_edmf_output=1 extra_output=edmf_a)
stage edmf-convective-off "$run_root/convective_init.nc" "$model" 00:05:00 "${edmf[@]}" transport=.false.
stage edmf-convective-on "$run_root/convective_init.nc" "$model" 00:05:00 "${edmf[@]}" transport=.true.
cases=(local-stable-off local-stable-on local-convective-off local-convective-on edmf-convective-off edmf-convective-on)
if [[ -n "$pre_p4_model" ]]; then
    stage local-stable-off-pre-p4 "$run_root/stable_init.nc" "$pre_p4_model" 00:05:00 mynn_edmf=0
    stage edmf-convective-off-pre-p4 "$run_root/convective_init.nc" "$pre_p4_model" 00:05:00 "${edmf[@]}"
    cases+=(local-stable-off-pre-p4 edmf-convective-off-pre-p4)
fi
ran=1
for c in "${cases[@]}"; do
    run_model "$c" || { fail "$c: the model exited with an error; see $run_root/$c"; ran=0; }
done
for c in local-stable-on local-convective-on edmf-convective-on; do
    has_log "$c" "PBL gas transport enabled = T" || fail "$c: PBL gas transport not logged as enabled"
    has_log "$c" "MYNN gas operator = local diffusion" || fail "$c: MYNN gas operator not logged"
done
has_log edmf-convective-on "+ updraft EDMF" || fail "edmf-convective-on: updraft EDMF not logged"
has_log local-convective-on "+ updraft EDMF" && fail "local-convective-on: updraft EDMF logged with edmf = 0"
has_log local-stable-off "PBL gas transport enabled = F" || fail "local-stable-off: transport logged as enabled"

if (( ran )); then
    echo "== local diffusion, stable surface layer"
    python "$here/check_pbl_ysu.py" mixing "$run_root/local-stable-on" "$run_root/local-stable-off" \
        --json "$run_root/local-stable.json" || status=1
    echo "== local diffusion, convective surface layer"
    python "$here/check_pbl_ysu.py" mixing "$run_root/local-convective-on" "$run_root/local-convective-off" \
        --json "$run_root/local-convective.json" || status=1
    python "$here/check_pbl_ysu.py" stronger "$run_root/local-convective.json" "$run_root/local-stable.json" \
        --gases qNO qO3 || status=1
    echo "== updraft EDMF, convective surface layer"
    python - "$run_root/edmf-convective-on" <<'PY' || status=1
import glob
import sys
import numpy as np
from netCDF4 import Dataset
peak, columns = 0.0, set()
for path in sorted(glob.glob(f"{sys.argv[1]}/output*.nc")):
    with Dataset(path) as ds:
        a = np.asarray(ds["edmf_a"][:])
        peak = max(peak, float(a.max()))
        columns |= set(np.nonzero((a > 0).any(axis=(0, 2)))[0].tolist())
print(f"plume area: peak {peak:.3g}, columns with plumes in some record: {len(columns)}")
if peak <= 0.0:
    print("FAIL edmf-convective-on: MYNN's updraft plumes never active")
    sys.exit(1)
print("PASS edmf-convective-on: MYNN's updraft plumes are active")
PY
    python "$here/check_pbl_ysu.py" mixing "$run_root/edmf-convective-on" "$run_root/edmf-convective-off" \
        --json "$run_root/edmf-convective.json" || status=1
    echo "== transport off is unchanged from the pre-P4 executable"
    if [[ -n "$pre_p4_model" ]]; then
        python "$here/check_pbl_ysu.py" compare "$run_root/local-stable-off" "$run_root/local-stable-off-pre-p4" \
            || status=1
        python "$here/check_pbl_ysu.py" compare "$run_root/edmf-convective-off" \
            "$run_root/edmf-convective-off-pre-p4" || status=1
    else
        echo "SKIP no PRE_P4_MODEL given"
    fi
fi

if (( status == 0 )); then
    echo "PASS: MYNN gas PBL qualification"
else
    echo "FAIL: MYNN gas PBL qualification"
fi
exit "$status"

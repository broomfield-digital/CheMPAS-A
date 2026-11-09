#!/usr/bin/env bash
# In-model qualification of MYNN PBL mixing of MICM gases (P4 of
# docs/plans/pbl-transport.md), on the supercell with MYNN and its surface
# layer (MPAS selects sf_mynn for bl_mynn), an inert mechanism and lightning
# off (pbl_cases.sh). The gases are passive.
#
# Usage: check_pbl_mynn.sh BASELINE_DIR RUN_ROOT [MODEL [PRE_P4_MODEL]]
#
# PRE_P4_MODEL is the executable before P4 with the MYNN snow-content fix
# (commit 0c21a646, "define MYNN's snow content when the microphysics has no
# snow") applied: without it MYNN reads an undefined array in cloudy columns
# of this Kessler case, and its results are not reproducible.
#
# Cases:
#   interval-runs              MYNN gas transport with a nondefault PBL interval
#                              (6 s) runs (P5; qualified in check_pbl_intervals.sh)
#   local-stable-on / -off     config_mynn_edmf = 0, baseline surface (288 K):
#                              with and without gas transport, 300 s
#   local-convective-on / -off config_mynn_edmf = 0, skintemp = sst = 305 K
#   edmf-convective-on / -off  config_mynn_edmf = 1 (updraft plumes, the MPAS
#                              default), 305 K
#   edmfdd-convective-on / -off config_mynn_edmf = config_mynn_edmf_dd = 1
#                              (updrafts and downdrafts), 305 K
#   edmfdd-convective-off-pre-p4 edmfdd-convective-off with PRE_P4_MODEL:
#                              bitwise
# On this grid MYNN's mass flux does not activate: MYNN switches its plumes off
# unless the first one rises a full model level, and with the PBL top at the
# first interface (299 m) the plumes stop in the 413 m layer above it (at
# 315 K with scale awareness off, the PBL stays below 432 m and plume-free for
# 15 simulated minutes). The EDMF cases are therefore integration evidence: the mass-flux
# code paths run with gas transport on, and plume area (edmf_a) and downdraft
# triggers are reported, not required. The column tests
# (unit/test_mynn_gas_operator.F90, unit/test_pbl_edmf.F90) qualify the
# mass-flux operator on MYNN's own plumes.
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
stage interval-runs "$run_root/stable_init.nc" "$model" 00:00:06 transport=.true. pbl_interval=00:00:06
if run_model interval-runs && has_log interval-runs "PBL gas transport enabled = T"; then
    echo "PASS interval-runs: MYNN gas transport runs at a nondefault PBL interval (P5; check_pbl_intervals.sh)"
else
    fail "interval-runs: transport at a nondefault PBL interval did not run; see $run_root/interval-runs"
fi

echo "== MYNN runs (300 s, output every 60 s)"
stage local-stable-off "$run_root/stable_init.nc" "$model" 00:05:00 mynn_edmf=0
stage local-stable-on "$run_root/stable_init.nc" "$model" 00:05:00 mynn_edmf=0 transport=.true.
stage local-convective-off "$run_root/convective_init.nc" "$model" 00:05:00 mynn_edmf=0 transport=.false.
stage local-convective-on "$run_root/convective_init.nc" "$model" 00:05:00 mynn_edmf=0 transport=.true.
edmf=(mynn_edmf=1 mynn_edmf_dd=0 mynn_edmf_output=1 extra_output=edmf_a)
stage edmf-convective-off "$run_root/convective_init.nc" "$model" 00:05:00 "${edmf[@]}" transport=.false.
stage edmf-convective-on "$run_root/convective_init.nc" "$model" 00:05:00 "${edmf[@]}" transport=.true.
edmfdd=(mynn_edmf=1 mynn_edmf_dd=1 mynn_edmf_output=1 extra_output=edmf_a)
stage edmfdd-convective-off "$run_root/convective_init.nc" "$model" 00:05:00 "${edmfdd[@]}" transport=.false.
stage edmfdd-convective-on "$run_root/convective_init.nc" "$model" 00:05:00 "${edmfdd[@]}" transport=.true.
cases=(local-stable-off local-stable-on local-convective-off local-convective-on edmf-convective-off edmf-convective-on
       edmfdd-convective-off edmfdd-convective-on)
if [[ -n "$pre_p4_model" ]]; then
    stage local-stable-off-pre-p4 "$run_root/stable_init.nc" "$pre_p4_model" 00:05:00 mynn_edmf=0
    stage edmfdd-convective-off-pre-p4 "$run_root/convective_init.nc" "$pre_p4_model" 00:05:00 "${edmfdd[@]}"
    cases+=(local-stable-off-pre-p4 edmfdd-convective-off-pre-p4)
fi
ran=1
# two eight-rank cases at a time
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
for c in local-stable-on local-convective-on edmf-convective-on edmfdd-convective-on; do
    has_log "$c" "PBL gas transport enabled = T" || fail "$c: PBL gas transport not logged as enabled"
    has_log "$c" "MYNN gas operator = local diffusion" || fail "$c: MYNN gas operator not logged"
done
has_log edmf-convective-on "+ updraft EDMF" || fail "edmf-convective-on: updraft EDMF not logged"
has_log edmfdd-convective-on "+ downdraft EDMF" || fail "edmfdd-convective-on: downdraft EDMF not logged"
has_log local-convective-on "+ updraft EDMF" && fail "local-convective-on: updraft EDMF logged with edmf = 0"
has_log edmf-convective-on "+ downdraft EDMF" && fail "edmf-convective-on: downdraft EDMF logged with edmf_dd = 0"
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
    echo "== updraft EDMF (MPAS default), convective surface layer"
    plumes() {
    python - "$1" <<'PY'
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
print(f"INFO {sys.argv[1].rsplit('/', 1)[-1]}: MYNN plume area peak {peak:.3g}, "
      f"columns with plumes in some record: {len(columns)}")
PY
    }
    plumes "$run_root/edmf-convective-on"
    python "$here/check_pbl_ysu.py" mixing "$run_root/edmf-convective-on" "$run_root/edmf-convective-off" \
        --json "$run_root/edmf-convective.json" || status=1
    echo "== updraft and downdraft EDMF, convective surface layer"
    plumes "$run_root/edmfdd-convective-on"
    triggers=$(grep -c "qstar=" "$run_root/edmfdd-convective-on/run.out" || true)
    echo "INFO edmfdd-convective-on: MYNN's downdraft scheme found a stratocumulus top $triggers times (all ranks)"
    python "$here/check_pbl_ysu.py" mixing "$run_root/edmfdd-convective-on" "$run_root/edmfdd-convective-off" \
        --json "$run_root/edmfdd-convective.json" || status=1
    echo "== transport off is unchanged from the pre-P4 executable"
    if [[ -n "$pre_p4_model" ]]; then
        python "$here/check_pbl_ysu.py" compare "$run_root/local-stable-off" "$run_root/local-stable-off-pre-p4" \
            || status=1
        python "$here/check_pbl_ysu.py" compare "$run_root/edmfdd-convective-off" \
            "$run_root/edmfdd-convective-off-pre-p4" || status=1
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
